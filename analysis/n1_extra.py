"""N1 protocol v0.6, items 4-5: F9 deep-dive and transfer tests for F5, F9, F10.

Run from sampling-imbalance-biogeography-rules:  .venv/bin/python analysis/n1_extra.py emp|ma1000|ma5000
Writes results/final_extra_<which>.json
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).resolve().parent))
import all_findings as A  # noqa: E402
import ma_findings as MA  # noqa: E402

OUT = A.OUT
N_BOOT = 200


def boot_slope(ds, idx, x, y, scheme, rng):
    """Slope estimate and 95% percentile CI with study bootstrap (200)."""
    g = {k: MA.codes(pd.Series(v[idx])) for k, v in ds.g.items()}
    nproj = g["proj"].max() + 1
    est = A.wslope(x, y, MA.weights(g, scheme))
    bs = []
    for _ in range(N_BOOT):
        k = np.bincount(rng.integers(0, nproj, nproj), minlength=nproj).astype(float)
        w = MA.weights(g, scheme, k)
        if w.sum() > 0:
            bs.append(A.wslope(x, y, w))
    lo, hi = np.quantile(bs, [.025, .975])
    return {"slope": float(est), "ci": [float(lo), float(hi)], "significant": bool(not (lo <= 0 <= hi))}


def acido_classes(ds, which):
    """Class label per OTU index for Acidobacteria OTUs (others None)."""
    if which == "emp":
        with h5py.File(A.EMP / "emp_deblur_90bp.qc_filtered.rare_5000.biom", "r") as f:
            cls = np.array([x.decode() for x in f["observation/metadata/taxonomy"][:, 2]])
        return np.where(ds.acido, cls, None)
    otu_ids = np.load(A.MAD / "derived" / "otu_ids.npy")
    tax = pd.read_csv(A.MAD / "taxonomic_mapping.tsv", sep="\t", usecols=[1, 2], dtype=str).drop_duplicates("OTU_ID_MAPv3").set_index("OTU_ID_MAPv3").Taxonomy
    t = pd.Series(otu_ids).map(tax).fillna("")
    cls = t.str.split(";").str[2].fillna("unclassified")
    return np.where(ds.acido, cls.to_numpy(), None)


def logit(s):
    s = np.clip(s, 1e-3, 1 - 1e-3)
    return np.log(s / (1 - s))


def f9_deep(ds, which, rng):
    m = ds.meta
    idx = np.where(m.soil.to_numpy() & m.ph.between(0, 14).to_numpy())[0]
    x = m.ph.to_numpy()[idx]
    counts = [ds.count[ds.indptr[i]:ds.indptr[i + 1]] for i in idx]
    rows = [ds.row(i) for i in idx]
    share = np.array([c[ds.acido[r]].sum() for c, r in zip(counts, rows)]) / ds.depth
    y = logit(share)
    out = {"n": int(len(idx)), "S0": boot_slope(ds, idx, x, y, "S0", rng), "S3": boot_slope(ds, idx, x, y, "S3", rng)}

    # influence: drop one country at a time
    country = m.country.to_numpy()[idx]
    infl = []
    for c, n in pd.Series(country).value_counts().items():
        if n < 30 or c == "no_coordinates":
            continue
        keep = country != c
        r = boot_slope(ds, idx[keep], x[keep], y[keep], "S0", rng)
        infl.append({"dropped": c, "n_dropped": int(n), **r})
    out["drop_one_country_S0"] = infl

    # finer groups: Acidobacteria classes with >= 1% mean relative abundance
    cls = acido_classes(ds, which)
    per = {}
    labels = sorted({v for v in cls if v is not None})
    shares = {lab: np.zeros(len(idx)) for lab in labels}
    for j, (c_, r_) in enumerate(zip(counts, rows)):
        lab = cls[r_]
        for l_ in set(v for v in lab if v is not None):
            shares[l_][j] = c_[lab == l_].sum() / ds.depth
    for lab, sh in shares.items():
        if sh.mean() >= 0.01:
            per[lab] = {"mean_share": float(sh.mean()), "S0": boot_slope(ds, idx, x, logit(sh), "S0", rng),
                        "S3": boot_slope(ds, idx, x, logit(sh), "S3", rng)}
    out["classes"] = per

    # mixed model: random intercept per study (unweighted)
    df = pd.DataFrame({"y": y, "ph": x, "study": m.study.to_numpy()[idx]})
    try:
        fit = smf.mixedlm("y ~ ph", df, groups=df["study"]).fit(method="lbfgs")
        ci = fit.conf_int().loc["ph"].tolist()
        out["mixed_model"] = {"slope": float(fit.params["ph"]), "ci": [float(ci[0]), float(ci[1])],
                              "significant": bool(not (ci[0] <= 0 <= ci[1]))}
    except Exception as e:
        out["mixed_model"] = {"error": str(e)}
    return out


def transfer(ds, idx, x, y, min_country):
    rmse = lambda a, b: float(np.sqrt(np.mean((a - b) ** 2)))
    rng = np.random.default_rng(A.SEED)
    fit_pred = lambda tr, te: np.polyval(np.polyfit(x[tr], y[tr], 1), x[te])
    n = len(idx)
    folds = rng.integers(0, 10, n)
    cv = np.empty(n)
    for k in range(10):
        cv[folds == k] = fit_pred(folds != k, folds == k)
    country = ds.meta.country.to_numpy()[idx]
    res = []
    for c, cnt in pd.Series(country).value_counts().items():
        if cnt < min_country or c == "no_coordinates":
            continue
        te = country == c
        res.append({"country": c, "n": int(cnt), "ratio": rmse(y[te], fit_pred(~te, te)) / rmse(y[te], cv[te])})
    fails = sum(r["ratio"] > 1.5 for r in res)
    return {"n_tested": len(res), "n_fail": fails, "transfer_fails": bool(res and fails >= len(res) / 3),
            "worst": sorted(res, key=lambda r: -r["ratio"])[:5]}


def transfers(ds, which):
    m = ds.meta
    mc = 30 if which == "emp" else 200
    out = {}
    soil = m.soil.to_numpy()
    lat = np.abs(m.lat.to_numpy())
    i5 = np.where(soil & ~np.isnan(lat))[0]
    out["F5"] = transfer(ds, i5, lat[i5], ds.rich[i5], mc)
    i10 = np.where(soil & ~np.isnan(m.elev.to_numpy()))[0]
    out["F10"] = transfer(ds, i10, m.elev.to_numpy()[i10], ds.rich[i10], mc)
    i9 = np.where(soil & m.ph.between(0, 14).to_numpy())[0]
    sh = np.array([ds.count[ds.indptr[i]:ds.indptr[i + 1]][ds.acido[ds.row(i)]].sum() for i in i9]) / ds.depth
    out["F9"] = transfer(ds, i9, m.ph.to_numpy()[i9], logit(sh), mc)
    return out


def main():
    which = sys.argv[1]
    ds = {"emp": A.load_emp, "ma1000": lambda: A.load_ma(1000), "ma5000": lambda: A.load_ma(5000)}[which]()
    rng = np.random.default_rng([A.SEED, 9, len(which)])
    res = {"f9_deep": f9_deep(ds, which, rng), "transfer": transfers(ds, which)}
    (OUT / f"final_extra_{which}.json").write_text(json.dumps(res, indent=1, default=float))
    d = res["f9_deep"]
    print(which, "F9 S0", d["S0"], "S3", d["S3"], "mixed", d["mixed_model"])
    print("  drop-one-country changes S0 significance:",
          [(r["dropped"], round(r["slope"], 3), r["significant"]) for r in d["drop_one_country_S0"] if r["significant"] != d["S0"]["significant"]])
    print("  classes:", {k: (round(v["mean_share"], 3), round(v["S0"]["slope"], 3), v["S0"]["significant"], round(v["S3"]["slope"], 3), v["S3"]["significant"]) for k, v in d["classes"].items()})
    print("  transfer:", {k: (v["n_fail"], v["n_tested"], v["transfer_fails"]) for k, v in res["transfer"].items()})


if __name__ == "__main__":
    main()
