"""N1: re-test EMP findings F1 (pH), F2 (temperature) and F3 (free-living vs host-associated)
under four weighting schemes, following N1_protocol.md v0.3.

Run from the N1_sampling_imbalance folder:
    .venv/bin/python analysis/run_findings.py
Writes results/findings_results.json and results/findings.png.
"""
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "emp_release1"
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
SEED = 20261121
N_BOOT = 1000
Q = 0.99
MIN_BIN_N = 10
MIN_COUNTRY = 30
SCHEMES = ["S0", "S1", "S2", "S3"]
COV = {  # finding -> covariate column, bin width, plausible range, shift threshold, paper n
    "F1": dict(col="ph", bin=0.5, lo=0, hi=14, thr=0.5, paper_n=3986, unit="pH"),
    "F2": dict(col="temperature_deg_c", bin=2.0, lo=-20, hi=120, thr=5.0, paper_n=6976, unit="°C"),
}


def richness(table):
    with h5py.File(DATA / table, "r") as f:
        ids = [x.decode() for x in f["sample/ids"][:]]
        indptr = f["sample/matrix/indptr"][:]
    return pd.Series(np.diff(indptr), index=ids, name="richness")


def load(table="emp_deblur_90bp.qc_filtered.rare_5000.biom"):
    m = pd.read_csv(DATA / "emp_qiime_mapping_qc_filtered.tsv", sep="\t", dtype=str,
                    encoding_errors="replace").set_index("#SampleID")
    for c in ["ph", "temperature_deg_c", "latitude_deg", "longitude_deg"]:
        m[c] = pd.to_numeric(m[c], errors="coerce")
    m["country"] = m["country"].str.replace("GAZ:", "", regex=False)
    d = m.join(richness(table), how="inner")
    d["cell"] = (np.floor(d.latitude_deg / 5).astype("Int64").astype(str) + "_" +
                 np.floor(d.longitude_deg / 5).astype("Int64").astype(str)).fillna("no_coordinates")
    d["richness"] = d["richness"].astype(float)
    return d


def weights(d, scheme):
    if scheme == "S0":
        return np.ones(len(d))
    key = {"S1": "study_id", "S2": "country", "S3": "cell"}[scheme]
    return (1.0 / d.groupby(key)[key].transform("size")).to_numpy()


def wquantile(x, w, q):
    o = np.argsort(x)
    x, w = x[o], w[o]
    c = np.cumsum(w) - 0.5 * w
    return np.interp(q * w.sum(), c, x)


def laplace(x, a, mu, b):
    return a * np.exp(-np.abs(x - mu) / b)


def covariate_estimate(d, w, cfg, min_bin_n=MIN_BIN_N):
    x, y = d[cfg["col"]].to_numpy(), d.richness.to_numpy()
    edges = np.arange(np.floor(x.min()), np.ceil(x.max()) + cfg["bin"], cfg["bin"])
    idx = np.digitize(x, edges)
    xs, ys = [], []
    for k in np.unique(idx):
        s = idx == k
        if s.sum() >= min_bin_n:
            xs.append(x[s].mean())
            ys.append(wquantile(y[s], w[s], Q))
    xs, ys = np.array(xs), np.array(ys)
    p, _ = curve_fit(laplace, xs, ys, p0=[ys.max(), xs[np.argmax(ys)], 3 * cfg["bin"]], maxfev=20000)
    X = np.column_stack([np.ones_like(x), x, x ** 2])
    sw = np.sqrt(w)
    b = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
    return {"ub_peak": float(p[1]), "ub_width": float(abs(p[2])),
            "mean_peak": float(-b[1] / (2 * b[2])), "mean_curvature": float(b[2])}


def f3_estimate(d, w):
    fl = d.empo_1.to_numpy() == "Free-living"
    ha = d.empo_1.to_numpy() == "Host-associated"
    y = d.richness.to_numpy()
    m_fl = np.average(y[fl], weights=w[fl])
    m_ha = np.average(y[ha], weights=w[ha])
    return {"log_ratio_fl_ha": float(np.log(m_fl / m_ha)), "mean_free_living": float(m_fl),
            "mean_host_associated": float(m_ha)}


def run_schemes(d, fn, rng):
    studies = d.study_id.unique()
    groups = {s: d[d.study_id == s] for s in studies}
    out = {}
    for s in SCHEMES:
        est = fn(d, weights(d, s))
        rows = []
        for _ in range(N_BOOT):
            b = pd.concat([groups[k] for k in rng.choice(studies, len(studies), replace=True)])
            try:
                rows.append(fn(b, weights(b, s)))
            except Exception:
                pass
        bt = pd.DataFrame(rows)
        out[s] = {"estimate": est, "n_boot_ok": len(bt),
                  "ci95": {k: [float(bt[k].quantile(.025)), float(bt[k].quantile(.975))] for k in bt}}
    return out


def transfer(d, col, rng):
    """Leave-one-country-out vs 10-fold CV: richness ~ habitat + covariate + covariate^2."""
    def design(df, habitats):
        H = np.column_stack([(df.empo_3.to_numpy() == h).astype(float) for h in habitats])
        if col is None:
            return H
        x = df[col].to_numpy()
        return np.column_stack([H, x, x ** 2])

    def fit_pred(tr, te):
        hab = sorted(tr.empo_3.unique())
        b = np.linalg.lstsq(design(tr, hab), tr.richness.to_numpy(), rcond=None)[0]
        Xt = design(te, hab)
        pred = Xt @ b
        unseen = Xt[:, :len(hab)].sum(axis=1) == 0
        pred[unseen] += b[:len(hab)].mean()  # unseen habitat: average habitat effect
        return pred

    rmse = lambda a, b: float(np.sqrt(np.mean((a - b) ** 2)))
    folds = rng.integers(0, 10, len(d))
    cv = np.empty(len(d))
    for k in range(10):
        cv[folds == k] = fit_pred(d[folds != k], d[folds == k])
    rows = []
    for c, g in d.groupby("country"):
        if len(g) < MIN_COUNTRY:
            continue
        pos = d.index.get_indexer(g.index)
        held = rmse(g.richness.to_numpy(), fit_pred(d[d.country != c], g))
        cvr = rmse(g.richness.to_numpy(), cv[pos])
        rows.append({"country": c, "n": int(len(g)), "loco_rmse": held, "cv_rmse": cvr, "ratio": held / cvr})
    fails = sum(r["ratio"] > 1.5 for r in rows)
    return {"countries": rows, "n_tested": len(rows), "n_fail": fails,
            "transfer_fails": bool(rows and fails >= len(rows) / 3)}


def decide_cov(res, thr):
    base = res["S0"]["estimate"]
    out = {}
    for s in SCHEMES[1:]:
        e, ci = res[s]["estimate"], res[s]["ci95"]
        claim_shift = e["ub_peak"] - base["ub_peak"]
        mean_shift = e["mean_peak"] - base["mean_peak"]
        curv_pct = 100 * (e["mean_curvature"] - base["mean_curvature"]) / abs(base["mean_curvature"])
        curv_zero = ci["mean_curvature"][0] <= 0 <= ci["mean_curvature"][1]
        out[s] = {"claim_peak_shift": claim_shift, "claim_changed": bool(abs(claim_shift) > thr),
                  "mean_peak_shift": mean_shift, "curvature_change_pct": curv_pct,
                  "curvature_ci_includes_zero": bool(curv_zero),
                  "average_changed": bool(abs(mean_shift) > thr or curv_zero or abs(curv_pct) > 50)}
    return out


def decide_f3(res):
    base = res["S0"]["estimate"]["log_ratio_fl_ha"]
    out = {}
    for s in SCHEMES[1:]:
        e = res[s]["estimate"]["log_ratio_fl_ha"]
        lo, hi = res[s]["ci95"]["log_ratio_fl_ha"]
        pct = 100 * (e - base) / abs(base)
        out[s] = {"log_ratio": e, "change_pct": pct, "sign_flip": bool(np.sign(e) != np.sign(base)),
                  "ci_includes_zero": bool(lo <= 0 <= hi)}
        out[s]["claim_changed"] = bool(out[s]["sign_flip"] or out[s]["ci_includes_zero"] or abs(pct) > 50)
    return out


def describe(d):
    return {"n": int(len(d)), "studies": int(d.study_id.nunique()), "countries": int(d.country.nunique()),
            "top_countries": {k: int(v) for k, v in d.country.value_counts().head(5).items()}}


def main():
    rng = np.random.default_rng(SEED)
    d = load()
    d150 = load("emp_deblur_150bp.qc_filtered.rare_5000.biom")
    results = {}
    for f, cfg in COV.items():
        col = cfg["col"]
        sub = d[d[col].between(cfg["lo"], cfg["hi"])]
        fn = lambda df, w, cfg=cfg: covariate_estimate(df, w, cfg)
        res = run_schemes(sub, fn, rng)
        fl = sub[sub.empo_1 == "Free-living"]
        s150 = d150[d150[col].between(cfg["lo"], cfg["hi"])]
        sens = {
            "free_living_only": {s: fn(fl, weights(fl, s)) for s in ["S0", "S3"]} | {"n": int(len(fl))},
            "150bp": {s: fn(s150, weights(s150, s)) for s in ["S0", "S3"]} | {"n": int(len(s150))},
        }
        # Post hoc (added after the first F2 run): bins need at least 20 samples
        sens["posthoc_min_bin_20"] = {s: covariate_estimate(sub, weights(sub, s), cfg, 20) for s in SCHEMES}
        if f == "F1":
            soil = sub[sub.empo_3 == "Soil (non-saline)"]
            sens["soil_only"] = {s: fn(soil, weights(soil, s)) for s in ["S0", "S3"]} | {"n": int(len(soil))}
        results[f] = {"data": describe(sub) | {"paper_n": cfg["paper_n"]}, "schemes": res,
                      "decisions": decide_cov(res, cfg["thr"]), "transfer": transfer(sub, col, rng),
                      "sensitivity": sens}
    f3d = d[d.empo_1.isin(["Free-living", "Host-associated"])]
    res = run_schemes(f3d, f3_estimate, rng)
    by_hab = {s: f3d.assign(w=weights(f3d, s)).groupby("empo_3").apply(
        lambda g: float(np.average(g.richness, weights=g.w)), include_groups=False).round(0).to_dict()
        for s in ["S0", "S3"]}
    results["F3"] = {"data": describe(f3d) | {"paper_n": 23828}, "schemes": res, "decisions": decide_f3(res),
                     "transfer": transfer(f3d, None, rng), "mean_richness_by_habitat": by_hab,
                     "sensitivity": {"150bp": {s: f3_estimate(x, weights(x, s)) for s in ["S0", "S3"]
                                               for x in [d150[d150.empo_1.isin(["Free-living", "Host-associated"])]]}}}
    (OUT / "findings_results.json").write_text(json.dumps(results, indent=2, ensure_ascii=False))
    plot(results)
    summary(results)


def summary(r):
    for f in ["F1", "F2"]:
        print(f, r[f]["data"])
        for s in SCHEMES:
            e, ci = r[f]["schemes"][s]["estimate"], r[f]["schemes"][s]["ci95"]
            print(f"  {s} ub_peak {e['ub_peak']:.2f} [{ci['ub_peak'][0]:.2f},{ci['ub_peak'][1]:.2f}]  "
                  f"mean_peak {e['mean_peak']:.2f} [{ci['mean_peak'][0]:.2f},{ci['mean_peak'][1]:.2f}]  "
                  f"curv {e['mean_curvature']:.3f} [{ci['mean_curvature'][0]:.3f},{ci['mean_curvature'][1]:.3f}]")
        for s, v in r[f]["decisions"].items():
            print("  ", s, {k: (round(x, 2) if isinstance(x, float) else x) for k, x in v.items()})
        t = r[f]["transfer"]
        print("   transfer:", t["n_fail"], "of", t["n_tested"], "fail;",
              [(c["country"], round(c["ratio"], 2)) for c in t["countries"]])
        print("   sensitivity:", json.dumps(r[f]["sensitivity"])[:600])
    print("F3", r["F3"]["data"])
    for s in SCHEMES:
        e, ci = r["F3"]["schemes"][s]["estimate"], r["F3"]["schemes"][s]["ci95"]
        print(f"  {s} log ratio {e['log_ratio_fl_ha']:.3f} [{ci['log_ratio_fl_ha'][0]:.3f},{ci['log_ratio_fl_ha'][1]:.3f}] "
              f"FL {e['mean_free_living']:.0f} HA {e['mean_host_associated']:.0f}")
    print("  ", r["F3"]["decisions"])
    t = r["F3"]["transfer"]
    print("   transfer:", t["n_fail"], "of", t["n_tested"], "fail;", [(c["country"], round(c["ratio"], 2)) for c in t["countries"]])
    print("   by habitat:", r["F3"]["mean_richness_by_habitat"])
    print("   150bp:", r["F3"]["sensitivity"])


def plot(r):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 4, figsize=(17, 4.2))
    names = ["S0 standard", "S1 study", "S2 country", "S3 grid"]
    for j, (f, lim, lab) in enumerate([("F1", (0, 14), "pH at peak"), ("F2", (-10, 75), "°C at peak")]):
        for i, s in enumerate(SCHEMES):
            for k, (key, mk, col) in enumerate([("ub_peak", "o", "#1f5f8b"), ("mean_peak", "s", "#b8661a")]):
                e, ci = r[f]["schemes"][s]["estimate"][key], r[f]["schemes"][s]["ci95"][key]
                ax[j].errorbar(e, i + (0.18 if k else -0.18), xerr=[[max(e - ci[0], 0)], [max(ci[1] - e, 0)]],
                               fmt=mk, color=col, capsize=3,
                               label=("upper bound (claim)" if k == 0 else "mean") if i == 0 else None)
        ax[j].set_xlim(*lim)
        ax[j].set_yticks(range(4), names)
        ax[j].invert_yaxis()
        ax[j].set_xlabel(f"{lab} (95% CI; axis clipped)")
        ax[j].set_title(f"{f}: {'richness vs pH' if f == 'F1' else 'richness vs temperature'}")
    ax[0].legend(frameon=False, fontsize=8, loc="lower left")
    for i, s in enumerate(SCHEMES):
        e, ci = r["F3"]["schemes"][s]["estimate"]["log_ratio_fl_ha"], r["F3"]["schemes"][s]["ci95"]["log_ratio_fl_ha"]
        ax[2].errorbar(e, i, xerr=[[e - ci[0]], [ci[1] - e]], fmt="o", color="#5d7a1f", capsize=3)
    ax[2].axvline(0, color="#888", lw=1)
    ax[2].set_yticks(range(4), names)
    ax[2].invert_yaxis()
    ax[2].set_xlabel("log(free-living ÷ host-associated richness)")
    ax[2].set_title("F3: free-living vs host")
    for k, (f, col) in enumerate([("F1", "#1f5f8b"), ("F2", "#b8661a"), ("F3", "#5d7a1f")]):
        t = r[f]["transfer"]["countries"]
        ax[3].scatter([c["ratio"] for c in t], [k] * len(t), color=col, alpha=.7)
    ax[3].axvline(1.5, color="#b3261e", ls="--", lw=1)
    ax[3].set_yticks(range(3), ["F1", "F2", "F3"])
    ax[3].invert_yaxis()
    ax[3].set_xlabel("Held-out ÷ CV error, per country")
    ax[3].set_title("Leave-one-country-out")
    fig.tight_layout()
    fig.savefig(OUT / "findings.png", dpi=150)


if __name__ == "__main__":
    main()
