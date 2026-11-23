"""N1 MicrobeAtlas extension, step 2: re-test F1-F3 under S0-S3 (protocol v0.4).

Needs data/microbeatlas/derived/meta.pkl (ma_prepare.py) and the decompressed OTU table
data/microbeatlas/samples-otus.97.metag.minfilter.biom (HDF5).
Run from sampling-imbalance-biogeography-rules:  .venv/bin/python analysis/ma_findings.py
Writes results/ma_findings.json, results/ma_findings.png, results/ma_richness.pkl (cache).
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "data" / "microbeatlas"
OUT = ROOT / "results"
SEED = 20261121
DEPTH = 1000
N_BOOT = 200
Q = 0.99
MIN_BIN_N, MIN_BIN_PROJ = 20, 3
MIN_COUNTRY = 200
SCHEMES = ["S0", "S1", "S2", "S3"]
COV = {
    "F1": dict(col="ph", bin=0.5, lo=0, hi=14, thr=0.5, emp=7.37),
    "F2": dict(col="temp", bin=2.0, lo=-20, hi=120, thr=5.0, emp=9.05),
}


# ---------- richness ----------
def rarefied_richness(ids_wanted):
    rng = np.random.default_rng(SEED)
    with h5py.File(D / "samples-otus.97.metag.minfilter.biom", "r") as f:
        ids = np.array([x.decode() for x in f["sample/ids"][:]])
        indptr = f["sample/matrix/indptr"][:]
        data = f["sample/matrix/data"][:]  # ~2.6 GB in memory; much faster than per-sample disk reads
        pos = {s: i for i, s in enumerate(ids)}
        out = {}
        todo = sorted(pos[s] for s in ids_wanted if s in pos)
        for n, i in enumerate(todo):
            c = data[indptr[i]:indptr[i + 1]].astype(np.int64)
            if c.sum() < DEPTH:
                continue
            out[ids[i]] = int((rng.multivariate_hypergeometric(c, DEPTH) > 0).sum())
            if n % 200000 == 0:
                print(f"  rarefied {n}/{len(todo)}", flush=True)
    return pd.Series(out, name="richness")


def load():
    meta = pd.read_pickle(D / "derived" / "meta.pkl")
    cache = OUT / "ma_richness.pkl"
    with h5py.File(D / "samples-otus.97.metag.minfilter.biom", "r") as f:
        table_ids = set(x.decode() for x in f["sample/ids"][:])
    key = "map_sid" if len(table_ids & set(meta.map_sid.head(200000))) else "sample"
    print("table id key:", key, "| table samples:", len(table_ids))
    if cache.exists():
        r = pd.read_pickle(cache)
    else:
        r = rarefied_richness(set(meta[key]) & table_ids)
        r.to_pickle(cache)
    d = meta.merge(r, left_on=key, right_index=True, how="inner")
    return d.drop_duplicates(key)


# ---------- weights with project-level bootstrap multiplicities ----------
def codes(s):
    return pd.factorize(s)[0]


def weights(g, scheme, k=None):
    """g: dict of integer code arrays (proj, country, cell); k: copies of each project (None = 1)."""
    p = g["proj"]
    kp = np.ones(p.max() + 1) if k is None else k
    kk = kp[p]
    if scheme == "S0":
        return kk
    if scheme == "S1":
        return (kk > 0) / np.bincount(p)[p]
    key = g["country"] if scheme == "S2" else g["cell"]
    tot = np.bincount(key, weights=kk)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(kk > 0, kk / tot[key], 0.0)


def wquantile(x, w, q):
    o = np.argsort(x)
    x, w = x[o], w[o]
    c = np.cumsum(w) - 0.5 * w
    return np.interp(q * w.sum(), c, x)


def laplace(x, a, mu, b):
    return a * np.exp(-np.abs(x - mu) / b)


def cov_est(x, y, proj, w, cfg):
    keep = w > 0
    x, y, proj, w = x[keep], y[keep], proj[keep], w[keep]
    edges = np.arange(np.floor(x.min()), np.ceil(x.max()) + cfg["bin"], cfg["bin"])
    idx = np.digitize(x, edges)
    xs, ys, y95 = [], [], []
    for b in np.unique(idx):
        s = idx == b
        if s.sum() >= MIN_BIN_N and len(np.unique(proj[s])) >= MIN_BIN_PROJ:
            xs.append(x[s].mean())
            ys.append(wquantile(y[s], w[s], Q))
            y95.append(wquantile(y[s], w[s], 0.95))
    xs, ys, y95 = np.array(xs), np.array(ys), np.array(y95)
    p, _ = curve_fit(laplace, xs, ys, p0=[ys.max(), xs[np.argmax(ys)], 3 * cfg["bin"]], maxfev=20000)
    X = np.column_stack([np.ones_like(x), x, x ** 2])
    sw = np.sqrt(w)
    b = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
    smooth = lambda v: np.convolve(v, np.ones(3) / 3, mode="same") if len(v) >= 3 else v  # v0.6 steadier peak
    return {"ub_peak": float(p[1]), "ub_width": float(abs(p[2])),
            "sm99_peak": float(xs[np.argmax(smooth(ys))]), "sm95_peak": float(xs[np.argmax(smooth(y95))]),
            "mean_peak": float(-b[1] / (2 * b[2])), "mean_curvature": float(b[2])}


def f3_est(y, ha, w):
    m_fl = np.sum(w[~ha] * y[~ha]) / np.sum(w[~ha])
    m_ha = np.sum(w[ha] * y[ha]) / np.sum(w[ha])
    return {"log_ratio_fl_ha": float(np.log(m_fl / m_ha)), "mean_free_living": float(m_fl),
            "mean_host_associated": float(m_ha)}


def run(d, fn, rng):
    g = {"proj": codes(d.project), "country": codes(d.country), "cell": codes(d.cell)}
    nproj = g["proj"].max() + 1
    res = {}
    for s in SCHEMES:
        est = fn(weights(g, s))
        rows = []
        for _ in range(N_BOOT):
            k = np.bincount(rng.integers(0, nproj, nproj), minlength=nproj).astype(float)
            try:
                rows.append(fn(weights(g, s, k)))
            except Exception:
                pass
        bt = pd.DataFrame(rows)
        res[s] = {"estimate": est, "n_boot_ok": len(bt),
                  "ci95": {c: [float(bt[c].quantile(.025)), float(bt[c].quantile(.975))] for c in bt}}
    return res


def decide_cov(res, thr):
    base = res["S0"]["estimate"]
    out = {}
    for s in SCHEMES[1:]:
        e, ci = res[s]["estimate"], res[s]["ci95"]
        cs = e["ub_peak"] - base["ub_peak"]
        ms = e["mean_peak"] - base["mean_peak"]
        cp = 100 * (e["mean_curvature"] - base["mean_curvature"]) / abs(base["mean_curvature"])
        cz = ci["mean_curvature"][0] <= 0 <= ci["mean_curvature"][1]
        out[s] = {"claim_peak_shift": cs, "claim_changed": bool(abs(cs) > thr), "mean_peak_shift": ms,
                  "curvature_change_pct": cp, "curvature_ci_includes_zero": bool(cz),
                  "average_changed": bool(abs(ms) > thr or cz or abs(cp) > 50)}
    return out


def decide_f3(res):
    base = res["S0"]["estimate"]["log_ratio_fl_ha"]
    out = {}
    for s in SCHEMES[1:]:
        e = res[s]["estimate"]["log_ratio_fl_ha"]
        lo, hi = res[s]["ci95"]["log_ratio_fl_ha"]
        pct = 100 * (e - base) / abs(base)
        out[s] = {"log_ratio": e, "change_pct": pct, "sign_flip": bool(np.sign(e) != np.sign(base)),
                  "ci_includes_zero": bool(lo <= 0 <= hi),
                  "claim_changed": bool(np.sign(e) != np.sign(base) or lo <= 0 <= hi or abs(pct) > 50)}
    return out


def transfer(d, col, rng):
    hab = d.env_2.where(d.env_2.isin(d.env_2.value_counts().head(40).index), "other")
    y = d.richness.to_numpy().astype(float)

    def fit_pred(tr, te):
        hm = pd.Series(y[tr]).groupby(hab.to_numpy()[tr]).mean()
        base_tr = hab.iloc[tr].map(hm).to_numpy()
        base_te = hab.iloc[te].map(hm).fillna(hm.mean()).to_numpy()
        if col is None:
            return base_te
        x = d[col].to_numpy()
        X = np.column_stack([x[tr], x[tr] ** 2, np.ones(len(tr))])
        b = np.linalg.lstsq(X, y[tr] - base_tr, rcond=None)[0]
        return base_te + np.column_stack([x[te], x[te] ** 2, np.ones(len(te))]) @ b

    rmse = lambda a, b: float(np.sqrt(np.mean((a - b) ** 2)))
    n = len(d)
    folds = rng.integers(0, 10, n)
    cv = np.empty(n)
    for k in range(10):
        cv[folds == k] = fit_pred(np.where(folds != k)[0], np.where(folds == k)[0])
    rows = []
    country = d.country.to_numpy()
    for c, cnt in d.country.value_counts().items():
        if cnt < MIN_COUNTRY or c == "no_coordinates":
            continue
        te = np.where(country == c)[0]
        tr = np.where(country != c)[0]
        held = rmse(y[te], fit_pred(tr, te))
        rows.append({"country": c, "n": int(cnt), "ratio": held / rmse(y[te], cv[te])})
    fails = sum(r["ratio"] > 1.5 for r in rows)
    return {"countries": rows, "n_tested": len(rows), "n_fail": fails,
            "transfer_fails": bool(rows and fails >= len(rows) / 3)}


def describe(d):
    vc = d.country.value_counts()
    coord = vc.drop("no_coordinates", errors="ignore")
    return {"n": int(len(d)), "projects": int(d.project.nunique()), "countries": int(len(coord)),
            "top_country_share": float(coord.iloc[0] / coord.sum()), "top3_share": float(coord.head(3).sum() / coord.sum()),
            "top_countries": {k: int(v) for k, v in coord.head(6).items()}}


def main():
    rng = np.random.default_rng(SEED)
    d = load()
    print("loaded", len(d), "samples with richness")
    amp = d[d.technology == "AMPLICON"]
    wgs = d[d.technology == "WGS"]
    results = {}
    for f, cfg in COV.items():
        sub = amp[amp[cfg["col"]].between(cfg["lo"], cfg["hi"])].reset_index(drop=True)
        x, y, pr = sub[cfg["col"]].to_numpy(), sub.richness.to_numpy().astype(float), codes(sub.project)
        fn = lambda w, x=x, y=y, pr=pr, cfg=cfg: cov_est(x, y, pr, w, cfg)
        res = run(sub, fn, rng)
        sens = {}
        for name, part in [("free_living_only", sub[~sub.host_associated]),
                           ("wgs", wgs[wgs[cfg["col"]].between(cfg["lo"], cfg["hi"])])]:
            part = part.reset_index(drop=True)
            if len(part) < 200:
                sens[name] = {"n": int(len(part))}
                continue
            gp = {"proj": codes(part.project), "country": codes(part.country), "cell": codes(part.cell)}
            sens[name] = {"n": int(len(part))}
            for s in ["S0", "S3"]:
                try:
                    sens[name][s] = cov_est(part[cfg["col"]].to_numpy(), part.richness.to_numpy().astype(float),
                                            gp["proj"], weights(gp, s), cfg)
                except Exception as e:
                    sens[name][s] = {"error": str(e)}
        results[f] = {"data": describe(sub), "emp_S0_ub_peak": cfg["emp"], "schemes": res,
                      "decisions": decide_cov(res, cfg["thr"]), "transfer": transfer(sub, cfg["col"], rng),
                      "sensitivity": sens}
        print(f, "done", flush=True)
    f3 = amp[amp.env_top != "unknown"].reset_index(drop=True)
    ha = f3.host_associated.to_numpy()
    yy = f3.richness.to_numpy().astype(float)
    res = run(f3, lambda w: f3_est(yy, ha, w), rng)
    results["F3"] = {"data": describe(f3), "schemes": res, "decisions": decide_f3(res),
                     "transfer": transfer(f3, None, rng),
                     "mean_by_env_top_S0": f3.groupby("env_top").richness.mean().round(0).to_dict()}
    (OUT / "ma_findings.json").write_text(json.dumps(results, indent=2, ensure_ascii=False))
    print(json.dumps({f: {"data": r["data"], "S": {s: r["schemes"][s]["estimate"] for s in SCHEMES},
                          "decisions": r["decisions"], "transfer_fail": [r["transfer"]["n_fail"], r["transfer"]["n_tested"]]}
                      for f, r in results.items()}, indent=1))


if __name__ == "__main__":
    main()
