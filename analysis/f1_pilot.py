"""N1 pilot: F1 (richness vs pH, non-host-associated samples) under four weighting schemes.

Follows N1_protocol.md v0.2. Run from the N1_sampling_imbalance folder:
    .venv/bin/python analysis/f1_pilot.py
Writes results/f1_pilot.json, results/f1_pilot.png and prints a summary.
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
RNG = np.random.default_rng(20261121)
N_BOOT = 1000
BIN = 0.5          # pH units per bin for the upper bound
MIN_BIN_N = 10     # samples needed in a bin
Q = 0.99           # upper-bound quantile
MIN_COUNTRY = 30   # samples needed for a held-out country


def richness(biom_path):
    with h5py.File(biom_path, "r") as f:
        ids = [x.decode() for x in f["sample/ids"][:]]
        indptr = f["sample/matrix/indptr"][:]
    return pd.Series(np.diff(indptr), index=ids, name="richness")


def load(table="emp_deblur_90bp.qc_filtered.rare_5000.biom"):
    m = pd.read_csv(DATA / "emp_qiime_mapping_qc_filtered.tsv", sep="\t", dtype=str, encoding_errors="replace")
    m = m.set_index("#SampleID")
    m["ph"] = pd.to_numeric(m["ph"], errors="coerce")
    m["lat"] = pd.to_numeric(m["latitude_deg"], errors="coerce")
    m["lon"] = pd.to_numeric(m["longitude_deg"], errors="coerce")
    m["country"] = m["country"].str.replace("GAZ:", "", regex=False)
    r = richness(DATA / table)
    d = m.join(r, how="inner")
    d = d[d["ph"].between(0, 14) & (d["empo_1"] == "Free-living")].copy()
    d["cell"] = (np.floor(d["lat"] / 5).astype("Int64").astype(str) + "_" +
                 np.floor(d["lon"] / 5).astype("Int64").astype(str))
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


def upper_bound_peak(d, w):
    edges = np.arange(np.floor(d.ph.min()), np.ceil(d.ph.max()) + BIN, BIN)
    idx = np.digitize(d.ph.to_numpy(), edges)
    xs, ys = [], []
    for k in np.unique(idx):
        s = idx == k
        if s.sum() >= MIN_BIN_N:
            xs.append(d.ph.to_numpy()[s].mean())
            ys.append(wquantile(d.richness.to_numpy()[s].astype(float), w[s], Q))
    xs, ys = np.array(xs), np.array(ys)
    p, _ = curve_fit(laplace, xs, ys, p0=[ys.max(), xs[np.argmax(ys)], 2.0], maxfev=20000)
    return {"peak": float(p[1]), "width": float(abs(p[2])), "max": float(p[0])}


def mean_quadratic(d, w):
    x, y = d.ph.to_numpy(), d.richness.to_numpy().astype(float)
    X = np.column_stack([np.ones_like(x), x, x ** 2])
    W = np.sqrt(w)
    b = np.linalg.lstsq(X * W[:, None], y * W, rcond=None)[0]
    return {"peak": float(-b[1] / (2 * b[2])), "curvature": float(b[2])}


def estimate(d, scheme):
    w = weights(d, scheme)
    return {**{f"ub_{k}": v for k, v in upper_bound_peak(d, w).items()},
            **{f"mean_{k}": v for k, v in mean_quadratic(d, w).items()}}


def bootstrap(d, scheme):
    studies = d.study_id.unique()
    groups = {s: d[d.study_id == s] for s in studies}
    rows = []
    for _ in range(N_BOOT):
        pick = RNG.choice(studies, size=len(studies), replace=True)
        b = pd.concat([groups[s] for s in pick], ignore_index=True)
        try:
            rows.append(estimate(b, scheme))
        except Exception:
            continue
    b = pd.DataFrame(rows)
    return {k: [float(b[k].quantile(0.025)), float(b[k].quantile(0.975))] for k in b}, len(b)


def transfer(d):
    """Leave-one-country-out vs 10-fold CV error for the S0 quadratic mean model."""
    def fit_pred(tr, te):
        x, y = tr.ph.to_numpy(), tr.richness.to_numpy().astype(float)
        b = np.polyfit(x, y, 2)
        return np.polyval(b, te.ph.to_numpy())
    rmse = lambda a, b: float(np.sqrt(np.mean((a - b) ** 2)))
    folds = RNG.integers(0, 10, len(d))
    pred = np.empty(len(d))
    for k in range(10):
        pred[folds == k] = fit_pred(d[folds != k], d[folds == k])
    cv_by_country = {}
    out = []
    for c, g in d.groupby("country"):
        if len(g) < MIN_COUNTRY:
            continue
        held = rmse(g.richness.to_numpy(), fit_pred(d[d.country != c], g))
        cv = rmse(g.richness.to_numpy(), pred[d.index.get_indexer(g.index)])
        out.append({"country": c, "n": len(g), "loco_rmse": held, "cv_rmse": cv, "ratio": held / cv})
    return out


def main():
    d = load()
    summary = {
        "n_samples": len(d), "n_studies": int(d.study_id.nunique()), "n_countries": int(d.country.nunique()),
        "top_countries": d.country.value_counts().head(6).to_dict(),
    }
    res = {}
    for s in ["S0", "S1", "S2", "S3"]:
        est = estimate(d, s)
        ci, nb = bootstrap(d, s)
        res[s] = {"estimate": est, "ci95": ci, "n_boot_ok": nb}
    base = res["S0"]["estimate"]
    decisions = {}
    for s in ["S1", "S2", "S3"]:
        e, ci = res[s]["estimate"], res[s]["ci95"]
        decisions[s] = {
            "ub_peak_shift": e["ub_peak"] - base["ub_peak"],
            "mean_peak_shift": e["mean_peak"] - base["mean_peak"],
            "curvature_change_pct": 100 * (e["mean_curvature"] - base["mean_curvature"]) / abs(base["mean_curvature"]),
            "curvature_ci_includes_zero": ci["mean_curvature"][0] <= 0 <= ci["mean_curvature"][1],
        }
        dd = decisions[s]
        dd["changed"] = bool(abs(dd["ub_peak_shift"]) > 0.5 or abs(dd["mean_peak_shift"]) > 0.5
                             or dd["curvature_ci_includes_zero"] or abs(dd["curvature_change_pct"]) > 50)
    tr = transfer(d)
    fails = sum(t["ratio"] > 1.5 for t in tr)
    sens = {}
    soil = d[d.empo_3 == "Soil (non-saline)"]
    sens["soil_only_S0"] = estimate(soil, "S0") | {"n": len(soil)}
    sens["soil_only_S1"] = estimate(soil, "S1")
    d150 = load("emp_deblur_150bp.qc_filtered.rare_5000.biom")
    sens["150bp_S0"] = estimate(d150, "S0") | {"n": len(d150)}
    sens["150bp_S1"] = estimate(d150, "S1")
    out = {"summary": summary, "schemes": res, "decisions": decisions,
           "transfer": {"countries": tr, "n_fail_ratio_gt_1.5": fails, "n_tested": len(tr),
                        "transfer_fails": bool(len(tr) and fails >= len(tr) / 3)},
           "sensitivity": sens}
    (OUT / "f1_pilot.json").write_text(json.dumps(out, indent=2))
    plot(out)
    print(json.dumps(out, indent=2))


def plot(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    labels = {"S0": "S0 standard", "S1": "S1 study-balanced", "S2": "S2 country-balanced", "S3": "S3 grid-balanced"}
    for j, (key, name) in enumerate([("ub_peak", "upper bound"), ("mean_peak", "mean")]):
        for i, s in enumerate(labels):
            e, ci = out["schemes"][s]["estimate"][key], out["schemes"][s]["ci95"][key]
            y = i + (0.18 if j else -0.18)
            ax[0].errorbar(e, y, xerr=[[e - ci[0]], [ci[1] - e]], fmt="o" if j == 0 else "s",
                           color="#1f5f8b" if j == 0 else "#b8661a", label=name if i == 0 else None, capsize=3)
    ax[0].set_yticks(range(4), list(labels.values()))
    ax[0].invert_yaxis()
    ax[0].set_xlim(0, 14)
    ax[0].set_xlabel("pH at peak richness (95% CI, study bootstrap; clipped to 0-14)")
    ax[0].legend(frameon=False)
    ax[0].set_title("F1: where richness peaks")
    t = sorted(out["transfer"]["countries"], key=lambda r: r["ratio"])
    ax[1].barh([f"{r['country']} (n={r['n']})" for r in t], [r["ratio"] for r in t], color="#5d7a1f")
    ax[1].axvline(1.5, color="#b3261e", ls="--", lw=1)
    ax[1].set_xlabel("Held-out error ÷ cross-validation error")
    ax[1].set_title("Leave-one-country-out")
    fig.tight_layout()
    fig.savefig(OUT / "f1_pilot.png", dpi=150)


if __name__ == "__main__":
    main()
