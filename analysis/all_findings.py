"""N1: F4-F10 in EMP and MicrobeAtlas, F1-F3 at 5,000 reads in MicrobeAtlas, and one combined
F1-F10 table (protocol v0.5).

Needs: EMP release 1 files; MicrobeAtlas derived files from ma_prepare.py and ma_rarefy.py.
Run from sampling-imbalance-biogeography-rules:  .venv/bin/python analysis/all_findings.py
Protocol v0.6: writes results/final_{emp,ma1000,ma5000}.csv/.json and, with `merge`, results/final_findings.csv.
(The v0.5 outputs all_findings*.csv were produced by the previous version of this script; see git history.)
"""
import json
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ma_findings as MA  # noqa: E402  (reuse weights, cov_est, f3_est, codes)

ROOT = Path(__file__).resolve().parents[1]
EMP = ROOT / "data" / "emp_release1"
MAD = ROOT / "data" / "microbeatlas"
OUT = ROOT / "results"
SEED = 20261121
N_BOOT = 1000  # v0.6 (v0.5 used 200)
N_PAIRS = 100_000
SCHEMES = ["S0", "S1", "S2", "S3"]

HERB = ["cattle", "cow", "sheep", "goat", "horse", "rabbit", "deer", "buffalo", "yak", "camel", "elephant",
        "kangaroo", "wombat", "sambar", "langur", "snub-nosed", "proboscis", "guereza", "colobus", "sloth"]
CARN = ["cat", "dog", "fox", "wolf", "seal", "lion", "tiger", "ferret", "mink", "anteater", "aardvark", "armadillo"]


def diet_of(name):
    s = str(name).lower()
    if any(re.search(rf"\b{w}s?\b", s) for w in HERB):
        return "herbivore"
    if any(re.search(rf"\b{w}s?\b", s) for w in CARN):
        return "carnivore"
    return np.nan


# ---------------- datasets ----------------
def csr_take(indptr, rows):
    """Positions of the nonzeros of `rows` in a CSR array, plus the new indptr (vectorised)."""
    starts, lens = indptr[rows], indptr[rows + 1] - indptr[rows]
    new_ip = np.concatenate([[0], np.cumsum(lens)])
    pos = np.repeat(starts - new_ip[:-1], lens) + np.arange(new_ip[-1])
    return pos, new_ip


class Dataset:
    """meta rows align with CSR rows (indptr/otu/count)."""

    def __init__(self, name, depth, meta, indptr, otu, count, acido):
        self.name, self.depth, self.meta = name, depth, meta.reset_index(drop=True)
        self.indptr, self.otu, self.count, self.acido = indptr, otu, count, acido
        self.rich = np.diff(indptr).astype(float)
        m = self.meta
        self.g = {"proj": MA.codes(m.study), "country": MA.codes(m.country), "cell": MA.codes(m.cell)}

    def row(self, i):
        return self.otu[self.indptr[i]:self.indptr[i + 1]]

    def sub(self, mask):
        idx = np.where(mask)[0]
        return idx


def load_emp():
    m = pd.read_csv(EMP / "emp_qiime_mapping_qc_filtered.tsv", sep="\t", dtype=str,
                    encoding_errors="replace").set_index("#SampleID")
    with h5py.File(EMP / "emp_deblur_90bp.qc_filtered.rare_5000.biom", "r") as f:
        ids = np.array([x.decode() for x in f["sample/ids"][:]])
        indptr = f["sample/matrix/indptr"][:]
        otu = f["sample/matrix/indices"][:]
        cnt = f["sample/matrix/data"][:]
        phylum = np.array([x.decode() for x in f["observation/metadata/taxonomy"][:, 1]])
    m = m.loc[ids]
    num = lambda c: pd.to_numeric(m[c], errors="coerce").to_numpy()
    lat, lon = num("latitude_deg"), num("longitude_deg")
    cell = np.array([f"{int(np.floor(a / 5))}_{int(np.floor(o / 5))}" if not (np.isnan(a) or np.isnan(o))
                     else "no_coordinates" for a, o in zip(lat, lon)])
    fl = (m.empo_1 == "Free-living").to_numpy()
    meta = pd.DataFrame({
        "study": m.study_id.to_numpy(), "country": m.country.str.replace("GAZ:", "", regex=False).to_numpy(),
        "cell": cell, "lat": lat, "lon": lon, "ph": num("ph"), "temp": num("temperature_deg_c"),
        "elev": num("elevation_m"), "habitat": m.empo_3.to_numpy(), "free_living": fl,
        "salinity": np.where(fl, m.empo_2.map({"Saline": "saline", "Non-saline": "nonsaline"}), np.nan),
        "soil": (m.empo_3 == "Soil (non-saline)").to_numpy(),
        "labelled": m.empo_1.isin(["Free-living", "Host-associated"]).to_numpy(),
        "diet": np.where((m.empo_3 == "Animal distal gut") & (m.host_class == "c__Mammalia"),
                         m.host_common_name.map({h: diet_of(h) for h in m.host_common_name.unique()}), np.nan),
    })
    return Dataset("EMP", 5000, meta, indptr, otu, cnt, phylum == "p__Acidobacteria")


def load_ma(depth):
    meta = pd.read_pickle(MAD / "derived" / "meta.pkl").drop_duplicates("map_sid").set_index("map_sid")
    z = np.load(MAD / "derived" / f"rare_{depth}.npz")
    ids = np.load(MAD / "derived" / "sample_ids.npy")[z["sample_row"]]
    otu_ids = np.load(MAD / "derived" / "otu_ids.npy")
    tax = pd.read_csv(MAD / "taxonomic_mapping.tsv", sep="\t", usecols=[1, 2], dtype=str)
    acido_ids = set(tax.loc[tax.Taxonomy.str.contains("Acidobacteri", na=False), "OTU_ID_MAPv3"])
    acido = pd.Series(otu_ids).isin(acido_ids).to_numpy()
    meas = pd.read_csv(MAD / "samples.measurements.tsv", sep="\t", usecols=[0, 1, 2])
    alt = meas[meas.Meta_variable == "altitude_m"].drop_duplicates("SampleID").set_index("SampleID").Value
    keep = pd.Index(meta.index).get_indexer(ids) >= 0  # hash lookup; np.isin on text is O(n*m)
    m = meta.loc[ids[keep]].reset_index()
    rows = np.where(keep)[0]
    sel, indptr = csr_take(z["indptr"], rows)
    otu, cnt = z["otu"][sel], z["count"][sel]
    amp = (m.technology == "AMPLICON").to_numpy()
    e2 = m.env_2.fillna("")
    sal = np.full(len(m), np.nan, dtype=object)
    sal[e2.isin(["aquatic;marine", "aquatic;sea", "aquatic;ocean"]).to_numpy()] = "saline"
    sal[(e2.isin(["aquatic;lake", "aquatic;river", "aquatic;groundwater", "aquatic;reservoir"]) |
         e2.str.startswith("soil")).to_numpy()] = "nonsaline"
    gut = m.body_site.fillna("").str.contains(r"gut|fec|faec|stool|intestin|rum|cecum|caecum|colon", case=False).to_numpy()
    host = e2.str.replace("animal;", "", regex=False)
    host_diet = {h: diet_of(h) for h in host.unique()}  # check each distinct label once
    diet = np.where(gut & (m.env_top == "animal").to_numpy(), host.map(host_diet), np.nan)
    hab_counts = e2.value_counts()
    meta2 = pd.DataFrame({
        "study": m.project.to_numpy(), "country": m.country.to_numpy(), "cell": m.cell.to_numpy(),
        "lat": m.lat.to_numpy(), "lon": m.lon.to_numpy(), "ph": m.ph.to_numpy(), "temp": m.temp.to_numpy(),
        "elev": pd.to_numeric(m["sample"].map(alt), errors="coerce").to_numpy(),
        "habitat": np.where(e2.map(hab_counts).to_numpy() >= 50, e2.to_numpy(), np.nan),
        "free_living": (~m.host_associated & (m.env_top != "unknown")).to_numpy(),
        "salinity": sal, "soil": (m.env_top == "soil").to_numpy(), "diet": diet,
        "labelled": (m.env_top != "unknown").to_numpy(),
    })
    # amplicon only (protocol v0.4)
    a = np.where(amp)[0]
    sel, ip_a = csr_take(indptr, a)
    return Dataset("MicrobeAtlas", depth, meta2.iloc[a], ip_a, otu[sel], cnt[sel], acido)


# ---------------- generic machinery ----------------
def bc_interval(b, est):
    """Bias-corrected percentile interval (v0.6 sensitivity check)."""
    from scipy.stats import norm
    if len(b) < 20 or not np.isfinite(est):
        return [np.nan, np.nan]
    p = np.clip(np.mean(b < est), 1 / len(b), 1 - 1 / len(b))
    z0 = norm.ppf(p)
    lo, hi = norm.cdf(2 * z0 + norm.ppf(0.025)), norm.cdf(2 * z0 + norm.ppf(0.975))
    return [float(np.quantile(b, lo)), float(np.quantile(b, hi))]


def bootstrap(ds, idx, fn, rng):
    """fn(w_full) -> dict; w_full is a weight vector over idx."""
    g = {k: MA.codes(pd.Series(v[idx])) for k, v in ds.g.items()}
    nproj = g["proj"].max() + 1
    out = {}
    for s in SCHEMES:
        est = fn(MA.weights(g, s))
        rows = []
        for _ in range(N_BOOT):
            k = np.bincount(rng.integers(0, nproj, nproj), minlength=nproj).astype(float)
            try:
                rows.append(fn(MA.weights(g, s, k)))
            except Exception:
                pass
        bt = pd.DataFrame(rows)
        out[s] = {"estimate": est, "n_boot_ok": len(bt),
                  "ci95": {c: [float(bt[c].quantile(.025)), float(bt[c].quantile(.975))] for c in bt},
                  "median": {c: float(bt[c].median()) for c in bt},
                  "bc95": {c: bc_interval(bt[c].dropna().to_numpy(), est.get(c, np.nan)) for c in bt}}
    m = ds.meta.iloc[idx]
    info = {"n": int(len(idx)), "studies": int(m.study.nunique()),
            "countries": int(m.country[~m.country.isin(["no_coordinates"])].nunique())}
    return out, info


def wslope(x, y, w):
    X = np.column_stack([np.ones_like(x), x])
    sw = np.sqrt(w)
    return float(np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0][1])


def wcorr(x, y, w):
    w = w / w.sum()
    mx, my = np.sum(w * x), np.sum(w * y)
    c = np.sum(w * (x - mx) * (y - my))
    return float(c / np.sqrt(np.sum(w * (x - mx) ** 2) * np.sum(w * (y - my) ** 2)))


def jaccard(ds, i, j):
    a, b = ds.row(i), ds.row(j)
    inter = np.intersect1d(a, b, assume_unique=True).size
    return 1.0 - inter / (a.size + b.size - inter)


# ---------------- findings ----------------
def two_class(x, lo, hi, a, b):
    out = pd.Series(None, index=x.index, dtype=object)
    out[x < lo] = a
    out[x > hi] = b
    return out


def f4(ds, rng):
    m = ds.meta
    splits = {
        "salinity": m.salinity.where(m.free_living),
        "ph": two_class(m.ph, 6.5, 7.5, "low", "high").where(m.free_living),
        "temp": two_class(m.temp, 15, 25, "cold", "warm").where(m.free_living),
    }
    pairs = {}
    for k, lab in splits.items():
        idx = np.where(lab.notna().to_numpy() & (lab.astype(str) != "nan").to_numpy())[0]
        if len(idx) < 100 or lab.iloc[idx].nunique() < 2:
            continue
        a, b = rng.choice(idx, N_PAIRS), rng.choice(idx, N_PAIRS)
        ok = a != b
        a, b = a[ok], b[ok]
        same = (lab.iloc[a].to_numpy() == lab.iloc[b].to_numpy())
        d = np.array([jaccard(ds, i, j) for i, j in zip(a, b)])
        pairs[k] = (a, b, same, d)
    all_idx = np.unique(np.concatenate([np.concatenate([p[0], p[1]]) for p in pairs.values()]))

    def fn(w):
        res = {}
        for k, (a, b, same, d) in pairs.items():
            pw = w[np.searchsorted(all_idx, a)] * w[np.searchsorted(all_idx, b)]
            res[f"D_{k}"] = float(np.sum(pw[~same] * d[~same]) / np.sum(pw[~same]) - np.sum(pw[same] * d[same]) / np.sum(pw[same]))
        if "D_salinity" in res:
            others = [v for k2, v in res.items() if k2 != "D_salinity"]
            res["salinity_over_max_other"] = res["D_salinity"] / max(others) if others else np.nan
        return res
    return bootstrap(ds, all_idx, fn, rng), "D_salinity", {k: int(len(p[0])) for k, p in pairs.items()}


def slope_finding(ds, mask, xcol, scale, rng, ylog=False):
    idx = np.where(mask & ~np.isnan(ds.meta[xcol].to_numpy()))[0]
    x = ds.meta[xcol].to_numpy()[idx] / scale
    y = ds.rich[idx]
    return bootstrap(ds, idx, lambda w: {"slope": wslope(x, y, w)}, rng)


def f5(ds, rng):
    m = ds.meta
    ds.meta["abslat"] = np.abs(m.lat)
    return slope_finding(ds, m.soil.to_numpy(), "abslat", 10.0, rng), "slope", {}


def f10(ds, rng):
    return slope_finding(ds, ds.meta.soil.to_numpy(), "elev", 1000.0, rng), "slope", {}


def f6(ds, rng):
    m = ds.meta
    ok = m.habitat.notna().to_numpy() & ~np.isnan(m.lat.to_numpy())
    idx = np.where(ok)[0]
    hab = m.habitat.to_numpy()[idx]
    groups = pd.Series(np.arange(len(idx))).groupby(hab).apply(np.array)
    groups = groups[groups.map(len) >= 10]
    sizes = groups.map(len).to_numpy().astype(float)
    p = sizes ** 2 / np.sum(sizes ** 2)
    hsel = rng.choice(len(groups), N_PAIRS, p=p)
    a = np.array([idx[rng.choice(groups.iloc[h])] for h in hsel])
    b = np.array([idx[rng.choice(groups.iloc[h])] for h in hsel])
    ok = a != b
    a, b = a[ok], b[ok]
    sim = np.array([1 - jaccard(ds, i, j) for i, j in zip(a, b)])
    lat1, lon1, lat2, lon2 = (np.radians(v) for v in (m.lat.to_numpy()[a], m.lon.to_numpy()[a], m.lat.to_numpy()[b], m.lon.to_numpy()[b]))
    h = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    km = 6371 * 2 * np.arcsin(np.sqrt(np.clip(h, 0, 1)))
    x = np.log10(km + 1)
    all_idx = np.unique(np.concatenate([a, b]))
    ia, ib = np.searchsorted(all_idx, a), np.searchsorted(all_idx, b)
    return bootstrap(ds, all_idx, lambda w: {"slope": wslope(x, sim, w[ia] * w[ib])}, rng), "slope", {"pairs": int(len(a))}


def f7(ds, rng):
    m = ds.meta
    idx = np.where(m.habitat.notna().to_numpy())[0]
    hab = m.habitat.to_numpy()[idx]
    meanprev = np.empty(len(idx))
    for hname, pos in pd.Series(np.arange(len(idx))).groupby(hab):
        rows = idx[pos.to_numpy()]
        allotu = np.concatenate([ds.row(r) for r in rows])
        prev = np.bincount(allotu, minlength=len(ds.acido)) / len(rows)
        for p_, r in zip(pos.to_numpy(), rows):
            meanprev[p_] = prev[ds.row(r)].mean()
    lr = np.log(ds.rich[idx])
    return bootstrap(ds, idx, lambda w: {"corr": wcorr(lr, meanprev, w)}, rng), "corr", {}


def f8(ds, rng):
    d = ds.meta.diet.astype(str).to_numpy()
    idx = np.where(np.isin(d, ["herbivore", "carnivore"]))[0]
    herb = d[idx] == "herbivore"
    y = ds.rich[idx]
    if herb.sum() < 10 or (~herb).sum() < 10:
        return ({}, {"n": int(len(idx)), "herbivore": int(herb.sum()), "carnivore": int((~herb).sum())}), "log_ratio", {}

    def fn(w):
        return {"log_ratio": float(np.log(np.sum(w[herb] * y[herb]) / np.sum(w[herb]) /
                                          (np.sum(w[~herb] * y[~herb]) / np.sum(w[~herb]))))}
    res = bootstrap(ds, idx, fn, rng)
    return res, "log_ratio", {"herbivore": int(herb.sum()), "carnivore": int((~herb).sum())}


def f9(ds, rng):
    m = ds.meta
    idx = np.where(m.soil.to_numpy() & ~np.isnan(m.ph.to_numpy()) & (m.ph.between(0, 14)).to_numpy())[0]
    share = np.array([ds.count[ds.indptr[i]:ds.indptr[i + 1]][ds.acido[ds.row(i)]].sum() for i in idx]) / ds.depth
    s = np.clip(share, 1e-3, 1 - 1e-3)
    y = np.log(s / (1 - s))
    x = m.ph.to_numpy()[idx]
    return bootstrap(ds, idx, lambda w: {"slope": wslope(x, y, w)}, rng), "slope", {"mean_share": float(share.mean())}


def f123(ds, rng):
    """F1-F3 at this depth with the MicrobeAtlas v0.4 estimators (used for the 5,000-read check)."""
    out = {}
    m = ds.meta
    for f, cfg in MA.COV.items():
        col = cfg["col"]
        idx = np.where(m[col].between(cfg["lo"], cfg["hi"]).to_numpy())[0]
        x, y = m[col].to_numpy()[idx], ds.rich[idx]
        pr = MA.codes(pd.Series(m.study.to_numpy()[idx]))
        out[f] = (bootstrap(ds, idx, lambda w, x=x, y=y, pr=pr, cfg=cfg: MA.cov_est(x, y, pr, w, cfg), rng), "ub_peak", {})
    idx = np.where(m.labelled.to_numpy())[0]
    ha = ~m.free_living.to_numpy()[idx]
    y = ds.rich[idx]
    out["F3"] = (bootstrap(ds, idx, lambda w: MA.f3_est(y, ha, w), rng), "log_ratio_fl_ha", {})
    return out


# ---------------- decisions and table ----------------
RULES = {"F1": ("shift", 0.5), "F2": ("shift", 5.0), "F3": ("ratio", None), "F4": ("f4", None),
         "F5": ("ci_status", None), "F6": ("ratio", None), "F7": ("ratio", None), "F8": ("ratio", None),
         "F9": ("ratio", None), "F10": ("ratio", None)}


def changed(f, s0, s, ci0, ci):
    kind, thr = RULES[f]
    if kind == "shift":
        return abs(s - s0) > thr
    inc0, inc = ci0[0] <= 0 <= ci0[1], ci[0] <= 0 <= ci[1]
    if kind == "ci_status":
        return inc0 != inc or np.sign(s) != np.sign(s0)
    return bool(np.sign(s) != np.sign(s0) or (inc and not inc0) or abs(s - s0) / abs(s0) > 0.5)


def rows_for(f, dataset, depth, res, key, extra):
    (sch, info) = res
    out = []
    if not sch:
        return [{"finding": f, "dataset": dataset, "depth": depth, "scheme": "all", "effect": key,
                 "note": f"not enough samples {info}"}]
    s0, ci0 = sch["S0"]["estimate"][key], sch["S0"]["ci95"][key]
    for s in SCHEMES:
        e, ci = sch[s]["estimate"][key], sch[s]["ci95"][key]
        r = {"finding": f, "dataset": dataset, "depth": depth, "scheme": s, "effect": key,
             "estimate": e, "ci_lo": ci[0], "ci_hi": ci[1], **info,
             "changed_vs_S0": "" if s == "S0" else bool(changed(f, s0, e, ci0, ci)), "extra": json.dumps(extra)}
        if "median" in sch[s]:
            r["boot_median"] = sch[s]["median"].get(key)
            r["bc_lo"], r["bc_hi"] = sch[s]["bc95"].get(key, [np.nan, np.nan])
        if f in ("F1", "F2"):
            est = sch[s]["estimate"]
            thr = RULES[f][1]
            peaks = [est["ub_peak"], est.get("sm99_peak", np.nan), est.get("sm95_peak", np.nan)]
            r["sm99_peak"], r["sm95_peak"] = peaks[1], peaks[2]
            r["peak_stable"] = bool(np.nanmax(peaks) - np.nanmin(peaks) <= thr)
        if f == "F4":
            o, o0 = sch[s]["estimate"].get("salinity_over_max_other"), sch["S0"]["estimate"].get("salinity_over_max_other")
            r["extra"] = json.dumps({**extra, "salinity_over_max_other": o})
            if s != "S0" and o is not None and o0 is not None:
                r["changed_vs_S0"] = bool((o > 1) != (o0 > 1) or abs(e - s0) / abs(s0) > 0.5)
        out.append(r)
    return out


def previous_f123():
    """F1-F3 for EMP (v0.3) and MicrobeAtlas at 1,000 reads (v0.4) from their saved results."""
    rows = []
    emp = json.loads((OUT / "findings_results.json").read_text())
    ma = json.loads((OUT / "ma_findings.json").read_text())
    for src, data, depth in [("EMP", emp, 5000), ("MicrobeAtlas", ma, 1000)]:
        for f, key in [("F1", "ub_peak"), ("F2", "ub_peak"), ("F3", "log_ratio_fl_ha")]:
            sch = data[f]["schemes"]
            info = {k: data[f]["data"].get(k) for k in ("n",)}
            info["studies"] = data[f]["data"].get("studies", data[f]["data"].get("projects"))
            info["countries"] = data[f]["data"].get("countries")
            rows += rows_for(f, src, depth, (sch, info), key, {"source": "earlier run"})
    return rows


FINDERS = {"F4": f4, "F5": f5, "F6": f6, "F7": f7, "F8": f8, "F9": f9, "F10": f10}
DS = None  # set before forking workers; shared copy-on-write


def task(name):
    """One finding on the loaded dataset, in its own process with its own seeded RNG."""
    rng = np.random.default_rng([SEED, DS.depth, sum(map(ord, name))])
    if name == "F123":
        out = []
        for f, (res, key, extra) in f123(DS, rng).items():
            out.append((f, res, key, extra))
        return out
    res, key, extra = FINDERS[name](DS, rng)
    return [(name, res, key, extra)]


def run_dataset(which):
    """which: emp | ma1000 | ma5000. Findings run in parallel (one process each)."""
    import multiprocessing as mp
    global DS
    DS = {"emp": load_emp, "ma1000": lambda: load_ma(1000), "ma5000": lambda: load_ma(5000)}[which]()
    label = "EMP" if which == "emp" else "MicrobeAtlas"
    print(f"== {label} depth {DS.depth}: {len(DS.meta)} samples", flush=True)
    names = list(FINDERS) + ["F123"]  # v0.6: F1-F3 re-run in every dataset
    rows, full = [], {}
    with mp.get_context("fork").Pool(len(names)) as pool:
        for out in pool.imap_unordered(task, names):
            for f, res, key, extra in out:
                rows += rows_for(f, label, DS.depth, res, key, extra)
                full[f"{f}_{label}_{DS.depth}"] = res[0]
                print("  ", f, "done", res[1], flush=True)
    pd.DataFrame(rows).to_csv(OUT / f"final_{which}.csv", index=False)
    (OUT / f"final_{which}.json").write_text(json.dumps(full, indent=1, default=float))


def merge():
    order = {f"F{i}": i for i in range(1, 11)}
    parts = [pd.read_csv(OUT / f"final_{w}.csv") for w in ("emp", "ma1000", "ma5000")]  # v0.6: one estimator everywhere
    df = pd.concat(parts, ignore_index=True)
    df = df.sort_values(["finding", "dataset", "depth", "scheme"], key=lambda c: c.map(order) if c.name == "finding" else c)
    df.to_csv(OUT / "final_findings.csv", index=False)
    print(df[["finding", "dataset", "depth", "scheme", "estimate", "ci_lo", "ci_hi", "n", "changed_vs_S0"]].to_string())


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "merge"
    if which == "merge":
        merge()
    else:
        run_dataset(which)


if __name__ == "__main__":
    main()
