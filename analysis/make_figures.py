"""N1 manuscript figures and the sample-imbalance table (reporting only; no new estimates).

Reads the frozen v0.6 outputs (results/final_findings.csv, results/final_extra_*.json) and, for the
descriptive imbalance table only, the sample metadata. Writes to results/figures/.

Run from sampling-imbalance-biogeography-rules:  .venv/bin/python analysis/make_figures.py

Figure 1  scorecard: rules x (dataset x weighting), cell = agreement with the published direction
Figure 2  F9: phylum-level Acidobacteria-pH slopes (all datasets, S0-S3 and mixed model) and EMP class split
Figure 3  transfer tests: held-out / cross-validation error ratios (F1-F3 all countries from the v0.3/v0.4 runs;\n          F5, F9, F10 five worst countries from the v0.6 runs)
Figure S1 forest plots of every estimate with 95% study/project-bootstrap intervals
Table 1   sample imbalance by country (results/figures/table1_imbalance.csv)
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = RES / "figures"
FIG.mkdir(exist_ok=True)

# palette (validated: dataviz validate_palette.js, light mode, all checks pass)
AGREE, DISAGREE, UNRES = "#2a78d6", "#eb6834", "#d6d5d0"
INK, INK2 = "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.edgecolor": INK2,
                     "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "savefig.dpi": 300, "savefig.bbox": "tight"})

RULES = {
    "F1": "F1 Richness peaks near pH 7",
    "F2": "F2 Richness peaks near 10 °C",
    "F3": "F3 Free-living > host-associated",
    "F4": "F4 Salinity is the strongest divider",
    "F5": "F5 No latitude gradient (soil)",
    "F6": "F6 Similarity decays with distance",
    "F7": "F7 Nestedness signature",
    "F8": "F8 Herbivore > carnivore gut",
    "F9": "F9 Acidobacteria rise as pH falls",
    "F10": "F10 Richness falls with elevation",
}
EXPECT = {"F3": 1, "F5": 0, "F6": -1, "F7": -1, "F8": 1, "F9": -1, "F10": -1}
DSETS = [("EMP", 5000, "EMP"), ("MicrobeAtlas", 1000, "MA 1k"), ("MicrobeAtlas", 5000, "MA 5k")]
SCHEMES = ["S0", "S1", "S2", "S3"]
SCHEME_LABEL = {"S0": "Sample", "S1": "Study", "S2": "Country", "S3": "Grid"}


def load():
    return pd.read_csv(RES / "final_findings.csv")


def status(r):
    """+1 agrees with the published rule, -1 contradicts it, 0 unresolved (CI includes zero).
    F1/F2: both smoothed peak measures within 1 pH unit of 7 / 5 °C of 10 -> +1, else -1.
    F4: salinity separation larger than both pH and temperature (ratio > 1) -> +1, else -1.
    F5 (a null claim): CI includes zero -> +1; CI excludes zero -> -1."""
    f = r.finding
    if f in ("F1", "F2"):
        c, tol = (7.0, 1.0) if f == "F1" else (10.0, 5.0)
        ok = abs(r.sm99_peak - c) <= tol and abs(r.sm95_peak - c) <= tol
        return 1 if ok else -1
    if f == "F4":
        return 1 if json.loads(r.extra)["salinity_over_max_other"] > 1 else -1
    inc = r.ci_lo <= 0 <= r.ci_hi
    if f == "F5":
        return 1 if inc else -1
    if inc:
        return 0
    return 1 if np.sign(r.estimate) == EXPECT[f] else -1


def cell_text(r):
    f = r.finding
    if f in ("F1", "F2"):
        lo, hi = sorted([r.sm99_peak, r.sm95_peak])  # range of the two smoothed peak measures
        return f"{lo:.1f}–{hi:.1f}" if f == "F1" else f"{lo:.0f}–{hi:.0f}"
    if f == "F4":
        return f"{json.loads(r.extra)['salinity_over_max_other']:.2f}"
    if f in ("F5", "F10"):
        return f"{r.estimate:.0f}" if abs(r.estimate) >= 10 else f"{r.estimate:.1f}"
    if f == "F6":
        return f"{r.estimate:.3f}"
    return f"{r.estimate:.2f}"


def fig1(df):
    order = list(RULES)
    cols = [(d, dep, lab, s) for d, dep, lab in DSETS for s in SCHEMES]
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    rows = []
    for i, f in enumerate(order):
        for j, (d, dep, lab, s) in enumerate(cols):
            r = df[(df.finding == f) & (df.dataset == d) & (df.depth == dep) & (df.scheme == s)].iloc[0]
            st = status(r)
            rows.append({"finding": f, "dataset": lab, "scheme": s, "status": st, "label": cell_text(r)})
            col = {1: AGREE, -1: DISAGREE, 0: UNRES}[st]
            ax.add_patch(plt.Rectangle((j + 0.04, i + 0.06), 0.92, 0.88, color=col, lw=0))
            sym = {1: "✓", -1: "✗", 0: "○"}[st]
            tc = "white" if st != 0 else INK
            ax.text(j + 0.5, i + 0.36, sym, ha="center", va="center", fontsize=7, color=tc)
            ax.text(j + 0.5, i + 0.72, cell_text(r), ha="center", va="center", fontsize=5.2, color=tc)
    ax.set_xlim(0, len(cols))
    ax.set_ylim(len(order), 0)
    ax.set_yticks(np.arange(len(order)) + 0.5, [RULES[f] for f in order])
    ax.set_xticks(np.arange(len(cols)) + 0.5, [SCHEME_LABEL[c[3]] for c in cols], rotation=90)
    for k, (_, _, lab) in enumerate(DSETS):
        ax.text(4 * k + 2, -0.35, lab, ha="center", va="bottom", fontsize=8, fontweight="bold", color=INK)
        if k:
            ax.axvline(4 * k, color="white", lw=3)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (AGREE, UNRES, DISAGREE)]
    ax.legend(handles, ["✓ agrees with published rule", "○ 95% CI includes zero", "✗ contradicts published rule"],
              loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, frameon=False, fontsize=7)
    fig.savefig(FIG / "fig1_scorecard.png")
    fig.savefig(FIG / "fig1_scorecard.pdf")
    plt.close(fig)
    pd.DataFrame(rows).to_csv(FIG / "fig1_scorecard_cells.csv", index=False)


def fig2(df):
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw={"width_ratios": [1.15, 1]})
    y, labels = 0, []
    for d, dep, lab in DSETS:
        ext = json.loads((RES / f"final_extra_{ {'EMP': 'emp', 'MA 1k': 'ma1000', 'MA 5k': 'ma5000'}[lab]}.json").read_text())
        for s in SCHEMES:
            r = df[(df.finding == "F9") & (df.dataset == d) & (df.depth == dep) & (df.scheme == s)].iloc[0]
            st = status(r)
            c = {1: AGREE, -1: DISAGREE, 0: INK2}[st]
            a.plot([r.ci_lo, r.ci_hi], [y, y], color=c, lw=1.5, solid_capstyle="round")
            a.plot(r.estimate, y, "o", color=c, ms=4.5, mec="white", mew=0.8)
            labels.append(f"{lab} · {SCHEME_LABEL[s]}")
            y += 1
        mm = ext["f9_deep"]["mixed_model"]
        a.plot(mm["ci"], [y, y], color=AGREE, lw=1.5)
        a.plot(mm["slope"], y, "D", color=AGREE, ms=4.5, mec="white", mew=0.8)
        labels.append(f"{lab} · mixed model")
        y += 1.6
    ypos = []
    yy = 0
    for _ in DSETS:
        ypos += [yy + k for k in range(5)]
        yy += 6.6 - 1
    a.set_yticks(ypos[:len(labels)], labels, fontsize=6.5)
    a.invert_yaxis()
    a.axvline(0, color=INK2, lw=0.8)
    a.set_xlabel("Slope of logit(Acidobacteria share) on pH")
    a.set_title("a  Phylum level", loc="left", fontsize=8, fontweight="bold")
    emp = json.loads((RES / "final_extra_emp.json").read_text())["f9_deep"]["classes"]
    names = sorted(emp, key=lambda k: emp[k]["S0"]["slope"])
    for i, k in enumerate(names):
        for off, s, mk in [(-0.15, "S0", "o"), (0.15, "S3", "s")]:
            v = emp[k][s]
            c = INK2 if not v["significant"] else (AGREE if v["slope"] < 0 else DISAGREE)
            b.plot(v["ci"], [i + off] * 2, color=c, lw=1.5)
            b.plot(v["slope"], i + off, mk, color=c, ms=4.5, mec="white", mew=0.8)
    b.set_yticks(range(len(names)), [f"{n.replace('c__', '')}\n({100 * emp[n]['mean_share']:.1f}% of reads)" for n in names], fontsize=6.5)
    b.axvline(0, color=INK2, lw=0.8)
    b.set_xlabel("Slope of logit(class share) on pH")
    b.set_title("b  EMP soil, Acidobacteria classes", loc="left", fontsize=8, fontweight="bold")
    b.plot([], [], "o", color=INK2, label="sample-weighted (S0)")
    b.plot([], [], "s", color=INK2, label="grid-weighted (S3)")
    b.legend(frameon=False, fontsize=6.5, loc="lower right")
    for ax in (a, b):
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "fig2_f9_acidobacteria.png")
    fig.savefig(FIG / "fig2_f9_acidobacteria.pdf")
    plt.close(fig)


def fig3():
    """F1-F3: every tested country (EMP v0.3 run, results/findings_results.json; MicrobeAtlas 1,000-read v0.4 run,
    results/ma_findings.json). F5, F9, F10: v0.6 runs, which stored only the five worst countries per test."""
    fig, ax = plt.subplots(figsize=(5.0, 4.4))
    y, ticks, ypos = 0, [], []
    early = {"EMP": json.loads((RES / "findings_results.json").read_text()),
             "MA 1k": json.loads((RES / "ma_findings.json").read_text())}
    for key, lab in [("emp", "EMP"), ("ma1000", "MA 1k"), ("ma5000", "MA 5k")]:
        tr = json.loads((RES / f"final_extra_{key}.json").read_text())["transfer"]
        tests = []
        if lab in early:
            for f in ("F1", "F2", "F3"):
                t = early[lab][f]["transfer"]
                tests.append((f, [c["ratio"] for c in t["countries"]], t["n_fail"], t["n_tested"], "all"))
        for f in ("F5", "F9", "F10"):
            t = tr[f]
            tests.append((f, [w["ratio"] for w in t["worst"]], t["n_fail"], t["n_tested"], "worst5"))
        for f, r, nf, nt, kind in tests:
            r = np.array(r)
            c = np.where(r > 1.5, DISAGREE, AGREE)
            ax.scatter(np.minimum(r, 2.95), [y] * len(r), s=14 if kind == "all" else 16,
                       c=c, marker="o" if kind == "all" else "D", edgecolors="white", linewidths=0.5, zorder=3)
            ax.text(3.05, y, f"{nf} / {nt} countries > 1.5", va="center", fontsize=6.2, color=INK2)
            ticks.append(f"{lab} · {f}")
            ypos.append(y)
            y += 1
        y += 0.7
    ax.set_yticks(ypos, ticks, fontsize=6.5)
    ax.invert_yaxis()
    ax.axvline(1.0, color=INK2, lw=0.8, ls=":")
    ax.axvline(1.5, color=DISAGREE, lw=1.0)
    ax.set_xlim(0.8, 3.0)
    ax.set_xlabel("Held-out country error ÷ cross-validation error")
    ax.scatter([], [], marker="o", c=INK2, s=14, label="every tested country (F1–F3)")
    ax.scatter([], [], marker="D", c=INK2, s=16, label="five worst countries stored (F5, F9, F10)")
    ax.legend(frameon=False, fontsize=6.2, loc="lower right", bbox_to_anchor=(1.0, -0.22), ncol=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(FIG / "fig3_transfer.png")
    fig.savefig(FIG / "fig3_transfer.pdf")
    plt.close(fig)


def figS1(df):
    fig, axes = plt.subplots(5, 2, figsize=(7.2, 9.0))
    units = {"F1": "Peak pH (Laplace mode)", "F2": "Peak °C (Laplace mode)", "F3": "log(FL ÷ HA richness)",
             "F4": "Salinity separation D", "F5": "Δ richness per 10° latitude", "F6": "Δ Jaccard similarity per log10 km",
             "F7": "Weighted correlation", "F8": "log(herbivore ÷ carnivore)", "F9": "Slope per pH unit",
             "F10": "Δ richness per 1,000 m"}
    for ax, f in zip(axes.flat, RULES):
        y, ticks = 0, []
        for d, dep, lab in DSETS:
            for s in SCHEMES:
                r = df[(df.finding == f) & (df.dataset == d) & (df.depth == dep) & (df.scheme == s)].iloc[0]
                c = {1: AGREE, -1: DISAGREE, 0: INK2}[status(r)]
                lo, hi = r.ci_lo, r.ci_hi
                if f in ("F1", "F2"):  # clip absurd Laplace intervals for display
                    lim = (0, 14) if f == "F1" else (-10, 60)
                    lo, hi = max(lo, lim[0]), min(hi, lim[1])
                ax.plot([lo, hi], [y, y], color=c, lw=1.2)
                ax.plot(r.estimate, y, "o", color=c, ms=3.5, mec="white", mew=0.5)
                ticks.append(f"{lab} {s}")
                y += 1
            y += 0.6
        ax.set_yticks([i + 0.6 * (i // 4) for i in range(12)], ticks, fontsize=5.5)
        ax.invert_yaxis()
        if f not in ("F1", "F2"):
            ax.axvline(0, color=INK2, lw=0.6)
        ax.set_title(RULES[f], loc="left", fontsize=7, fontweight="bold")
        ax.set_xlabel(units[f], fontsize=6.5)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "figS1_forest_all.png")
    fig.savefig(FIG / "figS1_forest_all.pdf")
    plt.close(fig)


def table1():
    """Descriptive country shares of the analysed samples (EMP QC-filtered; MicrobeAtlas amplicon at 1,000 reads)."""
    emp = pd.read_csv(ROOT / "data/emp_release1/emp_qiime_mapping_qc_filtered.tsv", sep="\t", dtype=str,
                      usecols=["#SampleID", "country", "study_id"], encoding_errors="replace")
    emp["country"] = emp.country.str.replace("GAZ:", "", regex=False)
    meta = pd.read_pickle(ROOT / "data/microbeatlas/derived/meta.pkl")
    z = np.load(ROOT / "data/microbeatlas/derived/rare_1000.npz")
    ids = np.load(ROOT / "data/microbeatlas/derived/sample_ids.npy")[z["sample_row"]]
    ma = meta[meta.map_sid.isin(set(ids)) & (meta.technology == "AMPLICON")].drop_duplicates("map_sid")
    out = []
    for name, d, study in [("EMP", emp, "study_id"), ("MicrobeAtlas (1,000 reads, amplicon)", ma, "project")]:
        c = d.country.fillna("no_coordinates")
        located = c[~c.isin(["no_coordinates", "ocean"])]
        vc = located.value_counts()
        out.append({"dataset": name, "samples": len(d), "studies": d[study].nunique(),
                    "samples_with_country": len(located), "countries": len(vc),
                    "top_country": vc.index[0], "top_country_share": round(vc.iloc[0] / len(located), 3),
                    "top3_countries": "; ".join(vc.index[:3]), "top3_share": round(vc.iloc[:3].sum() / len(located), 3),
                    "ocean_samples": int((c == "ocean").sum()), "no_coordinates": int((c == "no_coordinates").sum())})
    pd.DataFrame(out).to_csv(FIG / "table1_imbalance.csv", index=False)
    print(pd.DataFrame(out).T.to_string())


if __name__ == "__main__":
    df = load()
    fig1(df)
    fig2(df)
    fig3()
    figS1(df)
    table1()
    print("written to", FIG)
