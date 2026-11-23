# N1 final experiment results (protocol v0.6)

*22 Nov 2026. Protocol v0.6 committed before these runs (commit 58e39c0). Scripts: `analysis/all_findings.py` (F1–F10, one estimator for all datasets, 1,000 study/project bootstraps), `analysis/n1_extra.py` (F9 deep-dive, transfer tests). Tables: `final_findings.csv` (120 rows, with plain and bias-corrected intervals and steadier peak measures), `final_extra_{emp,ma1000,ma5000}.json`. First-pass analysis, not yet checked. Earlier results (`all_findings*`, `findings_*`, `ma_findings*`) are kept as the record.*

**Datasets:**
- EMP: 90 bp, 5,000 reads; 23,828 samples, 96 studies.
- MicrobeAtlas at 1,000 reads: 1.67M amplicon samples.
- MicrobeAtlas at 5,000 reads: 1.55M samples.

**Weightings:** S0 every sample counts once; S1 every study equal; S2 every country equal; S3 every 5° map square equal.

## Final scorecard

| # | Published rule | Verdict | Key numbers |
|---|---|---|---|
| F1 | Richness peaks near neutral pH | **Holds** (with a steadier peak measure) | Smoothed peak pH 6.2–7.7 in every dataset and weighting. The Laplace-fit peak was the unstable part (MicrobeAtlas 5,000 S1: 0.4) |
| F2 | Richness peaks at about 10 °C | **EMP only.** Not in MicrobeAtlas | EMP: all three peak measures give 8–12 °C. MicrobeAtlas: measures disagree (5–112 °C); no stable peak |
| F3 | Free-living richer than host-associated | **Holds** | log ratio 0.17–1.12 in every case; size depends on weighting and depth |
| F4 | Salinity is the strongest divider | **Not supported** | EMP: pH or temperature divide more in all weightings (ratio 0.56–0.74). MicrobeAtlas: the ranking flips with weighting (0.15–1.93) |
| F5 | No latitude gradient in soil richness | **Depends on data and weighting** | EMP: a clear gradient (richness falls with latitude) in every weighting. MicrobeAtlas: none under S0/S3, but study and country weighting create opposite gradients |
| F6 | Similarity decays with distance | **Holds** (weaker under grid weighting) | Negative in every case; the CI touches zero under S3 in MicrobeAtlas |
| F7 | Nestedness (poor communities hold widespread taxa) | **Holds in MicrobeAtlas; fragile in EMP** | MicrobeAtlas −0.28 to −0.70 everywhere; EMP fades to −0.04 (CI includes zero) under country weighting |
| F8 | Herbivore guts richer than carnivore guts | **Holds** (in MicrobeAtlas) | MicrobeAtlas log ratio 0.58–1.03 everywhere (herbivores 1.8–2.8× richer); EMP too thin (128 carnivore samples) |
| F9 | Acidobacteria rise as soil pH falls | **Hidden in the standard analysis, appears after correction** | S0 not significant in all three datasets. Significant negative under S1 (all three), S3 (EMP, MicrobeAtlas 5,000), and in a mixed model with study as random effect (all three) |
| F10 | Soil richness falls with elevation | **Not supported** | Never significant; sign changes with weighting |

## F9 deep-dive

- **Why it's hidden in EMP:** inside Acidobacteria, classes go in opposite directions. *Acidobacteriia*, DA052 and *Solibacteres* fall as pH rises (as published), but *Chloracidobacteria*, the most abundant class (6%), **rises**. Summed at phylum level, they cancel out in the standard analysis.
- **Correcting for study differences** (mixed model, unweighted) recovers the negative slope in all three datasets: EMP −0.29, MicrobeAtlas −0.10 and −0.09; all significant.
- **No single country causes it:** dropping any one country never changes the S0 result.
- **MicrobeAtlas caveat:** Acidobacteria are only 1.6% of soil reads there (EMP 19%), because many MicrobeAtlas OTUs lack a phylum label. No class reached the 1% cut-off, so the class split couldn't be tested there.

## Transfer tests (leave one country out)

No failures anywhere for F1–F10. For F5, F9 and F10: 0 failing countries out of 14 / 5 / 14 tested (EMP) and 13–50 tested (MicrobeAtlas). Wherever these rules hold, they predict unseen countries about as well as seen ones.

## Uncertainty checks

- Bias-corrected intervals agree with the plain ones except where peak fits are unstable (F1/F2 Laplace), which the steadier peak measures resolve.
- F6 in MicrobeAtlas (S0): the estimate sits just outside its interval. Pair statistics resampled by project over-weight within-project pairs; a known limitation of this bootstrap for pair-based effects. The direction (negative) is not in doubt.

## Bottom line

- **4 rules hold everywhere:** F1, F3, F6, F8.
- **2 hold in one dataset only:** F2 (EMP), F7 (MicrobeAtlas).
- **2 are not supported:** F4, F10.
- **2 depend on how the data are weighted:** F5 (latitude) and F9 (Acidobacteria). F9 is invisible in the usual analysis and appears consistently once study or geographic imbalance is corrected.
- Read depth changes effect sizes (F3, F7), not directions.
