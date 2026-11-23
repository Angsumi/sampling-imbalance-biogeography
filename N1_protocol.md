# N1 protocol: Does geographic sampling imbalance change what we conclude about microbial ecology?

*Version 0.6, 22 Nov 2026. Amendment 0.6 (final experiments) made after the v0.5 results and before running anything below.*

**Amendment 0.6 (22 Nov 2026): final experiments.** Motivated by the v0.5 caveats (unstable peak fits, odd bootstrap intervals, F9 appearing only after rebalancing). All of these are checks on the v0.5 results; the v0.5 decision rules stay the primary result.

1. **One estimator for all datasets.** F1–F10 are re-run in EMP, MicrobeAtlas 1,000 and MicrobeAtlas 5,000 with the same code. F1/F2 use the v0.4 estimator (bins ≥ 20 samples from ≥ 3 studies) everywhere. **1,000** study/project bootstraps (compute now allows it). The earlier EMP v0.3 and MicrobeAtlas v0.4 results are kept as the record.
2. **Steadier peak measure (F1, F2).** Besides the Laplace upper-bound peak, report the argmax of the 3-bin moving average of the weighted 99th and 95th percentiles. A peak counts as **stable** if all three measures fall within 0.5 pH / 5 °C of each other under S0.
3. **Better uncertainty ranges.** Report the bootstrap median and **bias-corrected (BC) percentile intervals** alongside the plain percentile intervals. Decisions use the plain intervals (as pre-set); BC intervals are a sensitivity check.
4. **F9 deep-dive** (soil Acidobacteria vs pH):
   - **Influence:** drop one country at a time (countries with ≥ 30 samples) and refit the S0 slope; report which removals change its sign or significance.
   - **Finer groups:** repeat the slope for each Acidobacteria class (EMP `c__` level; MicrobeAtlas third taxonomy level) with ≥ 1% mean relative abundance.
   - **Mixed model:** logit share ~ pH with a random intercept per study (unweighted), as an alternative to reweighting.
5. **Transfer tests for F5, F9, F10** (leave-one-country-out; model = the finding's predictor only; same failure rule as F1–F3).

*Version 0.5 note:*

**Amendment 0.5 (21 Nov 2026): seven more findings, and a 5,000-read check of all ten.**

**A. New findings.** Each is tested in EMP (90 bp, rarefied to 5,000) and MicrobeAtlas (97% OTUs, rarefied to 1,000), under S0–S3, with 200 study/project bootstraps. "Changed" = the sign flips, or the 95% CI includes zero (when S0's does not), or the effect changes by more than 50%, unless a different rule is given below.

| ID | Published rule | Samples | Effect measured |
|---|---|---|---|
| **F4** | Salinity is the strongest divider of free-living communities, stronger than pH or temperature (Lozupone & Knight 2007, *PNAS*) | Free-living. EMP: EMPO-2 Saline vs Non-saline. MicrobeAtlas: saline = `aquatic;marine/sea/ocean`; non-saline = `aquatic;lake/river/groundwater/reservoir` and `soil*` | Separation D = mean between-class Jaccard distance − mean within-class distance, from 100,000 random pairs. Same D for pH (< 6.5 vs > 7.5) and temperature (< 15 vs > 25 °C). Claim holds if D_salinity > both; **changed if the ordering flips** or D_salinity changes > 50% |
| **F5** | Soil bacterial richness shows no latitude gradient (Fierer & Jackson 2006, *PNAS*) | Soil samples with coordinates | Weighted slope of richness on absolute latitude (per 10°). **Changed if** whether the CI includes zero differs from S0, or the sign flips |
| **F6** | Community similarity decays with geographic distance | Pairs within the same habitat (EMPO-3 / second-level label), both with coordinates; 100,000 pairs | Slope of Jaccard similarity on log10 distance (km) |
| **F7** | Communities are nested: poor communities hold mainly widespread taxa (Thompson et al. 2017) | All samples with a habitat label | Weighted correlation between log richness and the mean within-habitat prevalence of a sample's taxa (expected negative) |
| **F8** | Herbivore guts are richer than carnivore guts (Ley et al. 2008, *Science*) | Mammal gut samples. Herbivores: cattle, sheep, goat, horse, rabbit, deer, buffalo, yak, camel, elephant, kangaroo, wombat, Sambar, colobine monkeys (langurs, snub-nosed, proboscis, guereza), sloth. Carnivores: cat, dog, fox, wolf, seal, lion, tiger, ferret, mink, anteater, aardvark, armadillo. Omnivores (human, pig, rodents, other primates) excluded. MicrobeAtlas gut = body-site field mentions gut/faeces/stool/intestine | log(herbivore ÷ carnivore mean richness) |
| **F9** | Acidobacteria become more common as soil pH falls | Soil samples with pH | Weighted slope of logit(Acidobacteria relative abundance) on pH (expected negative) |
| **F10** | Soil bacterial richness falls with elevation (Bryant et al. 2008, *PNAS*) | Soil samples with elevation/altitude | Weighted slope of richness on elevation (per 1,000 m) |

Transfer tests are not run for F4–F10 (compute limit). Pair-based findings (F4, F6) are reweighted by the product of the two samples' weights.

**B. 5,000-read check.** All of F1–F10 are re-run in MicrobeAtlas using only samples with ≥ 5,000 reads, rarefied to 5,000, so depth matches EMP. EMP results are already at 5,000 reads. A finding is **depth-sensitive** if its S0 conclusion differs between 1,000 and 5,000 reads.

**C. Output:** one table of F1–F10 × dataset × scheme in `results/all_findings.csv`, plus a summary in `results/all_findings_report.md`.

*Version 0.4 note:*

**Amendment 0.4 (21 Nov 2026): extension to MicrobeAtlas (section 6, step 3).** EMP results (`results/findings_report.md`) showed the headline claims are largely robust. As planned, F1–F3 are now re-tested in MicrobeAtlas, which is about 14 times larger and more unevenly sampled.

1. **Data** (downloaded 21 Nov 2026 from `https://microbeatlas.org/downloads/`, MAPv3 / MAPseq v1.2.6):
   - `samples-otus.97.mapped.metag.minfilter.refilt.biom.gz`: 97% OTUs; non-transcriptomic samples with ≥ 1,000 SSU reads and ≥ 20 OTUs.
   - `samples.measurements.tsv` (pH, temperature), `samples.info.latlon.parsed.tsv` (coordinates), `samples.env.info.tsv` (habitat labels, technology, project), `runs.info.tsv.gz`.
2. **Richness:** observed 97% OTUs after rarefying every sample to **1,000 reads** (fixed seed). **Amplicon samples only** in the main analysis; shotgun (WGS) samples are a sensitivity analysis.
3. **Sample sets:**
   - **F1:** samples with pH 0–14.
   - **F2:** samples with temperature −20 to 120 °C.
   - **F3:** host-associated (top-level environment `animal` or `plant`) vs free-living (`soil`, `aquatic` and the remaining non-host environments).
   - Free-living-only versions of F1 and F2 are sensitivity analyses.
4. **Study = MicrobeAtlas project.** **Country** is assigned from coordinates using Natural Earth admin-0 boundaries (1:110m). Points outside any country are grouped as `ocean`; samples without coordinates are grouped as `no_coordinates`.
5. **Upper-bound bins** need **≥ 20 samples from ≥ 3 projects**. This is set prospectively, based on the EMP F2 lesson.
6. **Bootstrap:** 200 project-level resamples (compute limit).
7. **Transfer:** leave-one-country-out for countries with ≥ 200 samples; model = second-level habitat label + covariate + covariate². Same failure rule.
8. **Decision rules:** unchanged from 0.3. In addition, the EMP and MicrobeAtlas S0 estimates are compared, as an independent replication of each claim.
9. **Imbalance is reported first:** share of samples from the top country and top 3 countries, and samples per million people by country.

*Version 0.3 note:*

**Amendment 0.3 (21 Nov 2026, after the F1 pilot, before F2 and F3):**

1. **Sample sets match the paper's counts exactly.** All EMP samples with pH (n = 3,986) and all with temperature (n = 6,976) are the main sets, because these equal the paper's stated n. The free-living-only subsets (pH 3,974) become sensitivity analyses.
2. **Two levels, judged separately.**
   - **Published claim:** F1/F2 upper-bound peak (99th-percentile Laplace fit); F3 log ratio of mean richness. The claim has **changed** if, under S1, S2 or S3: the upper-bound peak moves by more than 0.5 pH or 5 °C; or the F3 ratio flips sign, its CI includes zero, or it changes by more than 50%.
   - **Average pattern (F1/F2 only):** quadratic mean model. It has **changed** if the mean peak moves by more than 0.5 pH or 5 °C, or the curvature's CI includes zero, or the curvature changes by more than 50%.
3. **Fairer transfer test.** The leave-one-country-out model is richness ~ EMPO-3 habitat + covariate + covariate² (F3: habitat only), fitted on the other countries. Habitats unseen in training are predicted from the covariate terms and the training mean. Failure rule unchanged (held-out error > 1.5× CV error in at least one-third of countries with ≥ 30 samples).
4. **F2 bins:** 2 °C bins with at least 10 samples.
5. **S2 caveat:** S2 is kept unchanged but reported as the noisiest scheme; S3 is the main rebalancing for interpretation. This was decided after the pilot and is reported as such.
6. The F1 pilot (v0.2, commit 875c549) stays in `results/` as the record of what was run first.

*Version 0.2 note:*

**Amendment 0.2 (21 Nov 2026, before any analysis):** the source paper's F1 used all non-host-associated samples with pH (n = 3,986), not soil only; the 90 bp Deblur table; and a 99th-percentile upper-bound fit. F1–F3 below now match the paper (see `findings.md`). Soil only becomes a sensitivity analysis. The 90 bp rarefied table was added to the data.

**Rule:** secondary data only. All data are public (Earth Microbiome Project release 1; MicrobeAtlas for extension). No lab work.

---

## 1. Question

Microbiome data come mostly from a few rich countries. When every study, country or region is given equal weight, or when a model is tested on a region it never saw, do well-known findings in microbial ecology still hold?

**Pre-specified hypotheses**

- **H1:** at least one published EMP finding changes materially (see the decision rules in section 5) under study-balanced or geography-balanced weighting.
- **H2:** models of richness fitted on over-sampled regions predict held-out regions worse than within-region cross-validation suggests.
- **H0 (reported either way):** if all findings are stable under every scheme, the paper reports this as evidence that these findings are robust to sampling imbalance.

---

## 2. Novelty check (stage 2)

| Closest work | What it shows | Gap N1 fills |
|---|---|---|
| Abdill et al. 2022, *PLOS Biology* | >71% of human microbiome samples with known origin come from Europe, USA and Canada; Central and South Asia hold 26% of people but 2% of samples | Documents the skew; does not test whether it changes conclusions |
| *Missing microbiomes* (2024, PMC11245146); *Microbiome* commentary (2026) | Under-representation limits who benefits; structural causes | Argument, not a re-analysis |
| Thompson et al. 2017, *Nature* (EMP) | Global patterns, e.g. soil richness vs pH | Findings estimated with sample-level weighting |
| Regional soil studies (e.g. Europe 715 sites, 2024; Australia 2025) | pH and climate effects vary by region | Single regions; no cross-region transfer test |

**Verdict:** passes provisionally. Before freezing, do one more targeted search for re-weighting or leave-region-out re-analyses of EMP or MicrobeAtlas.

---

## 3. Data

### 3.1 Earth Microbiome Project release 1 (downloaded 21 Nov 2026)

Source: `ftp://ftp.microbio.me/emp/release1/` (kept locally in `data/emp_release1/`, not in Git).

| File | Use |
|---|---|
| `emp_qiime_mapping_qc_filtered.tsv` | Sample metadata: study, country, latitude/longitude, EMPO habitat, pH, temperature |
| `emp_deblur_90bp.qc_filtered.rare_5000.biom` | 16S V4 Deblur 90 bp tag sequences, rarefied to 5,000 reads (the table the paper used); richness recomputed from this |
| `emp_deblur_150bp.qc_filtered.rare_5000.biom` | 150 bp version, for sensitivity checks |
| `emp_deblur_150bp.qc_filtered.biom` | Unrarefied table for sensitivity checks |

**What is in it (QC-filtered set):** 23,828 samples, 96 studies, 43 countries. **50% of samples are from the USA**; Australia 16%, Germany 5%.

### 3.2 Data for each finding

| Finding | Samples | Studies | Countries | US share |
|---|---:|---:|---:|---:|
| Soil samples with pH | 935 | 22 | 12 | 46% |
| Any sample with pH | 3,986 | 34 | 18 | 53% |
| Any sample with temperature | 6,976 | 22 | 11 | 45% |
| Animal distal gut | 4,158 | 18 | 17 | 25% (Australia 46%) |

### 3.3 Extension set (only if the pilot passes)

MicrobeAtlas processed tables (2M+ samples) to repeat the main analyses at larger scale.

---

## 4. Findings to re-test

Each is taken as stated in Thompson et al. 2017. The exact published wording and estimate are copied into `findings.md` before any analysis.

| ID | Finding | Effect measured |
|---|---|---|
| **F1** | Richness of non-host-associated communities peaks near neutral pH (n = 3,986) | pH at peak of (a) the 99th-percentile upper bound, as published, and (b) the mean (quadratic fit); curvature |
| **F2** | Richness peaks at cool temperatures, about 10 °C (n = 6,976) | Same two measures for temperature |
| **F3** | Free-living communities are richer than host-associated ones | Log ratio of mean richness (EMPO level 1) |

Richness = observed tag sequences at 5,000 reads in the 90 bp table, recomputed from the BIOM table. Upper bound: weighted 99th percentile in 0.5-pH-unit bins (bins with at least 10 samples), Laplace fit. Soil-only and 150 bp versions are sensitivity analyses. Shannon and Faith's PD are secondary metrics.

---

## 5. Analysis

### 5.1 Four ways to estimate each finding

| Scheme | Weight per sample |
|---|---|
| **S0 Standard** | 1 (as published) |
| **S1 Study-balanced** | 1 / samples in its study |
| **S2 Country-balanced** | 1 / samples in its country |
| **S3 Grid-balanced** | 1 / samples in its 5° × 5° grid cell |

Mixed models with study as a random effect are run alongside as a check.

### 5.2 Transfer test

- **Leave-one-country-out** and **leave-one-continent-out**: fit on the rest, predict richness in the held-out region.
- Compare the held-out error with ordinary 10-fold cross-validation error.

### 5.3 Decision rules (fixed now)

A finding has **changed materially** if, under S1–S3 compared with S0, any of these holds:

1. **F1 and F2:** the peak moves by more than 0.5 pH units or 5 °C, or the curvature's 95% CI includes zero.
2. **F3:** the sign of the effect flips, or its 95% CI includes zero.
3. **Any finding:** the effect size changes by more than 50%.

**Transfer fails (H2)** if the held-out error is more than 1.5 times the cross-validation error in at least one-third of held-out regions.

### 5.4 Robustness

- Bootstrap by study (1,000 resamples) for all CIs.
- 100 bp and 90 bp Deblur tables.
- Unrarefied table with coverage-based richness.
- Drop the largest study per finding.

### 5.5 Main figure (fixed now)

- **Panel A:** world map of EMP sample density.
- **Panel B:** forest plot of each finding's effect under S0–S3.
- **Panel C:** held-out vs cross-validation prediction error by region.

---

## 6. Pilot and go/no-go

1. Run F1 (pH, all non-host-associated samples) only, under S0–S3 and leave-one-country-out.
2. **Go** if F1 changes materially or transfer fails, then add F2 and F3.
3. **If F1 is stable:** still run F2 and F3. If all are stable, either write the robustness result or extend to MicrobeAtlas, where imbalance is larger.

---

## 7. Outputs and resources

- Paper: *Nature Ecology & Evolution* or *ISME J*; robustness-only result → *Global Ecology & Biogeography* or *mSystems*.
- Code and derived tables on Zenodo; raw EMP data cited, not redistributed.
- **Compute:** local only (tables are ~300 MB). Cloud ₹0 unless MicrobeAtlas needs it.

---

## 8. Risks

| Risk | Plan |
|---|---|
| Findings are stable (weaker paper) | Reported honestly; extend to MicrobeAtlas; decided at the pilot |
| Few countries per finding (12 for soil pH) | Grid-balanced and study-balanced schemes; state the limit |
| Study and country are confounded | Study random effects; study-level bootstrap |
| Metadata errors (pH, coordinates) | Range checks; drop impossible values; report how many |
