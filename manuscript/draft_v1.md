# Uneven sampling hides and creates rules of microbial biogeography

*Draft v1, 23 Nov 2026. Target: Nature Communications (Article) or ISME Journal (Original Article). First full draft, written by an AI agent (Claude) from the frozen N1 outputs; **provisional until a person has checked every number, quote and reference.***

**Traceability convention.** Every number is followed by its source in the form ⟦file · row/key⟧. `FF` = `results/final_findings.csv`, rows identified as finding/dataset/depth/scheme (EMP rows are at 5,000 reads; "MA1k"/"MA5k" = MicrobeAtlas at 1,000/5,000 reads). `XE`, `XM1`, `XM5` = `results/final_extra_{emp,ma1000,ma5000}.json`. `T1` = `results/figures/table1_imbalance.csv`. `EMPv3` = `results/findings_results.json` (EMP, protocol v0.3), `MAv4` = `results/ma_findings.json` (MicrobeAtlas, v0.4). `SC` = `results/figures/fig1_scorecard_cells.csv`. Source quotes are in `findings.md`. These tags are removed in the submission version.

---

**Authors:** [PLACEHOLDER — author list, affiliations, ORCID, corresponding author]

---

## Abstract

Global microbiome surveys draw most of their samples from a few countries and a few large studies, yet the rules of microbial biogeography derived from them are estimated as if every sample were independent. We re-tested ten widely cited rules — four from the Earth Microbiome Project (EMP) flagship paper — in the EMP (23,828 samples, 96 studies) and in 1.67 million amplicon samples from MicrobeAtlas, weighting every sample, every study, every country or every 5° map cell equally, and predicting held-out countries. Four rules held in every dataset and weighting: richness peaks near neutral pH, free-living communities are richer than host-associated ones, similarity decays with distance and herbivore guts are richer than carnivore guts. Two held in only one dataset, two were not supported, and two depended on weighting. The decline of soil Acidobacteria with rising pH was invisible in the standard analysis of all three data sets but appeared after correcting for study imbalance, because acid-loving and alkaline-tolerant Acidobacteria classes cancel out at phylum level. Where rules held, they transferred to unseen countries. Reporting study-balanced estimates alongside sample-weighted ones should become routine in microbial macroecology.

*(~180 words; Nature Communications limit 150 — to be shortened.)*

---

## Introduction

Large compilations of microbial community data have made it possible to look for general rules of microbial biogeography: how richness varies with pH, temperature, latitude and elevation, how composition separates saline from non-saline habitats, and how similarity decays with distance [1–10]. The Earth Microbiome Project (EMP) catalogue was built to make such cross-study comparisons possible with standardised protocols [1, 11], and its flagship analysis reported that richness peaks near neutral pH and at relatively cool temperatures, that host-associated communities are less rich than free-living ones and that communities are nested [1].

These data sets are not evenly sampled. In the EMP release analysed here, 49.7% of samples come from the United States and 70.4% from three countries ⟦T1 · EMP · top_country_share, top3_share⟧; in the 1.67 million MicrobeAtlas amplicon samples we analysed, 39.3% of samples with a known country come from the United States ⟦T1 · MicrobeAtlas · top_country_share⟧. The same skew has been documented for human microbiome data, in which more than 71% of samples with known origin come from Europe, the United States and Canada [12], and has been argued to limit who benefits from microbiome research [13, 14]. Whether it also changes the ecological conclusions drawn from these data has not, to our knowledge, been tested. The EMP authors themselves cautioned that their pH and temperature patterns "should therefore be interpreted with caution" [1] and called for sampling that spans gradients of geography and chemistry more evenly [1].

When samples are clustered within studies and countries, a sample-weighted estimate is dominated by the largest studies. A rule can then appear general when it reflects a few well-sampled settings, or be hidden when large studies with unusual conditions swamp a consistent pattern. We therefore asked, for ten published rules (Table 2), whether the conclusion changes when every study, every country or every 5° map cell contributes equally, whether the result replicates in an independent, much larger compilation (MicrobeAtlas [15]), and whether models fitted on over-sampled countries predict held-out countries. We pre-specified two hypotheses: that at least one published finding changes materially under study- or geography-balanced weighting (H1), and that models fitted on over-sampled regions predict held-out regions worse than cross-validation suggests (H2).

## Results

### Two data sets with strong geographic imbalance

The EMP quality-filtered release contains 23,828 samples from 96 studies and 43 countries ⟦T1 · EMP · samples, studies, countries⟧; the United States, Australia and Germany supply 70.4% ⟦T1 · EMP · top3_share⟧. After rarefaction to 1,000 reads, MicrobeAtlas provided 1,665,984 amplicon samples from 23,742 projects ⟦T1 · MicrobeAtlas · samples, studies⟧; 956,833 could be placed in one of 166 countries, 154,612 fell in the ocean and 554,539 had no coordinates ⟦T1 · MicrobeAtlas · samples_with_country, countries, ocean_samples, no_coordinates⟧. Among located samples, the United States, China and Canada supply 57.9% ⟦T1 · MicrobeAtlas · top3_share⟧. At 5,000 reads, 1,545,275 MicrobeAtlas samples remained ⟦`results/final_run_ma5000.log` · line 1⟧. Sample sizes per rule range from 935 soil samples with pH in 22 EMP studies ⟦FF · F9/EMP/S0 · n, studies⟧ to 1.67 million MicrobeAtlas samples ⟦FF · F7/MA1k/S0 · n⟧ (Supplementary Table S1).

### A scorecard of ten rules

Across the 90 comparisons of a re-weighted estimate with its sample-weighted counterpart (10 rules × 3 data sets × 3 weightings), 39 met a pre-specified criterion for a material change ⟦FF · changed_vs_S0 = True, counted over S1–S3 rows⟧, so H1 is supported. Fig. 1 summarises each estimate as agreeing with the published rule, contradicting it, or unresolved (95% study-bootstrap interval includes zero). The rules fell into four groups.

**Rules that hold everywhere (F1, F3, F6, F8).** Richness peaked near neutral pH in every data set and weighting: the two smoothed upper-bound peak measures lay between pH 6.2 and 7.7 in all 12 cases ⟦FF · F1/all · sm99_peak, sm95_peak⟧. The Laplace upper-bound fit used in the original paper was unstable — in MicrobeAtlas at 5,000 reads with study weighting its mode fell at pH 0.4 ⟦FF · F1/MA5k/S1 · estimate⟧ — which is why the smoothed measures added in protocol v0.6 were used for interpretation. Free-living communities were richer than host-associated ones in all 12 cases, with log richness ratios from 0.16 to 1.12 ⟦FF · F3/MA1k/S1 · estimate; F3/EMP/S2 · estimate⟧; the size of the difference depended on weighting and read depth but its direction did not. Community similarity decayed with log distance in every case (slopes −0.008 to −0.043 per log10 km) ⟦FF · F6/MA1k/S2; F6/MA5k/S3 · estimate⟧, although under grid weighting in MicrobeAtlas the interval reached zero ⟦FF · F6/MA1k/S3, F6/MA5k/S3 · ci_hi⟧. Herbivore guts were richer than carnivore guts in MicrobeAtlas under every weighting (log ratio 0.58–1.03, i.e. 1.8–2.8-fold) ⟦FF · F8/MA1k/S2; F8/MA5k/S0 · estimate⟧; the EMP has too few carnivore samples (128) to resolve this ⟦FF · F8/EMP/S0 · extra.carnivore⟧, and its interval included zero in three of four weightings ⟦SC · F8/EMP⟧.

**Rules that hold in one data set only (F2, F7).** In the EMP, all three peak measures placed maximal richness at 6–12 °C under every weighting ⟦FF · F2/EMP/all · estimate, sm99_peak, sm95_peak⟧, reproducing the published optimum of about 10 °C. In MicrobeAtlas the measures disagreed widely (5–112 °C) and no stable peak emerged ⟦FF · F2/MA1k, F2/MA5k · sm99_peak, sm95_peak⟧; only one of eight MicrobeAtlas cells met our agreement criterion ⟦SC · F2/MA5k/S3⟧. Part of this difference may reflect the sample set: 24% of MicrobeAtlas amplicon samples with a temperature value are host-associated ⟦`findings.md` · F2 row⟧, whereas the published claim concerns non-host-associated samples. The prevalence signature of nestedness — poorer communities holding mainly widespread taxa — was strong in MicrobeAtlas under all weightings (weighted correlation −0.28 to −0.70) ⟦FF · F7/MA5k/S3; F7/MA1k/S2 · estimate⟧, but in the EMP it weakened from −0.40 to −0.04, with an interval including zero, under country weighting ⟦FF · F7/EMP/S0, F7/EMP/S2⟧.

**Rules not supported (F4, F10).** Salinity did not separate EMP communities more strongly than our pH or temperature contrasts under any weighting: its separation was 0.56–0.74 times the larger of the other two ⟦FF · F4/EMP/all · extra.salinity_over_max_other⟧. In MicrobeAtlas the ranking flipped with weighting, from 0.15 to 1.93 ⟦FF · F4/MA5k/S2, F4/MA5k/S1 · extra.salinity_over_max_other⟧. Because our pH and temperature contrasts were moderate splits rather than the extremes compared by Lozupone and Knight [2], this tests whether salinity is the strongest of three pre-defined binary splits, not their original claim. Soil richness showed no consistent elevation gradient: no estimate in any data set or weighting excluded zero, and the sign changed with weighting ⟦FF · F10/all · ci_lo, ci_hi, estimate⟧, in line with later work that found no elevational gradient in bacterial diversity [6] rather than with the original report [5].

**Rules that depend on weighting (F5, F9).** Fierer and Jackson reported that soil bacterial diversity "was unrelated to … latitude" [3]. In the EMP, soil richness fell with absolute latitude under every weighting (−68 tag sequences per 10°, 95% CI −133 to −11, sample-weighted) ⟦FF · F5/EMP/S0⟧, consistent with the weak latitudinal trend reported in the EMP paper itself [1]. In MicrobeAtlas there was no gradient under sample or grid weighting, but study weighting and country weighting produced significant gradients of opposite sign at 5,000 reads (−2.7 and +7.9 OTUs per 10°) ⟦FF · F5/MA5k/S1, F5/MA5k/S2⟧. The Acidobacteria–pH relationship is examined next.

### A soil rule hidden by study imbalance

The relative abundance of Acidobacteria is a textbook correlate of soil pH (R = −0.80 in a continental survey [7]). In the standard, sample-weighted analysis this slope was not significant in any of the three data sets; in the EMP it was even slightly positive (0.12, 95% CI −0.32 to 0.26) ⟦FF · F9/EMP/S0; F9/MA1k/S0; F9/MA5k/S0⟧ (Fig. 2a). With study weighting it became significantly negative in all three (EMP −0.23, MicrobeAtlas −0.16 and −0.18) ⟦FF · F9/EMP/S1, F9/MA1k/S1, F9/MA5k/S1⟧, and with grid weighting in the EMP and in MicrobeAtlas at 5,000 reads ⟦FF · F9/EMP/S3, F9/MA5k/S3⟧. A mixed model with a random intercept per study, without any weighting, recovered a negative slope in all three data sets (EMP −0.29, 95% CI −0.38 to −0.21; MicrobeAtlas −0.10 and −0.09) ⟦XE, XM1, XM5 · f9_deep.mixed_model⟧. No single country drove the sample-weighted null: dropping any one country never made the slope significant ⟦XE, XM1, XM5 · f9_deep.drop_one_country_S0⟧ (dropping Antarctica in the EMP turned the point estimate slightly negative, −0.006, but not significant ⟦XE · drop_one_country_S0 · Antarctica⟧).

The cancellation is visible below phylum level (Fig. 2b). In EMP soils, *Acidobacteriia*, DA052 and *Solibacteres* declined as pH rose under both sample and grid weighting, whereas *Chloracidobacteria* — the most abundant class, 6.2% of reads — increased (slope 0.90, 95% CI 0.60 to 1.07, sample-weighted) ⟦XE · f9_deep.classes⟧. Summed at phylum level, these opposing responses cancel when the sample set is dominated by a few large studies, and the expected negative relationship appears only once study differences are accounted for. This class-level test could not be repeated in MicrobeAtlas, where Acidobacteria made up only 1.6–1.7% of soil reads (EMP 19.2%) ⟦FF · F9/MA1k, F9/MA5k, F9/EMP · extra.mean_share⟧ because many MicrobeAtlas OTUs lack a phylum label, and no class reached the 1% abundance threshold ⟦XM1, XM5 · f9_deep.classes = {}⟧.

### Rules that hold transfer to unseen countries

Leave-one-country-out prediction never met the pre-specified failure rule (held-out error more than 1.5 times the cross-validation error in at least one-third of countries) (Fig. 3), so H2 is not supported. For F5, F9 and F10, no tested country exceeded 1.5 in any data set (0 of 5–50 countries per test) ⟦XE, XM1, XM5 · transfer.*.n_fail, n_tested⟧. For F1–F3, isolated countries did: in the EMP temperature model, Australia (1,559 samples) was predicted 2.8 times worse than cross-validation suggested ⟦EMPv3 · F2.transfer.countries · Australia⟧, and Spain (pH) and Antarctica (temperature) exceeded the threshold in MicrobeAtlas ⟦MAv4 · F1, F2 .transfer.countries; `results/ma_findings_report.md` lines 29, 44⟧. These are single-country failures, not the systematic failure H2 anticipated.

## Discussion

Sampling imbalance did not overturn most of what the EMP taught about microbial biogeography. The neutral-pH richness optimum, the richness advantage of free-living over host-associated communities, the decay of similarity with distance and the richness of herbivore guts survive study, country and grid weighting in two independent compilations that differ in sequencing processing, taxonomic resolution and size by two orders of magnitude. That is useful: it identifies which rules can be taught as general.

The weighting-dependent rules carry the more important message. The Acidobacteria–pH relationship, among the most reproducible single-taxon relationships in soil microbiology [7, 8], was undetectable in the standard analysis of three data sets and recovered once each study counted equally or study was modelled as a random effect. The mechanism we found — opposing class-level responses cancelling at phylum level, with large studies tipping the balance — is a form of aggregation bias that standard global analyses would not flag. The latitude result shows the converse risk: in MicrobeAtlas, study and country weighting produced significant gradients of opposite sign from data in which the sample-weighted analysis found none. A rule whose sign depends on how studies are weighted should not be reported as a global pattern in either direction.

Three rules did not replicate as stated. For temperature, the EMP's cool optimum was stable across weightings within the EMP but did not reappear in MicrobeAtlas; whether this reflects host-associated samples in the MicrobeAtlas temperature set, different habitat composition, or a pattern specific to the EMP's studies remains to be tested. Salinity was not the strongest of our three contrasts, but our contrasts differ from the original comparison of salinity with extremes of pH and temperature [2]. The elevation rule had little support in either data set, consistent with the mixed record in the literature [5, 6].

Our analysis builds on balancing steps the EMP authors already used: their nestedness test used a subset "with even representation across environments and studies" and their latitude analysis included study as a random effect [1]. We extend that logic to all rules, to a second data set and to transfer tests. MicrobeAtlas [15] makes this kind of replication possible at a scale the EMP could not reach, but it inherits a stronger geographic skew and coarser metadata.

**Practical recommendation.** Global microbiome analyses should report, alongside the sample-weighted estimate, a study-balanced estimate (or a model with study as a random effect) and a geography-balanced estimate, and should flag rules whose sign or significance depends on the choice.

### Limitations

1. **First-pass analysis.** All analyses were run once by an AI agent under a pre-registered, amended protocol and have not been independently re-run or code-reviewed by a person.
2. **MicrobeAtlas phylum-label gap.** Acidobacteria account for 1.6–1.7% of MicrobeAtlas soil reads versus 19% in the EMP ⟦FF · F9 · extra.mean_share⟧, because many OTUs lack a phylum assignment. F9 in MicrobeAtlas therefore rests on the labelled fraction, and the class-level explanation could be tested only in the EMP.
3. **Bootstrap limitation for pair statistics.** For pair-based effects (F4, F6) we resampled projects, which over-weights within-project pairs; in MicrobeAtlas the sample-weighted F6 estimate (−0.0104 at 1,000 reads) lies just outside its own interval (−0.0148 to −0.0111) ⟦FF · F6/MA1k/S0⟧. Directions are not in doubt, but intervals for F4 and F6 should be read as approximate.
4. **Proxy tests.** F7 uses a prevalence–richness correlation at tag-sequence/OTU level, not the NODF statistic at phylum to order level used in the original paper [1]; F4 uses moderate rather than extreme pH and temperature contrasts; F10's source [5] may concern Acidobacteria rather than all bacteria (to be confirmed) (`findings.md`).
5. **Sample sets.** The MicrobeAtlas F1/F2 sets include host-associated samples, whereas the published claim concerns non-host-associated samples; the free-living-only sensitivity analysis planned in protocol v0.4 is not in the final tables.
6. **Transfer tests are weak for single-predictor models.** For F5, F9 and F10 the transfer model contains only the rule's predictor; when the predictor explains little variance, held-out and cross-validated errors are both dominated by residual variance and will be similar. F1–F3 transfer used earlier protocol versions (v0.3/v0.4) with habitat terms. Only the five worst countries were stored for the v0.6 tests.
7. **Countries from coordinates.** MicrobeAtlas countries were assigned from coordinates at 1:110m resolution, with ocean samples pooled as one group and 33% of samples lacking coordinates ⟦T1⟧; country counts in Supplementary Table S1 include the ocean group.
8. **Rarefaction and taxonomic units differ** between data sets (Deblur 90-bp tag sequences vs 97% OTUs), so effect sizes are comparable in direction, not magnitude.

## Methods

### Protocol

The question, hypotheses, weighting schemes, decision rules and main figure were fixed in a version-controlled protocol (`N1_protocol.md`, original version) before the first analysis. Five dated amendments (v0.2–v0.6), each committed before the analyses it governs, matched the sample sets to the published counts (v0.2–v0.3), added MicrobeAtlas (v0.4), added F4–F10 and a read-depth check (v0.5), and added steadier peak measures, bias-corrected intervals, the F9 deep-dive and transfer tests for F5, F9 and F10 (v0.6). The v0.6 results (`results/final_*`) are the primary results; earlier runs are kept in `results/` as the record. Secondary data only: no samples were collected.

### Data

**EMP.** Release 1 (downloaded 21 Nov 2026 from `ftp://ftp.microbio.me/emp/release1/`): QC-filtered mapping file and the Deblur [16] 90-bp table rarefied to 5,000 reads per sample (23,828 samples), the table used for the published diversity analyses [1]. Richness is the number of observed tag sequences.

**MicrobeAtlas.** MAPv3 downloads (21 Nov 2026, `https://microbeatlas.org/downloads/`; MAPseq [17] v1.2.6): 97% OTU table of non-transcriptomic samples with ≥ 1,000 SSU reads and ≥ 20 OTUs, sample measurements (pH, temperature, altitude), coordinates, environment labels, technology and project [15]. Every sample was rarefied to 1,000 reads, and samples with ≥ 5,000 reads also to 5,000 reads, with a fixed seed (`analysis/ma_rarefy.py`). The main analysis used amplicon samples only. "Study" is the MicrobeAtlas project. Countries were assigned from coordinates using Natural Earth admin-0 boundaries at 1:110m (`analysis/ma_prepare.py`); points outside any country were grouped as "ocean".

### Rules and effect measures

Each rule's source wording and location are given in `findings.md` and Table 2. F1, F2: upper bound on richness against pH (bins of 0.5 pH units) or temperature (2 °C), using the weighted 99th percentile in bins with ≥ 20 samples from ≥ 3 studies; peak = mode of a Laplace fit, plus the argmax of the three-bin moving average of the weighted 99th and 95th percentiles (v0.6). F3: log ratio of mean richness, free-living to host-associated (EMPO level 1; MicrobeAtlas top-level environment). F4: separation D = mean between-class minus mean within-class Jaccard distance from 100,000 random pairs of free-living samples, for saline vs non-saline, pH < 6.5 vs > 7.5 and < 15 vs > 25 °C. F5, F10: weighted slope of soil richness on absolute latitude (per 10°) or elevation (per 1,000 m). F6: weighted slope of Jaccard similarity on log10(distance + 1 km) for 100,000 within-habitat pairs. F7: weighted correlation between log richness and the mean within-habitat prevalence of a sample's taxa. F8: log ratio of mean richness, herbivore to carnivore mammal gut (omnivores excluded; species lists in `analysis/all_findings.py`). F9: weighted slope of logit(Acidobacteria relative abundance) on soil pH.

### Weighting schemes

S0, every sample weight 1 (as published); S1, 1/(samples in its study); S2, 1/(samples in its country); S3, 1/(samples in its 5° × 5° cell). Pair-based effects used the product of the two samples' weights.

### Uncertainty and decision rules

95% intervals are percentiles of 1,000 study (project) bootstrap resamples, with bias-corrected percentile intervals [18] as a sensitivity check (`analysis/all_findings.py`). A re-weighted estimate counted as materially changed if it flipped sign, if its interval included zero when the S0 interval did not, or if it changed by more than 50% (F1/F2: peak shift > 0.5 pH or 5 °C; F4: the salinity ranking flipped; F5: whether the interval includes zero changed). For Fig. 1, each estimate was additionally classified as agreeing with, contradicting or not resolving the published rule (rules in the docstring of `analysis/make_figures.py:status`).

### F9 deep-dive

Drop-one-country refits of the S0 slope (countries with ≥ 30 samples), class-level slopes for Acidobacteria classes with ≥ 1% mean relative abundance, and a linear mixed model of logit share on pH with a random intercept per study, fitted with statsmodels [19] (`analysis/n1_extra.py`; 200 study bootstraps for the deep-dive intervals).

### Transfer tests

Leave-one-country-out prediction compared with 10-fold cross-validation [20], for countries with ≥ 30 (EMP) or ≥ 200 (MicrobeAtlas) samples. Transfer fails if held-out error exceeds 1.5 times cross-validation error in at least one-third of countries. F1–F3 (v0.3/v0.4) used richness ~ habitat + covariate + covariate²; F5, F9 and F10 (v0.6) used the rule's predictor only.

### Software

Python 3 with NumPy, pandas, SciPy, statsmodels [19], h5py, shapely and matplotlib; exact versions in `requirements.txt`. Every result can be regenerated with the commands in the docstrings of `analysis/ma_prepare.py`, `ma_rarefy.py`, `ma_findings.py`, `all_findings.py`, `n1_extra.py` and `make_figures.py`.

## Data availability

EMP release 1 and MicrobeAtlas MAPv3 are public (URLs above) and are cited, not redistributed. Derived per-sample tables are too large for the repository and will be deposited on Zenodo [PLACEHOLDER — DOI].

## Code availability

All scripts are in `projects/sampling-imbalance-biogeography-rules/analysis/` of the project repository [PLACEHOLDER — public repository URL and Zenodo DOI].

## Acknowledgements

[PLACEHOLDER]. We thank the EMP and MicrobeAtlas teams and the hundreds of investigators who made their data public.

## Author contributions

[PLACEHOLDER]. An AI agent (Claude, Anthropic) ran the analyses and drafted the manuscript under the authors' direction; [PLACEHOLDER — who checked what].

## Funding

[PLACEHOLDER]

## Competing interests

[PLACEHOLDER]

---

## Tables

**Table 1. Sample imbalance by country.** ⟦T1⟧

| Data set | Samples | Studies | Located in a country | Countries | Top country (share) | Top three (share) |
|---|---:|---:|---:|---:|---|---|
| EMP (QC-filtered, 5,000 reads) | 23,828 | 96 | 23,828 | 43 | United States (49.7%) | United States, Australia, Germany (70.4%) |
| MicrobeAtlas (amplicon, 1,000 reads) | 1,665,984 | 23,742 | 956,833 | 166 | United States (39.3%) | United States, China, Canada (57.9%) |

Shares are of samples located in a country; MicrobeAtlas also has 154,612 ocean samples and 554,539 without coordinates.

**Table 2. The ten rules re-tested.** Full verbatim quotes and locations in `findings.md`.

| ID | Rule | Source | Effect measured | Verdict |
|---|---|---|---|---|
| F1 | Richness peaks near neutral pH | Thompson 2017 [1] | Peak of upper-bound richness | Holds everywhere |
| F2 | Richness peaks at ~10 °C | Thompson 2017 [1] | Peak of upper-bound richness | EMP only |
| F3 | Host-associated less rich than free-living | Thompson 2017 [1] | log ratio of mean richness | Holds everywhere |
| F4 | Salinity is the main divider | Lozupone & Knight 2007 [2]; Thompson 2017 [1] | Jaccard separation vs pH and temperature splits | Not supported (as tested) |
| F5 | No latitude gradient in soil richness | Fierer & Jackson 2006 [3] | Slope on absolute latitude | Depends on data set and weighting |
| F6 | Similarity decays with distance | Hanson 2012 [10]; Martiny 2006 [9] | Slope of similarity on log distance | Holds (weaker under grid weighting) |
| F7 | Communities are nested | Thompson 2017 [1] | Prevalence–richness correlation (proxy) | MicrobeAtlas only |
| F8 | Herbivore guts richer than carnivore guts | Ley 2008 [4] | log ratio of mean richness | Holds (MicrobeAtlas; EMP too small) |
| F9 | Acidobacteria rise as soil pH falls | Jones 2009 [7]; Lauber 2009 [8] | Slope of logit share on pH | Hidden under sample weighting; appears after correction |
| F10 | Soil richness falls with elevation | Bryant 2008 [5] | Slope on elevation | Not supported |

## Figure legends

**Figure 1. Scorecard of ten biogeography rules across data sets and weightings.** Rows, rules; columns, data set (EMP; MicrobeAtlas at 1,000 and 5,000 reads) × weighting (every sample, study, country or 5° grid cell equal). Blue ✓, estimate agrees with the published rule; grey ○, 95% study-bootstrap interval includes zero; orange ✗, contradicts the rule. F1/F2 agree if both smoothed peak measures lie within 1 pH unit of 7 or 5 °C of 10 (cell text: range of the two peaks); F4 agrees if salinity separates communities more than both pH and temperature (cell text: ratio); F5, a null claim, agrees if the interval includes zero. Other cell text: point estimate. Source: `results/final_findings.csv`; cell values in `results/figures/fig1_scorecard_cells.csv`. Script: `analysis/make_figures.py`.

**Figure 2. The soil Acidobacteria–pH relationship is hidden by study imbalance.** **a**, Slope of logit(Acidobacteria relative abundance) on soil pH (95% study-bootstrap intervals) under four weightings and in a mixed model with a random intercept per study (diamonds), in three data sets. Blue, significant in the published (negative) direction; grey, interval includes zero. **b**, The same slope for Acidobacteria classes in EMP soils (≥ 1% mean abundance) under sample (circles) and grid (squares) weighting; orange, significant positive. Sources: `final_findings.csv` (a, weightings), `final_extra_*.json` (a, mixed models; b).

**Figure 3. Rules transfer to held-out countries.** Ratio of held-out-country prediction error to 10-fold cross-validation error. Orange line, failure threshold (1.5); a rule fails if one-third of countries exceed it. Circles, every tested country (F1–F3, protocol v0.3/v0.4 runs); diamonds, the five worst countries (F5, F9, F10, v0.6 runs; only these were stored). Right, countries above threshold / countries tested. Values above 3 are drawn at the axis edge. Sources: `results/findings_results.json`, `results/ma_findings.json`, `results/final_extra_*.json`.

**Supplementary Figure S1. All estimates with 95% intervals** (`results/figures/figS1_forest_all.png`). F1/F2 show the Laplace-mode peak, whose intervals are clipped for display; F4 shows the salinity separation D, coloured by the salinity ranking.

**Supplementary Table S1.** Samples, studies and countries per rule and data set (`results/final_findings.csv`, columns n, studies, countries; S0 rows). [To be generated as a formatted table.]

---

## References

*Checked against Europe PMC or Crossref on 23 Nov 2026 (`manuscript/references_verified.json`). Numbering follows first citation.*

1. Thompson LR, Sanders JG, McDonald D, Amir A, Ladau J, Locey KJ, et al. A communal catalogue reveals Earth's multiscale microbial diversity. *Nature* 2017;551:457–463. doi:10.1038/nature24621
2. Lozupone CA, Knight R. Global patterns in bacterial diversity. *Proc Natl Acad Sci USA* 2007;104:11436–11440. doi:10.1073/pnas.0611525104
3. Fierer N, Jackson RB. The diversity and biogeography of soil bacterial communities. *Proc Natl Acad Sci USA* 2006;103:626–631. doi:10.1073/pnas.0507535103
4. Ley RE, Hamady M, Lozupone C, Turnbaugh PJ, Ramey RR, Bircher JS, et al. Evolution of mammals and their gut microbes. *Science* 2008;320:1647–1651. doi:10.1126/science.1155725
5. Bryant JA, Lamanna C, Morlon H, Kerkhoff AJ, Enquist BJ, Green JL. Microbes on mountainsides: contrasting elevational patterns of bacterial and plant diversity. *Proc Natl Acad Sci USA* 2008;105 Suppl 1:11505–11511. doi:10.1073/pnas.0801920105
6. Fierer N, McCain CM, Meir P, Zimmermann M, Rapp JM, Silman MR, et al. Microbes do not follow the elevational diversity patterns of plants and animals. *Ecology* 2011;92:797–804. doi:10.1890/10-1170.1
7. Jones RT, Robeson MS, Lauber CL, Hamady M, Knight R, Fierer N. A comprehensive survey of soil acidobacterial diversity using pyrosequencing and clone library analyses. *ISME J* 2009;3:442–453. doi:10.1038/ismej.2008.127
8. Lauber CL, Hamady M, Knight R, Fierer N. Pyrosequencing-based assessment of soil pH as a predictor of soil bacterial community structure at the continental scale. *Appl Environ Microbiol* 2009;75:5111–5120. doi:10.1128/AEM.00335-09
9. Martiny JBH, Bohannan BJM, Brown JH, Colwell RK, Fuhrman JA, Green JL, et al. Microbial biogeography: putting microorganisms on the map. *Nat Rev Microbiol* 2006;4:102–112. doi:10.1038/nrmicro1341
10. Hanson CA, Fuhrman JA, Horner-Devine MC, Martiny JBH. Beyond biogeographic patterns: processes shaping the microbial landscape. *Nat Rev Microbiol* 2012;10:497–506. doi:10.1038/nrmicro2795
11. Gilbert JA, Jansson JK, Knight R. The Earth Microbiome project: successes and aspirations. *BMC Biol* 2014;12:69. doi:10.1186/s12915-014-0069-1
12. Abdill RJ, Adamowicz EM, Blekhman R. Public human microbiome data are dominated by highly developed countries. *PLoS Biol* 2022;20:e3001536. doi:10.1371/journal.pbio.3001536
13. Blake KS. Missing microbiomes: global underrepresentation restricts who research will benefit. *J Clin Invest* 2024;134:e183884. doi:10.1172/JCI183884
14. Van Daele W, Tito Tadeo RY, Perera J, Tshokey T, et al. Toward ethical human microbiome research: improving health through radical interdisciplinary and intercultural co-laboration. *Microbiome* 2026;14:107. doi:10.1186/s40168-026-02340-w
15. Matias Rodrigues JF, Tackmann J, Malfertheiner L, Patsch D, Perez-Molphe-Montoya E, Näpflin N, et al. The MicrobeAtlas database: global trends and insights into Earth's microbial ecosystems. *Cell* 2026;189:2092–2107.e17. doi:10.1016/j.cell.2026.01.021
16. Amir A, McDonald D, Navas-Molina JA, Kopylova E, Morton JT, Zech Xu Z, et al. Deblur rapidly resolves single-nucleotide community sequence patterns. *mSystems* 2017;2:e00191-16. doi:10.1128/mSystems.00191-16
17. Matias Rodrigues JF, Schmidt TSB, Tackmann J, von Mering C. MAPseq: highly efficient k-mer search with confidence estimates, for rRNA sequence analysis. *Bioinformatics* 2017;33:3808–3810. doi:10.1093/bioinformatics/btx517
18. Efron B. Better bootstrap confidence intervals. *J Am Stat Assoc* 1987;82:171–185. doi:10.1080/01621459.1987.10478410
19. Seabold S, Perktold J. Statsmodels: econometric and statistical modeling with Python. *Proc 9th Python in Science Conf* 2010:92–96. doi:10.25080/Majora-92bf1922-011
20. Roberts DR, Bahn V, Ciuti S, Boyce MS, Elith J, Guillera-Arroita G, et al. Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. *Ecography* 2017;40:913–929. doi:10.1111/ecog.02881

*Verified but not yet cited (candidates for the Discussion): Delgado-Baquerizo et al. 2018 Science; Bahram et al. 2018 Nature; Shaffer et al. 2022 Nat Microbiol. Natural Earth (data source) still needs a version/date.*
