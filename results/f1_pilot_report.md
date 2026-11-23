# N1 pilot report: F1 (richness vs pH)

*21 Nov 2026. Protocol v0.2 (commit 5274007, made before this run). Script: `analysis/f1_pilot.py`. Raw output: `f1_pilot.json`, figure: `f1_pilot.png`. First-pass analysis, not yet checked.*

## Data

- **3,974** free-living EMP samples with pH (the paper used 3,986; the 12-sample gap is not yet explained), from **34 studies** in **18 countries**.
- **USA 53%** (2,117) and **Germany 26%** (1,040) together make up 80% of the samples.
- Richness = observed 90 bp Deblur tag sequences at 5,000 reads.

## Results (95% CIs from 1,000 study-level bootstraps)

| Scheme | Upper-bound peak pH (as published) | Mean-model peak pH | Mean-model curvature |
|---|---|---|---|
| S0 standard | **7.36** (6.89–8.12) | 5.94 (1.38–6.88) | −26.7 (−67.6 to −3.0) |
| S1 study-balanced | **7.44** (6.64–8.42) | 5.89 (3.53–6.79) | −20.7 (−64.8 to −10.4) |
| S2 country-balanced | **7.21** (6.58–7.98) | no clear peak | −4.5 (−53.8 to **+3.7**) |
| S3 grid-balanced | **7.18** (6.45–8.59) | 4.37 (−10.5 to 6.6) | −14.0 (−46.1 to −0.3) |

**Transfer (leave-one-country-out):** held-out error was 1.00–1.18 times cross-validation error in all 11 countries tested (Germany highest at 1.18). No country passed the 1.5 failure line.

**Sensitivity:** soil only (n = 935): upper-bound peak 6.55 (S0) and 6.55 (S1). 150 bp table (n = 2,080): 7.04 (S0) and 7.04 (S1).

## Against the pre-set decision rules

| Rule | S1 | S2 | S3 |
|---|---|---|---|
| Upper-bound peak moves > 0.5 pH | No (+0.08) | No (−0.16) | No (−0.18) |
| Mean peak moves > 0.5 pH | No | **Yes** (no peak) | **Yes** (−1.6) |
| Curvature CI includes zero | No | **Yes** | No |
| Curvature changes > 50% | No (22%) | **Yes** (83%) | No (48%) |
| **Changed?** | No | **Yes** | **Yes** |
| Transfer fails (H2) | No | | |

**By the protocol's rules, F1 counts as changed under S2 and S3, so the pilot is a GO** for F2 and F3.

## What it means

1. **The published claim holds.** The upper bound of richness peaks near pH 7.2–7.4 whatever the weighting. Thompson et al.'s finding is robust to sampling imbalance.
2. **The average pattern does not.** The hump in *typical* richness against pH comes mostly from the over-sampled USA and Germany. With countries weighted equally, it flattens until it is no longer distinguishable from no relationship. This is the kind of "conclusion changes with rebalancing" result N1 is looking for, but it concerns a pattern the paper did not explicitly claim.
3. **The transfer test is weak as specified.** pH alone explains little of the variation in richness (RMSE about 400 tag sequences), so a near-flat model transfers equally badly everywhere. A model with habitat plus pH would be a fairer test. This needs a dated protocol amendment before F2/F3 are run.

## Caveats

- Country-balanced weights give single samples from small countries the same total weight as the USA's 2,117 samples. This makes S2 noisy by design; S3 (grid) is the more stable rebalancing.
- The decision rules did not say whether they apply to the upper bound, the mean, or both. Both were applied here; the next protocol amendment should state this explicitly.
- The 12-sample gap with the paper's n, and the exact published peak value, still need checking against the PDF.
