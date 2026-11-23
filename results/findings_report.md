# N1 results: F1–F3 under four weighting schemes

*21 Nov 2026. Protocol v0.3 (commit 794a150, made before F2/F3 were run). Script: `analysis/run_findings.py`. Output: `findings_results.json`, `findings_stdout.txt`, figure `findings.png`. First-pass analysis, not yet checked. 95% CIs from 1,000 study-level bootstraps.*

**Sample sets match the paper exactly:** F1 n = 3,986 (34 studies, 18 countries); F2 n = 6,976 (22 studies, 11 countries); F3 n = 23,828 (96 studies, 43 countries).

## F1: richness vs pH

| Scheme | Upper-bound peak pH (published claim) | Mean-model curvature |
|---|---|---|
| S0 standard | **7.37** (6.89–8.12) | −26.5 (−66.6 to −2.8) |
| S1 study | **7.44** (6.64–8.40) | −19.6 (−62.4 to −9.5) |
| S2 country | **7.21** (6.58–7.99) | −4.5 (−54.0 to **+3.8**) |
| S3 grid | **7.19** (6.45–8.58) | −13.5 (−45.5 to **+0.3**) |

- **Published claim: unchanged** under all schemes (largest shift 0.18 pH).
- **Average pattern: changed** under S2 and S3. The hump in typical richness is no longer distinguishable from no relationship once countries or grid cells are weighted equally.
- **Transfer:** fails in 1 of 11 countries (England 1.55) → no failure by the rule (needs one-third).
- **Sensitivity:** free-living only 7.36 / 7.18 (S0 / S3); 150 bp 7.04 / 6.90; soil only 6.55 / 6.53.

## F2: richness vs temperature

| Scheme | Upper-bound peak °C (published claim) |
|---|---|
| S0 standard | **9.1** (6.2–71.1) |
| S1 study | **9.9** (8.1–70.2) |
| S2 country | **9.1** (5.6–70.1) |
| S3 grid | **69.3** (9.0–70.6) |

- The S0 estimate (about 9 °C) reproduces the paper's "about 10 °C".
- **Published claim: changed under S3 by the pre-set rule.** But this rests on **one 10-sample bin at about 69 °C** whose 99th-percentile richness (1,846) beats the cool bins once grid-balancing lowers them (from about 2,500 to about 1,500). The wide upper CIs in every scheme show the same instability.
- **Post hoc check (bins ≥ 20 samples, added after the first run):** S0 9.0, S1 9.8, S2 9.2, **S3 11.7 °C**: no change. 150 bp table: 13.0 / 13.6 °C: no change.
- **Average pattern:** there is no hump in average richness against temperature even under S0 (curvature CI includes zero in every scheme), so "average changed" is not meaningful for F2.
- **Transfer:** fails in 2 of 10 countries (**Australia 2.78**, England 1.71) → no failure by the rule, but Australia (1,559 samples) is predicted very badly from the rest.

## F3: free-living vs host-associated richness

| Scheme | log(free-living ÷ host) | Free-living mean | Host mean |
|---|---|---:|---:|
| S0 standard | **0.79** (0.43–1.26) | 648 | 294 |
| S1 study | **1.00** (0.78–1.29) | 649 | 239 |
| S2 country | **1.12** (0.78–1.54) | 674 | 221 |
| S3 grid | **1.05** (0.86–1.31) | 643 | 224 |

- **Published claim: unchanged** (largest change 41%, below 50%). Rebalancing makes the gap **bigger**: free-living communities go from about 2.2× to about 2.9× richer than host-associated ones.
- The rhizosphere exception holds (mean 1,397 S0, 1,248 S3).
- **Transfer:** fails in 1 of 32 countries (Norway 1.52) → no failure.

## Overall against the hypotheses

| | F1 | F2 | F3 |
|---|---|---|---|
| Published claim changed (pre-set rule) | No | **Yes, S3 only** (fragile: one 10-sample bin; gone with ≥ 20-sample bins) | No (gap grows) |
| Average pattern changed | **Yes** (S2, S3) | Not meaningful (no average pattern) | – |
| Transfer fails (H2) | No | No (but Australia 2.8×) | No |

- **H1 (at least one finding changes materially):** formally yes (F2 under S3; F1's average pattern), but the published headline claims are **largely robust**.
- **H2 (transfer fails):** not supported by the rule; isolated large failures (Australia for temperature).

## What this means for the paper

The honest headline so far: **the EMP's headline biogeographic claims mostly survive rebalancing.** What does not survive is (a) the average-level pH pattern, and (b) temperature predictions for under-represented regions such as Australia. One effect (free-living vs host richness) is *understated* by the current imbalance.

That is a solid robustness paper (*Global Ecology & Biogeography*, *mSystems*, *ISME Communications*), not yet a *Nature Ecology & Evolution* paper. Protocol section 6 says that in this case we either write up the robustness result or extend to MicrobeAtlas, where imbalance is larger.

## Caveats

- S2 gives single samples from small countries very large weight; S3 is the main rebalancing for interpretation (amendment 0.3).
- The bins ≥ 20 check is post hoc and is reported as such.
- Findings are still taken from a summary of the paper; check the exact published wording and peak values against the PDF.
