# N1 MicrobeAtlas extension: F1–F3 in 1.9 million samples

*21 Nov 2026. Protocol v0.4 (commit 8ec1924, made before this analysis). Scripts: `analysis/ma_prepare.py`, `analysis/ma_findings.py`. Output: `ma_findings.json`, `ma_findings_stdout.txt`. First-pass analysis, not yet checked. 95% CIs from 200 project-level bootstraps.*

## Data

- MicrobeAtlas MAPv3, 97% OTUs, filtered table: **1,884,129 samples** × 102,824 OTUs. Every sample rarefied to **1,000 reads**; amplicon samples only in the main analysis.
- Countries assigned from coordinates (Natural Earth, 1:110m). Points outside any country count as `ocean`.

| Finding | Samples | Projects | Countries | Top country | Top 3 share | EMP had |
|---|---:|---:|---:|---|---:|---|
| F1 pH | 43,890 | 571 | 82 | USA 39% | 57% | 3,986 / 34 studies |
| F2 temperature | 54,021 | 778 | 80 | USA 31% (ocean 24%) | 66% | 6,976 / 22 studies |
| F3 free-living vs host | 1,342,615 | 18,618 | 157 | USA 35% | 61% | 23,828 / 96 studies |

## F1: richness vs pH: replicated; the peak position depends on weighting

| Scheme | Upper-bound peak pH (claim) | Mean-model curvature |
|---|---|---|
| S0 standard | **6.42** (6.01–7.07) | −2.42 (−3.35 to −1.70) |
| S1 project | **6.26** (5.98–7.08) | −1.46 (−1.90 to −1.19) |
| S2 country | **6.76** (6.38–7.18) | −1.85 (−2.38 to −1.23) |
| S3 grid | **7.20** (6.50–7.49) | −1.09 (−1.63 to −0.74) |

- **Replication:** the near-neutral peak is confirmed in 11× more samples, but the standard estimate (6.42) sits about 1 pH unit lower than EMP's (7.37).
- **Claim changed under S3 by the pre-set rule** (+0.78 pH). Grid-balancing moves the peak from slightly acidic **back to EMP's value**. Where the peak appears to sit depends on which regions dominate the data.
- **Average pattern:** the hump survives every scheme (curvature CI never includes zero) but halves under S3 (−55%) → "changed" under S3 by rule.
- Sensitivity: free-living only 6.37 / 7.15 (S0 / S3); shotgun samples 6.92 / 6.97 (n = 1,095).
- Transfer: 1 of 25 countries fails (Spain 2.05×) → no failure by rule.

## F2: richness vs temperature: the EMP peak does not replicate

| Scheme | Upper-bound peak °C (claim) |
|---|---|
| S0 standard | **28.9** (8.8–39.0) |
| S1 project | **16.8** (5.0–20.9) |
| S2 country | **13.7** (−19.1–20.9) |
| S3 grid | **18.3** (10.2–20.9) |

- **Claim changed under S1, S2 and S3** (shifts of 11–15 °C, threshold 5 °C).
- But the bin profile shows why: the 99th-percentile richness is **flat at about 140–150 from 3 °C to 40 °C** under S0, with noisy single-bin spikes. There is no sharp optimum, so the fitted "peak" lands wherever the noise is highest, and moves with weighting. Under S3, richness clearly falls above about 40 °C.
- **Interpretation:** EMP's "peak at about 10 °C" is **not reproduced** in the larger data. What holds is a broad plateau from cold to warm, with a decline at hot temperatures.
- Average pattern: no curvature under S0 (CI includes zero), so "average changed" is not meaningful.
- Transfer: 1 of 28 fails (Antarctica 1.85×) → no failure by rule.

## F3: free-living vs host-associated: replicated, size depends on weighting

| Scheme | log(free-living ÷ host) | Free-living | Host |
|---|---|---:|---:|
| S0 | **0.31** (0.28–0.35) | 66 | 48 |
| S1 project | **0.16** (0.15–0.18) | 65 | 55 |
| S2 country | **0.42** (0.29–0.45) | 65 | 43 |
| S3 grid | **0.26** (0.21–0.31) | 61 | 47 |

- **Direction replicated** in every scheme; never changed by rule. But the size swings from 1.18× (S1) to 1.51× (S2); S1 comes close to the 50% line (−47%).
- Soil is the richest top-level habitat, animals the poorest (S0 means 71 vs 47 OTUs at 1,000 reads).
- Transfer: 0 of 87 countries fail.

## Overall

| | EMP (v0.3) | MicrobeAtlas (v0.4) |
|---|---|---|
| F1 claim (pH peak) | Robust (7.2–7.4) | Replicated; **peak moves 6.3 → 7.2 with weighting** (S3 changed) |
| F2 claim (temp peak) | Robust at about 9 °C (one fragile bin) | **Not replicated**: no clear peak; estimate moves 14–29 °C with weighting |
| F3 claim | Robust; gap grows with rebalancing | Replicated; **size varies about 2.5-fold across schemes** |
| Transfer (H2) | No failure | No failure |

**The emerging paper story:** in a dataset 14× larger and more unevenly sampled, rebalancing matters more than in EMP. The *direction* of each finding holds, but the *numbers* people quote (where richness peaks, how much richer free-living communities are) move by amounts that exceed the pre-set thresholds. EMP's temperature optimum does not replicate at all.

## Caveats

- **Rarefaction depth differs:** 1,000 reads here vs 5,000 in EMP. Lower depth compresses richness (maximum about 150–180 here), which flattens curves, especially for temperature. The protocol fixed 1,000 because the filtered table's minimum is 1,000 reads. A 5,000-read subset is a needed sensitivity check.
- MicrobeAtlas mixes 16S regions and primers. Project weighting (S1) absorbs some of this, but region is not modelled.
- `ocean` (24% of temperature samples) is a single "country" in S2.
- 200 bootstraps (compute limit), not 1,000.
- No figure yet for this extension.
