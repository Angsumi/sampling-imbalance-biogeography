# N1: all ten findings (F1–F10), both datasets, both read depths

*22 Nov 2026. Protocol v0.5 (commit 861668d, made before F4–F10 and the 5,000-read check). Scripts: `analysis/all_findings.py` (with `ma_rarefy.py`, `ma_prepare.py`, `ma_findings.py`, `run_findings.py`). Full table: `all_findings.csv` (120 rows: finding × dataset × depth × scheme, with estimate, 95% CI, n, studies, countries, changed vs S0). First-pass analysis, not yet checked.*

- **Datasets:** EMP (90 bp, 5,000 reads, 1,000 bootstraps for F1–F3 and 200 for F4–F10); MicrobeAtlas at 1,000 reads and at 5,000 reads (200 project bootstraps).
- **Weightings:** S0 every sample counts once; S1 every study equal; S2 every country equal; S3 every 5° map square equal.

## One-line verdicts

| # | Published rule | Holds at S0 (as usually analysed)? | Survives rebalancing? | Depth-sensitive (MicrobeAtlas 1,000 vs 5,000)? |
|---|---|---|---|---|
| F1 | Richness peaks near neutral pH | Yes: EMP 7.4; MicrobeAtlas 6.4 / 5.9 | **Peak moves**: 6.3 → 7.2 (1,000) and 5.9 → 7.0 (5,000) with grid weighting | Borderline (6.4 vs 5.9) |
| F2 | Richness peaks at about 10 °C | EMP yes (9 °C); **MicrobeAtlas no** (no clear peak; 29 / 38 °C) | **No**: estimate jumps 13–38 °C between schemes | No real peak at either depth |
| F3 | Free-living richer than host-associated | Yes, everywhere | **Yes** (direction never changes); size varies about 2× | Size larger at 5,000 (0.31 → 0.51) |
| F4 | Salinity is the strongest divider | **No in EMP** (pH or temperature divide more: ratio 0.56); MicrobeAtlas borderline (0.84 / 1.03) | **No**: the ordering flips with weighting in MicrobeAtlas | Yes (ratio crosses 1) |
| F5 | No latitude gradient in soil richness | **No in EMP** (richness falls away from the equator); yes in MicrobeAtlas (CI includes zero) | **Mixed**: study and country weighting create a gradient in MicrobeAtlas (in opposite directions) | No |
| F6 | Similarity decays with distance | Yes, everywhere | **Mostly**: CI touches zero under grid weighting in MicrobeAtlas | No |
| F7 | Poor communities hold widespread taxa (nestedness) | Yes, everywhere | **EMP no** (fades to zero under country weighting); **MicrobeAtlas yes** | Weaker at 5,000 (−0.64 → −0.32) |
| F8 | Herbivore guts richer than carnivore guts | EMP unclear (128 carnivore samples, 7 studies); **MicrobeAtlas yes** (about 2–2.8× richer) | **Yes** in MicrobeAtlas, every scheme | No |
| F9 | Acidobacteria rise as soil pH falls | **No** in any dataset (CI includes zero) | **Appears only after rebalancing**: significant negative slope under S1 / S3 in all three datasets | No |
| F10 | Soil richness falls with elevation | **No** in any dataset | **No**: sign unstable, never significant | No |

## What this says

1. **Four rules hold up well:** F3 (free-living richer), F6 (distance decay), F8 (herbivore guts richer, in the larger data) and, in MicrobeAtlas, F7 (nestedness).
2. **Three published claims are not supported by these data at all**, before or after rebalancing: F2's 10 °C optimum (in the larger data), F4's "salinity is the strongest divider", and F10's elevation decline.
3. **Three findings depend on which places dominate the data:**
   - F1: where the pH peak sits.
   - F5: whether a latitude gradient appears.
   - F9: the Acidobacteria–pH link is *hidden* by the standard analysis and appears once over-sampled regions are down-weighted, consistently in all three datasets. This is the clearest example of imbalance masking a real pattern.
4. **Read depth matters for size, not usually direction:** F3 and F7 change size between 1,000 and 5,000 reads; only F4's ordering is depth-sensitive.

**The paper story:** "Sampling imbalance and analysis choices decide several textbook rules of microbial biogeography." The strongest single result is F9 (a pattern visible only after rebalancing, in three datasets). The broadest is the scorecard itself: 10 rules, 3 datasets, 4 weightings.

## Caveats

- Upper-bound peak fits (F1, F2) are unstable where richness is flat. Some bootstrap CIs are very wide, and a few percentile CIs do not contain the point estimate (e.g. F1 MicrobeAtlas 5,000 S0, F6 MicrobeAtlas 1,000 S0); this needs a closer look.
- F4 uses Jaccard distances, which are close to 1 for most pairs at fine resolution, so separations are small numbers; the ordering between factors is what matters.
- F8 diet groups come from a fixed host list; MicrobeAtlas gut samples are found by keywords in the body-site field.
- MicrobeAtlas mixes 16S regions and primers; S1 absorbs some of this, but region is not modelled.
- Transfer tests were not run for F4–F10 (protocol v0.5).
- Published claims are taken from summaries; check wording and numbers against each paper before writing.
