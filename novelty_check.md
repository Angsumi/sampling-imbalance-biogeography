# N1 focused novelty check

*23 Nov 2026. Provisional, AI-run; needs a person's check. Searches: Europe PMC REST API (queries below), run 23 Nov 2026. Builds on the stage-2 check in `N1_protocol.md` §2.*

## Question

Has anyone already tested whether published microbial biogeography rules (in particular the EMP findings of Thompson et al. 2017) change when samples are re-weighted by study, country or map cell, or when models are tested on held-out countries?

## Searches run (Europe PMC, 23 Nov 2026)

| # | Query (abridged) | Hits | Relevant |
|---|---|---:|---|
| 1 | "Earth Microbiome Project" AND (bias OR imbalance OR reweight* OR geographic), 2019–2026 | 902 | None that re-weight EMP; top hits use EMP as a reference set |
| 2 | "MicrobeAtlas" | 29 | Matias Rodrigues et al. 2026 *Cell* (database paper) — see below |
| 3 | "Missing microbiomes" | 4 | Blake 2024 *J Clin Invest* (commentary) |
| 4 | (microbiome OR microbial) AND ("sampling bias" OR "geographic bias" OR "sampling imbalance") AND (biogeograph* OR "global patterns"), 2020–2026 | 249 | None (pathogen genomics, plankton atlases) |
| 5 | ("leave-one-country-out" OR "leave-one-region-out" OR "spatial cross-validation") AND microbial, 2018–2026 | 70 | None for community-ecology rules; used in AMR mapping and clinical prediction |
| 6 | EMP AND richness AND (pH OR temperature) AND ("upper bound" OR Laplace OR re-analysis), 2018–2026 | 35 | None |
| 7 | (microbial) AND ("study effect" OR "batch effect" OR Simpson) AND (biogeograph* OR macroecolog*), 2018–2026 | 1,014 | None in top results; no re-test of published rules |
| 8 | TITLE "Earth Microbiome Project", 2019–2026 | 2 | Oral saliva re-analyses; not relevant |
| 9 | ABSTRACT nestedness AND microbial AND global, 2018–2026 | 11 | Shaffer et al. 2022 *Nat Microbiol* (EMP500 multi-omics; new data, not a re-weighting) |
| 10 | Microbiome-journal titles on under-representation, 2025–2026 | 14 | Van Daele et al. 2026; Passmore et al. 2026 (ethics/equity commentaries) |

## Closest work

| Work | What it does | Overlap with N1 |
|---|---|---|
| Abdill, Adamowicz & Blekhman 2022, *PLoS Biol* 20:e3001536 | Maps the country of origin of public human microbiome samples: ">71% of samples with a known origin come from Europe, the United States, and Canada" (abstract) | Documents the skew (human samples only); does not test whether it changes ecological conclusions |
| Blake 2024, *J Clin Invest* 134:e183884 ("Missing microbiomes") | Commentary on global under-representation | Argument; no re-analysis |
| Van Daele et al. 2026; Passmore et al. 2026, *Microbiome* | Equity and ethical-design commentaries | Argument; no re-analysis |
| Matias Rodrigues et al. 2026, *Cell* 189:2092 (MicrobeAtlas) | Database of > 2 million harmonised samples; reports "global trends", "geographic preferences" and "taxonomy-specific generalism trends" (abstract) | **Main risk.** Abstract does not mention re-testing EMP rules or re-weighting. **Full text could not be retrieved (Europe PMC 403, bioRxiv 429, Elsevier JS redirect); a person must read it** and check whether it reports richness vs pH/temperature/latitude or any study/country weighting |
| Thompson et al. 2017, *Nature* | Original findings. Its nestedness test already used "a 2,000-sample subset with even representation across environments and studies"; its latitude analysis used a mixed model with study as random effect | N1 extends this idea to all ten rules, three weightings, a second dataset and transfer tests. The manuscript must credit these existing balancing steps |
| Regional soil studies (Europe, Australia; protocol §2) | Show that pH/climate effects vary by region | Single regions; no cross-region transfer test of published rules |

## Verdict

**Passes, provisionally.** No published study was found that re-estimates a set of published microbial biogeography rules under study-, country- and grid-balanced weighting in two global datasets and tests transfer to held-out countries. The specific result that the soil Acidobacteria–pH relationship is masked at phylum level by opposing class-level trends and recovered after correcting for study imbalance was not found in the literature searched; class-level differences in pH response within Acidobacteria are known (Jones et al. 2009), so the manuscript should present the novelty as *the masking by sampling imbalance*, not the class differences themselves.

**Open before submission (person to do):**
1. Read the MicrobeAtlas *Cell* 2026 full text (and supplement) for any richness–environment or weighting analysis.
2. A Google Scholar / Web of Science pass on "re-weighting" + "Earth Microbiome Project" and on citing articles of Thompson 2017 from 2024–2026 (Europe PMC covers citing articles only partly).
