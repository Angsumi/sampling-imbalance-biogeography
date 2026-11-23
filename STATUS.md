# N1 status

**Stage 6 (writing), 23 Nov 2026.** First full manuscript draft written. All AI work is provisional until a person has checked it.

## Where things are

| Item | File | State |
|---|---|---|
| Protocol | `N1_protocol.md` (v0.6) | Frozen; no new analyses since v0.6 |
| Final results | `results/final_report.md`, `results/final_findings.csv`, `results/final_extra_*.json` | Complete (first pass) |
| Earlier runs (record) | `results/findings_*`, `results/ma_findings*`, `results/all_findings*`, `results/f1_pilot*` | Restored 23 Nov 2026 from git history (dropped by the repo reorganisation, commit f03b35f) |
| Source wording F1–F10 | `findings.md` | Verified: Thompson 2017 word for word (PMC full text); F4–F6, F8–F10 against source abstracts |
| Novelty | `novelty_check.md` | Passes provisionally; MicrobeAtlas *Cell* 2026 full text not read |
| Figures | `results/figures/` (script `analysis/make_figures.py`) | Fig. 1 scorecard, Fig. 2 F9, Fig. 3 transfer, Fig. S1 forest, Table 1 imbalance |
| Manuscript | `manuscript/draft_v1.md` | Draft v1 with source tags on every number |
| References | `manuscript/references_verified.json` | 23 checked against Europe PMC / Crossref; Natural Earth version missing |

## Needs the user's input

- Funding, competing interests, author contributions and acknowledgements; authors and affiliations have been added to the draft and title page.
- Target journal (draft follows Nature Communications; abstract ~180 words, limit 150).
- Archival DOI for data and code; the public repository URL is included in the draft.

## Needs a person's check

1. Every number in `draft_v1.md` against its ⟦source⟧ tag; every quote in `findings.md` against the papers.
2. **Flags in `findings.md`:** F7 is a proxy (not NODF); F4 contrasts are moderate, not extremes; F10's source (Bryant 2008) may concern Acidobacteria only — confirm from the PDF; MicrobeAtlas F1/F2 include host-associated samples (24% of temperature samples).
3. Read MicrobeAtlas *Cell* 2026 (doi:10.1016/j.cell.2026.01.021) for overlapping analyses.
4. Decide whether to run, under a dated amendment, the free-living-only MicrobeAtlas F1/F2 sensitivity analysis listed in v0.4 but absent from the final tables.
5. Code review of `analysis/all_findings.py` and `analysis/n1_extra.py` (first-pass analysis).
