# Geographic sampling imbalance and microbial biogeography

Manuscript, analysis code, protocol, supporting results and figures for a reanalysis of Earth Microbiome Project (EMP) and MicrobeAtlas data. The manuscript is a working draft; see `STATUS.md` for outstanding verification and submission details.

## Contents

- `manuscript/`: manuscript draft, reference metadata, supplementary table and `submission/` drafts (title page, cover letter, declarations and checklist).
- `N1_protocol.md`, `findings.md`, `novelty_check.md`: protocol and source checks.
- `analysis/`: data preparation, analyses and figure generation.
- `results/`: recorded analysis outputs and manuscript figures.
- `data/`: smaller metadata and Natural Earth boundaries. Large source data and generated caches are excluded from GitHub; the original local data are retained separately.

## Data setup

The original inputs were retrieved on 21 November 2026. Obtain EMP release 1 from <https://ftp.microbio.me/emp/release1/> and MicrobeAtlas MAPv3/MAPseq v1.2.6 from <https://microbeatlas.org/downloads/>. Sources and dataset selection are detailed in `N1_protocol.md`. Place the following files at the indicated paths (the MicrobeAtlas BIOM must be decompressed and renamed as shown):

```text
data/emp_release1/emp_deblur_90bp.qc_filtered.rare_5000.biom
data/emp_release1/emp_deblur_150bp.qc_filtered.rare_5000.biom
data/emp_release1/emp_deblur_150bp.qc_filtered.biom
data/microbeatlas/samples-otus.97.metag.minfilter.biom
data/microbeatlas/samples.env.info.tsv
data/microbeatlas/samples.info.latlon.parsed.tsv
data/microbeatlas/taxonomic_mapping.tsv
```

The MicrobeAtlas source BIOM is `samples-otus.97.mapped.metag.minfilter.refilt.biom.gz`. The metadata and boundaries already included here are also inputs. Some original auxiliary data are listed in the protocol. Large data are not archived in this repository; obtain the matching source release, or use the original local copies. No archival DOI or checksums have yet been established.

## Reproduction

Use Python 3.12 and install `requirements.txt` in a virtual environment. These are computationally intensive analyses; the MicrobeAtlas preparation loads large arrays and uses 16 worker processes. The original run produced the committed results; independent reproduction and code review remain outstanding.

```sh
python -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
python analysis/ma_prepare.py
python analysis/ma_rarefy.py
python analysis/run_findings.py
python analysis/ma_findings.py
python analysis/all_findings.py emp
python analysis/all_findings.py ma1000
python analysis/all_findings.py ma5000
python analysis/all_findings.py merge
python analysis/n1_extra.py emp
python analysis/n1_extra.py ma1000
python analysis/n1_extra.py ma5000
python analysis/make_figures.py
```

The earlier F1 pilot can be run with `python analysis/f1_pilot.py`. Scripts overwrite corresponding results and figures. Inspect the protocol and scripts for parameters and assumptions before running.
