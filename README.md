# sc-rna-atac-regulon-benchmark

A reproducible single-cell multiome project extending an RNA-only immune cell-state benchmark into paired RNA+ATAC integration, regulatory program analysis, and optional pseudotime exploration.

This repository uses the 10x PBMC Multiome dataset to demonstrate RNA/ATAC quality control, Seurat/Signac preprocessing, WNN integration, RNA-derived TF activity proxy analysis, ATAC motif validation with chromVAR, and monocyte pseudotime analysis.

## Why this project matters

Single-cell RNA-seq identifies transcriptional cell states, while ATAC-seq adds chromatin accessibility context. This project compares RNA-defined immune programs with ATAC-derived regulatory signals in a shared-cell multiome workflow.

The emphasis is on reproducible computational biology, transparent QC, and careful interpretation rather than overclaiming biological discovery.

## What this repo demonstrates

- Paired RNA+ATAC preprocessing and multimodal QC
- WNN integration for joint RNA+ATAC cell-state recovery
- Comparison against an RNA-only single-cell benchmark
- RNA-derived TF activity proxy analysis across immune populations
- ATAC motif enrichment using JASPAR/chromVAR
- Cross-modal RNA-vs-ATAC regulatory program comparison
- Optional monocyte pseudotime analysis using Slingshot
- Python/R package structure with tests and environment files

## Key results

- Recovered 7 of 8 broad RNA-defined immune populations after WNN label harmonization.
- Identified expected immune regulatory patterns:
  - antigen presentation in dendritic cells and B cells
  - cytotoxic programs in NK and CD8 T cells
  - inflammatory and monocyte-activation programs in monocytes
- Added traceable multimodal QC outputs with one failure flag per active QC criterion.
- Used ATAC motif deviations to independently evaluate RNA-derived TF activity proxies.
- Built pseudotime-ready monocyte outputs for linking chromatin accessibility with regulatory dynamics.

## Workflow overview

| Phase | Goal | Main scripts |
|---|---|---|
| 1 | Import RNA-only baseline outputs | scripts/01_import_rna_baseline.py |
| 2 | Preprocess RNA+ATAC multiome data | scripts/03_preprocess_multiome.R |
| 2 | Plot QC and optional doublet/fragment QC | scripts/03b_plot_qc.R, scripts/03c_doublet_qc.R, scripts/03d_plot_atac_fragment_qc.R |
| 2 | Run WNN integration | scripts/04_integrate_rna_atac_wnn.R |
| 3 | Compare RNA and WNN cell labels | scripts/05_compare_rna_vs_wnn.py |
| 3 | Compute and plot TF activity proxies | scripts/06_tf_activity.py, scripts/07_plot_tf_activity.py |
| 4 | Run ATAC motif enrichment | scripts/08_motif_enrichment_atac.R |
| 4 | Compare RNA and ATAC TF signals | scripts/09_compare_rna_atac_tf.py |
| 5 | Run monocyte pseudotime analysis | scripts/10_pseudotime_subset.R, scripts/11_plot_pseudotime.py |

## Quick start

Create the environment:

    mamba env create -f environment.yml
    conda activate sc-rna-atac

Install the Python package in editable mode:

    pip install -e .

Run tests:

    pytest

Run the main workflow:

    python scripts/01_import_rna_baseline.py
    Rscript scripts/02_download_or_prepare_multiome.R --download
    Rscript scripts/03_preprocess_multiome.R
    Rscript scripts/03b_plot_qc.R
    Rscript scripts/04_integrate_rna_atac_wnn.R --cluster-annotations results/integration/cluster_annotations.tsv
    python scripts/05_compare_rna_vs_wnn.py --label-map results/comparison/wnn_to_broad_labels.tsv
    python scripts/06_tf_activity.py
    python scripts/07_plot_tf_activity.py
    Rscript scripts/08_motif_enrichment_atac.R
    python scripts/09_compare_rna_atac_tf.py
    Rscript scripts/10_pseudotime_subset.R
    python scripts/11_plot_pseudotime.py

Optional QC scripts:

    Rscript scripts/03c_doublet_qc.R
    Rscript scripts/03d_plot_atac_fragment_qc.R

## Data

This project uses:

- RNA-only baseline outputs from the companion sc-cell-state-benchmark repository
- 10x Genomics PBMC Multiome RNA+ATAC data

Expected input locations:

| Path | Purpose |
|---|---|
| data/external/ | RNA-only baseline score tables |
| data/raw/ | raw 10x multiome files |
| data/processed/ | processed Seurat objects |

Large raw and processed data files are not intended to be committed.

## Main outputs

| Directory | Contents |
|---|---|
| results/baseline/ | standardized RNA-only baseline outputs |
| results/integration/ | QC summaries, WNN metadata, integration outputs |
| results/regulons/ | TF activity and ATAC motif summaries |
| results/comparison/ | RNA-vs-WNN and RNA-vs-ATAC comparison tables |
| results/dynamics/ | monocyte pseudotime outputs |
| results/figures/ | QC, TF activity, motif, and pseudotime plots |

## Method summary

RNA preprocessing uses SCTransform with mitochondrial percentage regression. ATAC preprocessing uses TF-IDF and LSI. Joint RNA+ATAC structure is modeled using WNN integration.

Regulatory analysis is intentionally lightweight and interpretable. RNA program scores are used as TF activity proxies, then compared against ATAC motif deviation scores from chromVAR. This provides a cross-modal regulatory consistency check without claiming full GRN reconstruction.

## Interpretation notes

This is a portfolio and methods-development repository, not a clinically validated biological study.

The TF analysis should be interpreted as regulatory program scoring and motif-based validation, not full causal regulon inference. Pseudotime analysis is exploratory and currently focused on monocyte-lineage cells.

## Environment notes

The full project environment is defined in environment.yml and includes Python 3.11, R 4.4, Seurat, Signac, chromVAR, Slingshot-related dependencies, and Bioconductor packages.

pyproject.toml defines the Python package metadata, dependencies, and test configuration.

renv.lock records the R package state used during development.

For detailed installation notes, troubleshooting, and R/Bioconductor caveats, see docs/setup.md.

## Repository structure

    sc-rna-atac-regulon-benchmark/
    ├── README.md
    ├── environment.yml
    ├── pyproject.toml
    ├── renv.lock
    ├── data/
    ├── docs/
    ├── results/
    ├── scripts/
    ├── src/sc_rna_atac_regulon_benchmark/
    └── tests/

## Related repository

This project builds on the RNA-only sc-cell-state-benchmark repository, which provides baseline immune cell-state outputs used here for multimodal comparison.

## License

MIT
