# sc-rna-atac-regulon-benchmark

Single-cell RNA + ATAC analysis of paired 10x Genomics PBMC Multiome data using Seurat and Signac.

The project focuses on understanding how transcriptomic and chromatin-accessibility information can be analyzed independently and combined using Weighted Nearest Neighbor (WNN) integration.

## Biological motivation

scRNA-seq measures the transcriptional state of a cell, while scATAC-seq measures regions of accessible chromatin that may support regulatory activity.

These modalities provide complementary information: chromatin accessibility can indicate regulatory potential, but accessibility alone does not imply transcription. Joint analysis can therefore provide additional context for interpreting cell identity and regulatory state.

## Analysis workflow

For the paired multiome data, RNA and ATAC are first processed separately:

```text
RNA counts                         ATAC peak counts
    │                                    │
SCTransform                              TF-IDF
    │                                    │
   PCA                                  LSI
    │                                    │
    └────────────── WNN ─────────────────┘
                       │
                 WNN graph
                       │
              clustering + UMAP
```

WNN combines the modality-specific cell representations by learning how informative RNA and ATAC are for the neighborhood of each cell, rather than concatenating genes and peaks into a single feature matrix.

For ATAC, the first LSI dimension is excluded from WNN because it is often strongly associated with sequencing depth.

## What this repository evaluates

The current analysis focuses on:

- RNA and ATAC quality control performed separately before integration
- RNA representation using PCA
- ATAC representation using TF-IDF and LSI
- WNN integration of the paired modalities
- comparison of RNA-only, ATAC-only, and WNN cell representations
- cell-specific RNA/ATAC modality weights
- preservation of known PBMC cell populations after integration

The main question is whether the joint representation adds useful information while preserving known biological structure.

## Downstream analyses

The repository also contains exploratory analyses extending beyond WNN:

- RNA-derived transcription-factor program scores
- ATAC motif accessibility using chromVAR
- cross-modal comparison of RNA regulatory programs and ATAC motif activity
- exploratory monocyte pseudotime analysis

These analyses are interpreted as associations between transcriptional and chromatin-accessibility signals; not evidence of causal regulation.

## Workflow

| Step | Analysis | Main script |
|---|---|---|
| 1 | Import RNA reference outputs | `scripts/01_import_rna_baseline.py` |
| 2 | Prepare multiome data | `scripts/02_download_or_prepare_multiome.R` |
| 3 | RNA + ATAC preprocessing and QC | `scripts/03_preprocess_multiome.R` |
| 4 | WNN integration | `scripts/04_integrate_rna_atac_wnn.R` |
| 5 | RNA vs WNN comparison | `scripts/05_compare_rna_vs_wnn.py` |
| 6 | RNA TF program analysis | `scripts/06_tf_activity.py` |
| 7 | ATAC motif analysis | `scripts/08_motif_enrichment_atac.R` |
| 8 | Cross-modal TF comparison | `scripts/09_compare_rna_atac_tf.py` |
| 9 | Monocyte pseudotime | `scripts/10_pseudotime_subset.R` |

## Quick start

```bash
mamba env create -f environment.yml
conda activate sc-rna-atac

Rscript scripts/02_download_or_prepare_multiome.R --download
Rscript scripts/03_preprocess_multiome.R
Rscript scripts/04_integrate_rna_atac_wnn.R
```

## Data

The project uses:

- 10x Genomics PBMC Multiome paired RNA + ATAC data
- RNA-only reference outputs from the companion `sc-cell-state-benchmark` repository

Large raw and processed data files are not committed.

## Interpretation

Integration quality is evaluated by checking whether known biological structure is preserved across RNA, ATAC, and WNN representations.

Important limitations include the sparsity of both modalities, sensitivity of ATAC analysis to sequencing depth and QC, and regulatory potential of open chromatin.

## Repository structure

```text
sc-rna-atac-regulon-benchmark/
├── data/
├── docs/
├── results/
├── scripts/
├── src/
├── tests/
├── environment.yml
├── pyproject.toml
└── renv.lock
```

## License

MIT
