# sc-rna-atac-multiome-analysis

Single-cell RNA + ATAC analysis of paired 10x Genomics PBMC Multiome data using Seurat and Signac.

The project focuses on understanding how transcriptomic and chromatin-accessibility information can be analyzed independently and combined using Weighted Nearest Neighbor (WNN) integration.

## Biological motivation

scRNA-seq measures the transcriptional state of a cell, while scATAC-seq measures regions of accessible chromatin associated with regulatory potential.

These modalities provide complementary information. Chromatin accessibility can indicate which regulatory regions are available for use, while RNA expression reflects the transcriptional state that is actually observed.

Joint analysis therefore provides a way to ask whether combining the two modalities improves representation of cell identity and cell state.

## Analysis workflow

RNA and ATAC are first processed independently before multimodal integration.

```text
RNA counts                         ATAC peak counts
    │                                    │
LogNormalize                            TF-IDF
    │                                    │
Variable features                       LSI
    │                                    │
Scale + PCA                             │
    │                                   │
    └────────────── WNN ─────────────────┘
                       │
                 WNN graph
                       │
              clustering + UMAP
                       │
             biological validation
                       │
               cell annotation
```

WNN does not simply concatenate RNA and ATAC features. Instead, it learns a cell-specific contribution from each modality and combines modality-specific neighborhood information.

## Current workflow

The main analysis is implemented in:

```text
scripts/01_multiome_core_analysis.R
```

The analysis workflow includes:

1. RNA + ATAC data loading and QC
2. RNA PCA and ATAC LSI dimensionality reduction
3. WNN integration, clustering, and RNA/ATAC/WNN visualization
4. integration validation using modality weights and technical QC metrics
5. RNA marker detection and SingleR-supported cell-type annotation

Downstream analyses are organized into:

- `02_gene_activity_analysis.R`
- `03_chromvar_motif_activity.R`
- `04_rna_atac_regulatory_comparison.R`
- `05_peak_to_gene_analysis.R`

## Validation and design decisions

Following steps were taken to evaluate whether multimodal integration captures biological structure:

- RNA-only, ATAC-only, and WNN representations are compared.
- Cell-specific RNA and ATAC modality weights are inspected.
- WNN clusters are evaluated against RNA depth, ATAC depth, mitochondrial fraction, TSS enrichment, and nucleosome signal.
- Broad cell identities are assigned using SingleR and validated using canonical RNA markers and cluster-specific differential expression.
- Cluster sizes and rare populations are inspected before accepting final annotations.
- Clustering resolution and dimensionality are treated as analysis choices and evaluated through sensitivity analyses.

### ATAC depth and LSI1

The first ATAC LSI dimension was strongly correlated with ATAC sequencing depth:

```text
Spearman correlation ≈ 0.924
```

Because LSI1 primarily captured library-depth variation, it was excluded from WNN integration.

The current multimodal representation uses:

```text
RNA:  PCA dimensions 1–40
ATAC: LSI dimensions 2–40
```

The upper dimensionality and clustering resolution are being evaluated through sensitivity analysis.

## Cell-type annotation

Cluster annotation uses two complementary sources of evidence.

### Reference-based annotation

SingleR is used with the Human Primary Cell Atlas reference to assign broad cell lineages at the cluster level.

### Marker-based annotation

SingleR labels are then validated and refined using:

- cluster-specific RNA markers
- canonical immune-cell markers
- marker gene expressions
- cluster size
- technical QC characteristics



## Downstream analyses

After finalizing WNN integration and cell-type annotation, the focus is on regulatory interpretation.

Planned steps include:

- clustering-resolution sensitivity analysis
- RNA / ATAC dimensionality sensitivity analysis
- final cell-type annotation
- gene-activity analysis from chromatin accessibility
- motif enrichment and chromVAR activity
- comparison of RNA expression and ATAC regulatory activity
- identification of candidate transcription-factor programs
- exploratory peak-to-gene relationships


## Data

The project uses the 10x Genomics PBMC Granulocyte Sorted 10k Multiome dataset containing paired RNA and ATAC measurements from the same cells.

The workflow downloads:

```text
pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi
```


## Quick start

Create and activate the environment:

```bash
mamba env create -f environment.yml
conda activate sc-rna-atac
```

Run the main analysis:

```bash
Rscript scripts/01_multiome_core_analysis.R
```

The script downloads the required 10x files if they are not already present.

## Repository structure

```text
sc-rna-atac-multiome-analysis/
├── data/
│   ├── raw/
│   └── processed/
├── docs/
├── results/
├── scripts/
│   ├── 01_multiome_core_analysis.R
│   ├── 02_gene_activity_analysis.R
│   ├── 03_chromvar_motif_activity.R
│   ├── 04_rna_atac_regulatory_comparison.R
│   └── 05_peak_to_gene_analysis.R
├── environment.yml
├── renv.lock
└── README.md
```


## License

MIT
