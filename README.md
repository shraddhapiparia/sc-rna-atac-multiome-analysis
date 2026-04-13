# sc-rna-atac-regulon-benchmark

A follow-up single-cell multiome repository extending the RNA-only
`sc-cell-state-benchmark` into paired RNA+ATAC integration, regulatory inference,
and optional dynamic analysis.

## What this repo demonstrates

- Import of RNA-only baseline outputs from the benchmark repo
- Paired RNA+ATAC integration using a shared-cell multiome dataset
- Joint clustering and embedding with multimodal neighborhood analysis
- TF/regulon inference and regulatory interpretation
- Optional pseudotime analysis of dynamic regulatory responses

## Relationship to the previous repo

This repository reuses RNA-only baseline outputs from `sc-cell-state-benchmark`
as the reference point for multimodal and regulatory comparisons.

## Repository structure

sc-rna-atac-regulon-benchmark/
├── README.md
├── environment.yml
├── requirements.txt
├── pyproject.toml
├── .gitignore
├── data/
│   ├── raw/
│   ├── processed/
│   └── external/
├── docs/
├── figures/
├── results/
│   ├── baseline/
│   ├── integration/
│   ├── regulons/
│   ├── dynamics/
│   └── project_summary.md
├── scripts/
├── src/sc_rna_atac_regulon_benchmark/
└── tests/

## Repo progression

1. RNA-only baseline import
2. RNA+ATAC integration
3. TF / regulon inference
4. Optional dynamic analysis

## Phase 1 — RNA baseline import

This phase reuses exported Kang et al. RNA baseline outputs from the upstream
`sc-cell-state-benchmark` repository. No re-processing is performed; the script
standardises column names and writes clean outputs for downstream comparison.

**Required inputs** (`data/external/`):

| File | Contents |
|------|----------|
| `kang_interferon_scores.csv` | Per-cell IFN scores (`cell_id`, `condition`, raw score columns) |
| `kang_program_scores.csv` | Per-cell program scores with `cell_type` annotations |

**Generated outputs** (`results/baseline/`):

| File | Contents |
|------|----------|
| `rna_cell_metadata.csv` | `cell_id`, `condition`, `cell_type` |
| `rna_ifn_scores.csv` | Per-cell IFN scores with standardised column names |
| `rna_program_scores.csv` | Per-cell program scores with standardised column names |
| `rna_baseline_summary.md` | Cell counts, cell-type counts, condition breakdown |

**Run:**

```bash
python scripts/01_import_rna_baseline.py
```

## Phase 2 — RNA+ATAC WNN integration

Paired RNA+ATAC integration using the 10x Genomics PBMC Multiome 10k dataset.
WNN clustering produces joint multimodal cell states; outputs feed into Phase 3
regulon inference and Phase 5 cross-dataset comparison.

### Broad cell-type recovery after WNN integration

The Phase 1 RNA-only baseline identified 8 broad immune cell classes. An exact
string comparison between RNA baseline labels and WNN cluster labels showed zero
overlap, because WNN produced finer-grained multimodal labels rather than
reusing the same label vocabulary.

After harmonising WNN labels to broad immune categories using
`results/comparison/wnn_to_broad_labels.tsv`, 7 of the 8 RNA-defined populations
were recovered:

- B cells
- CD14+ Monocytes
- CD4 T cells
- CD8 T cells
- Dendritic cells
- FCGR3A+ Monocytes
- NK cells

The only RNA-defined label not recovered was **Megakaryocytes**, likely absent
or extremely rare in the multiome dataset.

Example refined WNN labels not present in the RNA baseline:

- FOXP3+ Treg
- FCGR3A+ Monocyte
- pDC-like
- NK subset
- Activated CD4 T state

WNN preserved the major biological structure from the RNA-only analysis while
refining it into more specific multimodal cell states.

## Phase 3 — TF activity proxy analysis

Rather than running pySCENIC or ArchR, Phase 3 uses the five RNA program score
columns from `results/baseline/rna_program_scores.csv` as lightweight proxies for
transcription-factor activity.  Each program maps to a canonical TF axis:

| Program | TF axis |
|---------|---------|
| `interferon_alpha_beta_response` | IRF3 / IRF7 / STAT1/2 |
| `inflammatory_response` | NF-κB / AP-1 |
| `antigen_presentation_mhc` | CIITA / IRF1 |
| `cytotoxic_nk_like` | EOMES / TBX21 |
| `monocyte_activation` | SPI1 (PU.1) / CEBPB |

`06_tf_activity.py` computes mean ± std per RNA cell type and writes flat summary
tables used by the plotting script.  `07_plot_tf_activity.py` produces z-scored
heatmaps and per-program bar charts.

### Biological findings

- **Antigen presentation (CIITA / IRF1):** Dendritic cells and B cells show the
  strongest activity, consistent with their professional APC roles.
- **Cytotoxic program (EOMES / TBX21):** NK cells and CD8 T cells rank highest,
  as expected for cytolytic populations.
- **Inflammatory & monocyte-activation programs (NF-κB / AP-1, PU.1 / CEBPB):**
  CD14+ and FCGR3A+ monocytes dominate both axes.
- **Interferon α/β response (IRF3 / IRF7 / STAT1/2):** FCGR3A+ monocytes show
  the strongest signal among all populations.
- **CD4 T cells** are consistently low across all innate-associated programs.

The TF-activity proxy layer validates that the WNN-aligned populations preserve
expected immune regulatory programs.

### Run

```bash
python scripts/06_tf_activity.py
python scripts/07_plot_tf_activity.py
```

### Generated outputs

| File | Location |
|------|----------|
| `tf_activity_rna_by_celltype.csv` | `results/regulons/` |
| `tf_activity_rna_by_broad.csv` | `results/regulons/` |
| `tf_activity_wnn_celltype_map.csv` | `results/regulons/` |
| `tf_activity_summary.md` | `results/regulons/` |
| `tf_heatmap_rna_celltype.png` | `results/figures/tf_activity/` |
| `tf_heatmap_rna_broad.png` | `results/figures/tf_activity/` |
| `tf_barplot_interferon_alpha_beta_response.png` | `results/figures/tf_activity/` |
| `tf_barplot_inflammatory_response.png` | `results/figures/tf_activity/` |
| `tf_barplot_antigen_presentation_mhc.png` | `results/figures/tf_activity/` |
| `tf_barplot_cytotoxic_nk_like.png` | `results/figures/tf_activity/` |
| `tf_barplot_monocyte_activation.png` | `results/figures/tf_activity/` |
