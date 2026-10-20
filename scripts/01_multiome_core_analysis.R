#!/usr/bin/env Rscript

# core paired RNA + ATAC seq multiome workflow
# Steps:
#   1. Download/validate the 10x PBMC Multiome inputs
#   2. Load the paired feature matrix
#   3. Qulaity control of RNA-Seq independently + PCA on RNA
#   4. Quality control of ATAC-Seq independently + TFIDF/LSI on ATAC
#   5. Weighted Nearest Neighbour integration
#   6. Technical Validation
#   7. Biological validation/annotation

suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(patchwork)
    library(Signac)
    library(EnsDb.Hsapiens.v86)
    library(GenomeInfoDb)
    library(future)
    library(dplyr)
    library(SingleR)
    library(celldex)
})

set.seed(42)

# Paths and setting

raw_dir <- "data/raw"
processed_dir <- "data/processed"
results_dir <- "results"

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

base_url <- paste0(
  "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/",
  "pbmc_granulocyte_sorted_10k/"
)

files_to_download <- c(
  "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5",
  "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
  "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi"
)

# RNA QC thresholds
min_features_rna <- 200
max_features_rna <- 7500
max_percent_mt <- 20

# download data first time
for (filename in files_to_download) {
  destination <- file.path(raw_dir, filename)

  if (!file.exists(destination)) {
    message("Downloading: ", filename)
    download.file(
      url      = paste0(base_url, filename),
      destfile = destination,
      mode     = "wb"
    )
  } else {
    message("Already present: ", filename)
  }
}

h5_file <- file.path(
  raw_dir,
  "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
)

fragment_file <- file.path(
  raw_dir,
  "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
)

stopifnot(
  file.exists(h5_file),
  file.exists(fragment_file),
  file.exists(paste0(fragment_file, ".tbi"))
)

# Read input files

counts <- Read10X_h5(h5_file)

rna_counts  <- counts[["Gene Expression"]]
atac_counts <- counts[["Peaks"]]

message("RNA matrix:  ", nrow(rna_counts), " genes x ",
        ncol(rna_counts), " cells")
message("ATAC matrix: ", nrow(atac_counts), " peaks x ",
        ncol(atac_counts), " cells")

shared_barcodes <- intersect(colnames(rna_counts), colnames(atac_counts))

message(
  length(shared_barcodes), " shared RNA+ATAC barcodes out of ",
  ncol(rna_counts), " RNA barcodes"
)

if (length(shared_barcodes) == 0) {
  stop("No shared RNA/ATAC barcodes were found.")
}

rna_counts  <- rna_counts[, shared_barcodes]
atac_counts <- atac_counts[, shared_barcodes]

rna_obj <- CreateSeuratObject(
  counts       = rna_counts,
  assay        = "RNA",
  project      = "PBMC_multiome",
  min.cells    = 10,
  min.features = 0
)

# RNA-seq QC:

rna_obj$percent_mt <- PercentageFeatureSet(rna_obj, pattern = "^MT-")

rna_obj$fail_rna_features <- (rna_obj$nFeature_RNA < min_features_rna | rna_obj$nFeature_RNA > max_features_rna)
rna_obj$fail_percent_mt <- rna_obj$percent_mt >= max_percent_mt
rna_obj$pass_rna_qc <- !(rna_obj$fail_rna_features | rna_obj$fail_percent_mt)

message(
  "RNA QC: ", sum(rna_obj$pass_rna_qc), " / ", ncol(rna_obj),
  " cells pass (", round(100 * sum(rna_obj$pass_rna_qc) / ncol(rna_obj), 1), "%)"
)

qc_plot_df <- rbind(
  data.frame(
    stage = "Before QC",
    nFeature_RNA = rna_obj$nFeature_RNA,
    nCount_RNA = rna_obj$nCount_RNA,
    percent_mt = rna_obj$percent_mt
  ),
  data.frame(
    stage = "After RNA QC",
    nFeature_RNA = rna_obj$nFeature_RNA[rna_obj$pass_rna_qc],
    nCount_RNA = rna_obj$nCount_RNA[rna_obj$pass_rna_qc],
    percent_mt = rna_obj$percent_mt[rna_obj$pass_rna_qc]
  )
)

qc_plot_df$stage <- factor(
  qc_plot_df$stage,
  levels = c("Before QC", "After RNA QC")
)

p_features <- ggplot(qc_plot_df, aes(stage, nFeature_RNA)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  geom_hline(
    yintercept = c(min_features_rna, max_features_rna),
    linetype = "dashed"
  ) +
  labs(
    title = "Genes detected",
    x = NULL,
    y = "nFeature_RNA"
  ) +
  theme_classic()

p_counts <- ggplot(qc_plot_df, aes(stage, nCount_RNA)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  scale_y_log10() +
  labs(
    title = "RNA UMI depth",
    x = NULL,
    y = "nCount_RNA"
  ) +
  theme_classic()

p_mt <- ggplot(qc_plot_df, aes(stage, percent_mt)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  geom_hline(
    yintercept = max_percent_mt,
    linetype = "dashed"
  ) +
  labs(
    title = "Mitochondrial fraction",
    x = NULL,
    y = "percent.mt"
  ) +
  theme_classic()

rna_qc_plot <- (p_features | p_counts | p_mt) +
  plot_annotation(
    title = "RNA quality control",
    subtitle = paste0(
      ncol(rna_obj), " cells before QC; ",
      sum(rna_obj$pass_rna_qc), " cells retained"
    )
  )

ggsave(
  file.path(results_dir, "rna_qc_before_after.png"),
  rna_qc_plot,
  width = 11,
  height = 4.5,
  dpi = 300
)

# ATAC-seq QC:

# ATAC QC thresholds
min_count_atac <- 1000
max_count_atac <- 100000
min_tss <- 2
max_nucleosome <- 2

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"

atac_assay <- CreateChromatinAssay(
    counts = atac_counts[, colnames(rna_obj)],
    sep = c(":", "-"),
    fragments = fragment_file,
    annotation = annotations,
    min.cells = 10,
    min.features = 0
)
rna_obj[["ATAC"]] <- atac_assay

# Calculate ATAC QC metrics
qc_checkpoint <- file.path(
  processed_dir,
  "multiome_qc_metrics.rds"
)

if (file.exists(qc_checkpoint)) {

  message("Loading saved ATAC QC object: ", qc_checkpoint)

  rna_obj <- readRDS(qc_checkpoint)

} else {

  DefaultAssay(rna_obj) <- "ATAC"

  rna_obj <- NucleosomeSignal(rna_obj)

  rna_obj <- TSSEnrichment(
    rna_obj,
    fast = FALSE
  )

  saveRDS(
    rna_obj,
    qc_checkpoint
  )

  message("Saved ATAC QC checkpoint: ", qc_checkpoint)
}

rna_obj$fail_atac_counts <- (
    rna_obj$nCount_ATAC < min_count_atac |
    rna_obj$nCount_ATAC > max_count_atac
)

rna_obj$fail_tss <- rna_obj$TSS.enrichment < min_tss
rna_obj$fail_nucleosome <- (rna_obj$nucleosome_signal >= max_nucleosome)

rna_obj$pass_atac_qc <- !(
  rna_obj$fail_atac_counts |
  rna_obj$fail_tss |
  rna_obj$fail_nucleosome
)

message(
  "ATAC QC: ", sum(rna_obj$pass_atac_qc), " / ", ncol(rna_obj),
  " cells pass (",
  round(100 * sum(rna_obj$pass_atac_qc) / ncol(rna_obj), 1),
  "%)"
)

# plot ATAC QC before and after filtering

atac_qc_plot_df <- rbind(
  data.frame(
    stage = "Before QC",
    nCount_ATAC = rna_obj$nCount_ATAC,
    nFeature_ATAC = rna_obj$nFeature_ATAC,
    TSS_enrichment = rna_obj$TSS.enrichment,
    nucleosome_signal = rna_obj$nucleosome_signal
  ),
  data.frame(
    stage = "After ATAC QC",
    nCount_ATAC = rna_obj$nCount_ATAC[rna_obj$pass_atac_qc],
    nFeature_ATAC = rna_obj$nFeature_ATAC[rna_obj$pass_atac_qc],
    TSS_enrichment = rna_obj$TSS.enrichment[rna_obj$pass_atac_qc],
    nucleosome_signal = rna_obj$nucleosome_signal[rna_obj$pass_atac_qc]
  )
)

atac_qc_plot_df$stage <- factor(
  atac_qc_plot_df$stage,
  levels = c("Before QC", "After ATAC QC")
)

p_atac_count <- ggplot(
  atac_qc_plot_df,
  aes(stage, nCount_ATAC)
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  geom_hline(
    yintercept = c(min_count_atac, max_count_atac),
    linetype = "dashed"
  ) +
  scale_y_log10() +
  labs(
    title = "ATAC fragment depth",
    x = NULL,
    y = "nCount_ATAC"
  ) +
  theme_classic()

p_atac_features <- ggplot(
  atac_qc_plot_df,
  aes(stage, nFeature_ATAC)
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  scale_y_log10() +
  labs(
    title = "Accessible peaks detected",
    x = NULL,
    y = "nFeature_ATAC"
  ) +
  theme_classic()

p_tss <- ggplot(
  atac_qc_plot_df,
  aes(stage, TSS_enrichment)
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  geom_hline(
    yintercept = min_tss,
    linetype = "dashed"
  ) +
  labs(
    title = "TSS enrichment",
    x = NULL,
    y = "TSS enrichment"
  ) +
  theme_classic()

p_nucleosome <- ggplot(
  atac_qc_plot_df,
  aes(stage, nucleosome_signal)
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  geom_hline(
    yintercept = max_nucleosome,
    linetype = "dashed"
  ) +
  labs(
    title = "Nucleosome signal",
    x = NULL,
    y = "nucleosome signal"
  ) +
  theme_classic()

atac_qc_plot <- (
  p_atac_count |
  p_atac_features |
  p_tss |
  p_nucleosome
) +
  plot_annotation(
    title = "ATAC quality control",
    subtitle = paste0(
      ncol(rna_obj), " cells before QC; ",
      sum(rna_obj$pass_atac_qc), " cells retained"
    )
  )

ggsave(
  file.path(results_dir, "atac_qc_before_after.png"),
  atac_qc_plot,
  width = 13,
  height = 4.5,
  dpi = 300
)

# Combined filter:
rna_obj$pass_qc <- (
  rna_obj$pass_rna_qc &
  rna_obj$pass_atac_qc
)

message(
  "Combined QC: ",
  sum(rna_obj$pass_qc), " / ", ncol(rna_obj),
  " cells retained"
)

combined_qc_obj <- subset(
  rna_obj,
  subset = pass_qc
)

#   5. Weighted Nearest Neighbour integration

# RNA: SCTransform + PCA
DefaultAssay(combined_qc_obj) <- "RNA"

combined_qc_obj <- NormalizeData(
  combined_qc_obj,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = FALSE
)

combined_qc_obj <- FindVariableFeatures(
  combined_qc_obj,
  nfeatures = 3000,
  verbose = FALSE
)

combined_qc_obj <- ScaleData(
  combined_qc_obj,
  vars.to.regress = "percent_mt",
  verbose = FALSE
)

combined_qc_obj <- RunPCA(
  combined_qc_obj,
  npcs = 50,
  reduction.name = "pca",
  verbose = FALSE
)

# combined_qc_obj <- RunPCA(
#   combined_qc_obj,
#   assay = "SCT",
#   npcs = 50,
#   reduction.name = "pca",
#   verbose = FALSE
# )

ElbowPlot(combined_qc_obj, reduction = "pca", ndims= 50)
ggsave(
  file.path(results_dir, "rna_pca_elbow.png"),
  width = 6,
  height = 4,
  dpi = 300
)

# ATAC: TF-IDF + LSI

DefaultAssay(combined_qc_obj) <- "ATAC"

combined_qc_obj <- RunTFIDF(combined_qc_obj)

combined_qc_obj <- FindTopFeatures(
  combined_qc_obj,
  min.cutoff = "q0"
)

combined_qc_obj <- RunSVD(
  combined_qc_obj,
  n = 50,
  reduction.name = "lsi"
)

lsi1 <- Embeddings(combined_qc_obj, "lsi")[, 1]

lsi_depth_df <- data.frame(
  nCount_ATAC = combined_qc_obj$nCount_ATAC,
  LSI1 = lsi1
)

p_lsi_depth <- ggplot(
  lsi_depth_df,
  aes(x = nCount_ATAC, y = LSI1)
) +
  geom_point(size = 0.4, alpha = 0.5) +
  scale_x_log10() +
  labs(
    title = "LSI1 vs ATAC depth",
    x = "nCount_ATAC",
    y = "LSI1"
  ) +
  theme_classic()

ggsave(
  file.path(results_dir, "atac_lsi_depth.png"),
  p_lsi_depth,
  width = 6,
  height = 4,
  dpi = 300
)

lsi_depth_cor <- cor(
  combined_qc_obj$nCount_ATAC,
  lsi1,
  method = "spearman"
)

message(
  "Spearman correlation LSI1 vs ATAC depth: ",
  round(lsi_depth_cor, 3)
)

# save LSI variance plot
lsi_stdev <- Stdev(combined_qc_obj, reduction = "lsi")

lsi_var <- lsi_stdev^2 / sum(lsi_stdev^2) * 100

lsi_df <- data.frame(
  dimension = seq_along(lsi_var),
  variance = lsi_var
)

p_lsi <- ggplot(lsi_df, aes(dimension, variance)) +
  geom_col() +
  labs(
    title = "ATAC LSI variance",
    x = "LSI dimension",
    y = "Variance explained (%)"
  ) +
  theme_classic()

ggsave(
  file.path(results_dir, "atac_lsi_variance.png"),
  p_lsi,
  width = 6,
  height = 4,
  dpi = 300
)

message("RNA PCA and ATAC LSI complete.")

# Weighted Nearest Neighbor integration

rna_dims <- 1:40
atac_dims <- 2:40   # exclude LSI1 because it often tracks sequencing depth

combined_qc_obj <- FindMultiModalNeighbors(
  combined_qc_obj,
  reduction.list = list("pca", "lsi"),
  dims.list = list(rna_dims, atac_dims)
)

# WNN UMAP

combined_qc_obj <- RunUMAP(
  combined_qc_obj,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_",
  seed.use = 42
)

# WNN clustering

combined_qc_obj <- FindClusters(combined_qc_obj, graph.name = "wsnn", resolution = 0.3)

combined_qc_obj$wnn_cluster <- combined_qc_obj$seurat_clusters

message(
  "WNN complete: ",
  length(unique(combined_qc_obj$wnn_cluster)),
  " clusters"
)

# RNA-only UMAP

combined_qc_obj <- RunUMAP(
  combined_qc_obj,
  reduction = "pca",
  dims = rna_dims,
  reduction.name = "rna.umap",
  reduction.key = "rnaUMAP_",
  seed.use = 42
)

# ATAC-only UMAP

combined_qc_obj <- RunUMAP(
  combined_qc_obj,
  reduction = "lsi",
  dims = atac_dims,
  reduction.name = "atac.umap",
  reduction.key = "atacUMAP_",
  seed.use = 42
)

# 6. Validation: modality weights

p_weight_umap <- FeaturePlot(
  combined_qc_obj,
  features = "RNA.weight",
  reduction = "wnn.umap"
) +
  labs(title = "RNA modality weight across WNN embedding")

ggsave(
  file.path(results_dir, "wnn_rna_modality_weight.png"),
  p_weight_umap,
  width = 6,
  height = 5,
  dpi = 300
)

weight_df <- data.frame(
  cluster = combined_qc_obj$wnn_cluster,
  RNA.weight = combined_qc_obj$RNA.weight,
  ATAC.weight = combined_qc_obj$ATAC.weight
)

p_weight_cluster <- ggplot(
  weight_df,
  aes(x = factor(cluster), y = RNA.weight)
) +
  geom_boxplot(outlier.size = 0.4) +
  labs(
    title = "RNA modality weight by WNN cluster",
    x = "WNN cluster",
    y = "RNA modality weight"
  ) +
  theme_classic()

ggsave(
  file.path(results_dir, "wnn_rna_weight_by_cluster.png"),
  p_weight_cluster,
  width = 8,
  height = 4.5,
  dpi = 300
)

write.csv(
  weight_df,
  file.path(results_dir, "wnn_modality_weights.csv"),
  row.names = FALSE
)
# So far:
# 1. QC is strong: 95.5% of paired cells survive combined RNA+ATAC QC; 
#  WNN is not being built on a heavily degraded subset.
# 2. ATAC representation looks sensible: LSI1 strongly tracks depth
#  (Spearman ≈ 0.924) hence excluded it from WNN.
# 3. Both modalities contribute: weights are not collapsing to ~0 
#  or ~1 globally; RNA and ATAC both meaningfully influence neighborhoods.
# 4. Weights vary by cell/cluster
# 5. The graph is well structured: WNN clustering produced clear 
#  communities with high modularity ~0.945.
# 6. No obvious pathological output: no sign of one modality completely 
#  dominating, fragmented graphs, bizarre loss of cells, or 
#  unstable-looking integration.
# Hence we validate biology first and fine tune WNN params if needed.

# 7. Cell-type annotation
# Compare RNA-only, ATAC-only and WNN embeddings

p_rna_umap <- DimPlot(
  combined_qc_obj,
  reduction = "rna.umap",
  group.by = "wnn_cluster",
  label = TRUE,
  repel = TRUE
) +
  labs(title = "RNA-only UMAP") +
  NoLegend()

p_atac_umap <- DimPlot(
  combined_qc_obj,
  reduction = "atac.umap",
  group.by = "wnn_cluster",
  label = TRUE,
  repel = TRUE
) +
  labs(title = "ATAC-only UMAP") +
  NoLegend()

p_wnn_umap <- DimPlot(
  combined_qc_obj,
  reduction = "wnn.umap",
  group.by = "wnn_cluster",
  label = TRUE,
  repel = TRUE
) +
  labs(title = "WNN UMAP") +
  NoLegend()

p_embedding_compare <- (
  p_rna_umap |
  p_atac_umap |
  p_wnn_umap
)

ggsave(
  file.path(results_dir, "rna_atac_wnn_umap_comparison.png"),
  p_embedding_compare,
  width = 15,
  height = 5,
  dpi = 300
)

# Find RNA markers for WNN clusters

DefaultAssay(combined_qc_obj) <- "RNA"

Idents(combined_qc_obj) <- combined_qc_obj$wnn_cluster

rna_markers <- FindAllMarkers(
  combined_qc_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  rna_markers,
  file.path(results_dir, "wnn_cluster_rna_markers.csv"),
  row.names = FALSE
)

top_rna_markers <- rna_markers |>
  dplyr::group_by(cluster) |>
  dplyr::slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()

write.csv(
  top_rna_markers,
  file.path(results_dir, "wnn_cluster_top10_rna_markers.csv"),
  row.names = FALSE
)

# Canonical RNA marker panel for broad PBMC / granulocyte annotation

marker_genes <- c(
  # T cells
  "CD3D", "CD3E", "IL7R", "CD8A",

  # Regulatory T cells
  "FOXP3", "IL2RA", "CTLA4",

  # NK cells
  "NKG7", "GNLY", "KLRF1", "NCR1",

  # B cells
  "MS4A1", "CD79A", "CD37",

  # Plasma cells / plasmablasts
  "SDC1", "JCHAIN", "MZB1",

  # Monocytes / myeloid
  "LYZ", "LST1", "CTSS", "FCN1", "VCAN",
  "FCGR3A", "S100A8", "S100A9",

  # Dendritic cells
  "FCER1A", "CD1C", "CLEC10A",

  # Plasmacytoid dendritic cells
  "CLEC4C", "GZMB",

  # Granulocytes / neutrophils
  "CSF3R", "FCGR3B", "TREM1",

  # Platelets
  "PPBP", "PF4"
)

marker_genes <- marker_genes[
  marker_genes %in% rownames(combined_qc_obj)
]

p_marker_dot <- DotPlot(
  combined_qc_obj,
  features = marker_genes,
  group.by = "wnn_cluster"
) +
  RotatedAxis() +
  labs(
    title = "RNA markers across WNN clusters",
    x = "Marker gene",
    y = "WNN cluster"
  )

ggsave(
  file.path(results_dir, "wnn_cluster_rna_marker_dotplot.png"),
  p_marker_dot,
  width = 14,
  height = 7,
  dpi = 300
)

# Check whether WNN structure is associated with technical variables

qc_features <- c(
  "nCount_RNA",
  "nFeature_RNA",
  "percent_mt",
  "nCount_ATAC",
  "nFeature_ATAC",
  "TSS.enrichment",
  "nucleosome_signal"
)

p_qc_wnn <- FeaturePlot(
  combined_qc_obj,
  features = qc_features,
  reduction = "wnn.umap",
  ncol = 3
)

ggsave(
  file.path(results_dir, "wnn_umap_qc_variables.png"),
  p_qc_wnn,
  width = 12,
  height = 9,
  dpi = 300
)

# Cluster sizes

cluster_sizes <- as.data.frame(
  table(combined_qc_obj$wnn_cluster)
)

colnames(cluster_sizes) <- c(
  "cluster",
  "n_cells"
)

cluster_sizes$cluster <- as.character(cluster_sizes$cluster)

write.csv(
  cluster_sizes,
  file.path(results_dir, "wnn_cluster_sizes.csv"),
  row.names = FALSE
)

# SingleR broad annotation

DefaultAssay(combined_qc_obj) <- "RNA"

ref <- celldex::HumanPrimaryCellAtlasData()

test_expr <- GetAssayData(
  combined_qc_obj,
  assay = "RNA",
  layer = "data"
)

singleR_pred <- SingleR(
  test = test_expr,
  ref = ref,
  labels = ref$label.main,
  clusters = combined_qc_obj$wnn_cluster
)

singleR_table <- data.frame(
  cluster = rownames(singleR_pred),
  singleR_label = singleR_pred$labels,
  singleR_pruned = singleR_pred$pruned.labels,
  stringsAsFactors = FALSE
)
top_marker_summary <- top_rna_markers |>
  dplyr::group_by(cluster) |>
  dplyr::summarise(
    top_markers = paste(gene, collapse = ", "),
    .groups = "drop"
  )

# Ensure matching cluster type

top_marker_summary$cluster <- as.character(
  top_marker_summary$cluster
)

singleR_table$cluster <- as.character(
  singleR_table$cluster
)

# Build annotation table

annotation_table <- cluster_sizes |>
  dplyr::left_join(
    top_marker_summary,
    by = "cluster"
  ) |>
  dplyr::left_join(
    singleR_table,
    by = "cluster"
  ) |>
  dplyr::mutate(
    proposed_identity = NA_character_,
    confidence = NA_character_,
    notes = NA_character_
  ) |>
  dplyr::arrange(
    as.numeric(cluster)
  )

write.csv(
  annotation_table,
  file.path(
    results_dir,
    "wnn_cluster_annotation_table.csv"
  ),
  row.names = FALSE
)

print(annotation_table)

# SingleR diagnostics

singleR_diagnostics <- data.frame(
  cluster = rownames(singleR_pred),
  label = singleR_pred$labels,
  pruned_label = singleR_pred$pruned.labels,
  delta_next = singleR_pred$delta.next,
  stringsAsFactors = FALSE
)

write.csv(
  singleR_diagnostics,
  file.path(
    results_dir,
    "singleR_cluster_diagnostics.csv"
  ),
  row.names = FALSE
)