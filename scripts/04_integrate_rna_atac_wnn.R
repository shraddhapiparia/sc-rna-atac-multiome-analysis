#!/usr/bin/env Rscript
#
# 04_integrate_rna_atac_wnn.R
#
# WNN integration, multi-modal and per-modality UMAP, graph-based clustering,
# and cluster marker export for the preprocessed 10x PBMC Multiome object.
#
# Reads:
#   data/processed/pbmc_multiome_preprocessed.rds   (from script 03)
#   [optional] --cluster-annotations <file.tsv>     (cluster -> cell_type mapping)
#
# Writes:
#   data/processed/multiome_wnn_object.rds
#   results/integration/multiome_cell_metadata.csv
#   results/integration/multiome_umap_coords.csv
#   results/integration/wnn_cluster_markers.csv
#
# Usage:
#   Rscript scripts/04_integrate_rna_atac_wnn.R [options]
#   Rscript scripts/04_integrate_rna_atac_wnn.R --help

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Signac)
  library(future)
})

# ---------------------------------------------------------------------------
# Resolve repo root from script location (works from any working directory)
# ---------------------------------------------------------------------------

.script_file <- normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)),
  mustWork = FALSE
)
.repo_root <- dirname(dirname(.script_file))

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------

option_list <- list(

  # I/O paths
  make_option("--input-rds",
    default = file.path(.repo_root, "data", "processed",
                        "pbmc_multiome_preprocessed.rds"),
    type = "character",
    help = "Path to preprocessed Seurat .rds from script 03 [default: %default]"),
  make_option("--processed-dir",
    default = file.path(.repo_root, "data", "processed"),
    type = "character",
    help = "Output directory for the WNN .rds object [default: %default]"),
  make_option("--results-dir",
    default = file.path(.repo_root, "results", "integration"),
    type = "character",
    help = "Output directory for CSV/TSV artefacts [default: %default]"),

  # Dimensionality
  make_option("--rna-dims", type = "integer", default = 40L,
    help = "Max PCA dim for RNA in WNN (uses 1:N) [default: %default]"),
  make_option("--atac-dims", type = "integer", default = 40L,
    help = "Max LSI dim for ATAC in WNN (uses 2:N; dim 1 excluded) [default: %default]"),

  # Clustering
  make_option("--resolution", type = "double", default = 0.3,
    help = "Clustering resolution for FindClusters [default: %default]"),
  make_option("--cluster-algorithm", type = "integer", default = 1L,
    help = paste("Clustering algorithm: 1=Louvain, 2=Louvain multilevel,",
                 "3=SLM, 4=Leiden [default: %default]")),

  # Optional cell-type annotation
  make_option("--cluster-annotations",
    default = "",
    type = "character",
    help = paste("Two-column TSV (cluster, cell_type) mapping WNN cluster IDs",
                 "to cell-type labels. If omitted, cell_type_wnn is set to NA.")),

  make_option("--seed", type = "integer", default = 42L,
    help = "Random seed [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

# FindAllMarkers uses future for parallelism; the WNN Seurat object on this
# ~10k-cell dataset can exceed future's default 500 MiB global export limit.
# Cap at 12 GiB (larger than script 03 because the object now carries
# multiple reductions and graphs) and force sequential execution to avoid it.
options(future.globals.maxSize = 12 * 1024^3)
future::plan("sequential")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

stop_clear <- function(...) {
  message("ERROR: ", ...)
  quit(status = 1, save = "no")
}

# ---------------------------------------------------------------------------
# Validate inputs and set up output directories
# ---------------------------------------------------------------------------

if (!file.exists(opt$`input-rds`)) {
  stop_clear("Input .rds not found: ", opt$`input-rds`,
             "\nRun scripts/03_preprocess_multiome.R first.")
}

annot_file <- opt$`cluster-annotations`
use_annot  <- nzchar(annot_file)
if (use_annot && !file.exists(annot_file)) {
  stop_clear("--cluster-annotations file not found: ", annot_file)
}

dir.create(opt$`processed-dir`, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$`results-dir`,   recursive = TRUE, showWarnings = FALSE)

RDS_OUT     <- file.path(opt$`processed-dir`, "multiome_wnn_object.rds")
META_CSV    <- file.path(opt$`results-dir`,   "multiome_cell_metadata.csv")
UMAP_CSV    <- file.path(opt$`results-dir`,   "multiome_umap_coords.csv")
MARKER_CSV  <- file.path(opt$`results-dir`,   "wnn_cluster_markers.csv")

rna_dims  <- seq_len(opt$`rna-dims`)          # 1:N  for RNA PCA
atac_dims <- seq(2, opt$`atac-dims`)          # 2:N  for ATAC LSI (dim 1 excluded)

# ---------------------------------------------------------------------------
# Load preprocessed object
# ---------------------------------------------------------------------------

message("Loading preprocessed object: ", opt$`input-rds`)
obj <- readRDS(opt$`input-rds`)

# Verify required reductions are present and requested dims do not exceed
# what was actually computed in script 03.
for (red in c("pca", "lsi")) {
  if (!red %in% Reductions(obj)) {
    stop_clear("Reduction '", red, "' not found in input object. ",
               "Ensure script 03 completed successfully.")
  }
}

n_pca_available  <- ncol(Embeddings(obj, "pca"))
n_lsi_available  <- ncol(Embeddings(obj, "lsi"))

if (opt$`rna-dims` > n_pca_available) {
  stop_clear("--rna-dims ", opt$`rna-dims`, " exceeds available PCA dims (",
             n_pca_available, "). Re-run script 03 with a higher --n-dims ",
             "or lower --rna-dims.")
}
if (opt$`atac-dims` > n_lsi_available) {
  stop_clear("--atac-dims ", opt$`atac-dims`, " exceeds available LSI dims (",
             n_lsi_available, "). Re-run script 03 with a higher --n-dims ",
             "or lower --atac-dims.")
}

# ---------------------------------------------------------------------------
# WNN: multi-modal neighbour graph
# ---------------------------------------------------------------------------

message("Building WNN graph (RNA dims 1:", opt$`rna-dims`,
        ", ATAC dims 2:", opt$`atac-dims`, ")...")

obj <- FindMultiModalNeighbors(
  obj,
  reduction.list = list("pca", "lsi"),
  dims.list      = list(rna_dims, atac_dims),
  verbose        = FALSE
)

# ---------------------------------------------------------------------------
# UMAP embeddings
# ---------------------------------------------------------------------------

message("Computing joint WNN UMAP...")
obj <- RunUMAP(
  obj,
  nn.name        = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key  = "wnnUMAP_",
  seed.use       = opt$seed,
  verbose        = FALSE
)

message("Computing RNA-only UMAP...")
obj <- RunUMAP(
  obj,
  reduction      = "pca",
  dims           = rna_dims,
  reduction.name = "rna.umap",
  reduction.key  = "rnaUMAP_",
  seed.use       = opt$seed,
  verbose        = FALSE
)

message("Computing ATAC-only UMAP...")
obj <- RunUMAP(
  obj,
  reduction      = "lsi",
  dims           = atac_dims,
  reduction.name = "atac.umap",
  reduction.key  = "atacUMAP_",
  seed.use       = opt$seed,
  verbose        = FALSE
)

# ---------------------------------------------------------------------------
# Clustering on the wsnn graph
# ---------------------------------------------------------------------------

message("Clustering (resolution=", opt$resolution,
        ", algorithm=", opt$`cluster-algorithm`, ")...")

obj <- FindClusters(
  obj,
  graph.name = "wsnn",
  algorithm  = opt$`cluster-algorithm`,
  resolution = opt$resolution,
  verbose    = FALSE
)

# Rename to a stable, script-agnostic column; Seurat writes to seurat_clusters
obj[["wnn_cluster"]] <- as.character(obj$seurat_clusters)
n_clusters <- length(unique(obj$wnn_cluster))
message(n_clusters, " clusters identified")

# ---------------------------------------------------------------------------
# Cell-type annotation (optional)
# ---------------------------------------------------------------------------

if (use_annot) {
  message("Reading cluster annotations: ", annot_file)
  annot_df <- read.table(annot_file, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)

  required_cols <- c("cluster", "cell_type")
  missing_cols  <- setdiff(required_cols, colnames(annot_df))
  if (length(missing_cols) > 0) {
    stop_clear("--cluster-annotations TSV is missing columns: ",
               paste(missing_cols, collapse = ", "),
               "\nExpected: cluster<TAB>cell_type")
  }

  dup_clusters <- annot_df$cluster[duplicated(as.character(annot_df$cluster))]
  if (length(dup_clusters) > 0) {
    stop_clear("--cluster-annotations TSV has duplicated cluster value(s): ",
               paste(sort(unique(dup_clusters)), collapse = ", "),
               "\nEach cluster must appear exactly once.")
  }

  annot_map <- setNames(annot_df$cell_type, as.character(annot_df$cluster))

  # Warn about clusters present in the object but absent from the TSV
  unmapped <- setdiff(obj$wnn_cluster, names(annot_map))
  if (length(unmapped) > 0) {
    message("  Warning: ", length(unmapped), " cluster(s) have no annotation",
            " — cell_type_wnn will be NA for these: ",
            paste(sort(unmapped), collapse = ", "))
  }

  cell_type_vec          <- unname(annot_map[as.character(obj$wnn_cluster)])
  obj[["cell_type_wnn"]] <- cell_type_vec   # NA for unmapped clusters

} else {
  message("No --cluster-annotations provided; cell_type_wnn set to NA")
  obj[["cell_type_wnn"]] <- NA_character_
}

# ---------------------------------------------------------------------------
# Cluster marker genes
# ---------------------------------------------------------------------------

message("Finding cluster markers (FindAllMarkers)...")

# Use SCT assay if available; fall back to RNA for the --no-sct path from
# script 03. PrepSCTFindMarkers() is intentionally skipped — it is only
# needed for multi-sample/batch SCT models, and this is a single-run object.
if ("SCT" %in% Assays(obj)) {
  DefaultAssay(obj) <- "SCT"
} else {
  DefaultAssay(obj) <- "RNA"
}

Idents(obj) <- "wnn_cluster"

markers <- FindAllMarkers(
  obj,
  only.pos       = TRUE,
  min.pct        = 0.25,
  logfc.threshold = 0.25,
  verbose        = FALSE
)

# Standardise output columns
marker_cols <- c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
marker_cols <- intersect(marker_cols, colnames(markers))   # guard against version diffs
markers_out <- markers[, marker_cols]
markers_out <- markers_out[order(markers_out$cluster, -markers_out$avg_log2FC), ]

# ---------------------------------------------------------------------------
# Build output data frames
# ---------------------------------------------------------------------------

# -- UMAP coordinates -------------------------------------------------------

wnn_coords  <- as.data.frame(Embeddings(obj, "wnn.umap"))
rna_coords  <- as.data.frame(Embeddings(obj, "rna.umap"))
atac_coords <- as.data.frame(Embeddings(obj, "atac.umap"))

umap_df <- data.frame(
  cell_id    = rownames(wnn_coords),
  wnn_umap_1 = wnn_coords[[1]],
  wnn_umap_2 = wnn_coords[[2]],
  rna_umap_1 = rna_coords[[1]],
  rna_umap_2 = rna_coords[[2]],
  atac_umap_1 = atac_coords[[1]],
  atac_umap_2 = atac_coords[[2]],
  stringsAsFactors = FALSE
)

# -- Cell metadata ----------------------------------------------------------

md <- obj[[]]

# Standardise TSS column name (Signac uses TSS.enrichment with a dot)
tss_col <- if ("TSS.enrichment" %in% colnames(md)) "TSS.enrichment" else
           if ("TSS_enrichment"  %in% colnames(md)) "TSS_enrichment" else NA

meta_base_cols <- c("nFeature_RNA", "nCount_RNA", "percent_mt",
                    "nCount_ATAC", tss_col, "nucleosome_signal",
                    "wnn_cluster", "cell_type_wnn")
meta_base_cols <- meta_base_cols[!is.na(meta_base_cols)]
meta_base_cols <- intersect(meta_base_cols, colnames(md))

meta_df        <- md[, meta_base_cols, drop = FALSE]
meta_df$cell_id <- rownames(meta_df)

# Rename TSS column to underscore convention for downstream Python phases
if (!is.na(tss_col) && tss_col == "TSS.enrichment") {
  names(meta_df)[names(meta_df) == "TSS.enrichment"] <- "TSS_enrichment"
}

# Merge UMAP coordinates into metadata
meta_df <- merge(meta_df, umap_df, by = "cell_id", sort = FALSE)

# Enforce column order
ordered_cols <- c("cell_id",
                  "nFeature_RNA", "nCount_RNA", "percent_mt",
                  "nCount_ATAC", "TSS_enrichment", "nucleosome_signal",
                  "wnn_cluster", "cell_type_wnn",
                  "wnn_umap_1", "wnn_umap_2",
                  "rna_umap_1", "rna_umap_2",
                  "atac_umap_1", "atac_umap_2")
ordered_cols <- intersect(ordered_cols, colnames(meta_df))
meta_df      <- meta_df[, ordered_cols]

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------

message("Writing cell metadata: ", META_CSV)
write.csv(meta_df, META_CSV, row.names = FALSE, quote = FALSE)

message("Writing UMAP coordinates: ", UMAP_CSV)
write.csv(umap_df, UMAP_CSV, row.names = FALSE, quote = FALSE)

message("Writing cluster markers: ", MARKER_CSV)
write.csv(markers_out, MARKER_CSV, row.names = FALSE, quote = FALSE)

message("Saving WNN object: ", RDS_OUT)
saveRDS(obj, RDS_OUT)

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------

message("\nDone.")
message("  RDS     : ", RDS_OUT)
message("  Metadata: ", META_CSV)
message("  UMAP    : ", UMAP_CSV)
message("  Markers : ", MARKER_CSV)
message("  Clusters: ", n_clusters,
        if (!is.na(obj$cell_type_wnn[[1]])) " (annotated)" else " (unannotated — run with --cluster-annotations to assign cell types)")
