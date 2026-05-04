#!/usr/bin/env Rscript
#
# 03c_doublet_qc.R
#
# Annotates doublets in a preprocessed Seurat multiome object using
# scDblFinder. Cells are annotated only — none are removed.
#
# Reads:
#   data/processed/pbmc_multiome_preprocessed.rds
#
# Writes:
#   data/processed/pbmc_multiome_preprocessed_doublet_annotated.rds
#   results/integration/doublet_qc_annotations.csv
#
# Usage:
#   Rscript scripts/03c_doublet_qc.R [options]
#   Rscript scripts/03c_doublet_qc.R --help

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
})

# ---------------------------------------------------------------------------
# Resolve repo root from script location
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

  make_option("--input-rds",
    type    = "character",
    default = file.path(.repo_root, "data", "processed",
                        "pbmc_multiome_preprocessed.rds"),
    help    = "Path to preprocessed Seurat RDS [default: %default]"),

  make_option("--output-rds",
    type    = "character",
    default = file.path(.repo_root, "data", "processed",
                        "pbmc_multiome_preprocessed_doublet_annotated.rds"),
    help    = "Path for doublet-annotated Seurat RDS [default: %default]"),

  make_option("--output-csv",
    type    = "character",
    default = file.path(.repo_root, "results", "integration",
                        "doublet_qc_annotations.csv"),
    help    = "Path for doublet annotation CSV [default: %default]"),

  make_option("--assay",
    type    = "character",
    default = "RNA",
    help    = "Assay to use for doublet detection [default: %default]"),

  make_option("--seed",
    type    = "integer",
    default = 123L,
    help    = "Random seed [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

stop_clear <- function(...) {
  message("ERROR: ", ...)
  quit(status = 1, save = "no")
}

# ---------------------------------------------------------------------------
# Check scDblFinder availability
# ---------------------------------------------------------------------------

if (!requireNamespace("scDblFinder", quietly = TRUE)) {
  stop_clear(
    "scDblFinder is not installed.\n",
    "Install it with:\n",
    "  if (!requireNamespace('BiocManager', quietly = TRUE))\n",
    "    install.packages('BiocManager')\n",
    "  BiocManager::install('scDblFinder')"
  )
}

if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  stop_clear(
    "SingleCellExperiment is not installed.\n",
    "Install it with:\n",
    "  BiocManager::install('SingleCellExperiment')"
  )
}

suppressPackageStartupMessages({
  library(scDblFinder)
  library(SingleCellExperiment)
})

# ---------------------------------------------------------------------------
# Validate inputs and create output directories
# ---------------------------------------------------------------------------

if (!file.exists(opt$`input-rds`)) {
  stop_clear("Input RDS not found: ", opt$`input-rds`,
             "\nRun scripts/03_preprocess_multiome.R first.")
}

dir.create(dirname(opt$`output-rds`), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$`output-csv`), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load Seurat object
# ---------------------------------------------------------------------------

message("Loading Seurat object: ", opt$`input-rds`)
obj <- readRDS(opt$`input-rds`)
message(format(ncol(obj), big.mark = ","), " cells loaded")

# ---------------------------------------------------------------------------
# Convert RNA counts to SingleCellExperiment and run scDblFinder
# ---------------------------------------------------------------------------

message("Converting to SingleCellExperiment (assay: ", opt$assay, ")...")

if (!opt$assay %in% names(obj@assays)) {
  stop_clear("Assay '", opt$assay, "' not found in Seurat object. ",
             "Available assays: ", paste(names(obj@assays), collapse = ", "))
}

# Use raw counts (slot "counts") for doublet detection
counts_mat <- tryCatch(
  {
    GetAssayData(obj, assay = opt$assay, layer = "counts")
  },
  error = function(e) {
    message("Could not use GetAssayData(..., layer = 'counts'); trying slot = 'counts' for older Seurat.")
    GetAssayData(obj, assay = opt$assay, slot = "counts")
  }
)

if (nrow(counts_mat) == 0 || ncol(counts_mat) == 0) {
  stop_clear("Counts matrix is empty for assay: ", opt$assay)
}

sce <- SingleCellExperiment(assays = list(counts = counts_mat))

message("Running scDblFinder (seed = ", opt$seed, ")...")
set.seed(opt$seed)
sce <- scDblFinder(sce)

# ---------------------------------------------------------------------------
# Transfer doublet annotations back to Seurat metadata
# ---------------------------------------------------------------------------

doublet_score <- sce$scDblFinder.score
doublet_class <- as.character(sce$scDblFinder.class)

obj[["doublet_score"]] <- doublet_score
obj[["doublet_class"]] <- doublet_class

n_doublet   <- sum(doublet_class == "doublet")
n_singlet   <- sum(doublet_class == "singlet")
pct_doublet <- 100 * n_doublet / ncol(obj)

message(sprintf(
  "scDblFinder complete: %s singlets, %s doublets (%.1f%%)",
  format(n_singlet,   big.mark = ","),
  format(n_doublet,   big.mark = ","),
  pct_doublet
))

# ---------------------------------------------------------------------------
# Save annotated RDS
# ---------------------------------------------------------------------------

message("Saving annotated object: ", opt$`output-rds`)
saveRDS(obj, opt$`output-rds`)

# ---------------------------------------------------------------------------
# Write annotation CSV
# ---------------------------------------------------------------------------

message("Writing annotation CSV: ", opt$`output-csv`)

# Always-present columns
csv_cols <- c("doublet_score", "doublet_class")

# Optional metadata columns — include when present
optional_cols <- c("nCount_RNA", "nFeature_RNA", "nCount_ATAC",
                   "TSS.enrichment", "nucleosome_signal", "pass_qc")
for (col in optional_cols) {
  if (col %in% colnames(obj[[]])) csv_cols <- c(csv_cols, col)
}

ann_df         <- obj[[csv_cols]]
ann_df$barcode <- rownames(ann_df)
ann_df         <- ann_df[, c("barcode", csv_cols)]

# Standardise to underscore convention
names(ann_df)[names(ann_df) == "TSS.enrichment"] <- "TSS_enrichment"

write.csv(ann_df, opt$`output-csv`, row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------

message("\nDone.")
message("  RDS : ", opt$`output-rds`)
message("  CSV : ", opt$`output-csv`)
