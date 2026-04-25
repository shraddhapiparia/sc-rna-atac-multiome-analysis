#!/usr/bin/env Rscript
#
# 08_motif_enrichment_atac.R
#
# Phase 4 — ATAC motif enrichment via chromVAR.
#
# Loads the Phase 2 WNN Seurat object, scans ATAC peaks for seven JASPAR2020
# motifs that correspond to the five Phase 3 TF-activity program axes, runs
# chromVAR to compute per-cell deviation z-scores, and writes flat-file
# summaries used by script 09 for cross-modal comparison.
#
# Inputs:
#   data/processed/multiome_wnn_object.rds          (from script 04)
#
# Outputs (results/regulons/):
#   atac_motif_deviations.csv     per-cell deviation z-scores (7 motifs)
#   atac_motif_by_celltype.csv    per-cell-type mean / std / n  (long form, 8 rows per cell type)
#   atac_motif_summary.md         run report
#
# Usage:
#   Rscript scripts/08_motif_enrichment_atac.R [options]
#   Rscript scripts/08_motif_enrichment_atac.R --help

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Signac)
  library(JASPAR2020)
  library(TFBSTools)
  library(chromVAR)
  library(BiocParallel)
  library(GenomeInfoDb)
})

# ---------------------------------------------------------------------------
# Path resolution (works from any working directory)
# ---------------------------------------------------------------------------

.script_file <- normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)),
  mustWork = FALSE
)
.repo_root <- dirname(dirname(.script_file))

# ---------------------------------------------------------------------------
# Motif mapping table
#
# 8 rows: MA0050.2 (IRF1) appears twice — once for interferon_alpha_beta_response
# and once for antigen_presentation_mhc (IRF1 is CIITA's DNA-binding partner).
# 7 distinct motif IDs are queried from JASPAR2020.
# ---------------------------------------------------------------------------

MOTIF_MAP <- data.frame(
  motif_id = c(
    "MA0050.2", "MA0137.3",   # interferon_alpha_beta_response
    "MA0107.1", "MA0099.3",   # inflammatory_response
    "MA0050.2",               # antigen_presentation_mhc (IRF1, shared)
    "MA0690.1",               # cytotoxic_nk_like
    "MA0080.4", "MA0466.2"    # monocyte_activation
  ),
  tf_name = c(
    "IRF1",    "STAT1",
    "RELA",    "FOS::JUN",
    "IRF1",
    "TBX21",
    "SPI1",    "CEBPB"
  ),
  program = c(
    "interferon_alpha_beta_response", "interferon_alpha_beta_response",
    "inflammatory_response",          "inflammatory_response",
    "antigen_presentation_mhc",
    "cytotoxic_nk_like",
    "monocyte_activation",            "monocyte_activation"
  ),
  stringsAsFactors = FALSE
)

DISTINCT_MOTIF_IDS <- unique(MOTIF_MAP$motif_id)   # 7 unique IDs

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------

option_list <- list(

  make_option("--input-rds",
    default = file.path(.repo_root, "data", "processed", "multiome_wnn_object.rds"),
    type    = "character",
    help    = "Path to WNN Seurat .rds from script 04 [default: %default]"),

  make_option("--results-dir",
    default = file.path(.repo_root, "results", "regulons"),
    type    = "character",
    help    = "Output directory for CSV and markdown outputs [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

stop_clear <- function(...) {
  message("ERROR: ", ...)
  quit(status = 1, save = "no")
}

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

if (!file.exists(opt$`input-rds`)) {
  stop_clear("Input .rds not found: ", opt$`input-rds`,
             "\nRun scripts/04_integrate_rna_atac_wnn.R first.")
}

# ---------------------------------------------------------------------------
# Load object and validate required structure
# ---------------------------------------------------------------------------

message("Loading WNN object: ", opt$`input-rds`)
obj <- readRDS(opt$`input-rds`)
message("  ", ncol(obj), " cells in object")

if (!"ATAC" %in% Assays(obj)) {
  stop_clear(
    "ATAC assay not found in object.\n",
    "Available assays: ", paste(Assays(obj), collapse = ", ")
  )
}

if (!"cell_type_wnn" %in% colnames(obj[[]])) {
  stop_clear(
    "'cell_type_wnn' column absent from object metadata.\n",
    "Available metadata columns: ", paste(colnames(obj[[]]), collapse = ", ")
  )
}

meta_df    <- obj[[]]
wnn_labels <- setNames(meta_df$cell_type_wnn, rownames(meta_df))

n_annotated <- sum(!is.na(wnn_labels))
if (n_annotated == 0) {
  stop_clear(
    "All cell_type_wnn values are NA — cannot compute per-celltype summaries.\n",
    "Re-run script 04 with --cluster-annotations to assign cell types."
  )
}
message("  ", n_annotated, " of ", length(wnn_labels),
        " cells have non-NA cell_type_wnn labels")

# ---------------------------------------------------------------------------
# Check and load BSgenome.Hsapiens.UCSC.hg38
# ---------------------------------------------------------------------------

if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
  stop_clear(
    "BSgenome.Hsapiens.UCSC.hg38 is not installed.\n",
    "Install with:\n",
    "  BiocManager::install(\"BSgenome.Hsapiens.UCSC.hg38\")"
  )
}
suppressPackageStartupMessages(
  library(BSgenome.Hsapiens.UCSC.hg38)
)
genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

# ---------------------------------------------------------------------------
# Load JASPAR2020 and validate all 7 motif IDs are present
#
# all_versions = TRUE is used (contract specifies FALSE, but FALSE returns only
# the latest revision per base ID and would silently exclude versioned IDs such
# as MA0050.2 if a newer revision exists in the package).
# ---------------------------------------------------------------------------

message("Loading JASPAR2020 motif database (Homo sapiens, all versions)...")
all_pwms <- TFBSTools::getMatrixSet( JASPAR2020::JASPAR2020, opts = list(species = 9606, collection = "CORE", all_versions = TRUE) )

available_ids <- sapply(all_pwms, TFBSTools::ID)
missing_ids   <- setdiff(DISTINCT_MOTIF_IDS, available_ids)

if (length(missing_ids) > 0) {
  stop_clear(
    "The following motif ID(s) were not found in the installed JASPAR2020 database:\n  ",
    paste(missing_ids, collapse = "\n  "),
    "\nVerify that the JASPAR2020 Bioconductor package is correctly installed."
  )
}
message("  All ", length(DISTINCT_MOTIF_IDS), " required motif IDs confirmed in JASPAR2020")

# Subset PFMatrixList to our 7 distinct IDs and rename elements to their motif IDs.
# Renaming ensures chromVAR stores row names as motif IDs (not TF display names),
# which makes downstream extraction unambiguous.
keep_idx   <- which(available_ids %in% DISTINCT_MOTIF_IDS)
pwm_subset <- all_pwms[keep_idx]
names(pwm_subset) <- sapply(pwm_subset, TFBSTools::ID)
message("  Prepared PFMatrixList with ", length(pwm_subset), " motifs")

# ---------------------------------------------------------------------------
# Add motif annotations to ATAC peaks
# Filter peaks to seqlevels present in BSgenome before scanning to avoid
# seqlevel mismatch errors (e.g. non-standard contigs absent from hg38).
# ---------------------------------------------------------------------------

DefaultAssay(obj) <- "ATAC"

atac_assay_full  <- obj[["ATAC"]]
atac_ranges_full <- granges(atac_assay_full)
valid_seqlevels  <- GenomeInfoDb::seqlevels(genome)

keep_mask <- as.character(GenomeInfoDb::seqnames(atac_ranges_full)) %in% valid_seqlevels
n_total   <- length(atac_ranges_full)
n_keep    <- sum(keep_mask)
n_dropped <- n_total - n_keep

if (n_dropped > 0L) {
  message("  Seqlevel filter: dropped ", n_dropped, " of ", n_total,
          " peaks whose seqnames are absent from BSgenome.Hsapiens.UCSC.hg38; ",
          n_keep, " peaks retained")
}

if (n_keep == 0L) {
  stop_clear(
    "No ATAC peaks remain after filtering to seqlevels present in ",
    "BSgenome.Hsapiens.UCSC.hg38.\n",
    "Check that the ATAC assay uses UCSC-style chromosome names (e.g. 'chr1')."
  )
}

retained_features <- rownames(atac_assay_full)[keep_mask]

if (length(retained_features) == 0L) {
  stop_clear(
    "No retained ATAC features remain after seqlevel filtering.\n",
    "Cannot continue motif scanning."
  )
}

# Pull matrices directly from the assay, not by subsetting the Seurat object
counts_mat <- GetAssayData(atac_assay_full, slot = "counts")
counts_mat <- counts_mat[retained_features, , drop = FALSE]

data_mat <- tryCatch(
  GetAssayData(atac_assay_full, slot = "data"),
  error = function(e) NULL
)
if (!is.null(data_mat) && nrow(data_mat) > 0L) {
  data_mat <- data_mat[retained_features, , drop = FALSE]
}

peak_ranges <- atac_ranges_full[keep_mask]
names(peak_ranges) <- retained_features

motif_assay <- Signac::CreateChromatinAssay(
  counts = counts_mat,
  ranges = peak_ranges,
  genome = "hg38",
  fragments = Fragments(atac_assay_full),
  annotation = Annotation(atac_assay_full)
)

# Restore normalized data if present and dimensions match
if (!is.null(data_mat) &&
    nrow(data_mat) == nrow(motif_assay) &&
    ncol(data_mat) == ncol(motif_assay)) {
  motif_assay <- SetAssayData(motif_assay, slot = "data", new.data = data_mat)
}

obj_motif <- obj
obj_motif[["ATAC"]] <- motif_assay
DefaultAssay(obj_motif) <- "ATAC"

message("  Filtered ATAC assay now contains ", nrow(obj_motif[["ATAC"]]), " peaks")

message("Adding motif annotations to ATAC peaks...")
obj_motif <- Signac::AddMotifs(
  object = obj_motif,
  genome = genome,
  pfm    = pwm_subset
)
message("  Motif annotations added")

# ---------------------------------------------------------------------------
# Run chromVAR
# ---------------------------------------------------------------------------

if (!"RunChromVAR" %in% getNamespaceExports("Signac")) {
  stop_clear(
    "The installed version of Signac does not export RunChromVAR.\n",
    "This occurs in newer Signac releases that removed the chromVAR wrapper.\n",
    "To resolve, either:\n",
    "  1. Install an older Signac version that includes RunChromVAR, e.g.:\n",
    "       remotes::install_version(\"Signac\", version = \"1.9.0\")\n",
    "  2. Refactor this script to call chromVAR::computeDeviations() directly,\n",
    "     using the peak-by-motif matrix from Signac::GetMotifMatrix()."
  )
}

message("Running chromVAR (per-cell TF motif deviation z-scores)...")
BiocParallel::register(BiocParallel::SerialParam())   # sequential for reproducibility

obj_motif <- Signac::RunChromVAR(
  object = obj_motif,
  genome = genome,
  assay  = "ATAC"
)
# Transfer the chromvar assay back onto the full object so that all downstream
# code (metadata access, cell-type labels, etc.) continues to use obj.
obj[["chromvar"]] <- obj_motif[["chromvar"]]
message("  chromVAR complete — 'chromvar' assay added to object")

# ---------------------------------------------------------------------------
# Extract deviation z-scores for the 7 distinct motif IDs
# ---------------------------------------------------------------------------

dev_mat <- tryCatch(
  GetAssayData(obj, assay = "chromvar", layer = "data"),   # Seurat v5
  error = function(e)
    GetAssayData(obj, assay = "chromvar", slot  = "data")  # Seurat v4
)

# Validate that expected IDs appear as row names.
# (After names(pwm_subset) <- motif IDs above, chromVAR should use IDs directly.)
missing_from_dev <- setdiff(DISTINCT_MOTIF_IDS, rownames(dev_mat))

if (length(missing_from_dev) > 0) {
  # Fallback: search for IDs as substrings of row names (e.g. "IRF1-MA0050.2")
  row_map <- vapply(DISTINCT_MOTIF_IDS, function(mid) {
    hit <- grep(mid, rownames(dev_mat), fixed = TRUE, value = TRUE)
    if (length(hit) == 0L) NA_character_ else hit[[1L]]
  }, character(1L))

  still_missing <- names(row_map)[is.na(row_map)]
  if (length(still_missing) > 0) {
    stop_clear(
      "chromVAR deviation matrix does not contain rows for: ",
      paste(still_missing, collapse = ", "),
      "\nAvailable row names (first 20): ",
      paste(head(rownames(dev_mat), 20L), collapse = ", ")
    )
  }
  # Rename rows to canonical motif IDs
  dev_sub <- dev_mat[row_map, , drop = FALSE]
  rownames(dev_sub) <- names(row_map)
} else {
  dev_sub <- dev_mat[DISTINCT_MOTIF_IDS, , drop = FALSE]
}

message("  Extracted deviation scores for ", nrow(dev_sub), " motifs × ",
        ncol(dev_sub), " cells")

# ---------------------------------------------------------------------------
# Per-cell CSV
# ---------------------------------------------------------------------------

message("Building per-cell deviation CSV...")
cell_ids  <- colnames(dev_sub)
ct_labels <- wnn_labels[cell_ids]   # preserves NA for unannotated cells

cell_df <- data.frame(
  cell_id       = cell_ids,
  cell_type_wnn = unname(ct_labels),
  t(as.matrix(dev_sub)),
  check.names      = FALSE,
  stringsAsFactors = FALSE,
  row.names        = NULL
)

dir.create(opt$`results-dir`, recursive = TRUE, showWarnings = FALSE)

out_cell <- file.path(opt$`results-dir`, "atac_motif_deviations.csv")
write.csv(cell_df, out_cell, row.names = FALSE, quote = FALSE)
message("  Wrote: ", out_cell)

# ---------------------------------------------------------------------------
# Per-cell-type long-form summary (8 rows × n_celltypes)
# ---------------------------------------------------------------------------

message("Computing per-cell-type motif summaries...")

annotated_mask <- !is.na(ct_labels)
dev_annotated  <- dev_sub[, annotated_mask, drop = FALSE]
ct_annotated   <- ct_labels[annotated_mask]
cell_types     <- sort(unique(ct_annotated))

# Warn for cell types with fewer than 10 cells
ct_counts <- table(ct_annotated)
low_ct    <- names(ct_counts)[ct_counts < 10L]
if (length(low_ct) > 0L) {
  message("  Warning: ", length(low_ct), " cell type(s) have fewer than 10 cells: ",
          paste(sort(low_ct), collapse = ", "))
}

# Build rows: one per (MOTIF_MAP row) × (cell type)
summary_rows <- vector("list", nrow(MOTIF_MAP) * length(cell_types))
idx <- 0L

for (i in seq_len(nrow(MOTIF_MAP))) {
  mid  <- MOTIF_MAP$motif_id[i]
  tf   <- MOTIF_MAP$tf_name[i]
  prog <- MOTIF_MAP$program[i]

  scores <- dev_annotated[mid, ]   # named numeric: cell -> deviation

  for (ct in cell_types) {
    ct_scores <- scores[ct_annotated == ct]
    n         <- length(ct_scores)
    m         <- mean(ct_scores, na.rm = TRUE)
    s         <- if (n > 1L) sd(ct_scores, na.rm = TRUE) else NA_real_

    idx <- idx + 1L
    summary_rows[[idx]] <- data.frame(
      cell_type = ct,
      motif_id  = mid,
      tf_name   = tf,
      program   = prog,
      mean      = m,
      std       = s,
      n_cells   = n,
      stringsAsFactors = FALSE
    )
  }
}

by_ct_df           <- do.call(rbind, summary_rows)
by_ct_df           <- by_ct_df[order(by_ct_df$cell_type, by_ct_df$motif_id), ]
rownames(by_ct_df) <- NULL

out_byct <- file.path(opt$`results-dir`, "atac_motif_by_celltype.csv")
write.csv(by_ct_df, out_byct, row.names = FALSE, quote = FALSE)
message("  Wrote: ", out_byct)

# ---------------------------------------------------------------------------
# Markdown summary
# ---------------------------------------------------------------------------

message("Writing summary report...")

motif_table_lines <- c(
  "| motif_id | tf_name | program |",
  "|----------|---------|---------|",
  vapply(seq_len(nrow(MOTIF_MAP)), function(i)
    sprintf("| %s | %s | %s |",
            MOTIF_MAP$motif_id[i], MOTIF_MAP$tf_name[i], MOTIF_MAP$program[i]),
    character(1L))
)

ct_count_lines <- c(
  "| cell_type | n_cells |",
  "|-----------|---------|",
  vapply(cell_types, function(ct)
    sprintf("| %s | %d |", ct, as.integer(ct_counts[ct])),
    character(1L))
)

low_block <- if (length(low_ct) > 0L){
  paste(paste0("  - ", sort(low_ct)), collapse = "\n")
}else{
  "  *(none)*"}

summary_md <- paste(
  "# ATAC Motif Enrichment Summary — Phase 4",
  "",
  "## Input",
  "",
  sprintf("- Object: `%s`", opt$`input-rds`),
  sprintf("- Total cells in object: %d", ncol(obj)),
  sprintf("- Cells with non-NA cell_type_wnn: %d", n_annotated),
  "",
  "## Motif → TF → program mapping (8 rows; MA0050.2 shared across two programs)",
  "",
  paste(motif_table_lines, collapse = "\n"),
  "",
  "## Per-cell-type cell counts",
  "",
  paste(ct_count_lines, collapse = "\n"),
  "",
  "## Warnings",
  "",
  "**Cell types with fewer than 10 cells:**",
  low_block,
  "",
  sep = "\n"
)

out_md <- file.path(opt$`results-dir`, "atac_motif_summary.md")
writeLines(summary_md, out_md)
message("  Wrote: ", out_md)

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------

message("\nDone.")
message("  Per-cell deviations  : ", out_cell)
message("  Per-cell-type summary: ", out_byct)
message("  Markdown summary     : ", out_md)
