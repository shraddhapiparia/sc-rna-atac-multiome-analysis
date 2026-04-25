#!/usr/bin/env Rscript
#
# 10_pseudotime_subset.R
#
# Phase 5 — Monocyte pseudotime trajectory via Slingshot.
#
# Subsets the Phase 2 WNN Seurat object to monocyte-lineage cells
# (CD14+ Monocytes and FCGR3A+ Monocytes), fits a Slingshot pseudotime
# trajectory on the existing WNN UMAP embedding, and writes flat-file outputs
# structured for downstream cell-level joins with Phase 4 ATAC motif deviations
# and population-level joins with Phase 3 RNA program scores.
#
# The RNA program scores in Phase 3 originate from the Kang et al. dataset and
# have no barcode overlap with the multiome cells here.  Output 2 therefore
# annotates each cell with the population-level mean / std for its broad label,
# making the population-level limitation explicit in the file structure.  If
# per-cell RNA program scores for multiome cells are generated in a future
# phase, they can be joined to Output 1 on cell_id.
#
# Inputs:
#   data/processed/multiome_wnn_object.rds              (from script 04)
#   results/comparison/wnn_to_broad_labels.tsv           (fine → broad label map)
#   results/regulons/tf_activity_rna_by_celltype.csv     (from script 06; optional)
#
# Outputs (results/dynamics/):
#   pseudotime_monocyte_cells.csv          per-cell pseudotime; join key for ATAC
#   pseudotime_monocyte_rna_population.csv per-cell × RNA population means (if available)
#   pseudotime_summary.md                  run report and join-status table
#
# Usage:
#   Rscript scripts/10_pseudotime_subset.R [options]
#   Rscript scripts/10_pseudotime_subset.R --help

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
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
# Constants
# ---------------------------------------------------------------------------

PROGRAM_COLS <- c(
  "interferon_alpha_beta_response",
  "inflammatory_response",
  "antigen_presentation_mhc",
  "cytotoxic_nk_like",
  "monocyte_activation"
)

MONOCYTE_BROAD_LABELS <- c("CD14+ Monocytes", "FCGR3A+ Monocytes")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

stop_clear <- function(...) {
  message("ERROR: ", ...)
  quit(status = 1, save = "no")
}

# ---------------------------------------------------------------------------
# Check and load Bioconductor dependencies
# ---------------------------------------------------------------------------

for (.pkg in c("SingleCellExperiment", "slingshot")) {
  if (!requireNamespace(.pkg, quietly = TRUE)) {
    stop_clear(
      .pkg, " is not installed.\n",
      "Install with: BiocManager::install(\"", .pkg, "\")"
    )
  }
}

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot)
})

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------

option_list <- list(

  make_option("--input-rds",
    default = file.path(.repo_root, "data", "processed", "multiome_wnn_object.rds"),
    type    = "character",
    help    = "Path to WNN Seurat .rds from script 04 [default: %default]"),

  make_option("--label-map",
    default = file.path(.repo_root, "results", "comparison", "wnn_to_broad_labels.tsv"),
    type    = "character",
    help    = "Two-column TSV (fine_label<TAB>broad_label) [default: %default]"),

  make_option("--rna-activity",
    default = file.path(.repo_root, "results", "regulons",
                        "tf_activity_rna_by_celltype.csv"),
    type    = "character",
    help    = paste("Path to tf_activity_rna_by_celltype.csv (script 06).",
                    "If absent, Output 2 is skipped with a warning. [default: %default]")),

  make_option("--results-dir",
    default = file.path(.repo_root, "results", "dynamics"),
    type    = "character",
    help    = "Output directory [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$`results-dir`, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Validate and load the WNN object
# ---------------------------------------------------------------------------

if (!file.exists(opt$`input-rds`)) {
  stop_clear(
    "Input .rds not found: ", opt$`input-rds`,
    "\nRun scripts/04_integrate_rna_atac_wnn.R first."
  )
}

message("Loading WNN object: ", opt$`input-rds`)
obj            <- readRDS(opt$`input-rds`)
n_total_cells  <- ncol(obj)
message("  ", n_total_cells, " cells in object")

if (!"wnn.umap" %in% Reductions(obj)) {
  stop_clear(
    "'wnn.umap' reduction not found in object.\n",
    "Available reductions: ", paste(Reductions(obj), collapse = ", ")
  )
}

if (!"cell_type_wnn" %in% colnames(obj[[]])) {
  stop_clear(
    "'cell_type_wnn' column not found in object metadata.\n",
    "Available columns: ", paste(colnames(obj[[]]), collapse = ", ")
  )
}

# ---------------------------------------------------------------------------
# 2. Load broad-label mapping
# ---------------------------------------------------------------------------

if (!file.exists(opt$`label-map`)) {
  stop_clear("Label map TSV not found: ", opt$`label-map`)
}

label_map_df <- read.table(
  opt$`label-map`, header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

if (!all(c("fine_label", "broad_label") %in% colnames(label_map_df))) {
  stop_clear(
    "Label map is missing required columns.\n",
    "Expected: fine_label<TAB>broad_label\n",
    "Found:    ", paste(colnames(label_map_df), collapse = ", ")
  )
}

fine_to_broad <- setNames(label_map_df$broad_label, label_map_df$fine_label)

# ---------------------------------------------------------------------------
# 3. Identify monocyte fine labels and subset to monocyte cells
# ---------------------------------------------------------------------------

message("Identifying monocyte-lineage cells...")

mono_fine_labels <- names(fine_to_broad)[fine_to_broad %in% MONOCYTE_BROAD_LABELS]

if (length(mono_fine_labels) == 0L) {
  stop_clear(
    "No fine labels in wnn_to_broad_labels.tsv map to monocyte target categories:\n  ",
    paste(MONOCYTE_BROAD_LABELS, collapse = "\n  "),
    "\nBroad labels present in the TSV: ",
    paste(sort(unique(label_map_df$broad_label)), collapse = ", ")
  )
}

message("  Monocyte fine label(s) identified: ",
        paste(mono_fine_labels, collapse = ", "))

wnn_labels_all <- obj[["cell_type_wnn"]][, 1]
mono_cell_mask <- !is.na(wnn_labels_all) & wnn_labels_all %in% mono_fine_labels
mono_cells     <- colnames(obj)[mono_cell_mask]
n_mono         <- length(mono_cells)
message("  ", n_mono, " monocyte cells found in object")

if (n_mono < 50L) {
  message("  Warning: fewer than 50 monocyte cells (", n_mono,
          ") — pseudotime estimates may be unreliable")
}

meta_sub <- obj[[]][mono_cells, , drop = FALSE]
cell_type_wnn_sub <- as.character(meta_sub$cell_type_wnn)
broad_vec         <- unname(fine_to_broad[cell_type_wnn_sub])

broad_labels_in_sub <- unique(broad_vec)
n_distinct_broad    <- length(broad_labels_in_sub)

if (n_distinct_broad < 2L) {
  stop_clear(
    "Only ", n_distinct_broad, " distinct broad label(s) present in the monocyte subset",
    " (", paste(broad_labels_in_sub, collapse = ", "), ").\n",
    "Slingshot requires at least 2 cluster labels to fit a lineage."
  )
}

message("  Broad labels in subset: ",
        paste(sort(broad_labels_in_sub), collapse = ", "))

# ---------------------------------------------------------------------------
# 4. Run Slingshot on the WNN UMAP embedding
# ---------------------------------------------------------------------------

message("Running Slingshot pseudotime on WNN UMAP embedding...")

cell_ids_sub <- mono_cells
umap_emb     <- Embeddings(obj, "wnn.umap")[cell_ids_sub, , drop = FALSE]

# Build a minimal SingleCellExperiment: UMAP in reducedDims, broad labels in colData
sce <- SingleCellExperiment(
  colData = data.frame(
    label     = broad_vec,
    row.names = cell_ids_sub,
    stringsAsFactors = FALSE
  )
)
reducedDim(sce, "UMAP") <- umap_emb

sds <- slingshot::slingshot(
  sce,
  clusterLabels = "label",
  reducedDim    = "UMAP",
  start.clus    = "CD14+ Monocytes",
  end.clus      = "FCGR3A+ Monocytes"
)

pt_matrix  <- slingshot::slingPseudotime(sds)
n_lineages <- ncol(pt_matrix)
message("  Slingshot identified ", n_lineages, " lineage(s)")

if (n_lineages > 1L) {
  message("  Warning: only lineage 1 is used; remaining lineages are ignored")
}

pseudotime <- pt_matrix[, 1]   # named numeric: cell_id -> raw pseudotime

# NA handling
n_na_pt <- sum(is.na(pseudotime))
if (n_na_pt > 0L) {
  message("  Warning: ", n_na_pt, " of ", length(pseudotime),
          " cells (", round(100 * n_na_pt / length(pseudotime), 1),
          "%) have NA pseudotime — included in outputs with NA values")
}

if (sum(!is.na(pseudotime)) == 0L) {
  stop_clear(
    "All pseudotime values are NA.\n",
    "Check that the WNN UMAP embedding separates CD14+ and FCGR3A+ monocytes ",
    "and that both broad labels are present in the subset."
  )
}

# Min-max scale to [0, 1]
pt_min            <- min(pseudotime, na.rm = TRUE)
pt_max            <- max(pseudotime, na.rm = TRUE)
pseudotime_scaled <- (pseudotime - pt_min) / (pt_max - pt_min)

message("  Pseudotime range (raw): [", round(pt_min, 4), ", ", round(pt_max, 4), "]")

# ---------------------------------------------------------------------------
# 5. Build Output 1: pseudotime_monocyte_cells.csv
#
# Canonical per-cell table.  cell_id is the join key for all downstream
# cell-level merges (currently ATAC motif deviations; future RNA cell-level).
# ---------------------------------------------------------------------------

message("Building Output 1 (per-cell pseudotime)...")

out1 <- data.frame(
  cell_id           = cell_ids_sub,
  cell_type_wnn     = cell_type_wnn_sub,
  broad_label       = broad_vec,
  wnn_umap_1        = umap_emb[cell_ids_sub, 1],
  wnn_umap_2        = umap_emb[cell_ids_sub, 2],
  pseudotime        = pseudotime[cell_ids_sub],
  pseudotime_scaled = pseudotime_scaled[cell_ids_sub],
  stringsAsFactors  = FALSE,
  row.names         = NULL
)

# Sort: broad_label ascending, pseudotime ascending, NAs last
out1 <- out1[order(out1$broad_label, out1$pseudotime, na.last = TRUE), ]
rownames(out1) <- NULL

out1_path <- file.path(opt$`results-dir`, "pseudotime_monocyte_cells.csv")
write.csv(out1, out1_path, row.names = FALSE, quote = FALSE)
message("  Wrote: ", out1_path)

# ---------------------------------------------------------------------------
# 6. Build Output 2: pseudotime_monocyte_rna_population.csv (conditional)
#
# Annotates each monocyte cell with the RNA baseline mean / std for its
# broad_label from tf_activity_rna_by_celltype.csv.  Every cell of the same
# broad_label receives identical RNA values — these are population statistics,
# NOT per-cell measurements for the multiome cells.
# ---------------------------------------------------------------------------

write_rna_output <- TRUE
rna_skip_reason  <- ""
rna_col_names    <- character(0L)
out2_path        <- file.path(opt$`results-dir`, "pseudotime_monocyte_rna_population.csv")

# --- Check file availability ---
if (!file.exists(opt$`rna-activity`)) {
  message("  Warning: RNA activity file not found: ", opt$`rna-activity`,
          " — skipping Output 2")
  write_rna_output <- FALSE
  rna_skip_reason  <- paste0("File not found: ", opt$`rna-activity`)
}

# --- Check required columns ---
if (write_rna_output) {
  rna_act          <- read.csv(opt$`rna-activity`, stringsAsFactors = FALSE)
  req_rna_cols     <- c("cell_type", "program", "mean", "std")
  missing_rna_cols <- setdiff(req_rna_cols, colnames(rna_act))

  if (length(missing_rna_cols) > 0L) {
    message("  Warning: RNA activity file missing column(s): ",
            paste(missing_rna_cols, collapse = ", "),
            " — skipping Output 2")
    write_rna_output <- FALSE
    rna_skip_reason  <- paste0("Missing column(s): ",
                                paste(missing_rna_cols, collapse = ", "))
  }
}

# --- Build and write Output 2 ---
if (write_rna_output) {
  rna_mono <- rna_act[
    rna_act$cell_type %in% MONOCYTE_BROAD_LABELS &
    rna_act$program   %in% PROGRAM_COLS,
  ]

  missing_rna_labels <- setdiff(MONOCYTE_BROAD_LABELS, rna_mono$cell_type)
  if (length(missing_rna_labels) > 0L) {
    message("  Warning: RNA activity file does not contain broad label(s): ",
            paste(missing_rna_labels, collapse = ", "),
            " — those cells will have NA in RNA columns")
  }

  # Wide lookup: one row per broad label, columns = <program>_rna_mean / _rna_std
  rna_wide <- data.frame(
    cell_type        = MONOCYTE_BROAD_LABELS,
    stringsAsFactors = FALSE
  )

  for (prog in PROGRAM_COLS) {
    prog_rows     <- rna_mono[rna_mono$program == prog, , drop = FALSE]
    mean_col_name <- paste0(prog, "_rna_mean")
    std_col_name  <- paste0(prog, "_rna_std")

    if (nrow(prog_rows) == 0L) {
      rna_wide[[mean_col_name]] <- NA_real_
      rna_wide[[std_col_name]]  <- NA_real_
    } else {
      mean_lookup <- setNames(prog_rows$mean, prog_rows$cell_type)
      std_lookup  <- setNames(prog_rows$std,  prog_rows$cell_type)
      # unname() prevents named-vector confusion during data.frame assignment
      rna_wide[[mean_col_name]] <- unname(mean_lookup[rna_wide$cell_type])
      rna_wide[[std_col_name]]  <- unname(std_lookup[rna_wide$cell_type])
    }
  }

  # Interleaved column order: mean then std for each program in PROGRAM_COLS order
  rna_col_names <- as.vector(rbind(
    paste0(PROGRAM_COLS, "_rna_mean"),
    paste0(PROGRAM_COLS, "_rna_std")
  ))

  # Join: all monocyte cells × RNA wide table on broad_label
  out2_base <- out1[, c("cell_id", "broad_label", "pseudotime_scaled")]

  out2 <- merge(
    out2_base,
    rna_wide,
    by.x  = "broad_label",
    by.y  = "cell_type",
    all.x = TRUE,   # retain all cells; NA for broad labels absent from RNA file
    sort  = FALSE
  )

  # Enforce column order per contract
  out2 <- out2[, c("cell_id", "broad_label", "pseudotime_scaled", rna_col_names)]

  # Sort: broad_label ascending, pseudotime_scaled ascending, NAs last
  out2 <- out2[order(out2$broad_label, out2$pseudotime_scaled, na.last = TRUE), ]
  rownames(out2) <- NULL

  write.csv(out2, out2_path, row.names = FALSE, quote = FALSE)
  message("  Wrote: ", out2_path)
}

# ---------------------------------------------------------------------------
# 7. Write pseudotime_summary.md
# ---------------------------------------------------------------------------

message("Writing summary report...")

# -- Per-broad-label cell counts --------------------------------------------
broad_count_lines <- c(
  "| broad_label | n_cells |",
  "|-------------|---------|",
  vapply(sort(MONOCYTE_BROAD_LABELS), function(bl) {
    sprintf("| %s | %d |", bl, sum(out1$broad_label == bl, na.rm = TRUE))
  }, character(1L))
)

# -- Per-cell_type_wnn counts -----------------------------------------------
ct_counts <- sort(table(out1$cell_type_wnn), decreasing = TRUE)
wnn_count_lines <- c(
  "| cell_type_wnn | broad_label | n_cells |",
  "|---------------|-------------|---------|",
  vapply(names(ct_counts), function(ct) {
    bls <- unique(out1$broad_label[out1$cell_type_wnn == ct])
    bl_str <- if (length(bls) == 1L) bls else paste(bls, collapse = " / ")
    sprintf("| %s | %s | %d |", ct, bl_str, as.integer(ct_counts[ct]))
  }, character(1L))
)

# -- Pseudotime range per broad label ---------------------------------------
pt_range_lines <- c(
  "| broad_label | min | mean | max |",
  "|-------------|-----|------|-----|",
  vapply(sort(MONOCYTE_BROAD_LABELS), function(bl) {
    sub_pt <- out1$pseudotime[out1$broad_label == bl & !is.na(out1$pseudotime)]
    if (length(sub_pt) == 0L)
      sprintf("| %s | NA | NA | NA |", bl)
    else
      sprintf("| %s | %.4f | %.4f | %.4f |",
              bl, min(sub_pt), mean(sub_pt), max(sub_pt))
  }, character(1L))
)

# -- Output 2 block ---------------------------------------------------------
if (write_rna_output) {
  out2_block <- paste(
    sprintf("- Path: `%s`", out2_path),
    "- Columns: cell_id, broad_label, pseudotime_scaled,",
    paste0("  ", paste(rna_col_names, collapse = ", ")),
    "- **Important:** `_rna_mean` / `_rna_std` columns are population-level statistics",
    "  from `tf_activity_rna_by_celltype.csv` (Phase 3, Kang et al. dataset).",
    "  Every cell with the same `broad_label` carries identical RNA values.",
    "  These are NOT per-cell measurements for multiome cells.",
    "  If per-cell RNA program scores for multiome cells are generated in a future",
    "  phase, join them to Output 1 on `cell_id` instead.",
    sep = "\n"
  )
} else {
  out2_block <- sprintf("- **Skipped** — %s", rna_skip_reason)
}

# -- Join status table ------------------------------------------------------
rna_join_status <- if (write_rna_output){
  "Ready (population means only)"
}else{
  paste0("Skipped — ", rna_skip_reason)}

rna_join_source <- if (write_rna_output){
  "`results/regulons/tf_activity_rna_by_celltype.csv`"
}else{
  "—"}

# -- Assemble markdown -------------------------------------------------------
summary_md <- paste(
  "# Pseudotime Analysis Summary — Phase 5",
  "",
  "## Input",
  "",
  sprintf("- Object: `%s`", opt$`input-rds`),
  sprintf("- Total cells in full object: %d", n_total_cells),
  "",
  "## Monocyte subset",
  "",
  paste(broad_count_lines, collapse = "\n"),
  "",
  sprintf("**Total monocyte cells:** %d", nrow(out1)),
  "",
  "### Cells by WNN fine label",
  "",
  paste(wnn_count_lines, collapse = "\n"),
  "",
  "## Slingshot pseudotime",
  "",
  sprintf("- Lineages identified: **%d** (using lineage 1 only%s)",
          n_lineages,
          if (n_lineages > 1L) "; additional lineages ignored" else ""),
  sprintf("- NA-pseudotime cells: **%d** (%.1f%% of subset)",
          n_na_pt, 100 * n_na_pt / nrow(out1)),
  "",
  "### Pseudotime range by broad label",
  "",
  paste(pt_range_lines, collapse = "\n"),
  "",
  "## Outputs",
  "",
  "### Output 1: `pseudotime_monocyte_cells.csv`",
  "",
  sprintf("- Path: `%s`", out1_path),
  "- Columns: cell_id, cell_type_wnn, broad_label, wnn_umap_1, wnn_umap_2,",
  "  pseudotime, pseudotime_scaled",
  "- Primary join key for cell-level downstream joins (`cell_id`).",
  "",
  "### Output 2: `pseudotime_monocyte_rna_population.csv`",
  "",
  out2_block,
  "",
  "## Downstream join status",
  "",
  "| Downstream join | Join key | Status | Source file |",
  "|-----------------|----------|--------|-------------|",
  paste0("| ATAC motif deviations (cell-level) | `cell_id` | Ready | ",
         "`results/regulons/atac_motif_deviations.csv` |"),
  paste0("| RNA program scores (population-level) | `broad_label` | ",
         rna_join_status, " | ", rna_join_source, " |"),
  paste0("| RNA program scores (cell-level) | `cell_id` | ",
         "Not yet available — requires per-cell RNA scoring of multiome cells | — |"),
  "",
  sep = "\n"
)

out_md <- file.path(opt$`results-dir`, "pseudotime_summary.md")
writeLines(summary_md, out_md)
message("  Wrote: ", out_md)

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------

message("\nDone.")
message("  Output 1 (per-cell pseudotime)  : ", out1_path)
if (write_rna_output) {
  message("  Output 2 (RNA population means) : ", out2_path)
} else {
  message("  Output 2 (RNA population means) : skipped — ", rna_skip_reason)
}
message("  Summary                         : ", out_md)
