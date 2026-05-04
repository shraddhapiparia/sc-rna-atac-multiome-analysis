#!/usr/bin/env Rscript
#
# 03d_plot_atac_fragment_qc.R
#
# ATAC fragment-length QC visualization using Signac::FragmentHistogram.
# Loads the preprocessed Seurat/Signac object, groups cells by nucleosome
# signal (if available), and saves a fragment-length histogram.
#
# Reads:
#   data/processed/pbmc_multiome_preprocessed.rds
#
# Writes:
#   results/figures/qc/atac_fragment_histogram.{fmt}
#
# Usage:
#   Rscript scripts/03d_plot_atac_fragment_qc.R [options]
#   Rscript scripts/03d_plot_atac_fragment_qc.R --help

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Signac)
  library(ggplot2)
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

  make_option("--outdir",
    type    = "character",
    default = file.path(.repo_root, "results", "figures", "qc"),
    help    = "Output directory for figures [default: %default]"),

  make_option("--nucleosome-threshold",
    type    = "double",
    default = 2.0,
    help    = "nucleosome_signal cutoff for Low/High grouping [default: %default]"),

  make_option("--dpi",
    type    = "integer",
    default = 300L,
    help    = "Figure resolution in DPI [default: %default]"),

  make_option("--fmt",
    type    = "character",
    default = "png",
    help    = "Output format: 'png' or 'pdf' [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!opt$fmt %in% c("png", "pdf")) {
  stop("--fmt must be 'png' or 'pdf'.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

stop_clear <- function(...) {
  message("ERROR: ", ...)
  quit(status = 1, save = "no")
}

# ---------------------------------------------------------------------------
# Validate input and create output directory
# ---------------------------------------------------------------------------

if (!file.exists(opt$`input-rds`)) {
  stop_clear("Input RDS not found: ", opt$`input-rds`,
             "\nRun scripts/03_preprocess_multiome.R first.")
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load Seurat object
# ---------------------------------------------------------------------------

message("Loading Seurat object: ", opt$`input-rds`)
obj <- readRDS(opt$`input-rds`)
message(format(ncol(obj), big.mark = ","), " cells loaded")

# ---------------------------------------------------------------------------
# Validate ATAC assay and fragment file
# ---------------------------------------------------------------------------

if (!"ATAC" %in% names(obj@assays)) {
  stop_clear("ATAC assay not found in the Seurat object.")
}

DefaultAssay(obj) <- "ATAC"

# Check that at least one fragment file is registered and accessible
frags <- Fragments(obj)
if (length(frags) == 0L) {
  stop_clear(
    "No fragment files are registered with the ATAC assay.\n",
    "Ensure the object was created with a valid fragments path."
  )
}

frag_path <- frags[[1]]@path
if (!file.exists(frag_path)) {
  stop_clear(
    "Fragment file not found at registered path:\n  ", frag_path,
    "\nMove the file back to that location or re-run script 03 ",
    "with the correct --raw-dir."
  )
}

message("Fragment file OK: ", frag_path)

# ---------------------------------------------------------------------------
# Add nucleosome_group grouping variable (if nucleosome_signal is present)
# ---------------------------------------------------------------------------

has_nucleosome <- "nucleosome_signal" %in% colnames(obj[[]])

if (has_nucleosome) {
  message(sprintf(
    "Grouping cells by nucleosome_signal (threshold = %.1f)...",
    opt$`nucleosome-threshold`
  ))
  obj[["nucleosome_group"]] <- ifelse(
    obj[["nucleosome_signal"]][[1]] < opt$`nucleosome-threshold`,
    "Low nucleosome",
    "High nucleosome"
  )
  group_by_col <- "nucleosome_group"
} else {
  message("nucleosome_signal not found in metadata; plotting without grouping.")
  group_by_col <- NULL
}

# ---------------------------------------------------------------------------
# Build fragment-length histogram using Signac
# ---------------------------------------------------------------------------

message("Building ATAC fragment-length histogram...")

if (!is.null(group_by_col)) {
  p <- FragmentHistogram(
    object   = obj,
    group.by = group_by_col
  ) +
    scale_color_manual(
      values = c("Low nucleosome"  = "#2166AC",
                 "High nucleosome" = "#D73027"),
      name   = "Nucleosome group"
    ) +
    labs(
      title    = "ATAC Fragment-Length Distribution",
      subtitle = sprintf(
        "Nucleosome groups: nucleosome_signal %s %.1f vs ≥ %.1f",
        "<", opt$`nucleosome-threshold`, opt$`nucleosome-threshold`
      ),
      x = "Fragment length (bp)",
      y = "Proportion"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "grey45"),
      legend.position = "top"
    )
} else {
  p <- FragmentHistogram(object = obj) +
    labs(
      title = "ATAC Fragment-Length Distribution",
      x     = "Fragment length (bp)",
      y     = "Proportion"
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12))
}

# ---------------------------------------------------------------------------
# Save figure
# ---------------------------------------------------------------------------

out_path <- file.path(opt$outdir, paste0("atac_fragment_histogram.", opt$fmt))

ggsave(
  filename = out_path,
  plot     = p,
  width    = 7,
  height   = 4.5,
  dpi      = opt$dpi,
  units    = "in"
)

message("Wrote: ", out_path)

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------

message("\nDone.")
