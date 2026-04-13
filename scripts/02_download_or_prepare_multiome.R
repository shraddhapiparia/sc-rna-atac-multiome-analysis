#!/usr/bin/env Rscript
#
# 02_download_or_prepare_multiome.R
#
# Validate (or optionally download) the three raw 10x PBMC Multiome files
# required for Phase 2, then write data/raw/manifest.txt.
#
# Modes:
#   validate (default) -- confirm all three files are present and compute md5s
#   download           -- fetch files from 10x CDN, then validate
#
# Usage:
#   Rscript scripts/02_download_or_prepare_multiome.R              # validate
#   Rscript scripts/02_download_or_prepare_multiome.R --download   # download then validate
#   Rscript scripts/02_download_or_prepare_multiome.R --raw-dir /absolute/path/to/data/raw

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Resolve repo root from script location so the script works regardless of the
# working directory from which it is invoked.
.script_file <- normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)),
  mustWork = FALSE
)
.repo_root <- dirname(dirname(.script_file))   # scripts/ is one level below repo root

RAW_DIR      <- file.path(.repo_root, "data", "raw")
MANIFEST     <- file.path(RAW_DIR, "manifest.txt")
BASE_URL     <- paste0(
  "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/",
  "pbmc_granulocyte_sorted_10k/"
)

TARGET_FILES <- c(
  "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5",
  "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
  "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi"
)

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------

args     <- commandArgs(trailingOnly = TRUE)
download <- "--download" %in% args

# --raw-dir <path> overrides the default derived from script location
raw_dir_flag <- grep("^--raw-dir$", args)
if (length(raw_dir_flag) > 0) {
  RAW_DIR  <- normalizePath(args[raw_dir_flag + 1], mustWork = FALSE)
  MANIFEST <- file.path(RAW_DIR, "manifest.txt")
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

stop_clear <- function(...) {
  message("ERROR: ", ...)
  quit(status = 1, save = "no")
}

# ---------------------------------------------------------------------------
# Optionally download missing files
# ---------------------------------------------------------------------------

if (download) {
  if (!dir.exists(RAW_DIR)) dir.create(RAW_DIR, recursive = TRUE)

  for (fname in TARGET_FILES) {
    dest <- file.path(RAW_DIR, fname)
    if (file.exists(dest)) {
      message("Already present, skipping: ", fname)
      next
    }
    url <- paste0(BASE_URL, fname)
    message("Downloading: ", fname)
    tryCatch(
      download.file(url, destfile = dest, mode = "wb", quiet = FALSE),
      error = function(e) stop_clear("Download failed for ", fname, ": ", conditionMessage(e))
    )
  }
}

# ---------------------------------------------------------------------------
# Validate presence
# ---------------------------------------------------------------------------

missing_files <- TARGET_FILES[
  !file.exists(file.path(RAW_DIR, TARGET_FILES))
]

if (length(missing_files) > 0) {
  stop_clear(
    length(missing_files), " required file(s) not found in ", RAW_DIR, ":\n",
    paste0("  ", missing_files, collapse = "\n"), "\n",
    "Re-run with --download to fetch them automatically."
  )
}

# ---------------------------------------------------------------------------
# Compute md5 checksums and write manifest
# ---------------------------------------------------------------------------

paths <- file.path(RAW_DIR, TARGET_FILES)
md5s  <- tools::md5sum(paths)   # named character vector: path -> md5

manifest_lines <- paste(TARGET_FILES, md5s, sep = "\t")
writeLines(manifest_lines, MANIFEST)

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

message("All files present. Manifest written to: ", MANIFEST)
for (line in manifest_lines) message("  ", line)
