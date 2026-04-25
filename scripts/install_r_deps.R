#!/usr/bin/env Rscript
#
# install_r_deps.R
#
# VERIFICATION AND FALLBACK script for R package dependencies.
#
# All packages listed here should already be installed by conda
# (via environment.yml).  This script checks that every package loads
# cleanly and, for any that are missing, attempts a BiocManager fallback
# install before reporting the final status.
#
# Run after activating the conda environment:
#   conda activate sc-rna-atac
#   Rscript scripts/install_r_deps.R
#
# Exit codes:
#   0 — all packages present
#   1 — one or more packages could not be loaded after fallback

message("=== sc-rna-atac R package verification ===")
message(sprintf("R %s | %s", R.version$version.string,
                if (.Machine$integer.max > 2^31) "64-bit" else "32-bit"))
message("")

# ---------------------------------------------------------------------------
# Expected packages — installed by conda (environment.yml).
# Grouped by where they come from so gaps are easy to diagnose.
# ---------------------------------------------------------------------------

# CRAN / conda-forge
cran_pkgs <- c(
  "optparse",      # CLI argument parsing       (all R scripts)
  "future",        # parallel backend            (scripts 03, 04)
  "ggplot2",       # plotting                    (script 03b)
  "patchwork",     # multi-panel figures         (script 03b)
  "Seurat",        # scRNA-seq analysis          (scripts 03, 04)
  "Signac"         # ATAC-seq analysis           (scripts 03, 04, 08)
)

# Bioconductor / bioconda
bioc_pkgs <- c(
  "GenomicRanges",                # GRanges intervals           (scripts 03, 08)
  "GenomeInfoDb",                 # seqlevels / genome metadata (scripts 03, 08)
  "BiocParallel",                 # parallel backend            (script 08)
  "chromVAR",                     # TF motif deviations         (script 08)
  "TFBSTools",                    # motif scanning              (script 08)
  "JASPAR2020",                   # motif database              (script 08)
  "SingleCellExperiment",         # Bioc container              (script 10)
  "slingshot",                    # pseudotime (script 10); expected MISSING from conda on
                                  # osx-arm64 — BiocManager fallback below installs it
  "EnsDb.Hsapiens.v86",           # hg38 gene annotations       (script 03)
  "BSgenome.Hsapiens.UCSC.hg38"   # hg38 genome sequence        (script 08)
)

all_pkgs <- c(cran_pkgs, bioc_pkgs)

# ---------------------------------------------------------------------------
# Check which packages are already loadable
# ---------------------------------------------------------------------------

message("── Checking installed packages ────────────────────────────────────")

status  <- vapply(all_pkgs, function(pkg) {
  ok  <- requireNamespace(pkg, quietly = TRUE)
  ver <- if (ok) as.character(packageVersion(pkg)) else "MISSING"
  message(sprintf("  %-42s %s", pkg, ver))
  ok
}, logical(1L))

missing_pkgs <- all_pkgs[!status]
message("")

# ---------------------------------------------------------------------------
# Fallback: attempt BiocManager install for anything still missing.
# This should only be needed if conda failed to install a package.
# ---------------------------------------------------------------------------

if (length(missing_pkgs) > 0L) {
  message("── Fallback: BiocManager install ──────────────────────────────────")
  message(sprintf("  %d package(s) not found via conda:", length(missing_pkgs)))
  for (pkg in missing_pkgs) message("    - ", pkg)
  message("")

  # Ensure BiocManager is available
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("  Installing BiocManager from CRAN...")
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    install.packages("BiocManager")
  }

  # Match Bioconductor release to this R version (auto-detected)
  r_minor <- as.integer(R.version$minor)
  r_major <- as.integer(R.version$major)
  bioc_ver <- if (r_major == 4L && r_minor >= 4L) "3.20" else
              if (r_major == 4L && r_minor == 3L) "3.18" else
              BiocManager::version()
  message(sprintf("  Using Bioconductor %s (R %s.%s detected)", bioc_ver,
                  r_major, r_minor))
  BiocManager::install(version = bioc_ver, ask = FALSE, update = FALSE)

  message("  Installing missing packages...")
  BiocManager::install(missing_pkgs, ask = FALSE, update = FALSE)

  # Re-check after fallback
  message("")
  message("── Re-checking after fallback ─────────────────────────────────────")
  status <- vapply(all_pkgs, function(pkg) {
    ok  <- requireNamespace(pkg, quietly = TRUE)
    ver <- if (ok) as.character(packageVersion(pkg)) else "STILL MISSING"
    message(sprintf("  %-42s %s", pkg, ver))
    ok
  }, logical(1L))

  missing_pkgs <- all_pkgs[!status]
  message("")
}

# ---------------------------------------------------------------------------
# Final result
# ---------------------------------------------------------------------------

if (length(missing_pkgs) > 0L) {
  message("ERROR: ", length(missing_pkgs), " package(s) could not be loaded:")
  for (pkg in missing_pkgs) message("  - ", pkg)
  message("")
  message("Possible causes:")
  message("  1. The conda environment was not activated before running this script.")
  message("     Run: conda activate sc-rna-atac")
  message("  2. env creation failed partway through (try: mamba env create -f environment.yml)")
  message("  3. A package requires a newer R version.")
  message("     R 4.4 is required; R 4.3 is not supported.")
  quit(status = 1L, save = "no")
}

n <- length(all_pkgs)
message(sprintf("All %d packages verified.", n))
message("")
message("Next: generate renv.lock to lock exact versions —")
message("  Rscript -e \"library(renv); renv::init(bare=TRUE); renv::snapshot()\"")
