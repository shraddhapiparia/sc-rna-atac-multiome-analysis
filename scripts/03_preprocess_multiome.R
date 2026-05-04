#!/usr/bin/env Rscript
#
# 03_preprocess_multiome.R
#
# QC, filtering, normalisation, and dimensionality reduction for the
# 10x PBMC Multiome 10k dataset.
#
# Reads:
#   data/raw/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
#   data/raw/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz[.tbi]
#
# Writes:
#   data/processed/pbmc_multiome_preprocessed.rds
#   results/integration/preprocessing_qc_stats.csv
#   results/integration/multiome_qc_summary.md
#
# No UMAP or WNN are computed here — handled in script 04.
#
# Usage:
#   Rscript scripts/03_preprocess_multiome.R [options]
#   Rscript scripts/03_preprocess_multiome.R --help

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Signac)
  library(future)
  library(EnsDb.Hsapiens.v86)
  library(GenomeInfoDb)
  library(GenomicRanges)
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

  # Directories
  make_option("--raw-dir",
    default = file.path(.repo_root, "data", "raw"),
    type = "character",
    help = "Directory containing raw 10x files [default: %default]"),
  make_option("--processed-dir",
    default = file.path(.repo_root, "data", "processed"),
    type = "character",
    help = "Output directory for the .rds object [default: %default]"),
  make_option("--results-dir",
    default = file.path(.repo_root, "results", "integration"),
    type = "character",
    help = "Output directory for QC artefacts [default: %default]"),

  # RNA QC thresholds
  make_option("--min-features-rna", type = "integer", default = 200L,
    help = "Minimum nFeature_RNA [default: %default]"),
  make_option("--max-features-rna", type = "integer", default = 7500L,
    help = "Maximum nFeature_RNA [default: %default]"),
  make_option("--min-count-rna", type = "integer", default = 0L,
    help = "Minimum nCount_RNA; 0 = no lower bound [default: %default]"),
  make_option("--max-count-rna", type = "integer", default = .Machine$integer.max,
    help = "Maximum nCount_RNA; no upper bound by default"),
  make_option("--max-pct-mt", type = "double", default = 20.0,
    help = "Maximum percent mitochondrial reads [default: %default]"),

  # ATAC QC thresholds
  make_option("--min-count-atac", type = "integer", default = 1000L,
    help = "Minimum nCount_ATAC [default: %default]"),
  make_option("--max-count-atac", type = "integer", default = 100000L,
    help = "Maximum nCount_ATAC [default: %default]"),
  make_option("--min-tss", type = "double", default = 2.0,
    help = "Minimum TSS enrichment score [default: %default]"),
  make_option("--max-nucleosome", type = "double", default = 2.0,
    help = "Maximum nucleosome signal [default: %default]"),

  # Optional filters
  make_option("--blacklist-filter", action = "store_true", default = FALSE,
    help = "Apply ENCODE hg38 blacklist ratio filter"),
  make_option("--max-blacklist-ratio", type = "double", default = 0.05,
    help = "Maximum blacklist_ratio when --blacklist-filter is set [default: %default]"),
  make_option("--blacklist-source", type = "character", default = "builtin",
    help = paste("Source for hg38 blacklist regions when --blacklist-filter is set.",
                 "Allowed: builtin (Signac built-in), annotationhub, none (skip).",
                 "[default: %default]")),

  # Normalisation
  make_option("--no-sct", action = "store_true", default = FALSE,
    help = "Use NormalizeData + FindVariableFeatures + ScaleData instead of SCTransform"),

  # Dimensionality reduction
  make_option("--n-dims", type = "integer", default = 50L,
    help = "Number of PCA / LSI dimensions [default: %default]"),

  make_option("--seed", type = "integer", default = 42L,
    help = "Random seed [default: %default]"),

  # FRiP-like ATAC metric from 10x per-barcode metadata
  make_option("--per-barcode-metrics", type = "character", default = NULL,
    help = paste("Optional path to 10x per-barcode metrics CSV",
                 "(columns: barcode, peak_region_fragments, passed_filters)",
                 "[default: %default]")),
  make_option("--min-pct-reads-in-peaks", type = "double", default = 15.0,
    help = "Minimum percent reads in peaks (FRiP-like) when --per-barcode-metrics is supplied [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!opt$`blacklist-source` %in% c("builtin", "annotationhub", "none")) {
  stop("--blacklist-source must be one of: builtin, annotationhub, none", call. = FALSE)
}

set.seed(opt$seed)

# SCTransform uses future for parallelism; on this ~10k-cell dataset the
# globals it tries to export can exceed future's default 500 MiB limit and
# cause an error. Cap at 3 GiB and force sequential execution to avoid it.
options(future.globals.maxSize = 3 * 1024^3)
future::plan("sequential")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

stop_clear <- function(...) {
  message("ERROR: ", ...)
  quit(status = 1, save = "no")
}

fmt_int <- function(x) formatC(as.integer(x), format = "d", big.mark = ",")

# ---------------------------------------------------------------------------
# Validate inputs and set up output directories
# ---------------------------------------------------------------------------

H5_FILE   <- file.path(opt$`raw-dir`,
  "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
FRAG_FILE <- file.path(opt$`raw-dir`,
  "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz")
RDS_OUT   <- file.path(opt$`processed-dir`, "pbmc_multiome_preprocessed.rds")
QC_CSV    <- file.path(opt$`results-dir`,   "preprocessing_qc_stats.csv")
QC_MD     <- file.path(opt$`results-dir`,   "multiome_qc_summary.md")

for (f in c(H5_FILE, FRAG_FILE, paste0(FRAG_FILE, ".tbi"))) {
  if (!file.exists(f)) stop_clear("Required input file not found: ", f)
}

dir.create(opt$`processed-dir`, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$`results-dir`,   recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load raw data
# ---------------------------------------------------------------------------

message("Loading h5 file: ", H5_FILE)
counts <- Read10X_h5(H5_FILE)

rna_counts  <- counts[["Gene Expression"]]
atac_counts <- counts[["Peaks"]]
n_cells_raw <- ncol(rna_counts)
message(fmt_int(n_cells_raw), " barcodes in raw data")

# ---------------------------------------------------------------------------
# Create Seurat object (RNA assay)
# ---------------------------------------------------------------------------

message("Creating Seurat object (RNA)...")
obj <- CreateSeuratObject(
  counts       = rna_counts,
  assay        = "RNA",
  min.cells    = 10,   # genes detected in < 10 cells are dropped
  min.features = 0     # per-cell filtering applied explicitly below
)

# ---------------------------------------------------------------------------
# Build genome annotation from EnsDb (hg38)
# ---------------------------------------------------------------------------

message("Building hg38 genome annotation from EnsDb.Hsapiens.v86...")
annotations <- suppressMessages(
  GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations)         <- "hg38"

# ---------------------------------------------------------------------------
# Add ATAC assay
# ---------------------------------------------------------------------------

message("Creating ChromatinAssay (ATAC)...")

# Align ATAC barcodes to those retained in the RNA object
shared_barcodes <- intersect(colnames(obj), colnames(atac_counts))
n_cells_shared  <- length(shared_barcodes)
if (n_cells_shared < ncol(obj)) {
  message("  Note: ", fmt_int(ncol(obj) - n_cells_shared),
          " RNA barcodes absent from ATAC peak matrix — these will be dropped")
  obj <- subset(obj, cells = shared_barcodes)
}
message(fmt_int(n_cells_shared), " shared RNA+ATAC barcodes entering QC")

atac_assay <- CreateChromatinAssay(
  counts      = atac_counts[, colnames(obj)],
  sep         = c(":", "-"),
  fragments   = FRAG_FILE,
  annotation  = annotations,
  min.cells   = 10,
  min.features = 0
)

obj[["ATAC"]] <- atac_assay

# ---------------------------------------------------------------------------
# RNA QC metrics
# ---------------------------------------------------------------------------

message("Computing RNA QC metrics...")
DefaultAssay(obj) <- "RNA"
# nFeature_RNA and nCount_RNA are set automatically by CreateSeuratObject
obj[["percent_mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-", assay = "RNA")

# ---------------------------------------------------------------------------
# ATAC QC metrics
# ---------------------------------------------------------------------------

message("Computing ATAC QC metrics...")
DefaultAssay(obj) <- "ATAC"

obj <- NucleosomeSignal(obj)

# TSSEnrichment: fast = FALSE computes per-cell scores (required for filtering)
obj <- TSSEnrichment(obj, fast = FALSE)

blacklist_available <- FALSE

if (!opt$`blacklist-filter`) {
  message("Blacklist filtering not requested (--blacklist-filter not set); skipping.")
} else if (opt$`blacklist-source` == "none") {
  message("Blacklist filtering skipped: --blacklist-source is 'none'.")
} else {
  bl_source <- opt$`blacklist-source`
  message("Computing blacklist ratio (source: ", bl_source, ")...")

  bl_regions <- if (bl_source == "builtin") {
    # Try the global environment first (Signac <= 1.3 exports it on attach),
    # then fall back to the package namespace (Signac >= 1.4 keeps it there).
    bl <- if (exists("blacklist_hg38_unified")) {
      blacklist_hg38_unified
    } else {
      tryCatch(
        get("blacklist_hg38_unified", envir = asNamespace("Signac")),
        error = function(e) NULL
      )
    }
    if (is.null(bl)) {
      stop_clear(
        "blacklist_hg38_unified not found in the Signac package.\n",
        "Options:\n",
        "  Upgrade Signac (>= 1.1):   install.packages('Signac')\n",
        "  Use AnnotationHub:          --blacklist-source annotationhub\n",
        "  Skip blacklist filtering:   --blacklist-source none"
      )
    }
    bl
  } else {
    # annotationhub — only loaded on demand; not a required dependency
    if (!requireNamespace("AnnotationHub", quietly = TRUE)) {
      stop_clear(
        "AnnotationHub is not installed.\n",
        "Install with:  BiocManager::install('AnnotationHub')\n",
        "Or use the Signac built-in object:  --blacklist-source builtin"
      )
    }
    tryCatch({
      ah  <- AnnotationHub::AnnotationHub()
      q   <- AnnotationHub::query(ah, c("hg38", "blacklist", "ENCODE"))
      ids <- names(q)
      if (length(ids) == 0L) {
        stop_clear(
          "No hg38 blacklist records found in AnnotationHub.\n",
          "Try --blacklist-source builtin instead."
        )
      }
      message("  Using AnnotationHub record: ", ids[1])
      bl <- ah[[ids[1]]]
      seqlevelsStyle(bl) <- "UCSC"
      bl
    }, error = function(e) {
      stop_clear(
        "Failed to retrieve hg38 blacklist from AnnotationHub:\n  ",
        conditionMessage(e), "\n",
        "Check your internet connection, or use --blacklist-source builtin."
      )
    })
  }

  obj[["blacklist_ratio"]] <- FractionCountsInRegion(
    obj,
    assay   = "ATAC",
    regions = bl_regions
  )
  blacklist_available <- TRUE
  message("  blacklist_ratio computed for ", fmt_int(ncol(obj)), " cells.")
}

DefaultAssay(obj) <- "RNA"

# ---------------------------------------------------------------------------
# Optional FRiP-like ATAC metric from 10x per-barcode metrics
# ---------------------------------------------------------------------------

frip_available <- FALSE

if (!is.null(opt$`per-barcode-metrics`)) {
  pbm_file <- opt$`per-barcode-metrics`
  if (!file.exists(pbm_file)) {
    stop_clear("--per-barcode-metrics file not found: ", pbm_file)
  }
  message("Reading per-barcode metrics: ", pbm_file)
  pbm <- read.csv(pbm_file, stringsAsFactors = FALSE)

  required_cols <- c("barcode", "peak_region_fragments", "passed_filters")
  missing_cols  <- setdiff(required_cols, colnames(pbm))
  if (length(missing_cols) > 0) {
    stop_clear(
      "per-barcode metrics CSV is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  rownames(pbm) <- pbm$barcode
  cells         <- colnames(obj)
  matched       <- pbm[cells, , drop = FALSE]   # NA rows for unmatched cells

  prf <- matched$peak_region_fragments
  pf  <- matched$passed_filters
  # Guard against division by zero and missing values
  pct <- ifelse(!is.na(prf) & !is.na(pf) & pf > 0, 100 * prf / pf, NA_real_)

  obj[["peak_region_fragments"]] <- prf
  obj[["passed_filters"]]        <- pf
  obj[["pct_reads_in_peaks"]]    <- pct

  n_frip_na <- sum(is.na(pct))
  message(
    "FRiP-like metric (pct_reads_in_peaks) computed. ",
    fmt_int(sum(!is.na(pct))), " cells with values; ",
    fmt_int(n_frip_na), " NAs (unmatched or zero passed_filters)."
  )
  frip_available <- TRUE
} else {
  message("--per-barcode-metrics not supplied; skipping FRiP/pct_reads_in_peaks.")
}

# ---------------------------------------------------------------------------
# Determine pass_qc — one fail_* column per criterion, then combine.
# Computed before subsetting so the full pre-filter object is auditable
# via the QC CSV.
# ---------------------------------------------------------------------------

message("Applying QC thresholds...")
md <- obj[[]]

# Individual per-criterion fail flags (TRUE = cell failed that criterion)
obj[["fail_rna_features"]] <-
  is.na(md$nFeature_RNA) |
  !(md$nFeature_RNA >= opt$`min-features-rna` &
    md$nFeature_RNA <= opt$`max-features-rna`)

obj[["fail_rna_counts"]] <-
  is.na(md$nCount_RNA) |
  !(md$nCount_RNA >= opt$`min-count-rna` &
    md$nCount_RNA <= opt$`max-count-rna`)
    
obj[["fail_percent_mt"]] <-
  is.na(md$percent_mt) |
  md$percent_mt >= opt$`max-pct-mt`

obj[["fail_atac_counts"]] <-
  is.na(md$nCount_ATAC) |
  !(md$nCount_ATAC >= opt$`min-count-atac` &
    md$nCount_ATAC <= opt$`max-count-atac`)

obj[["fail_tss_enrichment"]] <-
  is.na(md$TSS.enrichment) |
  md$TSS.enrichment < opt$`min-tss`

obj[["fail_nucleosome_signal"]] <-
  is.na(md$nucleosome_signal) |
  md$nucleosome_signal >= opt$`max-nucleosome`

if (blacklist_available) {
  obj[["fail_blacklist_ratio"]] <-
    is.na(md$blacklist_ratio) |
    md$blacklist_ratio >= opt$`max-blacklist-ratio`
}

if (frip_available) {
  # Cells with NA pct_reads_in_peaks also fail (unmatched / zero denominator)
  obj[["fail_pct_reads_in_peaks"]] <-
    is.na(md$pct_reads_in_peaks) |
    md$pct_reads_in_peaks < opt$`min-pct-reads-in-peaks`
}

# Refresh md after adding fail_* columns
md <- obj[[]]

# pass_qc is the inverse union of all active fail columns
fail_cols <- grep("^fail_", colnames(md), value = TRUE)
pass <- rowSums(md[, fail_cols, drop = FALSE]) == 0

obj[["pass_qc"]] <- pass
n_pass <- sum(pass)
message(fmt_int(n_pass), " / ", fmt_int(n_cells_shared), " shared RNA+ATAC cells pass QC")

# Per-criterion failure counts for reporting
fail_counts <- colSums(md[, fail_cols, drop = FALSE])

# ---------------------------------------------------------------------------
# Write QC stats CSV (all pre-filter cells, for audit and downstream plots)
# ---------------------------------------------------------------------------

message("Writing QC stats CSV: ", QC_CSV)

qc_cols <- c("nFeature_RNA", "nCount_RNA", "percent_mt",
             "nCount_ATAC", "nFeature_ATAC",
             "TSS.enrichment", "nucleosome_signal",
             "pass_qc",
             "fail_rna_features", "fail_rna_counts", "fail_percent_mt",
             "fail_atac_counts", "fail_tss_enrichment", "fail_nucleosome_signal")
if (blacklist_available)    qc_cols <- c(qc_cols, "blacklist_ratio", "fail_blacklist_ratio")
if (frip_available)         qc_cols <- c(qc_cols, "peak_region_fragments",
                                         "passed_filters", "pct_reads_in_peaks",
                                         "fail_pct_reads_in_peaks")

qc_df         <- obj[[qc_cols]]
qc_df$cell_id <- rownames(qc_df)
qc_df         <- qc_df[, c("cell_id", qc_cols)]

# Standardise column name to underscore convention for downstream Python phases
names(qc_df)[names(qc_df) == "TSS.enrichment"] <- "TSS_enrichment"

write.csv(qc_df, QC_CSV, row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# Subset to passing cells
# ---------------------------------------------------------------------------

obj <- subset(obj, subset = pass_qc == TRUE)

# Store cell_id as explicit metadata column
obj[["cell_id"]] <- colnames(obj)

# ---------------------------------------------------------------------------
# RNA normalisation and PCA
# ---------------------------------------------------------------------------

DefaultAssay(obj) <- "RNA"

if (!opt$`no-sct`) {
  message("Normalising RNA with SCTransform...")
  obj <- SCTransform(obj, vars.to.regress = "percent_mt", verbose = FALSE)
  obj <- RunPCA(obj, assay = "SCT", npcs = opt$`n-dims`, verbose = FALSE,
                reduction.name = "pca")
} else {
  message("Normalising RNA with NormalizeData + ScaleData (--no-sct)...")
  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = 1e4, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 3000, verbose = FALSE)
  obj <- ScaleData(obj, vars.to.regress = "percent_mt", verbose = FALSE)
  obj <- RunPCA(obj, npcs = opt$`n-dims`, verbose = FALSE,
                reduction.name = "pca")
}

# ---------------------------------------------------------------------------
# ATAC: TF-IDF + LSI
# Note: dim 1 of LSI typically correlates with sequencing depth and is
# dropped when building the WNN graph in script 04.
# ---------------------------------------------------------------------------

message("Computing LSI for ATAC...")
DefaultAssay(obj) <- "ATAC"

obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = "q0")   # all peaks retained post-QC
obj <- RunSVD(obj, n = opt$`n-dims`, reduction.name = "lsi")

# ---------------------------------------------------------------------------
# Save preprocessed object
# ---------------------------------------------------------------------------

message("Saving preprocessed object: ", RDS_OUT)
saveRDS(obj, RDS_OUT)

# ---------------------------------------------------------------------------
# Write QC summary markdown
# ---------------------------------------------------------------------------

message("Writing QC summary markdown: ", QC_MD)

norm_mode <- if (!opt$`no-sct`) "SCTransform" else
  "NormalizeData + FindVariableFeatures + ScaleData"

bl_line <- if (blacklist_available) {
  paste0("| blacklist_ratio (source: ", opt$`blacklist-source`, ") | < ",
         opt$`max-blacklist-ratio`, " |\n")
} else if (opt$`blacklist-filter`) {
  paste0("| blacklist_ratio | not computed (--blacklist-source ",
         opt$`blacklist-source`, ") |\n")
} else {
  "| blacklist_ratio filter | not requested |\n"
}

frip_line <- if (frip_available) {
  paste0("| pct_reads_in_peaks (FRiP) | >= ", opt$`min-pct-reads-in-peaks`, " |\n")
} else {
  "| pct_reads_in_peaks (FRiP) | not computed (--per-barcode-metrics not supplied) |\n"
}

# Per-criterion failure rows for the markdown table
fail_rows <- paste0(
  "| ", names(fail_counts), " | ", fail_counts, " |\n",
  collapse = ""
)

qc_md <- paste0(
  "# Multiome QC Summary\n\n",

  "## Cell counts\n\n",
  "| Stage | Cells |\n",
  "|-------|-------|\n",
  "| Raw RNA barcodes | ", fmt_int(n_cells_raw), " |\n",
  "| Shared RNA+ATAC barcodes (pre-QC) | ", fmt_int(n_cells_shared), " |\n",
  "| Pass QC | ", fmt_int(n_pass), " |\n",
  "| Fail QC (removed) | ", fmt_int(n_cells_shared - n_pass), " |\n\n",

  "## Cells failing each QC criterion\n\n",
  "| Criterion | Cells failing |\n",
  "|-----------|---------------|\n",
  fail_rows, "\n",

  "## QC thresholds applied\n\n",
  "| Metric | Threshold |\n",
  "|--------|-----------|\n",
  "| nFeature_RNA | ", opt$`min-features-rna`, " – ", opt$`max-features-rna`, " |\n",
  "| nCount_RNA | ", opt$`min-count-rna`, " – ", opt$`max-count-rna`, " |\n",
  "| percent_mt | < ", opt$`max-pct-mt`, " |\n",
  "| nCount_ATAC | ", opt$`min-count-atac`, " – ", opt$`max-count-atac`, " |\n",
  "| TSS_enrichment | >= ", opt$`min-tss`, " |\n",
  "| nucleosome_signal | < ", opt$`max-nucleosome`, " |\n",
  bl_line,
  frip_line, "\n",

  "## Normalisation\n\n",
  "| Modality | Method |\n",
  "|----------|--------|\n",
  "| RNA | ", norm_mode, " |\n",
  "| ATAC | TF-IDF + LSI (", opt$`n-dims`, " dims) |\n\n",

  "## Reductions stored in object\n\n",
  "| Reduction | Dims | Notes |\n",
  "|-----------|------|-------|\n",
  "| `pca` | ", opt$`n-dims`, " | RNA |\n",
  "| `lsi` | ", opt$`n-dims`, " | ATAC; dim 1 dropped in script 04 |\n"
)

writeLines(qc_md, QC_MD)

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------

message("\nDone.")
message("  RDS : ", RDS_OUT)
message("  CSV : ", QC_CSV)
message("  MD  : ", QC_MD)
