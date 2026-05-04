#!/usr/bin/env Rscript
#
# 04b_integration_diagnostics.R
#
# Lightweight QC/diagnostic script comparing before vs after QC filtering
# and summarising integration structure. Reads existing CSV and RDS outputs —
# does NOT rerun preprocessing or WNN.
#
# Reads (required):
#   results/integration/preprocessing_qc_stats.csv
#   results/integration/multiome_cell_metadata.csv
#
# Reads (optional — skipped gracefully if absent):
#   data/processed/pbmc_multiome_preprocessed.rds   (PCA elbow, LSI variance,
#                                                    SCT assay info)
#   data/processed/multiome_wnn_object.rds           (RNA/ATAC modality weights)
#
# Writes:
#   results/figures/integration_diagnostics/
#     rna_qc_before_after.{fmt}
#     atac_qc_filtered.{fmt}
#     pca_elbow.{fmt}                  (if preprocessed RDS found)
#     lsi_variance.{fmt}               (if preprocessed RDS found)
#     wnn_umap_cluster.{fmt}
#     modality_umap_comparison.{fmt}
#     wnn_umap_rna_weight.{fmt}        (if RNA.weight found in WNN object)
#     cluster_sizes.{fmt}
#   results/integration/integration_diagnostics_summary.md
#
# Usage:
#   Rscript scripts/04b_integration_diagnostics.R
#   Rscript scripts/04b_integration_diagnostics.R --dpi 150 --fmt pdf
#   Rscript scripts/04b_integration_diagnostics.R --help

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(patchwork)
})

# ---------------------------------------------------------------------------
# Resolve repo root from working directory (same convention as 03b_plot_qc.R)
# ---------------------------------------------------------------------------

.repo_root <- normalizePath(getwd(), mustWork = TRUE)

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------

option_list <- list(

  make_option("--qc-csv",
    default = file.path(.repo_root, "results", "integration",
                        "preprocessing_qc_stats.csv"),
    type = "character",
    help = "preprocessing_qc_stats.csv from script 03 [default: %default]"),

  make_option("--meta-csv",
    default = file.path(.repo_root, "results", "integration",
                        "multiome_cell_metadata.csv"),
    type = "character",
    help = "multiome_cell_metadata.csv from script 04 [default: %default]"),

  make_option("--preprocessed-rds",
    default = file.path(.repo_root, "data", "processed",
                        "pbmc_multiome_preprocessed.rds"),
    type = "character",
    help = "Preprocessed Seurat .rds from script 03 (optional) [default: %default]"),

  make_option("--wnn-rds",
    default = file.path(.repo_root, "data", "processed",
                        "multiome_wnn_object.rds"),
    type = "character",
    help = "WNN Seurat .rds from script 04 (optional) [default: %default]"),

  make_option("--outdir",
    default = file.path(.repo_root, "results", "figures",
                        "integration_diagnostics"),
    type = "character",
    help = "Output directory for figures [default: %default]"),

  make_option("--summary-md",
    default = file.path(.repo_root, "results", "integration",
                        "integration_diagnostics_summary.md"),
    type = "character",
    help = "Path for the summary markdown file [default: %default]"),

  make_option("--dpi",
    default = 300L, type = "integer",
    help = "Figure resolution in DPI [default: %default]"),

  make_option("--fmt",
    default = "png", type = "character",
    help = "Output format: 'png' or 'pdf' [default: %default]"),

  make_option("--skip-rds",
    action = "store_true", default = FALSE,
    help = "Skip loading RDS files (CSV-only mode)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!opt$fmt %in% c("png", "pdf")) {
  stop("--fmt must be 'png' or 'pdf'.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

msg <- function(...) message("[04b] ", ...)
skip <- function(...) message("[04b] SKIP: ", ...)

save_fig <- function(plot, stem, width, height) {
  path <- file.path(opt$outdir, paste0(stem, ".", opt$fmt))
  ggsave(filename = path, plot = plot,
         width = width, height = height,
         dpi = opt$dpi, units = "in")
  msg("  Wrote: ", path)
  invisible(path)
}

# Track what was generated and what was skipped, for the summary MD
generated <- character(0)
skipped   <- character(0)

note_generated <- function(label) { generated <<- c(generated, label) }
note_skipped   <- function(label) { skipped   <<- c(skipped,   label) }

# ---------------------------------------------------------------------------
# Design constants (consistent with 03b_plot_qc.R)
# ---------------------------------------------------------------------------

COL_PASS   <- "#2166AC"
COL_FAIL   <- "#D73027"
COL_THRESH <- "#E6550D"
COL_BOX    <- "grey25"

CLUSTER_PALETTE <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
  "#D37295", "#FABFD2", "#A0CBE8", "#86BCB6", "#8CD17D"
)

theme_diag <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title             = element_text(size = 10),
      axis.text              = element_text(size = 9),
      panel.grid.major.y     = element_line(color = "grey92", linewidth = 0.3),
      legend.title           = element_text(size = 9, face = "bold"),
      legend.text            = element_text(size = 8),
      legend.key.size        = unit(0.55, "lines"),
      plot.title             = element_text(face = "bold", size = 11),
      plot.subtitle          = element_text(size = 9, color = "grey45"),
      plot.caption           = element_text(size = 7.5, color = COL_THRESH,
                                            hjust = 0),
      plot.caption.position  = "plot"
    )
}

fill_pf <- scale_fill_manual(
  name   = "QC status",
  values = c("Pass" = COL_PASS, "Fail" = COL_FAIL)
)

# Violin helper for a single QC metric
violin_panel <- function(df, yvar, ylab, title,
                         thresh_min = NULL, thresh_max = NULL,
                         log_y = FALSE) {
  if (!yvar %in% colnames(df)) {
    return(ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0(yvar, "\nnot available"),
               size = 4, color = "grey60") +
      theme_void())
  }
  hlines <- list()
  if (!is.null(thresh_min) && thresh_min > 0)
    hlines <- c(hlines, list(
      geom_hline(yintercept = thresh_min, linetype = "dashed",
                 color = COL_THRESH, linewidth = 0.6)))
  if (!is.null(thresh_max) && is.finite(thresh_max))
    hlines <- c(hlines, list(
      geom_hline(yintercept = thresh_max, linetype = "dashed",
                 color = COL_THRESH, linewidth = 0.6)))

  p <- ggplot(df, aes(x = qc_group, y = .data[[yvar]], fill = qc_group)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.75, color = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA,
                 fill = "white", color = COL_BOX, linewidth = 0.35) +
    hlines +
    fill_pf +
    labs(title = title, x = NULL, y = ylab) +
    theme_diag() +
    theme(legend.position = "none")
  if (log_y)
    p <- p + scale_y_log10(labels = function(x) prettyNum(x, big.mark = ","))
  p
}

# UMAP scatter helper
umap_scatter <- function(df, colour_var, title,
                         x_col = "wnn_umap_1", y_col = "wnn_umap_2",
                         palette = NULL, continuous = FALSE) {
  if (!all(c(x_col, y_col, colour_var) %in% colnames(df))) return(NULL)
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                      colour = .data[[colour_var]])) +
    geom_point(size = 0.4, alpha = 0.5, stroke = 0) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    guides(colour = guide_legend(
      override.aes = list(size = 2.5, alpha = 1),
      ncol = 2
    )) +
    theme_diag() +
    theme(axis.ticks = element_blank(),
          axis.text  = element_blank(),
          axis.line  = element_blank())
  if (continuous) {
    p <- p + scale_colour_viridis_c(option = "plasma", name = colour_var)
  } else if (!is.null(palette)) {
    n_levels <- length(unique(df[[colour_var]]))
    pal      <- rep_len(palette, n_levels)
    p        <- p + scale_colour_manual(name = colour_var, values = pal)
  }
  p
}

dir.create(opt$outdir,              recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$`summary-md`), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Load required CSVs
# ---------------------------------------------------------------------------

for (f in c(opt$`qc-csv`, opt$`meta-csv`)) {
  if (!file.exists(f)) stop("Required file not found: ", f, call. = FALSE)
}

msg("Loading QC stats: ", opt$`qc-csv`)
qc <- read.csv(opt$`qc-csv`, stringsAsFactors = FALSE)
qc$pass_qc  <- as.logical(qc$pass_qc)
qc$qc_group <- factor(ifelse(qc$pass_qc, "Pass", "Fail"),
                      levels = c("Fail", "Pass"))

msg("Loading cell metadata: ", opt$`meta-csv`)
meta <- read.csv(opt$`meta-csv`, stringsAsFactors = FALSE)
meta$wnn_cluster <- as.character(meta$wnn_cluster)

n_qc_total <- nrow(qc)
n_qc_pass  <- sum(qc$pass_qc, na.rm = TRUE)
n_qc_fail  <- n_qc_total - n_qc_pass
n_cells_wnn <- nrow(meta)
n_clusters  <- length(unique(meta$wnn_cluster))

msg(sprintf("QC CSV: %s total  |  %s pass  |  %s fail",
            format(n_qc_total, big.mark = ","),
            format(n_qc_pass,  big.mark = ","),
            format(n_qc_fail,  big.mark = ",")))
msg(sprintf("Metadata: %s cells  |  %s WNN clusters",
            format(n_cells_wnn, big.mark = ","), n_clusters))

# ---------------------------------------------------------------------------
# 2. Load optional RDS files
# ---------------------------------------------------------------------------

preproc_obj <- NULL
wnn_obj     <- NULL

if (!opt$`skip-rds`) {

  if (file.exists(opt$`preprocessed-rds`)) {
    msg("Loading preprocessed RDS: ", opt$`preprocessed-rds`)
    msg("  (this may take a moment...)")
    preproc_obj <- tryCatch(
      readRDS(opt$`preprocessed-rds`),
      error = function(e) {
        skip("Could not load preprocessed RDS: ", conditionMessage(e))
        NULL
      }
    )
  } else {
    skip("Preprocessed RDS not found: ", opt$`preprocessed-rds`)
  }

  if (file.exists(opt$`wnn-rds`)) {
    msg("Loading WNN RDS: ", opt$`wnn-rds`)
    msg("  (this may take a moment...)")
    wnn_obj <- tryCatch(
      readRDS(opt$`wnn-rds`),
      error = function(e) {
        skip("Could not load WNN RDS: ", conditionMessage(e))
        NULL
      }
    )
  } else {
    skip("WNN RDS not found: ", opt$`wnn-rds`)
  }

} else {
  skip("RDS loading disabled via --skip-rds")
}

# Helper: extract metadata column from Seurat object safely
seurat_meta_col <- function(obj, col) {
  if (is.null(obj)) return(NULL)
  md <- tryCatch(obj[[]], error = function(e) NULL)
  if (is.null(md) || !col %in% colnames(md)) return(NULL)
  md[[col]]
}

# Helper: extract reduction stdev safely
reduction_stdev <- function(obj, reduc_name) {
  if (is.null(obj)) return(NULL)
  tryCatch({
    if (!reduc_name %in% names(obj@reductions)) return(NULL)
    obj@reductions[[reduc_name]]@stdev
  }, error = function(e) NULL)
}

# Helper: list assays safely
seurat_assays <- function(obj) {
  if (is.null(obj)) return(character(0))
  tryCatch(names(obj@assays), error = function(e) character(0))
}

# ---------------------------------------------------------------------------
# 3. RNA before/after filtering diagnostics
# ---------------------------------------------------------------------------

msg("Building Figure 1: RNA QC before/after filtering...")

subtitle_all <- sprintf(
  "%s total  ·  %s pass (%.1f%%)  ·  %s fail",
  format(n_qc_total, big.mark = ","),
  format(n_qc_pass,  big.mark = ","),
  100 * n_qc_pass / n_qc_total,
  format(n_qc_fail,  big.mark = ",")
)

# Panel row 1: all cells (coloured by pass/fail)
p_feat_all  <- violin_panel(qc, "nFeature_RNA", "Genes detected",
                             "nFeature_RNA",
                             thresh_min = 200, thresh_max = 7500)
p_count_all <- violin_panel(qc, "nCount_RNA", "UMI count (log₁₀)",
                             "nCount_RNA", log_y = TRUE)
p_mt_all    <- violin_panel(qc, "percent_mt", "Mitochondrial reads (%)",
                             "percent_mt", thresh_max = 20)

# Panel row 2: passing cells only (uniform fill)
qc_pass <- qc[qc$pass_qc, ]
qc_pass$qc_group <- factor("Pass", levels = c("Fail", "Pass"))

p_feat_pass  <- violin_panel(qc_pass, "nFeature_RNA", "Genes detected",
                              "nFeature_RNA (pass only)")
p_count_pass <- violin_panel(qc_pass, "nCount_RNA", "UMI count (log₁₀)",
                              "nCount_RNA (pass only)", log_y = TRUE)
p_mt_pass    <- violin_panel(qc_pass, "percent_mt",
                              "Mitochondrial reads (%)",
                              "percent_mt (pass only)")

fig_rna <- (p_feat_all | p_count_all | p_mt_all) /
           (p_feat_pass | p_count_pass | p_mt_pass) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "RNA QC: all cells (top) vs post-filter passing cells (bottom)",
    subtitle = subtitle_all,
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "grey45"),
      plot.margin   = margin(8, 8, 4, 8)
    )
  ) & theme(legend.position = "bottom")

save_fig(fig_rna, "rna_qc_before_after", width = 11, height = 8)
note_generated("rna_qc_before_after")

# ---------------------------------------------------------------------------
# 4. ATAC QC diagnostics (passing cells only)
# ---------------------------------------------------------------------------

msg("Building Figure 2: ATAC QC (filtered cells)...")

# Use metadata CSV for passing-cell ATAC metrics; QC CSV has nFeature_ATAC too
atac_df <- meta
atac_df$qc_group <- factor("Pass", levels = c("Fail", "Pass"))

p_atac_count <- violin_panel(atac_df, "nCount_ATAC",
                              "Peak count (log₁₀)", "nCount_ATAC",
                              thresh_min = 1000, thresh_max = 100000,
                              log_y = TRUE)
p_tss        <- violin_panel(atac_df, "TSS_enrichment",
                              "TSS enrichment score", "TSS enrichment",
                              thresh_min = 2)
p_nuc        <- violin_panel(atac_df, "nucleosome_signal",
                              "Nucleosome signal", "nucleosome_signal",
                              thresh_max = 2)

# nFeature_ATAC is in QC CSV but not metadata — add if available
qc_pass_atac <- qc[qc$pass_qc, ]
qc_pass_atac$qc_group <- factor("Pass", levels = c("Fail", "Pass"))

if ("nFeature_ATAC" %in% colnames(qc_pass_atac)) {
  p_nfeat_atac <- violin_panel(qc_pass_atac, "nFeature_ATAC",
                                "Peaks detected", "nFeature_ATAC",
                                log_y = TRUE)
  fig_atac <- (p_atac_count | p_tss | p_nuc | p_nfeat_atac) +
    plot_layout(guides = "collect", ncol = 4)
} else {
  fig_atac <- (p_atac_count | p_tss | p_nuc) +
    plot_layout(guides = "collect")
}

fig_atac <- fig_atac +
  plot_annotation(
    title    = "ATAC QC Metrics — post-filter passing cells",
    subtitle = sprintf("%s cells", format(n_cells_wnn, big.mark = ",")),
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "grey45"),
      plot.margin   = margin(8, 8, 4, 8)
    )
  ) & theme(legend.position = "none")

save_fig(fig_atac, "atac_qc_filtered", width = 11, height = 4.5)
note_generated("atac_qc_filtered")

# ---------------------------------------------------------------------------
# 5. PCA elbow plot (from preprocessed RDS)
# ---------------------------------------------------------------------------

pca_stdev <- reduction_stdev(preproc_obj, "pca")

if (!is.null(pca_stdev) && length(pca_stdev) > 0) {
  msg("Building Figure 3: PCA elbow plot...")
  pca_var <- pca_stdev^2 / sum(pca_stdev^2) * 100
  pca_df  <- data.frame(dim = seq_along(pca_var), pct_var = pca_var)

  fig_pca <- ggplot(pca_df, aes(x = dim, y = pct_var)) +
    geom_point(size = 1.8, colour = COL_PASS) +
    geom_line(colour = COL_PASS, linewidth = 0.5) +
    scale_x_continuous(breaks = pretty(pca_df$dim, n = 10)) +
    labs(
      title    = "PCA elbow — RNA (post-SCTransform / normalisation)",
      subtitle = sprintf("%d PCs computed", nrow(pca_df)),
      x        = "Principal Component",
      y        = "Variance explained (%)"
    ) +
    theme_diag()

  save_fig(fig_pca, "pca_elbow", width = 7, height = 4.5)
  note_generated("pca_elbow")
} else {
  skip("PCA elbow: preprocessed RDS not available or 'pca' reduction not found")
  note_skipped("pca_elbow (preprocessed RDS unavailable)")
}

# ---------------------------------------------------------------------------
# 6. LSI variance bar chart (from preprocessed RDS)
# ---------------------------------------------------------------------------

lsi_stdev <- reduction_stdev(preproc_obj, "lsi")

if (!is.null(lsi_stdev) && length(lsi_stdev) > 1) {
  msg("Building Figure 4: LSI variance bar chart...")
  lsi_var <- lsi_stdev^2 / sum(lsi_stdev^2) * 100
  lsi_df  <- data.frame(
    dim     = seq_along(lsi_var),
    pct_var = lsi_var,
    used    = ifelse(seq_along(lsi_var) == 1, "Dim 1 (excluded)", "Used in WNN")
  )

  fig_lsi <- ggplot(lsi_df, aes(x = dim, y = pct_var, fill = used)) +
    geom_col(width = 0.8) +
    scale_fill_manual(
      name   = NULL,
      values = c("Dim 1 (excluded)" = COL_FAIL, "Used in WNN" = COL_PASS)
    ) +
    scale_x_continuous(breaks = pretty(lsi_df$dim, n = 10)) +
    labs(
      title    = "LSI variance by dimension — ATAC (TF-IDF + SVD)",
      subtitle = "Dim 1 typically correlates with sequencing depth and is excluded from WNN",
      x        = "LSI Dimension",
      y        = "Variance explained (%)"
    ) +
    theme_diag() +
    theme(legend.position = "top")

  save_fig(fig_lsi, "lsi_variance", width = 8, height = 4.5)
  note_generated("lsi_variance")
} else {
  skip("LSI variance: preprocessed RDS not available or 'lsi' reduction not found")
  note_skipped("lsi_variance (preprocessed RDS unavailable)")
}

# ---------------------------------------------------------------------------
# 7. WNN UMAP coloured by cluster
# ---------------------------------------------------------------------------

msg("Building Figure 5: WNN UMAP by cluster...")

meta$wnn_cluster_fct <- factor(meta$wnn_cluster,
                                levels = sort(unique(meta$wnn_cluster)))

fig_umap_cluster <- umap_scatter(
  meta, "wnn_cluster_fct",
  title   = "WNN UMAP — coloured by cluster",
  palette = CLUSTER_PALETTE
)

if (!is.null(fig_umap_cluster)) {
  fig_umap_cluster <- fig_umap_cluster +
    labs(subtitle = sprintf("%s cells  ·  %d clusters",
                            format(n_cells_wnn, big.mark = ","), n_clusters))
  save_fig(fig_umap_cluster, "wnn_umap_cluster", width = 6.5, height = 5.5)
  note_generated("wnn_umap_cluster")
} else {
  skip("WNN UMAP cluster: required columns not found in metadata CSV")
  note_skipped("wnn_umap_cluster")
}

# ---------------------------------------------------------------------------
# 8. Three-panel modality UMAP comparison
#    (RNA UMAP | ATAC UMAP | WNN UMAP, all coloured by cluster)
# ---------------------------------------------------------------------------

msg("Building Figure 6: Modality UMAP comparison...")

p_rna_umap  <- umap_scatter(meta, "wnn_cluster_fct",
                             "RNA-only UMAP",
                             x_col = "rna_umap_1", y_col = "rna_umap_2",
                             palette = CLUSTER_PALETTE)
p_atac_umap <- umap_scatter(meta, "wnn_cluster_fct",
                             "ATAC-only UMAP",
                             x_col = "atac_umap_1", y_col = "atac_umap_2",
                             palette = CLUSTER_PALETTE)
p_wnn_umap2 <- umap_scatter(meta, "wnn_cluster_fct",
                             "WNN UMAP",
                             palette = CLUSTER_PALETTE)

modality_panels <- Filter(Negate(is.null),
                           list(p_rna_umap, p_atac_umap, p_wnn_umap2))

if (length(modality_panels) >= 2) {
  fig_modality <- wrap_plots(modality_panels, ncol = length(modality_panels)) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = "Per-modality vs WNN joint UMAP (coloured by WNN cluster)",
      subtitle = "RNA (PCA), ATAC (LSI), and joint WNN embeddings shown side-by-side",
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey45"),
        plot.margin   = margin(8, 8, 4, 8)
      )
    ) & theme(legend.position = "bottom")

  save_fig(fig_modality, "modality_umap_comparison",
           width = 5.5 * length(modality_panels), height = 5.5)
  note_generated("modality_umap_comparison")
} else {
  skip("Modality UMAP comparison: fewer than 2 UMAP embeddings available")
  note_skipped("modality_umap_comparison")
}

# ---------------------------------------------------------------------------
# 9. WNN UMAP coloured by RNA modality weight (from WNN RDS)
# ---------------------------------------------------------------------------

rna_weights <- seurat_meta_col(wnn_obj, "RNA.weight")

if (!is.null(rna_weights) && length(rna_weights) == n_cells_wnn) {
  msg("Building Figure 7: WNN UMAP by RNA modality weight...")
  meta$RNA.weight <- rna_weights

  fig_rna_weight <- umap_scatter(
    meta, "RNA.weight",
    title      = "WNN UMAP — RNA modality weight",
    continuous = TRUE
  )
  if (!is.null(fig_rna_weight)) {
    fig_rna_weight <- fig_rna_weight +
      labs(
        subtitle = sprintf(
          "RNA weight: median %.2f  ·  range [%.2f, %.2f]",
          median(rna_weights, na.rm = TRUE),
          min(rna_weights, na.rm = TRUE),
          max(rna_weights, na.rm = TRUE)
        )
      )
    save_fig(fig_rna_weight, "wnn_umap_rna_weight", width = 6.5, height = 5.5)
    note_generated("wnn_umap_rna_weight")
  }
} else {
  skip("RNA.weight not found in WNN object (WNN RDS absent or column missing)")
  note_skipped("wnn_umap_rna_weight (RNA.weight unavailable)")
}

# ---------------------------------------------------------------------------
# 10. Cluster size bar chart + table
# ---------------------------------------------------------------------------

msg("Building Figure 8: Cluster size bar chart...")

cluster_sizes <- as.data.frame(table(wnn_cluster = meta$wnn_cluster),
                                stringsAsFactors = FALSE)
colnames(cluster_sizes) <- c("wnn_cluster", "n_cells")
cluster_sizes$wnn_cluster <- factor(cluster_sizes$wnn_cluster,
                                    levels = cluster_sizes$wnn_cluster[
                                      order(as.numeric(cluster_sizes$wnn_cluster))
                                    ])
cluster_sizes$pct <- round(100 * cluster_sizes$n_cells / n_cells_wnn, 1)

pal_clust <- setNames(
  rep_len(CLUSTER_PALETTE, nlevels(cluster_sizes$wnn_cluster)),
  levels(cluster_sizes$wnn_cluster)
)

fig_clust_size <- ggplot(cluster_sizes,
                          aes(x = wnn_cluster, y = n_cells,
                              fill = wnn_cluster, label = n_cells)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(vjust = -0.4, size = 3, color = "grey30") +
  scale_fill_manual(values = pal_clust) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = "WNN cluster sizes",
    subtitle = sprintf("%s cells  ·  %d clusters", format(n_cells_wnn, big.mark = ","),
                       n_clusters),
    x        = "WNN cluster",
    y        = "Number of cells"
  ) +
  theme_diag()

save_fig(fig_clust_size, "cluster_sizes", width = 6, height = 4.5)
note_generated("cluster_sizes")

# ---------------------------------------------------------------------------
# 11. SCT assay summary from preprocessed RDS (informational only)
# ---------------------------------------------------------------------------

assays_preproc <- seurat_assays(preproc_obj)
assays_wnn     <- seurat_assays(wnn_obj)

sct_info <- list(found = FALSE)

if ("SCT" %in% assays_preproc) {
  sct_info$found  <- TRUE
  sct_info$source <- "preprocessed RDS"
  n_sct_genes <- tryCatch(
    nrow(preproc_obj@assays[["SCT"]]),
    error = function(e) NA_integer_
  )
  sct_info$n_genes <- n_sct_genes
  msg(sprintf("SCT assay found in preprocessed RDS (%s genes)",
              ifelse(is.na(n_sct_genes), "?", format(n_sct_genes, big.mark = ","))))
} else if ("SCT" %in% assays_wnn) {
  sct_info$found  <- TRUE
  sct_info$source <- "WNN RDS"
  n_sct_genes <- tryCatch(
    nrow(wnn_obj@assays[["SCT"]]),
    error = function(e) NA_integer_
  )
  sct_info$n_genes <- n_sct_genes
  msg(sprintf("SCT assay found in WNN RDS (%s genes)",
              ifelse(is.na(n_sct_genes), "?", format(n_sct_genes, big.mark = ","))))
} else if (length(assays_preproc) > 0 || length(assays_wnn) > 0) {
  all_assays <- union(assays_preproc, assays_wnn)
  msg("No SCT assay found; assay(s) in objects: ", paste(all_assays, collapse = ", "))
  skip("SCT assay summary (SCT assay not present — LogNormalize path used)")
  note_skipped("SCT assay info (LogNormalize was used instead of SCTransform)")
} else {
  skip("SCT assay summary (no RDS objects loaded)")
  note_skipped("SCT assay info (RDS not loaded)")
}

# ---------------------------------------------------------------------------
# 12. Write summary markdown
# ---------------------------------------------------------------------------

msg("Writing summary markdown: ", opt$`summary-md`)

norm_mode_str <- if (sct_info$found) {
  paste0("SCTransform (SCT assay in ", sct_info$source, "; ",
         ifelse(is.na(sct_info$n_genes), "?",
                format(sct_info$n_genes, big.mark = ",")), " genes)")
} else if (length(assays_preproc) > 0) {
  paste0("LogNormalize + ScaleData (assays: ",
         paste(assays_preproc, collapse = ", "), ")")
} else {
  "Unknown (RDS not loaded; check script 03 --no-sct flag)"
}

pca_str <- if (!is.null(pca_stdev) && length(pca_stdev) > 0) {
  sprintf("%d PCs; first PC explains %.1f%% of variance",
          length(pca_stdev), pca_stdev[1]^2 / sum(pca_stdev^2) * 100)
} else { "Not available (preprocessed RDS not loaded)" }

lsi_str <- if (!is.null(lsi_stdev) && length(lsi_stdev) > 1) {
  sprintf("%d LSI dims; dim 1 (excluded from WNN) = %.1f%% variance",
          length(lsi_stdev), lsi_stdev[1]^2 / sum(lsi_stdev^2) * 100)
} else { "Not available (preprocessed RDS not loaded)" }

rna_w_str <- if (!is.null(rna_weights)) {
  sprintf("median = %.2f  (range [%.2f, %.2f])",
          median(rna_weights, na.rm = TRUE),
          min(rna_weights, na.rm = TRUE),
          max(rna_weights, na.rm = TRUE))
} else { "Not available (WNN RDS not loaded)" }

cluster_md_rows <- paste0(
  apply(cluster_sizes, 1, function(r)
    sprintf("| %s | %s | %.1f%% |",
            r["wnn_cluster"],
            format(as.integer(r["n_cells"]), big.mark = ","),
            as.numeric(r["pct"]))),
  collapse = "\n"
)

generated_md <- if (length(generated) > 0) {
  paste0(paste0("- `", generated, "`"), collapse = "\n")
} else { "- (none)" }

skipped_md <- if (length(skipped) > 0) {
  paste0(paste0("- ", skipped), collapse = "\n")
} else { "- (none)" }

summary_text <- paste0(
  "# Integration Diagnostics Summary\n\n",
  "_Generated by `scripts/04b_integration_diagnostics.R`_\n\n",
  "> **Note:** These diagnostics assess QC filtering, dimensional reductions,",
  " and WNN integration structure. They do not constitute biological",
  " validation.\n\n",
  "---\n\n",

  "## Cell counts\n\n",
  "| Stage | Cells |\n",
  "|-------|-------|\n",
  sprintf("| Total barcodes in QC CSV | %s |\n", format(n_qc_total, big.mark = ",")),
  sprintf("| Pass QC | %s (%.1f%%) |\n", format(n_qc_pass, big.mark = ","), 100 * n_qc_pass / n_qc_total),
  sprintf("| Fail QC | %s |\n", format(n_qc_fail, big.mark = ",")),
  sprintf("| Cells in WNN object | %s |\n\n", format(n_cells_wnn, big.mark = ",")),

  "## Normalisation\n\n",
  "| Modality | Method |\n",
  "|----------|--------|\n",
  sprintf("| RNA | %s |\n", norm_mode_str),
  "| ATAC | TF-IDF + SVD (LSI) |\n\n",

  "## Dimensionality reduction\n\n",
  "| Reduction | Summary |\n",
  "|-----------|----------|\n",
  sprintf("| PCA (RNA) | %s |\n", pca_str),
  sprintf("| LSI (ATAC) | %s |\n\n", lsi_str),

  "## WNN integration\n\n",
  "| Metric | Value |\n",
  "|--------|-------|\n",
  sprintf("| WNN clusters | %d |\n", n_clusters),
  sprintf("| RNA modality weight | %s |\n\n", rna_w_str),

  "## Cluster sizes\n\n",
  "| Cluster | Cells | % of total |\n",
  "|---------|-------|------------|\n",
  cluster_md_rows, "\n\n",

  "## Generated figures\n\n",
  generated_md, "\n\n",

  "## Skipped (inputs unavailable)\n\n",
  skipped_md, "\n"
)

writeLines(summary_text, opt$`summary-md`)
msg("  Wrote: ", opt$`summary-md`)

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------

msg(sprintf(
  "\nDone.  %d figure(s) written to:  %s",
  length(generated), opt$outdir
))
msg("Summary markdown: ", opt$`summary-md`)

if (length(skipped) > 0) {
  msg(sprintf("Skipped %d item(s) (inputs unavailable; see summary for details)",
              length(skipped)))
}
