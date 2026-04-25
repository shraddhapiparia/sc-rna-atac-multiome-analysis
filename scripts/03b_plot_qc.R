#!/usr/bin/env Rscript
#
# 03b_plot_qc.R
#
# Publication-quality QC visualisation for the 10x PBMC Multiome dataset.
# Reads preprocessing_qc_stats.csv written by script 03 and produces four
# figures covering RNA and ATAC QC metrics.
#
# This script is standalone — it reads a flat CSV and requires only
# optparse, ggplot2, and patchwork (no Seurat/Signac object needed).
#
# Reads:
#   results/integration/preprocessing_qc_stats.csv
#
# Writes (results/figures/qc/):
#   rna_qc_violin.{fmt}   nFeature_RNA / nCount_RNA / percent_mt  (3-panel)
#   rna_qc_scatter.{fmt}  nFeature_RNA vs percent_mt              (scatter)
#   atac_qc_violin.{fmt}  nCount_ATAC / TSS_enrichment / nucleosome_signal (3-panel)
#   atac_qc_scatter.{fmt} TSS_enrichment vs nucleosome_signal     (scatter)
#
# Usage:
#   Rscript scripts/03b_plot_qc.R
#   Rscript scripts/03b_plot_qc.R --dpi 150 --fmt pdf
#   Rscript scripts/03b_plot_qc.R --help

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(patchwork)   # install.packages("patchwork") if absent
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
# Arguments — thresholds default to script 03 values
# ---------------------------------------------------------------------------

option_list <- list(

  make_option("--qc-csv",
    default = file.path(.repo_root, "results", "integration",
                        "preprocessing_qc_stats.csv"),
    type    = "character",
    help    = "Path to preprocessing_qc_stats.csv from script 03 [default: %default]"),

  make_option("--outdir",
    default = file.path(.repo_root, "results", "figures", "qc"),
    type    = "character",
    help    = "Output directory for figures [default: %default]"),

  make_option("--dpi",
    default = 300L,
    type    = "integer",
    help    = "Figure resolution in DPI [default: %default]"),

  make_option("--fmt",
    default = "png",
    type    = "character",
    help    = "Output format: 'png' or 'pdf' [default: %default]"),

  # ---- RNA thresholds (mirror script 03 defaults) ----
  make_option("--min-features-rna",
    type = "integer", default = 200L,
    help = "nFeature_RNA lower bound [default: %default]"),

  make_option("--max-features-rna",
    type = "integer", default = 7500L,
    help = "nFeature_RNA upper bound [default: %default]"),

  make_option("--max-pct-mt",
    type = "double", default = 20.0,
    help = "percent_mt upper bound [default: %default]"),

  make_option("--min-count-rna",
    type = "integer", default = 0L,
    help = "nCount_RNA lower bound; 0 = no bound drawn [default: %default]"),

  # ---- ATAC thresholds ----
  make_option("--min-count-atac",
    type = "integer", default = 1000L,
    help = "nCount_ATAC lower bound [default: %default]"),

  make_option("--max-count-atac",
    type = "integer", default = 100000L,
    help = "nCount_ATAC upper bound [default: %default]"),

  make_option("--min-tss",
    type = "double", default = 2.0,
    help = "TSS enrichment lower bound [default: %default]"),

  make_option("--max-nucleosome",
    type = "double", default = 2.0,
    help = "nucleosome_signal upper bound [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!opt$fmt %in% c("png", "pdf")) {
  stop("--fmt must be 'png' or 'pdf'.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Input validation + output directory
# ---------------------------------------------------------------------------

if (!file.exists(opt$`qc-csv`)) {
  stop(
    "QC CSV not found: ", opt$`qc-csv`,
    "\nRun scripts/03_preprocess_multiome.R first.",
    call. = FALSE
  )
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

message("Loading QC stats: ", opt$`qc-csv`)
qc <- read.csv(opt$`qc-csv`, stringsAsFactors = FALSE)
qc$pass_qc <- as.logical(qc$pass_qc)

n_total <- nrow(qc)
n_pass  <- sum(qc$pass_qc, na.rm = TRUE)
n_fail  <- n_total - n_pass

message(sprintf(
  "  %s cells total  |  %s pass (%.1f%%)  |  %s fail",
  format(n_total, big.mark = ","),
  format(n_pass,  big.mark = ","),
  100 * n_pass / n_total,
  format(n_fail,  big.mark = ",")
))

# Ordered factor: Fail left, Pass right
qc$qc_group <- factor(
  ifelse(qc$pass_qc, "Pass", "Fail"),
  levels = c("Fail", "Pass")
)

# ---------------------------------------------------------------------------
# Design constants
# ---------------------------------------------------------------------------

COL_PASS   <- "#2166AC"   # blue
COL_FAIL   <- "#D73027"   # red
COL_THRESH <- "#E6550D"   # orange — threshold lines and captions
COL_BOX    <- "grey25"    # boxplot outline

# Subtitle reused across all four figures
SUBTITLE <- sprintf(
  "%s cells total  ·  %s pass (%.1f%%)  ·  %s fail",
  format(n_total, big.mark = ","),
  format(n_pass,  big.mark = ","),
  100 * n_pass / n_total,
  format(n_fail,  big.mark = ",")
)

# ---------------------------------------------------------------------------
# Base theme
# ---------------------------------------------------------------------------

theme_qc <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title             = element_text(size = 10),
      axis.text              = element_text(size = 9),
      panel.grid.major.y     = element_line(color = "grey92", linewidth = 0.3),
      legend.title           = element_text(size = 9, face = "bold"),
      legend.text            = element_text(size = 8),
      legend.key.size        = unit(0.55, "lines"),
      plot.title             = element_text(face = "bold", size = 11),
      plot.caption           = element_text(size = 7.5, color = COL_THRESH,
                                            hjust = 0),
      plot.caption.position  = "plot"
    )
}

# Shared fill and colour scales for pass/fail
fill_scale <- scale_fill_manual(
  name   = "QC status",
  values = c("Pass" = COL_PASS, "Fail" = COL_FAIL)
)

colour_scale <- scale_colour_manual(
  name   = "QC status",
  values = c("Pass" = COL_PASS, "Fail" = COL_FAIL)
)

# ---------------------------------------------------------------------------
# Helper: build threshold caption string
# ---------------------------------------------------------------------------

make_caption <- function(...) {
  # Each argument is a named character pair: c(label = value)
  # Returns NULL if no thresholds are provided
  parts <- list(...)
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L) return(NULL)
  paste0("Threshold: ", paste(unlist(parts), collapse = "  ·  "))
}

# ---------------------------------------------------------------------------
# Helper: violin panel for one QC metric
#
# Renders two violins side-by-side (Fail / Pass), a slim boxplot overlay,
# and dashed threshold lines. Threshold values appear in the panel caption.
# ---------------------------------------------------------------------------

violin_panel <- function(df, yvar, ylab, title,
                         thresh_min = NULL,
                         thresh_max = NULL,
                         log_y      = FALSE) {

  # Build threshold lines and caption text
  hlines  <- list()
  cap_parts <- list()

  if (!is.null(thresh_min) && thresh_min > 0) {
    hlines    <- c(hlines, list(
      geom_hline(yintercept = thresh_min, linetype = "dashed",
                 color = COL_THRESH, linewidth = 0.65)
    ))
    cap_parts <- c(cap_parts,
                   list(paste0("min = ", prettyNum(thresh_min, big.mark = ","))))
  }

  if (!is.null(thresh_max) && is.finite(thresh_max)) {
    hlines    <- c(hlines, list(
      geom_hline(yintercept = thresh_max, linetype = "dashed",
                 color = COL_THRESH, linewidth = 0.65)
    ))
    cap_parts <- c(cap_parts,
                   list(paste0("max = ", prettyNum(thresh_max, big.mark = ","))))
  }

  caption_str <- if (length(cap_parts) > 0L)
    make_caption(unlist(cap_parts))
  else
    NULL

  p <- ggplot(df, aes(x = qc_group, y = .data[[yvar]], fill = qc_group)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.75, color = NA) +
    geom_boxplot(
      width         = 0.15,
      outlier.shape = NA,
      fill          = "white",
      color         = COL_BOX,
      linewidth     = 0.35
    ) +
    hlines +
    fill_scale +
    labs(title = title, x = NULL, y = ylab, caption = caption_str) +
    theme_qc() +
    theme(legend.position = "none")

  if (log_y) {
    p <- p + scale_y_log10(
      labels = function(x) prettyNum(x, big.mark = ",")
    )
  }

  p
}

# ---------------------------------------------------------------------------
# Helper: save a figure (ggplot or patchwork)
# ---------------------------------------------------------------------------

save_fig <- function(plot, stem, width, height) {
  path <- file.path(opt$outdir, paste0(stem, ".", opt$fmt))
  ggsave(
    filename = path,
    plot     = plot,
    width    = width,
    height   = height,
    dpi      = opt$dpi,
    units    = "in"
  )
  message("  Wrote: ", path)
}

# ---------------------------------------------------------------------------
# Figure 1 — RNA QC violin (3-panel)
#
# nFeature_RNA  |  nCount_RNA (log10)  |  percent_mt
# ---------------------------------------------------------------------------

message("Building Figure 1: RNA QC violins...")

p_feat <- violin_panel(
  df         = qc,
  yvar       = "nFeature_RNA",
  ylab       = "Genes detected",
  title      = "nFeature_RNA",
  thresh_min = opt$`min-features-rna`,
  thresh_max = opt$`max-features-rna`
)

# nCount_RNA: lower bound is 0 (no filter) and upper bound is .Machine$integer.max
# Neither threshold is biologically meaningful to draw, so none are annotated.
p_count <- violin_panel(
  df         = qc,
  yvar       = "nCount_RNA",
  ylab       = "UMI count (log₁₀ scale)",
  title      = "nCount_RNA",
  thresh_min = if (opt$`min-count-rna` > 0L) opt$`min-count-rna` else NULL,
  thresh_max = NULL,
  log_y      = TRUE
)

p_mt <- violin_panel(
  df         = qc,
  yvar       = "percent_mt",
  ylab       = "Mitochondrial reads (%)",
  title      = "percent_mt",
  thresh_min = NULL,
  thresh_max = opt$`max-pct-mt`
)

fig_rna_violin <- (p_feat | p_count | p_mt) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "RNA QC Metrics — 10x PBMC Multiome",
    subtitle = SUBTITLE,
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "grey45"),
      plot.margin   = margin(8, 8, 4, 8)
    )
  ) &
  theme(legend.position = "bottom")

save_fig(fig_rna_violin, "rna_qc_violin", width = 9.5, height = 5.5)

# ---------------------------------------------------------------------------
# Figure 2 — RNA QC scatter: nFeature_RNA vs percent_mt
#
# Points coloured by pass/fail. Three threshold lines annotated with their
# values using y = Inf / x = thresh positioning (no clip = "off" needed).
# ---------------------------------------------------------------------------

message("Building Figure 2: RNA QC scatter...")

fig_rna_scatter <- ggplot(
    qc,
    aes(x = nFeature_RNA, y = percent_mt, colour = qc_group)
  ) +
  geom_point(size = 0.55, alpha = 0.35, stroke = 0) +
  # Threshold lines
  geom_vline(xintercept = opt$`min-features-rna`,
             linetype = "dashed", color = COL_THRESH, linewidth = 0.65) +
  geom_vline(xintercept = opt$`max-features-rna`,
             linetype = "dashed", color = COL_THRESH, linewidth = 0.65) +
  geom_hline(yintercept = opt$`max-pct-mt`,
             linetype = "dashed", color = COL_THRESH, linewidth = 0.65) +
  # Threshold annotations — positioned at panel edge using Inf / vjust trick
  annotate("text",
    x = opt$`min-features-rna`, y = Inf,
    label = paste0("min = ", opt$`min-features-rna`),
    hjust = -0.08, vjust = 1.5,
    size = 2.8, color = COL_THRESH
  ) +
  annotate("text",
    x = opt$`max-features-rna`, y = Inf,
    label = paste0("max = ", format(opt$`max-features-rna`, big.mark = ",")),
    hjust = -0.08, vjust = 1.5,
    size = 2.8, color = COL_THRESH
  ) +
  annotate("text",
    x = Inf, y = opt$`max-pct-mt`,
    label = paste0("max = ", opt$`max-pct-mt`, "%"),
    hjust = 1.03, vjust = -0.45,
    size = 2.8, color = COL_THRESH
  ) +
  colour_scale +
  guides(colour = guide_legend(
    override.aes = list(size = 2.5, alpha = 1),
    title        = "QC status"
  )) +
  labs(
    title    = "RNA QC — nFeature_RNA vs percent_mt",
    subtitle = SUBTITLE,
    x        = "Genes detected (nFeature_RNA)",
    y        = "Mitochondrial reads (%)"
  ) +
  theme_qc()

save_fig(fig_rna_scatter, "rna_qc_scatter", width = 6.5, height = 5)

# ---------------------------------------------------------------------------
# Figure 3 — ATAC QC violin (3-panel)
#
# nCount_ATAC (log10)  |  TSS_enrichment  |  nucleosome_signal
# ---------------------------------------------------------------------------

message("Building Figure 3: ATAC QC violins...")

p_atac_count <- violin_panel(
  df         = qc,
  yvar       = "nCount_ATAC",
  ylab       = "Peak count (log₁₀ scale)",
  title      = "nCount_ATAC",
  thresh_min = opt$`min-count-atac`,
  thresh_max = opt$`max-count-atac`,
  log_y      = TRUE
)

p_tss <- violin_panel(
  df         = qc,
  yvar       = "TSS_enrichment",
  ylab       = "TSS enrichment score",
  title      = "TSS enrichment",
  thresh_min = opt$`min-tss`,
  thresh_max = NULL
)

p_nuc <- violin_panel(
  df         = qc,
  yvar       = "nucleosome_signal",
  ylab       = "Nucleosome signal",
  title      = "nucleosome_signal",
  thresh_min = NULL,
  thresh_max = opt$`max-nucleosome`
)

fig_atac_violin <- (p_atac_count | p_tss | p_nuc) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "ATAC QC Metrics — 10x PBMC Multiome",
    subtitle = SUBTITLE,
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "grey45"),
      plot.margin   = margin(8, 8, 4, 8)
    )
  ) &
  theme(legend.position = "bottom")

save_fig(fig_atac_violin, "atac_qc_violin", width = 9.5, height = 5.5)

# ---------------------------------------------------------------------------
# Figure 4 — ATAC QC scatter: TSS_enrichment vs nucleosome_signal
#
# The two threshold lines define a pass region (top-left quadrant excluded).
# ---------------------------------------------------------------------------

message("Building Figure 4: ATAC QC scatter...")

fig_atac_scatter <- ggplot(
    qc,
    aes(x = TSS_enrichment, y = nucleosome_signal, colour = qc_group)
  ) +
  geom_point(size = 0.55, alpha = 0.35, stroke = 0) +
  # Threshold lines
  geom_vline(xintercept = opt$`min-tss`,
             linetype = "dashed", color = COL_THRESH, linewidth = 0.65) +
  geom_hline(yintercept = opt$`max-nucleosome`,
             linetype = "dashed", color = COL_THRESH, linewidth = 0.65) +
  # Threshold annotations
  annotate("text",
    x = opt$`min-tss`, y = Inf,
    label = paste0("min = ", opt$`min-tss`),
    hjust = -0.08, vjust = 1.5,
    size = 2.8, color = COL_THRESH
  ) +
  annotate("text",
    x = Inf, y = opt$`max-nucleosome`,
    label = paste0("max = ", opt$`max-nucleosome`),
    hjust = 1.03, vjust = -0.45,
    size = 2.8, color = COL_THRESH
  ) +
  # Pass region shading (light blue, drawn under points)
  annotate("rect",
    xmin = opt$`min-tss`, xmax = Inf,
    ymin = -Inf,          ymax = opt$`max-nucleosome`,
    fill = COL_PASS, alpha = 0.04
  ) +
  colour_scale +
  guides(colour = guide_legend(
    override.aes = list(size = 2.5, alpha = 1),
    title        = "QC status"
  )) +
  labs(
    title    = "ATAC QC — TSS enrichment vs nucleosome signal",
    subtitle = SUBTITLE,
    x        = "TSS enrichment score",
    y        = "Nucleosome signal"
  ) +
  theme_qc()

save_fig(fig_atac_scatter, "atac_qc_scatter", width = 6.5, height = 5)

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------

message(sprintf("\nDone. %d figures written to: %s", 4L, opt$outdir))
