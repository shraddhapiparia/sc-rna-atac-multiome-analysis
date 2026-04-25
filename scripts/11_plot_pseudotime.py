"""
11_plot_pseudotime.py

Phase 5 — Monocyte pseudotime visualisation.

Reads the flat CSV outputs produced by 10_pseudotime_subset.R and the Phase 4
ATAC motif deviation summary, then writes up to four figures:

  Figure 1 (always):
    pseudotime_umap_monocyte.<fmt>
    UMAP coloured by pseudotime (left panel) and broad label (right panel).

  Figure 2 (always — stops if ATAC file absent or columns missing):
    pseudotime_atac_motif_trends.<fmt>
    1 × 4 scatter: pseudotime (x) vs ATAC chromVAR deviation z-score (y) for
    each motif in MOTIF_INFO, with a LOWESS smooth coloured by broad label.

  Figure 3 (conditional on results/dynamics/pseudotime_monocyte_rna_population.csv):
    pseudotime_rna_population_reference.<fmt>
    Bar chart of population-level RNA program scores by broad label.
    NOTE: these are population-level means from the Kang et al. reference dataset
    (different cells from the multiome ATAC cells).  No pseudotime axis is shown
    because the RNA values are uniform within each broad label.

  Figure 4 (only if --rna-cell-scores is supplied and passes validation):
    pseudotime_rna_cell_trends.<fmt>
    1 × 5 scatter: pseudotime (x) vs per-cell RNA program score (y) for each
    of the five TF-activity programs, with a LOWESS smooth.

Three-layer data separation is strictly enforced:
  - Layer 1 (cell-level ATAC): per-cell chromVAR deviation z-scores from script 08.
  - Layer 2 (population-level RNA reference): Kang et al. means from script 06.
  - Layer 3 (per-cell RNA): only used when --rna-cell-scores is provided.

Inputs (default paths, all overridable):
  results/dynamics/pseudotime_monocyte_cells.csv
  results/dynamics/pseudotime_monocyte_rna_population.csv   [optional]
  results/regulons/atac_motif_deviations.csv                [required for Figure 2]
  <--rna-cell-scores>                                       [optional, Figure 4]

Outputs (results/figures/dynamics/):
  pseudotime_umap_monocyte.<fmt>
  pseudotime_atac_motif_trends.<fmt>
  pseudotime_rna_population_reference.<fmt>
  pseudotime_rna_cell_trends.<fmt>

Usage:
  python scripts/11_plot_pseudotime.py
  python scripts/11_plot_pseudotime.py --help
  python scripts/11_plot_pseudotime.py --rna-cell-scores path/to/rna_cell_scores.csv
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PROGRAM_COLS = [
    "interferon_alpha_beta_response",
    "inflammatory_response",
    "antigen_presentation_mhc",
    "cytotoxic_nk_like",
    "monocyte_activation",
]

PROGRAM_LABELS = {
    "interferon_alpha_beta_response": "Interferon \u03b1/\u03b2",
    "inflammatory_response":          "Inflammatory response",
    "antigen_presentation_mhc":       "Antigen presentation",
    "cytotoxic_nk_like":              "Cytotoxic / NK-like",
    "monocyte_activation":            "Monocyte activation",
}

# Ordered list of (motif_id, TF_name) for Figure 2.
MOTIF_INFO = [
    ("MA0107.1", "RELA"),
    ("MA0099.3", "FOS::JUN"),
    ("MA0080.4", "SPI1"),
    ("MA0466.2", "CEBPB"),
]

VALID_FMTS = {"png", "pdf"}

_MUTED = sns.color_palette("muted")
BROAD_LABEL_COLORS = {
    "CD14+ Monocytes":    _MUTED[0],
    "FCGR3A+ Monocytes":  _MUTED[1],
}

# Minimum cells per broad label for Figure 2 ATAC trend warnings.
_MIN_CELLS_PER_LABEL = 10

# Minimum cells per broad label required to draw Figure 4 RNA trends.
# Both broad labels must meet this threshold; if either fails, Figure 4 is skipped entirely.
_MIN_RNA_CELLS_PER_LABEL = 25

# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT   = os.path.dirname(_SCRIPT_DIR)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot Phase 5 monocyte pseudotime figures."
    )
    parser.add_argument(
        "--dynamics-dir",
        default=os.path.join(_REPO_ROOT, "results", "dynamics"),
        help="Directory containing pseudotime CSV outputs from script 10 "
             "[default: %(default)s]",
    )
    parser.add_argument(
        "--regulons-dir",
        default=os.path.join(_REPO_ROOT, "results", "regulons"),
        help="Directory containing ATAC motif CSV outputs from script 08 "
             "[default: %(default)s]",
    )
    parser.add_argument(
        "--outdir",
        default=os.path.join(_REPO_ROOT, "results", "figures", "dynamics"),
        help="Output directory for figures [default: %(default)s]",
    )
    parser.add_argument(
        "--rna-cell-scores",
        default=None,
        metavar="PATH",
        help="Optional CSV with per-cell RNA program scores for multiome cells "
             "(enables Figure 4).  Must contain columns: cell_id plus one column "
             "per PROGRAM_COLS entry.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Figure resolution in DPI [default: %(default)s]",
    )
    parser.add_argument(
        "--fmt",
        default="png",
        choices=sorted(VALID_FMTS),
        help="Output format: png or pdf [default: %(default)s]",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def stop(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def warn(msg: str) -> None:
    print(f"WARNING: {msg}", file=sys.stderr)


def load_csv(path: str, required_cols: list, label: str) -> pd.DataFrame:
    if not os.path.isfile(path):
        stop(f"{label} not found: {path}")
    df = pd.read_csv(path)
    missing = set(required_cols) - set(df.columns)
    if missing:
        stop(f"{label} is missing columns: {sorted(missing)}")
    return df


def save_fig(fig: plt.Figure, path: str, dpi: int) -> None:
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Wrote: {path}")


def _broad_color_series(labels: pd.Series) -> list:
    """Return a list of colours for each row based on broad_label."""
    return [BROAD_LABEL_COLORS.get(l, "#888888") for l in labels]


def _add_broad_legend(ax: plt.Axes) -> None:
    """Add a broad-label colour legend to an axes."""
    handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=c, markersize=7, label=l)
        for l, c in BROAD_LABEL_COLORS.items()
    ]
    ax.legend(handles=handles, title="Broad label", fontsize=7,
              title_fontsize=7, loc="best", framealpha=0.7)


def _pearsonr_np(x: np.ndarray, y: np.ndarray) -> float:
    """Pearson r using numpy only (no scipy dependency)."""
    xm = x - x.mean()
    ym = y - y.mean()
    denom = np.sqrt((xm ** 2).sum() * (ym ** 2).sum())
    if denom == 0.0:
        return float("nan")
    return float(np.dot(xm, ym) / denom)


def _lowess_line(ax: plt.Axes, x: np.ndarray, y: np.ndarray,
                 color: str = "black", lw: float = 1.5) -> None:
    """Overlay a LOWESS smooth on an axes."""
    from statsmodels.nonparametric.smoothers_lowess import lowess
    if len(x) < 4:
        return
    order = np.argsort(x)
    xs, ys = x[order], y[order]
    try:
        smooth = lowess(ys, xs, frac=0.4, return_sorted=True)
        ax.plot(smooth[:, 0], smooth[:, 1], color=color, lw=lw, zorder=5)
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Figure 1 — UMAP coloured by pseudotime and broad label
# ---------------------------------------------------------------------------

def plot_umap(cells: pd.DataFrame, out_path: str, dpi: int) -> None:
    """
    Two-panel UMAP:
      Left:  scatter coloured by scaled pseudotime (viridis); NA cells in grey first.
      Right: scatter coloured by broad label (categorical colours).
    No axis ticks or axis labels on either panel.
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

    x = cells["wnn_umap_1"].values
    y = cells["wnn_umap_2"].values

    # --- Left: pseudotime ---
    pt = cells["pseudotime_scaled"].values
    na_mask  = np.isnan(pt)
    val_mask = ~na_mask

    # NA cells plotted first in light grey so valid cells overlay them
    if na_mask.any():
        axes[0].scatter(x[na_mask], y[na_mask],
                        c="#cccccc", s=6, linewidths=0, alpha=0.6, zorder=2)

    sc = axes[0].scatter(x[val_mask], y[val_mask],
                         c=pt[val_mask], cmap="viridis",
                         s=6, linewidths=0, alpha=0.8, zorder=3)
    cb = fig.colorbar(sc, ax=axes[0], fraction=0.046, pad=0.04)
    cb.set_label("Pseudotime (scaled)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    axes[0].set_title("Monocyte pseudotime (WNN UMAP)", fontsize=10)
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    axes[0].set_xlabel("")
    axes[0].set_ylabel("")
    sns.despine(ax=axes[0], left=True, bottom=True)

    # --- Right: broad label ---
    colors = _broad_color_series(cells["broad_label"])
    axes[1].scatter(x, y, c=colors, s=6, linewidths=0, alpha=0.8)
    _add_broad_legend(axes[1])
    axes[1].set_title("Cell type", fontsize=10)
    axes[1].set_xticks([])
    axes[1].set_yticks([])
    axes[1].set_xlabel("")
    axes[1].set_ylabel("")
    sns.despine(ax=axes[1], left=True, bottom=True)

    fig.tight_layout()
    save_fig(fig, out_path, dpi)


# ---------------------------------------------------------------------------
# Figure 2 — ATAC motif deviation trends along pseudotime
# ---------------------------------------------------------------------------

def plot_atac_trends(cells: pd.DataFrame, atac_wide: pd.DataFrame,
                     out_path: str, dpi: int) -> None:
    """
    1 × 4 panel — one sub-panel per motif in MOTIF_INFO.
    x: pseudotime_scaled, y: chromVAR deviation z-score.
    Points coloured by broad_label; one LOWESS smooth per broad_label.

    atac_wide must contain: cell_id, MA0107.1, MA0099.3, MA0080.4, MA0466.2.
    Join is inner on cell_id; NA-pseudotime rows are dropped after join.
    Stops if join yields 0 rows or fewer than 10 rows remain after NA removal.
    Warns if fewer than 10 rows per broad label.
    Shared legend placed below all panels.
    """
    # Join pseudotime + broad_label onto wide ATAC
    merged = cells[["cell_id", "pseudotime_scaled", "broad_label"]].merge(
        atac_wide, on="cell_id", how="inner",
    )

    if merged.empty:
        stop("Figure 2: join of pseudotime and ATAC deviations yielded 0 rows. "
             "Check that cell_id values match between "
             "pseudotime_monocyte_cells.csv and atac_motif_deviations.csv.")

    # Drop NA-pseudotime rows (Slingshot outliers)
    merged = merged.dropna(subset=["pseudotime_scaled"]).reset_index(drop=True)

    if len(merged) < 10:
        stop(f"Figure 2: only {len(merged)} rows remain after removing "
             "NA-pseudotime cells. Cannot produce reliable trend plots.")

    # Warn about sparse broad labels
    for broad, grp in merged.groupby("broad_label"):
        if len(grp) < _MIN_CELLS_PER_LABEL:
            warn(f"Figure 2: broad label '{broad}' has only {len(grp)} cells "
                 f"after join (< {_MIN_CELLS_PER_LABEL}).")

    n_motifs = len(MOTIF_INFO)
    fig, axes = plt.subplots(1, n_motifs, figsize=(14, 3.8))

    legend_handles = {}   # collect one handle per broad label for shared legend

    for ax, (motif_id, tf_name) in zip(axes, MOTIF_INFO):
        if motif_id not in merged.columns:
            ax.text(0.5, 0.5, f"{motif_id}\n{tf_name}\n(column absent)",
                    ha="center", va="center", transform=ax.transAxes, fontsize=8)
            ax.set_title(f"{tf_name}\n({motif_id})", fontsize=9)
            sns.despine(ax=ax)
            continue

        for broad, grp in merged.groupby("broad_label"):
            c = BROAD_LABEL_COLORS.get(broad, "#888888")
            sc = ax.scatter(grp["pseudotime_scaled"], grp[motif_id],
                            c=[c], s=5, alpha=0.5, linewidths=0)
            if broad not in legend_handles:
                legend_handles[broad] = plt.Line2D(
                    [0], [0], marker="o", color="w",
                    markerfacecolor=c, markersize=7, label=broad,
                )

        for broad, grp in merged.groupby("broad_label"):
            c = BROAD_LABEL_COLORS.get(broad, "#888888")
            _lowess_line(ax, grp["pseudotime_scaled"].values,
                         grp[motif_id].values, color=c, lw=1.8)

        ax.axhline(0, color="black", lw=0.6, linestyle="--", alpha=0.5)
        ax.set_title(f"{tf_name}\n({motif_id})", fontsize=9)
        ax.set_xlabel("Pseudotime (scaled)", fontsize=8)
        ax.set_ylabel("chromVAR z-score", fontsize=8)
        ax.tick_params(labelsize=7)
        sns.despine(ax=ax)

    fig.suptitle(
        "ATAC motif deviation along monocyte pseudotime (cell-level)",
        fontsize=10, y=1.01,
    )
    fig.tight_layout()

    # Shared legend below all panels
    if legend_handles:
        fig.legend(
            handles=list(legend_handles.values()),
            title="Broad label",
            fontsize=7,
            title_fontsize=7,
            loc="lower center",
            ncol=len(legend_handles),
            bbox_to_anchor=(0.5, -0.08),
            framealpha=0.7,
        )

    save_fig(fig, out_path, dpi)


# ---------------------------------------------------------------------------
# Figure 3 — RNA population reference (NO pseudotime axis)
# ---------------------------------------------------------------------------

def plot_rna_population_reference(rna_pop: pd.DataFrame,
                                   out_path: str, dpi: int) -> None:
    """
    Bar chart of population-level RNA program scores by broad label.

    rna_pop has columns: cell_id, broad_label, <program>_rna_mean,
    <program>_rna_std (and optionally pseudotime columns — ignored here).

    One grouped bar per program, two bars per group (one per broad label).
    No pseudotime axis — these values are uniform within each broad label.
    """
    # Collapse to one row per broad_label (values are uniform within label)
    mean_cols = [f"{p}_rna_mean" for p in PROGRAM_COLS]
    std_cols  = [f"{p}_rna_std"  for p in PROGRAM_COLS]

    available_mean = [c for c in mean_cols if c in rna_pop.columns]
    if not available_mean:
        warn("RNA population file has no <program>_rna_mean columns; skipping Figure 3.")
        return

    ref = (
        rna_pop[["broad_label"] + available_mean + std_cols]
        .drop_duplicates(subset="broad_label")
        .set_index("broad_label")
    )

    broad_labels = list(BROAD_LABEL_COLORS.keys())
    ref = ref.reindex([l for l in broad_labels if l in ref.index])

    if ref.empty:
        warn("No recognised broad labels in RNA population file; skipping Figure 3.")
        return

    # Build bar chart
    programs_present = [p for p in PROGRAM_COLS if f"{p}_rna_mean" in ref.columns]
    n_prog = len(programs_present)
    n_broad = len(ref)

    bar_width = 0.35
    x_pos = np.arange(n_prog)

    fig, ax = plt.subplots(figsize=(8, 4.5))

    for i, broad in enumerate(ref.index):
        means = [ref.loc[broad, f"{p}_rna_mean"] for p in programs_present]
        stds  = [ref.loc[broad, f"{p}_rna_std"]
                 if f"{p}_rna_std" in ref.columns else 0.0
                 for p in programs_present]
        offset = (i - (n_broad - 1) / 2) * bar_width
        color  = BROAD_LABEL_COLORS.get(broad, "#888888")
        ax.bar(
            x_pos + offset, means, bar_width,
            yerr=stds, capsize=4,
            color=color, alpha=0.85, label=broad,
            error_kw={"elinewidth": 1.0, "ecolor": "black"},
        )

    ax.set_xticks(x_pos)
    ax.set_xticklabels(
        [PROGRAM_LABELS.get(p, p) for p in programs_present],
        rotation=30, ha="right", fontsize=8,
    )
    ax.set_ylabel("Mean RNA program score (Phase 3 baseline)", fontsize=8)
    ax.set_title(
        "RNA baseline program scores \u2014 monocyte populations\n"
        "(population reference: Kang et al. dataset;"
        " values are uniform within each broad label)",
        fontsize=9,
    )
    ax.axhline(0, color="black", lw=0.6, linestyle="--", alpha=0.5)
    ax.legend(title="Broad label", fontsize=8, title_fontsize=8)
    ax.tick_params(axis="y", labelsize=7)
    sns.despine(ax=ax)

    fig.tight_layout()
    save_fig(fig, out_path, dpi)


# ---------------------------------------------------------------------------
# Figure 4 — per-cell RNA program trends (optional)
# ---------------------------------------------------------------------------

def plot_rna_cell_trends(cells: pd.DataFrame, rna_cell: pd.DataFrame,
                          out_path: str, dpi: int) -> None:
    """
    1 × 5 panel — one sub-panel per program in PROGRAM_COLS.
    x: pseudotime_scaled, y: per-cell RNA program score.
    Points coloured by broad_label; one LOWESS smooth per broad_label.

    rna_cell must contain: cell_id plus one column per PROGRAM_COLS entry.

    Both monocyte broad labels must have >= _MIN_RNA_CELLS_PER_LABEL cells after
    the join; if either fails, Figure 4 is skipped entirely (no partial plot).
    """
    # Join pseudotime onto per-cell RNA scores
    merged = cells[["cell_id", "pseudotime_scaled", "broad_label"]].merge(
        rna_cell, on="cell_id", how="inner",
    )

    if len(merged) < 10:
        warn(
            f"Figure 4: only {len(merged)} total cells after join "
            "(need >= 10); skipping Figure 4."
        )
        return

    # Both broad labels must meet _MIN_RNA_CELLS_PER_LABEL; skip entirely if either fails
    counts  = merged["broad_label"].value_counts()
    failing = [l for l in BROAD_LABEL_COLORS
               if counts.get(l, 0) < _MIN_RNA_CELLS_PER_LABEL]
    if failing:
        warn(
            f"Figure 4: broad label(s) {failing} have fewer than "
            f"{_MIN_RNA_CELLS_PER_LABEL} cells after join; "
            "skipping Figure 4 entirely."
        )
        return

    # --- Console reporting ---
    print(f"  Figure 4 joined cells: {len(merged):,}")
    label_counts = {l: int(counts.get(l, 0)) for l in BROAD_LABEL_COLORS}
    print(f"  Figure 4 cells by broad label: {label_counts}")
    for prog in PROGRAM_COLS:
        if prog not in merged.columns:
            continue
        for broad, grp in merged.groupby("broad_label"):
            sub = grp[["pseudotime_scaled", prog]].dropna()
            if len(sub) >= 3:
                r = _pearsonr_np(sub["pseudotime_scaled"].values, sub[prog].values)
                print(f"  {prog} | {broad} | r={r:.2f} | n={len(sub)}")
            else:
                print(f"  {prog} | {broad} | n<3")

    # --- Plot ---
    n_prog = len(PROGRAM_COLS)
    fig, axes = plt.subplots(1, n_prog, figsize=(17, 3.8))

    legend_handles = {}

    for ax, prog in zip(axes, PROGRAM_COLS):
        if prog not in merged.columns:
            ax.text(0.5, 0.5, f"{PROGRAM_LABELS.get(prog, prog)}\n(column absent)",
                    ha="center", va="center", transform=ax.transAxes, fontsize=8)
            ax.set_title(PROGRAM_LABELS.get(prog, prog), fontsize=9)
            sns.despine(ax=ax)
            continue

        for broad, grp in merged.groupby("broad_label"):
            c = BROAD_LABEL_COLORS.get(broad, "#888888")
            ax.scatter(grp["pseudotime_scaled"], grp[prog],
                       c=[c], s=5, alpha=0.5, linewidths=0)
            if broad not in legend_handles:
                legend_handles[broad] = plt.Line2D(
                    [0], [0], marker="o", color="w",
                    markerfacecolor=c, markersize=7, label=broad,
                )

        for broad, grp in merged.groupby("broad_label"):
            c = BROAD_LABEL_COLORS.get(broad, "#888888")
            _lowess_line(ax, grp["pseudotime_scaled"].values,
                         grp[prog].values, color=c, lw=1.8)

        ax.axhline(0, color="black", lw=0.6, linestyle="--", alpha=0.5)
        ax.set_title(PROGRAM_LABELS.get(prog, prog), fontsize=9)
        ax.set_xlabel("Pseudotime (scaled)", fontsize=8)
        ax.set_ylabel("RNA program score", fontsize=8)
        ax.tick_params(labelsize=7)
        sns.despine(ax=ax)

    fig.suptitle(
        "RNA program score along monocyte pseudotime (cell-level multiome data)",
        fontsize=10, y=1.01,
    )
    fig.tight_layout()

    # Shared legend below all panels — same style as Figure 2
    if legend_handles:
        fig.legend(
            handles=list(legend_handles.values()),
            title="Broad label",
            fontsize=7,
            title_fontsize=7,
            loc="lower center",
            ncol=len(legend_handles),
            bbox_to_anchor=(0.5, -0.08),
            framealpha=0.7,
        )

    save_fig(fig, out_path, dpi)


# ---------------------------------------------------------------------------
# ATAC cell-level deviation loader
# ---------------------------------------------------------------------------

_ATAC_MOTIF_COLS = [mid for mid, _ in MOTIF_INFO]


def load_atac_cell_deviations(regulons_dir: str) -> pd.DataFrame:
    """
    Load per-cell ATAC motif deviation z-scores (wide format).

    Expects atac_motif_deviations.csv with columns:
      cell_id, MA0107.1, MA0099.3, MA0080.4, MA0466.2

    Stops (sys.exit) if the file is absent or any required column is missing.
    """
    path = os.path.join(regulons_dir, "atac_motif_deviations.csv")
    if not os.path.isfile(path):
        stop(
            f"ATAC deviation file not found: {path}\n"
            "Run scripts/08_motif_enrichment_atac.R first."
        )
    df = pd.read_csv(path)
    required = {"cell_id"} | set(_ATAC_MOTIF_COLS)
    missing  = required - set(df.columns)
    if missing:
        stop(
            f"atac_motif_deviations.csv is missing required columns: "
            f"{sorted(missing)}\n"
            "Expected: cell_id, " + ", ".join(_ATAC_MOTIF_COLS)
        )
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    sns.set_theme(style="white")
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # -----------------------------------------------------------------------
    # Load required input: per-cell pseudotime
    # -----------------------------------------------------------------------

    cells_path = os.path.join(args.dynamics_dir, "pseudotime_monocyte_cells.csv")
    print("Loading pseudotime_monocyte_cells.csv...")
    cells = load_csv(
        cells_path,
        required_cols=["cell_id", "cell_type_wnn", "broad_label",
                        "wnn_umap_1", "wnn_umap_2",
                        "pseudotime", "pseudotime_scaled"],
        label="pseudotime_monocyte_cells.csv",
    )
    print(f"  {len(cells):,} cells loaded.")

    if cells.empty:
        stop("pseudotime_monocyte_cells.csv is empty; cannot produce any figures.")

    # -----------------------------------------------------------------------
    # Figure 1 — UMAP
    # -----------------------------------------------------------------------

    print("Plotting Figure 1: UMAP...")
    plot_umap(
        cells=cells,
        out_path=os.path.join(args.outdir, f"pseudotime_umap_monocyte.{args.fmt}"),
        dpi=args.dpi,
    )

    # -----------------------------------------------------------------------
    # Figure 2 — ATAC motif trends (required)
    # -----------------------------------------------------------------------

    print("Loading ATAC cell-level deviations for Figure 2...")
    atac_wide = load_atac_cell_deviations(args.regulons_dir)
    print(f"  {len(atac_wide):,} rows loaded.")
    print("Plotting Figure 2: ATAC motif trends...")
    plot_atac_trends(
        cells=cells,
        atac_wide=atac_wide,
        out_path=os.path.join(args.outdir, f"pseudotime_atac_motif_trends.{args.fmt}"),
        dpi=args.dpi,
    )

    # -----------------------------------------------------------------------
    # Figure 3 — RNA population reference
    # -----------------------------------------------------------------------

    rna_pop_path = os.path.join(args.dynamics_dir,
                                "pseudotime_monocyte_rna_population.csv")
    if os.path.isfile(rna_pop_path):
        print("Loading pseudotime_monocyte_rna_population.csv for Figure 3...")
        rna_pop = pd.read_csv(rna_pop_path)
        print(f"  {len(rna_pop):,} rows loaded.")
        print("Plotting Figure 3: RNA population reference...")
        plot_rna_population_reference(
            rna_pop=rna_pop,
            out_path=os.path.join(args.outdir,
                                  f"pseudotime_rna_population_reference.{args.fmt}"),
            dpi=args.dpi,
        )
    else:
        print(
            "  pseudotime_monocyte_rna_population.csv not found; "
            "skipping Figure 3.\n"
            "  Run 10_pseudotime_subset.R with the RNA program scores file present."
        )

    # -----------------------------------------------------------------------
    # Figure 4 — per-cell RNA trends (optional)
    # -----------------------------------------------------------------------

    if args.rna_cell_scores is not None:
        print(f"Loading per-cell RNA scores from: {args.rna_cell_scores}")
        if not os.path.isfile(args.rna_cell_scores):
            warn(f"--rna-cell-scores file not found: {args.rna_cell_scores}; "
                 "skipping Figure 4.")
        else:
            rna_cell = pd.read_csv(args.rna_cell_scores)
            missing_prog = [p for p in PROGRAM_COLS if p not in rna_cell.columns]
            if "cell_id" not in rna_cell.columns:
                warn("--rna-cell-scores file has no 'cell_id' column; "
                     "skipping Figure 4.")
            elif missing_prog:
                warn(
                    f"--rna-cell-scores file is missing program columns: "
                    f"{missing_prog}; skipping Figure 4."
                )
            else:
                print(f"  {len(rna_cell):,} rows loaded.")
                print("Plotting Figure 4: per-cell RNA trends...")
                plot_rna_cell_trends(
                    cells=cells,
                    rna_cell=rna_cell,
                    out_path=os.path.join(args.outdir,
                                          f"pseudotime_rna_cell_trends.{args.fmt}"),
                    dpi=args.dpi,
                )
    else:
        print(
            "  --rna-cell-scores not provided; Figure 4 skipped.\n"
            "  Supply per-cell RNA program scores for multiome cells to enable it."
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
