"""
07_plot_tf_activity.py

Phase 3 — TF activity proxy plots.

Reads the flat summary CSVs produced by 06_tf_activity.py and writes:
  - two z-scored heatmaps (cell-type and broad-label)
  - five per-program bar charts with error bars

Inputs  (results/regulons/):
  tf_activity_rna_by_celltype.csv
  tf_activity_rna_by_broad.csv

Outputs (results/figures/tf_activity/):
  tf_heatmap_rna_celltype.<fmt>
  tf_heatmap_rna_broad.<fmt>
  tf_barplot_interferon_alpha_beta_response.<fmt>
  tf_barplot_inflammatory_response.<fmt>
  tf_barplot_antigen_presentation_mhc.<fmt>
  tf_barplot_cytotoxic_nk_like.<fmt>
  tf_barplot_monocyte_activation.<fmt>

Usage:
  python scripts/07_plot_tf_activity.py
  python scripts/07_plot_tf_activity.py --help
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")           # non-interactive backend; must precede pyplot import
import matplotlib.pyplot as plt
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

# Shorter axis labels for figures
PROGRAM_LABELS = {
    "interferon_alpha_beta_response": "IFN α/β\n(IRF3/7·STAT1/2)",
    "inflammatory_response":          "Inflammation\n(NF-κB/AP-1)",
    "antigen_presentation_mhc":       "Antigen pres.\n(CIITA/IRF1)",
    "cytotoxic_nk_like":              "Cytotoxic\n(EOMES/TBX21)",
    "monocyte_activation":            "Mono. activ.\n(PU.1/CEBPB)",
}

VALID_FMTS = {"png", "pdf"}

# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT  = os.path.dirname(_SCRIPT_DIR)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot TF activity proxy summaries from Phase 3 CSV outputs."
    )
    parser.add_argument(
        "--results-dir",
        default=os.path.join(_REPO_ROOT, "results", "regulons"),
        help="Directory containing tf_activity_*.csv files [default: %(default)s]",
    )
    parser.add_argument(
        "--outdir",
        default=os.path.join(_REPO_ROOT, "results", "figures", "tf_activity"),
        help="Output directory for figures [default: %(default)s]",
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


def load_csv(path: str, required_cols: list) -> pd.DataFrame:
    if not os.path.isfile(path):
        stop(f"Input file not found: {path}\nRun 06_tf_activity.py first.")
    df = pd.read_csv(path)
    missing = set(required_cols) - set(df.columns)
    if missing:
        stop(f"{path} is missing columns: {sorted(missing)}")
    return df


def pivot_to_matrix(long_df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    """
    Pivot long-form (group, program, mean) to a matrix with groups as rows
    and programs as columns.  Only the 'mean' values are used.
    """
    return (
        long_df[long_df["program"].isin(PROGRAM_COLS)]
        .pivot_table(index=group_col, columns="program", values="mean")
        .reindex(columns=PROGRAM_COLS)       # enforce consistent column order
    )


def zscore_columns(mat: pd.DataFrame) -> pd.DataFrame:
    """
    Z-score each program (column) independently across the rows of mat.
    Programs with zero variance are left as zero.
    """
    result = mat.copy().astype(float)
    for col in result.columns:
        sigma = result[col].std(ddof=0)
        mu    = result[col].mean()
        result[col] = (result[col] - mu) / sigma if sigma > 0 else 0.0
    return result


def save_fig(fig: plt.Figure, path: str, dpi: int) -> None:
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Wrote: {path}")


# ---------------------------------------------------------------------------
# Plot functions
# ---------------------------------------------------------------------------

def plot_heatmap(
    long_df: pd.DataFrame,
    group_col: str,
    title: str,
    out_path: str,
    dpi: int,
) -> None:
    """
    Z-scored heatmap: rows = groups, columns = programs.

    Z-scoring is applied per-program across the groups present in this
    specific matrix, so each column has mean 0 and std 1 within the plot.
    This makes cross-program visual comparison meaningful regardless of
    absolute score magnitudes.
    """
    mat    = pivot_to_matrix(long_df, group_col)
    zmat   = zscore_columns(mat)

    # Replace program column names with short display labels
    zmat.columns = [PROGRAM_LABELS.get(c, c) for c in zmat.columns]

    n_rows = len(zmat)
    fig_h  = max(4, 0.45 * n_rows + 1.5)
    fig, ax = plt.subplots(figsize=(8, fig_h))

    sns.heatmap(
        zmat,
        ax=ax,
        cmap="RdBu_r",
        center=0,
        linewidths=0.4,
        linecolor="white",
        annot=True,
        fmt=".2f",
        annot_kws={"size": 7},
        cbar_kws={"label": "z-score", "shrink": 0.6},
    )

    ax.set_title(title, fontsize=11, pad=10)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelsize=8, rotation=30)
    ax.tick_params(axis="y", labelsize=8, rotation=0)

    save_fig(fig, out_path, dpi)


def plot_barplot(
    long_df: pd.DataFrame,
    program: str,
    out_path: str,
    dpi: int,
) -> None:
    """
    Grouped bar chart for one program.

    Cell types are sorted descending by mean score so the most active
    population is leftmost, making the figure easy to scan without a legend.
    """
    sub = (
        long_df[long_df["program"] == program]
        .sort_values("mean", ascending=False)
        .reset_index(drop=True)
    )

    if sub.empty:
        print(f"  Warning: no data for program '{program}', skipping bar chart")
        return

    n_bars = len(sub)
    fig_w  = max(6, 0.55 * n_bars + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, 4.5))

    x_pos = range(n_bars)
    ax.bar(
        x_pos,
        sub["mean"],
        yerr=sub["std"].fillna(0),
        capsize=4,
        color=sns.color_palette("muted", n_bars),
        edgecolor="white",
        linewidth=0.5,
        error_kw={"elinewidth": 1.0, "ecolor": "black"},
    )

    ax.set_xticks(list(x_pos))
    ax.set_xticklabels(sub["cell_type"], rotation=40, ha="right", fontsize=8)
    ax.set_ylabel("Mean program score", fontsize=9)
    ax.set_title(
        f"{PROGRAM_LABELS.get(program, program).replace(chr(10), ' ')}\n"
        f"sorted by mean, error bars = std",
        fontsize=10,
    )
    ax.axhline(0, color="black", linewidth=0.6, linestyle="--")
    ax.tick_params(axis="y", labelsize=8)
    sns.despine(ax=ax)

    save_fig(fig, out_path, dpi)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # -----------------------------------------------------------------------
    # Load inputs
    # -----------------------------------------------------------------------

    celltype_csv = os.path.join(args.results_dir, "tf_activity_rna_by_celltype.csv")
    broad_csv    = os.path.join(args.results_dir, "tf_activity_rna_by_broad.csv")

    print("Loading tf_activity_rna_by_celltype.csv...")
    celltype_df = load_csv(celltype_csv, ["cell_type", "program", "mean", "std", "n_cells"])

    print("Loading tf_activity_rna_by_broad.csv...")
    broad_df = load_csv(broad_csv, ["broad_label", "program", "mean", "std", "n_cells"])

    # -----------------------------------------------------------------------
    # Heatmap — RNA by cell type
    # -----------------------------------------------------------------------

    print("Plotting heatmap: RNA by cell type...")
    plot_heatmap(
        long_df   = celltype_df,
        group_col = "cell_type",
        title     = "TF activity proxy — RNA baseline by cell type\n(z-scored per program)",
        out_path  = os.path.join(args.outdir, f"tf_heatmap_rna_celltype.{args.fmt}"),
        dpi       = args.dpi,
    )

    # -----------------------------------------------------------------------
    # Heatmap — RNA by broad label
    # -----------------------------------------------------------------------

    print("Plotting heatmap: RNA by broad label...")
    plot_heatmap(
        long_df   = broad_df,
        group_col = "broad_label",
        title     = "TF activity proxy — RNA baseline by broad label\n(z-scored per program)",
        out_path  = os.path.join(args.outdir, f"tf_heatmap_rna_broad.{args.fmt}"),
        dpi       = args.dpi,
    )

    # -----------------------------------------------------------------------
    # Bar charts — one per program
    # -----------------------------------------------------------------------

    print("Plotting per-program bar charts...")
    for prog in PROGRAM_COLS:
        plot_barplot(
            long_df  = celltype_df,
            program  = prog,
            out_path = os.path.join(args.outdir, f"tf_barplot_{prog}.{args.fmt}"),
            dpi      = args.dpi,
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
