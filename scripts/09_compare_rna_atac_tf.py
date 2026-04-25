"""
09_compare_rna_atac_tf.py

Phase 4 — Cross-modal comparison of RNA TF-activity proxies vs ATAC motif deviations.

For each of the five Phase 3 TF-activity program axes, collapses ATAC chromVAR
motif deviation scores to program level (averaging across motifs within the
same program), then correlates per-cell-type RNA proxy scores against per-cell-type
ATAC deviation scores.  Writes a correlation table and six figures.

Inputs (results/regulons/):
  tf_activity_rna_by_celltype.csv     from script 06
  atac_motif_by_celltype.csv          from script 08

Outputs:
  results/comparison/
    tf_rna_atac_correlation.csv

  results/figures/tf_activity/
    tf_rna_atac_scatter_interferon_alpha_beta_response.<fmt>
    tf_rna_atac_scatter_inflammatory_response.<fmt>
    tf_rna_atac_scatter_antigen_presentation_mhc.<fmt>
    tf_rna_atac_scatter_cytotoxic_nk_like.<fmt>
    tf_rna_atac_scatter_monocyte_activation.<fmt>
    tf_heatmap_atac_celltype.<fmt>

Usage:
  python scripts/09_compare_rna_atac_tf.py
  python scripts/09_compare_rna_atac_tf.py --help
"""

import argparse
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")           # non-interactive backend; must precede pyplot import
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats

sns.set_theme(style="white")

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

TF_AXES = {
    "interferon_alpha_beta_response": "IRF3 / IRF7 / STAT1/2",
    "inflammatory_response":          "NF-κB / AP-1",
    "antigen_presentation_mhc":       "CIITA / IRF1",
    "cytotoxic_nk_like":              "EOMES / TBX21",
    "monocyte_activation":            "SPI1 (PU.1) / CEBPB",
}

# Short display labels (match script 07 exactly)
PROGRAM_LABELS = {
    "interferon_alpha_beta_response": "IFN α/β\n(IRF3/7·STAT1/2)",
    "inflammatory_response":          "Inflammation\n(NF-κB/AP-1)",
    "antigen_presentation_mhc":       "Antigen pres.\n(CIITA/IRF1)",
    "cytotoxic_nk_like":              "Cytotoxic\n(EOMES/TBX21)",
    "monocyte_activation":            "Mono. activ.\n(PU.1/CEBPB)",
}

ATAC_TO_BROAD = {
    "Activated_CD4_T_GATA3_high": "CD4 T cells",
    "CD4_T_naive_CCR7_IL7R": "CD4 T cells",
    "Treg_FOXP3": "CD4 T cells",
    "Proliferating_or_activated_T": "CD4 T cells",

    "CD8_T_naive_memory": "CD8 T cells",
    "T_cell_like_8": "CD8 T cells",

    "B_cell_like": "B cells",
    "Naive_B_cells_IGHD_TCL1A": "B cells",
    "Plasma_cell": "B cells",

    "CD14_Monocytes_inflammatory": "CD14+ Monocytes",
    "FCGR3A_Monocytes": "FCGR3A+ Monocytes",

    "CD1C_DC": "Dendritic cells",
    "pDC_like": "Dendritic cells",

    "Cytotoxic_NK_TRGC2": "NK cells",
    "NK_cells_NCAM1_KLRF1": "NK cells",
    "RORC_like_innate_lymphoid": "NK cells",
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
        description="Cross-modal comparison: RNA TF-activity proxies vs ATAC motif deviations."
    )
    parser.add_argument(
        "--rna-activity",
        default=os.path.join(_REPO_ROOT, "results", "regulons",
                             "tf_activity_rna_by_celltype.csv"),
        help="Path to tf_activity_rna_by_celltype.csv (script 06) [default: %(default)s]",
    )
    parser.add_argument(
        "--atac-motif",
        default=os.path.join(_REPO_ROOT, "results", "regulons",
                             "atac_motif_by_celltype.csv"),
        help="Path to atac_motif_by_celltype.csv (script 08) [default: %(default)s]",
    )
    parser.add_argument(
        "--outdir",
        default=os.path.join(_REPO_ROOT, "results", "figures", "tf_activity"),
        help="Output directory for figures [default: %(default)s]",
    )
    parser.add_argument(
        "--corr-outdir",
        default=os.path.join(_REPO_ROOT, "results", "comparison"),
        help="Output directory for correlation CSV [default: %(default)s]",
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


def load_csv(path: str, required_cols: list, run_hint: str) -> pd.DataFrame:
    if not os.path.isfile(path):
        stop(f"Input file not found: {path}\n{run_hint}")
    df = pd.read_csv(path)
    missing = set(required_cols) - set(df.columns)
    if missing:
        stop(f"{path} is missing columns: {sorted(missing)}")
    return df


def save_fig(fig: plt.Figure, path: str, dpi: int) -> None:
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Wrote: {path}")


def zscore_columns(mat: pd.DataFrame) -> pd.DataFrame:
    """
    Z-score each column independently.  Zero-variance columns → 0.0.
    Matches the zscore_columns logic used in script 07.
    """
    result = mat.copy().astype(float)
    for col in result.columns:
        sigma = result[col].std(ddof=0)
        mu    = result[col].mean()
        result[col] = (result[col] - mu) / sigma if sigma > 0 else 0.0
    return result


def map_atac_to_broad_labels(atac_df: pd.DataFrame) -> pd.DataFrame:
    out = atac_df.copy()
    out["cell_type_broad"] = out["cell_type"].map(ATAC_TO_BROAD)
    return out

# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def build_atac_program_scores(atac_df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse per-motif ATAC deviation means to one score per (cell_type, program).

    For programs with multiple motifs (e.g. monocyte_activation has SPI1 and
    CEBPB), the mean of the per-motif mean deviations is used.  For programs
    with a single motif, this is a no-op average.

    Returns columns: cell_type, program, atac_mean.
    """
    return (
        atac_df
        .groupby(["cell_type", "program"])["mean"]
        .mean()
        .reset_index()
        .rename(columns={"mean": "atac_mean"})
    )


def compute_correlations(merged_df: pd.DataFrame) -> pd.DataFrame:
    """
    Pearson r and two-tailed p-value across cell types for each program.

    Returns columns: program, tf_axis, n_celltypes, pearson_r, pearson_p.
    NaN is written for r and p if scipy raises (e.g. zero-variance input).
    """
    rows = []
    for prog in PROGRAM_COLS:
        sub = (
            merged_df[merged_df["program"] == prog]
            .dropna(subset=["rna_mean", "atac_mean"])
        )
        n = len(sub)
        try:
            r, p = stats.pearsonr(sub["rna_mean"], sub["atac_mean"])
        except Exception:
            r, p = float("nan"), float("nan")
        rows.append({
            "program":     prog,
            "tf_axis":     TF_AXES[prog],
            "n_celltypes": n,
            "pearson_r":   r,
            "pearson_p":   p,
        })
    return pd.DataFrame(
        rows,
        columns=["program", "tf_axis", "n_celltypes", "pearson_r", "pearson_p"],
    )


# ---------------------------------------------------------------------------
# Plot functions
# ---------------------------------------------------------------------------

def plot_scatter(
    sub_df:   pd.DataFrame,
    program:  str,
    r_val:    float,
    p_val:    float,
    out_path: str,
    dpi:      int,
) -> None:
    """
    Scatter: one point per cell type.  x = RNA program score, y = ATAC deviation.
    Cell types colored by alphabetical order using sns.color_palette("muted").
    A linear regression fit line (seaborn.regplot) is drawn over the points.
    """
    cell_types = sorted(sub_df["cell_type"].unique())
    palette    = sns.color_palette("muted", len(cell_types))
    color_map  = dict(zip(cell_types, palette))

    fig, ax = plt.subplots(figsize=(5.5, 4.5))

    for _, row in sub_df.iterrows():
        ax.scatter(
            row["rna_mean"],
            row["atac_mean"],
            s=80,
            color=color_map[row["cell_type"]],
            edgecolors="white",
            linewidths=0.5,
            zorder=3,
        )
        ax.annotate(
            row["cell_type"],
            (row["rna_mean"], row["atac_mean"]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7,
        )

    # Linear regression fit line
    sns.regplot(
        x="rna_mean",
        y="atac_mean",
        data=sub_df,
        ax=ax,
        scatter=False,
        color="black",
        line_kws={"linewidth": 1.0, "alpha": 0.8},
    )

    # Pearson annotation (upper-left corner)
    if not (math.isnan(r_val) or math.isnan(p_val)):
        annot = f"r = {r_val:.2f}\np = {p_val:.3f}"
    else:
        annot = "r = n/a"
    ax.text(
        0.05, 0.95, annot,
        transform=ax.transAxes,
        va="top", ha="left",
        fontsize=8,
    )

    display_label = PROGRAM_LABELS.get(program, program).replace("\n", " ")
    ax.set_title(f"{display_label}\nRNA vs ATAC", fontsize=10)
    ax.set_xlabel("RNA program score (mean)", fontsize=9)
    ax.set_ylabel("ATAC motif deviation (mean z-score)", fontsize=9)
    ax.tick_params(labelsize=8)
    sns.despine(ax=ax)

    save_fig(fig, out_path, dpi)


def plot_atac_heatmap(
    atac_program_df: pd.DataFrame,
    out_path:        str,
    dpi:             int,
) -> None:
    """
    Z-scored heatmap of ATAC motif deviation by cell type.
    Layout mirrors the Phase 3 RNA heatmap (script 07 plot_heatmap).
    """
    mat = (
        atac_program_df
        .pivot_table(index="cell_type", columns="program", values="atac_mean")
        .reindex(columns=PROGRAM_COLS)
    )
    zmat = zscore_columns(mat)
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
        cbar_kws={"label": "z-score (ATAC motif deviation)", "shrink": 0.6},
    )

    ax.set_title(
        "TF activity proxy — ATAC motif deviation by cell type\n(z-scored per program)",
        fontsize=11,
        pad=10,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelsize=8, rotation=30)
    ax.tick_params(axis="y", labelsize=8, rotation=0)

    save_fig(fig, out_path, dpi)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir,      exist_ok=True)
    os.makedirs(args.corr_outdir, exist_ok=True)

    # -----------------------------------------------------------------------
    # Load inputs
    # -----------------------------------------------------------------------

    print("Loading RNA TF activity by cell type...")
    rna_df = load_csv(
        args.rna_activity,
        required_cols=["cell_type", "program", "mean"],
        run_hint="Run scripts/06_tf_activity.py first.",
    )

    print("Loading ATAC motif deviation by cell type...")
    atac_df = load_csv(
        args.atac_motif,
        required_cols=["cell_type", "motif_id", "tf_name", "program", "mean"],
        run_hint="Run scripts/08_motif_enrichment_atac.R first.",
    )

    # -----------------------------------------------------------------------
    # Collapse motif means to program-level ATAC scores
    # -----------------------------------------------------------------------

    print("Mapping ATAC fine labels to broad RNA-compatible labels...")
    atac_df = map_atac_to_broad_labels(atac_df)

    unmapped = sorted(atac_df.loc[atac_df["cell_type_broad"].isna(), "cell_type"].unique())
    if unmapped:
        print(f"  Warning: dropping unmapped ATAC cell types: {unmapped}", file=sys.stderr)

    atac_df = atac_df.dropna(subset=["cell_type_broad"]).copy()
    atac_df["cell_type"] = atac_df["cell_type_broad"]

    print("Collapsing per-motif ATAC means to broad program-level scores...")
    atac_program_df = build_atac_program_scores(atac_df)

    # Identify programs absent from the ATAC data (warn; skip those programs later)
    atac_programs    = set(atac_program_df["program"].unique())
    missing_programs = [p for p in PROGRAM_COLS if p not in atac_programs]
    for prog in missing_programs:
        print(
            f"  Warning: program '{prog}' has no rows in ATAC data — "
            f"scatter and heatmap column will be skipped",
            file=sys.stderr,
        )

    # -----------------------------------------------------------------------
    # Inner join: RNA program scores × ATAC program scores
    # -----------------------------------------------------------------------

    rna_prog = (
        rna_df[["cell_type", "program", "mean"]]
        .rename(columns={"mean": "rna_mean"})
    )
    merged = pd.merge(
        rna_prog,
        atac_program_df,
        on=["cell_type", "program"],
        how="inner",
    )
    # Validate per-program cell-type counts after join
    for prog in PROGRAM_COLS:
        if prog in missing_programs:
            continue
        sub = merged[merged["program"] == prog]
        n   = len(sub)
        if n < 3:
            stop(
                f"Only {n} cell type(s) in the inner join for program '{prog}' "
                f"(minimum 3 required).\n"
                f"Check that cell_type labels overlap between "
                f"tf_activity_rna_by_celltype.csv and atac_motif_by_celltype.csv."
            )
        if n < 5:
            print(
                f"  Warning: only {n} cell types for program '{prog}' "
                f"(fewer than 5 — correlation estimate may be unreliable)",
                file=sys.stderr,
            )

    # -----------------------------------------------------------------------
    # Compute Pearson correlations
    # -----------------------------------------------------------------------

    print("Computing Pearson correlations (RNA proxy vs ATAC deviation)...")
    corr_df = compute_correlations(merged)

    corr_path = os.path.join(args.corr_outdir, "tf_rna_atac_correlation.csv")
    corr_df.to_csv(corr_path, index=False)
    print(f"  Wrote: {corr_path}")

    for _, row in corr_df.iterrows():
        r_str = f"{row['pearson_r']:.2f}" if not math.isnan(row["pearson_r"]) else "n/a"
        p_str = f"{row['pearson_p']:.3f}" if not math.isnan(row["pearson_p"]) else "n/a"
        print(f"    {row['program']}: r={r_str}, p={p_str}, n={row['n_celltypes']}")

    # -----------------------------------------------------------------------
    # Scatter plots — one per program
    # -----------------------------------------------------------------------

    print("Plotting RNA vs ATAC scatter plots...")
    for prog in PROGRAM_COLS:
        if prog in missing_programs:
            continue
        sub   = merged[merged["program"] == prog].copy()
        r_val = corr_df.loc[corr_df["program"] == prog, "pearson_r"].iloc[0]
        p_val = corr_df.loc[corr_df["program"] == prog, "pearson_p"].iloc[0]
        plot_scatter(
            sub_df   = sub,
            program  = prog,
            r_val    = r_val,
            p_val    = p_val,
            out_path = os.path.join(
                args.outdir, f"tf_rna_atac_scatter_{prog}.{args.fmt}"
            ),
            dpi=args.dpi,
        )

    # -----------------------------------------------------------------------
    # ATAC heatmap
    # -----------------------------------------------------------------------

    print("Plotting ATAC heatmap...")
    plot_atac_heatmap(
        atac_program_df=atac_program_df,
        out_path=os.path.join(args.outdir, f"tf_heatmap_atac_celltype.{args.fmt}"),
        dpi=args.dpi,
    )

    print("\nDone.")


if __name__ == "__main__":
    main()
