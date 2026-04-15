"""
06_tf_activity.py

Phase 3 — TF activity proxy summaries from RNA program scores.

Computes per-cell-type and per-broad-label mean/std/n summaries for the five
program score columns from the Phase 1 RNA baseline, and writes a deduplicated
WNN fine-to-broad label mapping table for downstream plotting.

Inputs:
  results/baseline/rna_program_scores.csv
  results/integration/multiome_cell_metadata.csv
  results/integration/cluster_annotations.tsv
  results/comparison/wnn_to_broad_labels.tsv

Outputs (results/regulons/):
  tf_activity_rna_by_celltype.csv
  tf_activity_rna_by_broad.csv
  tf_activity_wnn_celltype_map.csv
  tf_activity_summary.md

Usage:
  python scripts/06_tf_activity.py
  python scripts/06_tf_activity.py --help
"""

import argparse
import os
import sys

import pandas as pd

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

# Human-readable TF-axis label for each program (used in summary only)
PROGRAM_TF_AXES = {
    "interferon_alpha_beta_response": "IRF3 / IRF7 / STAT1/2",
    "inflammatory_response":          "NF-κB / AP-1",
    "antigen_presentation_mhc":       "CIITA / IRF1",
    "cytotoxic_nk_like":              "EOMES / TBX21",
    "monocyte_activation":            "SPI1 (PU.1) / CEBPB",
}

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
        description="Compute TF activity proxy summaries from RNA program scores."
    )
    parser.add_argument(
        "--rna-program-scores",
        default=os.path.join(_REPO_ROOT, "results", "baseline", "rna_program_scores.csv"),
        help="Path to rna_program_scores.csv [default: %(default)s]",
    )
    parser.add_argument(
        "--wnn-meta",
        default=os.path.join(_REPO_ROOT, "results", "integration", "multiome_cell_metadata.csv"),
        help="Path to multiome_cell_metadata.csv [default: %(default)s]",
    )
    parser.add_argument(
        "--cluster-annotations",
        default=os.path.join(_REPO_ROOT, "results", "integration", "cluster_annotations.tsv"),
        help="Path to cluster_annotations.tsv [default: %(default)s]",
    )
    parser.add_argument(
        "--label-map",
        default=os.path.join(_REPO_ROOT, "results", "comparison", "wnn_to_broad_labels.tsv"),
        help="Path to wnn_to_broad_labels.tsv [default: %(default)s]",
    )
    parser.add_argument(
        "--outdir",
        default=os.path.join(_REPO_ROOT, "results", "regulons"),
        help="Output directory [default: %(default)s]",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def stop(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def load_csv(path: str, required_cols: list, sep: str = ",") -> pd.DataFrame:
    """Load a delimited file and validate that required columns are present."""
    if not os.path.isfile(path):
        stop(f"Input file not found: {path}")
    df = pd.read_csv(path, sep=sep)
    missing = set(required_cols) - set(df.columns)
    if missing:
        stop(f"{path} is missing required columns: {sorted(missing)}")
    return df


def melt_summary(grouped: pd.DataFrame, group_col: str) -> pd.DataFrame:
    """
    Convert a wide per-group mean/std/n table to long form.

    grouped has columns: [group_col, program_mean..., program_std..., program_n...]
    Returns long form: [group_col, program, mean, std, n_cells]
    """
    rows = []
    for prog in PROGRAM_COLS:
        tmp = grouped[[group_col, f"{prog}_mean", f"{prog}_std", f"{prog}_n"]].copy()
        tmp.columns = [group_col, "mean", "std", "n_cells"]
        tmp.insert(1, "program", prog)
        rows.append(tmp)
    return pd.concat(rows, ignore_index=True).sort_values([group_col, "program"])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # -----------------------------------------------------------------------
    # 1–2. Load and validate RNA program scores
    # -----------------------------------------------------------------------

    print("Loading RNA program scores...")
    prog_df = load_csv(
        args.rna_program_scores,
        required_cols=["cell_id", "cell_type"] + PROGRAM_COLS,
    )
    print(f"  {len(prog_df):,} cells loaded from RNA program scores")

    # -----------------------------------------------------------------------
    # 3. Compute RNA per-cell_type summaries (mean, std, n_cells per program)
    # -----------------------------------------------------------------------

    print("Computing per-cell_type TF activity summaries...")

    agg_dict = {}
    for prog in PROGRAM_COLS:
        agg_dict[f"{prog}_mean"] = (prog, "mean")
        agg_dict[f"{prog}_std"]  = (prog, "std")
        agg_dict[f"{prog}_n"]    = (prog, "count")

    celltype_wide = (
        prog_df
        .groupby("cell_type", sort=True)
        .agg(**agg_dict)
        .reset_index()
    )

    celltype_long = melt_summary(celltype_wide, "cell_type")

    out_celltype = os.path.join(args.outdir, "tf_activity_rna_by_celltype.csv")
    celltype_long.to_csv(out_celltype, index=False)
    print(f"  Wrote: {out_celltype}")

    # -----------------------------------------------------------------------
    # 4–6. Broad-label table: RNA cell_type IS the broad label for this phase
    # -----------------------------------------------------------------------

    broad_long = celltype_long.rename(columns={"cell_type": "broad_label"}).copy()

    out_broad = os.path.join(args.outdir, "tf_activity_rna_by_broad.csv")
    broad_long.to_csv(out_broad, index=False)
    print(f"  Wrote: {out_broad}")

    # -----------------------------------------------------------------------
    # 7–8. Load and validate WNN cell metadata
    # -----------------------------------------------------------------------

    print("Loading WNN cell metadata...")
    wnn_meta = load_csv(
        args.wnn_meta,
        required_cols=["cell_id", "cell_type_wnn"],
    )
    print(f"  {len(wnn_meta):,} cells loaded from WNN metadata")

    # # -----------------------------------------------------------------------
    # # 9–10. Load cluster annotations and validate coverage
    # # -----------------------------------------------------------------------

    # print("Loading cluster annotations...")
    # cluster_annots = load_csv(
    #     args.cluster_annotations,
    #     required_cols=["cluster", "cell_type"],
    #     sep="\t",
    # )

    # # Coerce cluster IDs to str for consistent comparison
    # wnn_meta["cluster"]        = wnn_meta["cluster"].astype(str)
    # cluster_annots["cluster"]  = cluster_annots["cluster"].astype(str)

    # wnn_clusters    = set(wnn_meta["cluster"].unique())
    # annot_clusters  = set(cluster_annots["cluster"].unique())
    # uncovered       = wnn_clusters - annot_clusters
    # if uncovered:
    #     stop(
    #         f"cluster_annotations.tsv is missing {len(uncovered)} cluster(s) "
    #         f"present in WNN metadata: {sorted(uncovered)}"
    #     )

    # -----------------------------------------------------------------------
    # 11–13. Load broad-label map and perform WNN fine-to-broad mapping
    # -----------------------------------------------------------------------

    print("Loading WNN-to-broad label map...")
    label_map_df = load_csv(
        args.label_map,
        required_cols=["fine_label", "broad_label"],
        sep="\t",
    )

    fine_to_broad = dict(zip(label_map_df["fine_label"], label_map_df["broad_label"]))

    wnn_fine_labels  = set(wnn_meta["cell_type_wnn"].dropna().unique())
    mapped_labels    = wnn_fine_labels & set(fine_to_broad.keys())
    unmapped_labels  = sorted(wnn_fine_labels - set(fine_to_broad.keys()))

    if unmapped_labels:
        print(
            f"  Warning: {len(unmapped_labels)} cell_type_wnn value(s) have no "
            f"match in wnn_to_broad_labels.tsv — assigned 'Unassigned': "
            + ", ".join(unmapped_labels)
        )

    wnn_meta["broad_label"] = (
        wnn_meta["cell_type_wnn"]
        .map(lambda v: fine_to_broad.get(v, "Unassigned") if pd.notna(v) else "Unassigned")
    )

    # -----------------------------------------------------------------------
    # 14. Write deduplicated WNN fine-to-broad mapping table
    # -----------------------------------------------------------------------

    wnn_map_df = (
        wnn_meta[["cell_type_wnn", "broad_label"]]
        .dropna(subset=["cell_type_wnn"])
        .drop_duplicates()
        .sort_values("cell_type_wnn")
        .reset_index(drop=True)
    )

    out_map = os.path.join(args.outdir, "tf_activity_wnn_celltype_map.csv")
    wnn_map_df.to_csv(out_map, index=False)
    print(f"  Wrote: {out_map}")

    # -----------------------------------------------------------------------
    # 15. Write markdown summary
    # -----------------------------------------------------------------------

    n_rna_celltypes  = celltype_long["cell_type"].nunique()
    n_broad_labels   = broad_long["broad_label"].nunique()
    n_wnn_fine       = len(wnn_fine_labels)
    n_mapped         = len(mapped_labels)

    unmapped_block = (
        "\n".join(f"  - {lbl}" for lbl in unmapped_labels)
        if unmapped_labels else "  *(none)*"
    )

    program_block = "\n".join(
        f"  - `{p}` → {PROGRAM_TF_AXES[p]}" for p in PROGRAM_COLS
    )

    summary = (
        "# TF Activity Summary — Phase 3\n\n"
        "## Input row counts\n\n"
        "| File | Rows |\n"
        "|------|------|\n"
        f"| rna_program_scores.csv | {len(prog_df):,} |\n"
        f"| multiome_cell_metadata.csv | {len(wnn_meta):,} |\n"
        # f"| cluster_annotations.tsv | {len(cluster_annots):,} |\n"
        f"| wnn_to_broad_labels.tsv | {len(label_map_df):,} |\n\n"
        "## Detected program score columns\n\n"
        f"{program_block}\n\n"
        "## RNA baseline summary\n\n"
        f"- Cell types summarised: **{n_rna_celltypes}**\n"
        f"- Broad labels in output: **{n_broad_labels}** "
        "(RNA cell_type used directly as broad label)\n\n"
        "## WNN label mapping\n\n"
        f"- Unique WNN fine labels: **{n_wnn_fine}**\n"
        f"- Mapped to broad labels: **{n_mapped}**\n"
        f"- Unmapped (assigned 'Unassigned'): **{len(unmapped_labels)}**\n\n"
        f"**Unmapped labels:**\n{unmapped_block}\n"
    )

    out_summary = os.path.join(args.outdir, "tf_activity_summary.md")
    with open(out_summary, "w") as fh:
        fh.write(summary)
    print(f"  Wrote: {out_summary}")

    print("\nDone.")


if __name__ == "__main__":
    main()
