"""
05_compare_rna_vs_wnn.py

Population-level comparison of Phase 1 RNA baseline cell-type labels versus
Phase 2 multiome WNN cell-type labels.

Inputs:
  results/baseline/rna_cell_metadata.csv         -- cell_id, condition, cell_type
  results/integration/multiome_cell_metadata.csv -- cell_id, ..., cell_type_wnn

Reserved for future phases (not used here):
  results/baseline/rna_ifn_scores.csv
  results/baseline/rna_program_scores.csv

Outputs:
  results/comparison/rna_vs_wnn_celltype_counts.csv
  results/comparison/rna_vs_wnn_celltype_fractions.csv
  results/comparison/rna_vs_wnn_summary.md

Usage:
  python scripts/05_compare_rna_vs_wnn.py
  python scripts/05_compare_rna_vs_wnn.py --help
"""

import argparse
import os
import sys

import pandas as pd

# ---------------------------------------------------------------------------
# Resolve repo root relative to this script (works from any working directory)
# ---------------------------------------------------------------------------

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT  = os.path.dirname(_SCRIPT_DIR)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Population-level comparison of RNA baseline vs. multiome WNN cell types."
    )
    parser.add_argument(
        "--rna-meta",
        default=os.path.join(_REPO_ROOT, "results", "baseline", "rna_cell_metadata.csv"),
        help="Path to rna_cell_metadata.csv from Phase 1 [default: %(default)s]",
    )
    parser.add_argument(
        "--wnn-meta",
        default=os.path.join(_REPO_ROOT, "results", "integration", "multiome_cell_metadata.csv"),
        help="Path to multiome_cell_metadata.csv from Phase 2 [default: %(default)s]",
    )
    # Reserved for future phases — accepted but not used in this script.
    parser.add_argument(
        "--rna-ifn-scores",
        default=os.path.join(_REPO_ROOT, "results", "baseline", "rna_ifn_scores.csv"),
        help="[Reserved for future use] Path to rna_ifn_scores.csv",
    )
    parser.add_argument(
        "--rna-program-scores",
        default=os.path.join(_REPO_ROOT, "results", "baseline", "rna_program_scores.csv"),
        help="[Reserved for future use] Path to rna_program_scores.csv",
    )
    parser.add_argument(
        "--outdir",
        default=os.path.join(_REPO_ROOT, "results", "comparison"),
        help="Output directory [default: %(default)s]",
    )
    parser.add_argument(
        "--label-map",
        default=None,
        help=(
            "Optional two-column TSV (fine_label, broad_label) mapping WNN "
            "cell_type_wnn values to broad labels for a harmonised comparison. "
            "Unmapped labels are kept unchanged."
        ),
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_and_validate(path: str, required_cols: list) -> pd.DataFrame:
    if not os.path.isfile(path):
        sys.exit(f"ERROR: Input file not found: {path}")
    df = pd.read_csv(path)
    missing = set(required_cols) - set(df.columns)
    if missing:
        sys.exit(f"ERROR: {path} is missing required columns: {sorted(missing)}")
    return df


def celltype_counts(df: pd.DataFrame, label_col: str, dataset_name: str) -> pd.DataFrame:
    """Return per-cell-type counts for one dataset."""
    counts = (
        df[label_col]
        .dropna()
        .value_counts()
        .reset_index()
    )
    counts.columns = ["cell_type", "count"]
    counts.insert(0, "dataset", dataset_name)
    return counts.sort_values("cell_type").reset_index(drop=True)


def add_fraction(counts: pd.DataFrame) -> pd.DataFrame:
    """Add a fraction column in-place; fraction is within-dataset."""
    counts = counts.copy()
    counts["fraction"] = counts["count"] / counts["count"].sum()
    return counts


def bullet_list(items: list) -> str:
    return "\n".join(f"  - {x}" for x in items) if items else "  *(none)*"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    # -- Load inputs ---------------------------------------------------------
    rna_meta = load_and_validate(args.rna_meta, ["cell_id", "cell_type"])
    wnn_meta = load_and_validate(args.wnn_meta, ["cell_id", "cell_type_wnn"])

    os.makedirs(args.outdir, exist_ok=True)

    # -- Cell-type counts per dataset ----------------------------------------
    rna_counts = celltype_counts(rna_meta, "cell_type",     "rna_baseline")
    wnn_counts = celltype_counts(wnn_meta, "cell_type_wnn", "multiome_wnn")

    counts_df = pd.concat([rna_counts, wnn_counts], ignore_index=True)
    counts_path = os.path.join(args.outdir, "rna_vs_wnn_celltype_counts.csv")
    counts_df.to_csv(counts_path, index=False)

    # -- Cell-type fractions per dataset -------------------------------------
    fractions_df = pd.concat(
        [add_fraction(rna_counts), add_fraction(wnn_counts)],
        ignore_index=True,
    )
    fractions_path = os.path.join(args.outdir, "rna_vs_wnn_celltype_fractions.csv")
    fractions_df.to_csv(fractions_path, index=False)

    # -- Label overlap statistics --------------------------------------------
    rna_labels = set(rna_meta["cell_type"].dropna().unique())
    wnn_labels = set(wnn_meta["cell_type_wnn"].dropna().unique())

    shared   = sorted(rna_labels & wnn_labels)
    rna_only = sorted(rna_labels - wnn_labels)
    wnn_only = sorted(wnn_labels - rna_labels)

    # -- Markdown summary ----------------------------------------------------
    summary = (
        "# RNA vs. WNN Cell-Type Comparison Summary\n\n"
        "## Dataset sizes\n\n"
        "| Dataset | Cells |\n"
        "|---------|-------|\n"
        f"| RNA baseline (`rna_cell_metadata`) | {len(rna_meta):,} |\n"
        f"| Multiome WNN (`multiome_cell_metadata`) | {len(wnn_meta):,} |\n\n"
        "## Cell-type label counts\n\n"
        "| Dataset | Unique labels |\n"
        "|---------|---------------|\n"
        f"| RNA baseline | {len(rna_labels)} |\n"
        f"| Multiome WNN | {len(wnn_labels)} |\n\n"
        "## Label overlap\n\n"
        f"**Shared labels ({len(shared)}):**\n"
        f"{bullet_list(shared)}\n\n"
        f"**RNA baseline only ({len(rna_only)}):**\n"
        f"{bullet_list(rna_only)}\n\n"
        f"**Multiome WNN only ({len(wnn_only)}):**\n"
        f"{bullet_list(wnn_only)}\n\n"
        "## Notes\n\n"
        "- This is a population-level label comparison across independent datasets.\n"
        "- Barcodes are not shared; no cell-level join is performed.\n"
        "- Exact label overlap may underestimate biological overlap when label "
        "granularity differs between datasets (e.g. a broad RNA label may map to "
        "multiple fine-grained WNN clusters).\n"
        "- IFN and program score comparisons are reserved for a future phase.\n"
    )

    summary_path = os.path.join(args.outdir, "rna_vs_wnn_summary.md")
    with open(summary_path, "w") as fh:
        fh.write(summary)

    # -- Broad-label comparison (optional) -----------------------------------
    if args.label_map is not None:
        if not os.path.isfile(args.label_map):
            sys.exit(f"ERROR: --label-map file not found: {args.label_map}")

        map_df = pd.read_csv(args.label_map, sep="\t")
        missing_map_cols = {"fine_label", "broad_label"} - set(map_df.columns)
        if missing_map_cols:
            sys.exit(
                f"ERROR: --label-map TSV is missing columns: {sorted(missing_map_cols)}"
                "\nExpected: fine_label<TAB>broad_label"
            )

        # Build mapping dict; unmapped labels fall back to the original value
        label_mapping = dict(zip(map_df["fine_label"], map_df["broad_label"]))
        broad_wnn = wnn_meta["cell_type_wnn"].map(
            lambda v: label_mapping.get(v, v) if pd.notna(v) else v
        )

        # Counts and fractions for the broad comparison
        wnn_meta_broad = wnn_meta.copy()
        wnn_meta_broad["broad_cell_type"] = broad_wnn

        rna_broad_counts = celltype_counts(rna_meta,      "cell_type",       "rna_baseline")
        wnn_broad_counts = celltype_counts(wnn_meta_broad, "broad_cell_type", "multiome_wnn_broad")

        broad_counts_df   = pd.concat([rna_broad_counts, wnn_broad_counts], ignore_index=True)
        broad_counts_path = os.path.join(args.outdir, "rna_vs_wnn_broad_counts.csv")
        broad_counts_df.to_csv(broad_counts_path, index=False)

        # Label overlap for broad comparison
        wnn_broad_labels = set(broad_wnn.dropna().unique())
        broad_shared   = sorted(rna_labels & wnn_broad_labels)
        broad_rna_only = sorted(rna_labels - wnn_broad_labels)
        broad_wnn_only = sorted(wnn_broad_labels - rna_labels)

        broad_summary = (
            "# RNA vs. WNN Broad Cell-Type Comparison Summary\n\n"
            f"Label mapping file: `{os.path.basename(args.label_map)}`\n\n"
            "## Dataset sizes\n\n"
            "| Dataset | Cells |\n"
            "|---------|-------|\n"
            f"| RNA baseline (`rna_cell_metadata`) | {len(rna_meta):,} |\n"
            f"| Multiome WNN — broad labels | {wnn_broad_labels and len(wnn_meta_broad['broad_cell_type'].dropna()):,} |\n\n"
            "## Cell-type label counts (broad)\n\n"
            "| Dataset | Unique broad labels |\n"
            "|---------|---------------------|\n"
            f"| RNA baseline | {len(rna_labels)} |\n"
            f"| Multiome WNN (broad) | {len(wnn_broad_labels)} |\n\n"
            "## Label overlap (broad)\n\n"
            f"**Shared labels ({len(broad_shared)}):**\n"
            f"{bullet_list(broad_shared)}\n\n"
            f"**RNA baseline only ({len(broad_rna_only)}):**\n"
            f"{bullet_list(broad_rna_only)}\n\n"
            f"**WNN broad only ({len(broad_wnn_only)}):**\n"
            f"{bullet_list(broad_wnn_only)}\n\n"
            "## Notes\n\n"
            "- WNN `cell_type_wnn` labels were harmonised to broad labels via "
            f"`{os.path.basename(args.label_map)}`.\n"
            "- Labels absent from the mapping file were kept unchanged.\n"
            "- Barcodes are not shared; no cell-level join is performed.\n"
        )

        broad_summary_path = os.path.join(args.outdir, "rna_vs_wnn_broad_summary.md")
        with open(broad_summary_path, "w") as fh:
            fh.write(broad_summary)

        print(f"Wrote: {broad_counts_path}")
        print(f"Wrote: {broad_summary_path}")

    # -- Report --------------------------------------------------------------
    print(f"Wrote: {counts_path}")
    print(f"Wrote: {fractions_path}")
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()
