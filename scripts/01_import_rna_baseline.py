"""
01_import_rna_baseline.py

Standardise the Kang et al. RNA baseline score files and write four output artefacts:
  - results/baseline/rna_cell_metadata.csv
  - results/baseline/rna_ifn_scores.csv
  - results/baseline/rna_program_scores.csv
  - results/baseline/rna_baseline_summary.md
"""

import argparse
import os
import sys

import pandas as pd


# ---------------------------------------------------------------------------
# Column mappings
# ---------------------------------------------------------------------------

IFN_COL_MAP = {
    "ifn_scanpy": "scanpy_ifn_score",
    "ifn_meanexpr": "mean_ifn_score",
    "ifn_rank": "rank_ifn_score",
}

PROG_COL_MAP = {
    "prog_interferon_alpha_beta_response": "interferon_alpha_beta_response",
    "prog_inflammatory_response": "inflammatory_response",
    "prog_antigen_presentation_mhc": "antigen_presentation_mhc",
    "prog_cytotoxic_nk_like": "cytotoxic_nk_like",
    "prog_monocyte_activation": "monocyte_activation",
}

IFN_REQUIRED = {"cell_id", "condition"} | set(IFN_COL_MAP.keys())
PROG_REQUIRED = {"cell_id", "condition", "cell_type"} | set(PROG_COL_MAP.keys())


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def validate_columns(df: pd.DataFrame, required: set, path: str) -> None:
    missing = required - set(df.columns)
    if missing:
        sys.exit(
            f"ERROR: {path} is missing required columns: {sorted(missing)}"
        )


def load_and_validate(path: str, required: set) -> pd.DataFrame:
    if not os.path.isfile(path):
        sys.exit(f"ERROR: input file not found: {path}")
    df = pd.read_csv(path)
    validate_columns(df, required, path)
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Import and standardise RNA baseline scores from Kang et al."
    )
    parser.add_argument(
        "--ifn-scores",
        default="data/external/kang_interferon_scores.csv",
        help="Path to kang_interferon_scores.csv",
    )
    parser.add_argument(
        "--program-scores",
        default="data/external/kang_program_scores.csv",
        help="Path to kang_program_scores.csv",
    )
    parser.add_argument(
        "--outdir",
        default="results/baseline",
        help="Output directory (created if absent)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # -- Load inputs ---------------------------------------------------------
    ifn_df = load_and_validate(args.ifn_scores, IFN_REQUIRED)
    prog_df = load_and_validate(args.program_scores, PROG_REQUIRED)

    # -- Safety checks (before any output is written) ------------------------

    # 1. Duplicated cell_id in program scores
    dup_mask = prog_df["cell_id"].duplicated()
    if dup_mask.any():
        n_dup = dup_mask.sum()
        sys.exit(
            f"ERROR: kang_program_scores.csv has {n_dup} duplicated cell_id value(s)."
        )

    # 3. Condition consistency across overlapping cell_ids
    overlap = set(ifn_df["cell_id"]) & set(prog_df["cell_id"])
    if overlap:
        ifn_cond = ifn_df.set_index("cell_id")["condition"]
        prog_cond = prog_df.set_index("cell_id")["condition"]
        overlap_list = sorted(overlap)
        mismatches = [
            cid for cid in overlap_list if ifn_cond[cid] != prog_cond[cid]
        ]
        if mismatches:
            sys.exit(
                f"ERROR: {len(mismatches)} cell_id(s) have mismatched condition "
                f"between IFN and program tables. First offender: '{mismatches[0]}' "
                f"(ifn={ifn_cond[mismatches[0]]!r}, prog={prog_cond[mismatches[0]]!r})."
            )

    os.makedirs(args.outdir, exist_ok=True)

    # -- rna_program_scores --------------------------------------------------
    prog_out = prog_df[["cell_id", "condition", "cell_type"] + list(PROG_COL_MAP.keys())].copy()
    prog_out.rename(columns=PROG_COL_MAP, inplace=True)
    prog_path = os.path.join(args.outdir, "rna_program_scores.csv")
    prog_out.to_csv(prog_path, index=False)

    # -- rna_cell_metadata ---------------------------------------------------
    meta_out = prog_df[["cell_id", "condition", "cell_type"]].copy()
    meta_path = os.path.join(args.outdir, "rna_cell_metadata.csv")
    meta_out.to_csv(meta_path, index=False)

    # -- rna_ifn_scores (join cell_type from program scores) -----------------
    cell_type_map = prog_df.set_index("cell_id")["cell_type"]
    ifn_out = ifn_df[["cell_id", "condition"] + list(IFN_COL_MAP.keys())].copy()
    ifn_out.rename(columns=IFN_COL_MAP, inplace=True)
    ifn_out.insert(2, "cell_type", ifn_out["cell_id"].map(cell_type_map))

    # 2. Missing cell_type after join
    n_missing = ifn_out["cell_type"].isna().sum()
    if n_missing:
        sys.exit(
            f"ERROR: {n_missing} cell(s) in IFN scores have no matching cell_type "
            f"in program scores."
        )

    ifn_path = os.path.join(args.outdir, "rna_ifn_scores.csv")
    ifn_out.to_csv(ifn_path, index=False)

    # -- rna_baseline_summary.md ---------------------------------------------
    n_cells = len(meta_out)
    n_cell_types = meta_out["cell_type"].nunique()
    condition_counts = meta_out["condition"].value_counts().sort_index()
    program_cols = list(PROG_COL_MAP.values())

    condition_lines = "\n".join(
        f"  - {cond}: {cnt}" for cond, cnt in condition_counts.items()
    )
    program_lines = "\n".join(f"  - {col}" for col in program_cols)

    summary = (
        "# RNA Baseline Summary\n\n"
        f"- **Number of cells:** {n_cells}\n"
        f"- **Number of cell types:** {n_cell_types}\n"
        "- **Counts by condition:**\n"
        f"{condition_lines}\n"
        "- **Program score columns:**\n"
        f"{program_lines}\n"
    )

    summary_path = os.path.join(args.outdir, "rna_baseline_summary.md")
    with open(summary_path, "w") as fh:
        fh.write(summary)

    print(f"Wrote: {meta_path}")
    print(f"Wrote: {ifn_path}")
    print(f"Wrote: {prog_path}")
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()
