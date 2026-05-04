# Setup notes

This document contains detailed environment setup and troubleshooting notes for the `sc-rna-atac-regulon-benchmark` repository.

The main README is intentionally concise. Use this file when setting up the full Python/R/Bioconductor environment or troubleshooting package installation issues.

## Environment strategy

This project uses both Python and R.

- `environment.yml` defines the full conda environment.
- `pyproject.toml` defines the Python package metadata, editable install behavior, and Python test configuration.
- `renv.lock` records the R package state used during development.

Use `environment.yml` as the primary setup file because the project depends on compiled R and Bioconductor packages, including Seurat, Signac, chromVAR, JASPAR2020, and genome annotation packages.

## Recommended setup

Install Miniconda or Anaconda first.

Install `mamba` in the base environment:

    conda install -n base -c conda-forge mamba

Create the project environment:

    mamba env create -f environment.yml

Activate it:

    conda activate sc-rna-atac

Install the Python package in editable mode:

    pip install -e .

Run tests:

    pytest

## Updating an existing environment

After pulling repository changes, update the environment with:

    mamba env update -f environment.yml --prune

Then reinstall the local Python package if needed:

    pip install -e .

## R version note

This project is configured for R 4.4 through conda.

R 4.3 is not recommended. Dependency conflicts were observed between chromVAR, sparseMatrixStats, BSgenome, and related Bioconductor packages under R 4.3.

Use the R version specified in `environment.yml`.

## Verifying R packages

Run:

    Rscript scripts/install_r_deps.R

The script checks whether the expected R packages are available. Any missing package should be resolved before running the full workflow.

If conda did not install a package successfully, first try:

    mamba env update -f environment.yml --prune

Then rerun:

    Rscript scripts/install_r_deps.R

## Slingshot note

`bioconductor-slingshot` may not resolve cleanly through conda on all platforms, especially macOS Apple Silicon, because some dependencies may be unavailable for `osx-arm64`.

If Slingshot is missing, `scripts/install_r_deps.R` can attempt installation through BiocManager.

This only affects the optional pseudotime phase. The core RNA+ATAC preprocessing, WNN integration, QC, and TF/motif workflows can still be run without pseudotime.

## renv usage

The repository includes `renv.lock` to record the R package state used during development.

After creating the conda environment, `renv` can be used to restore the R package snapshot:

    Rscript -e "renv::restore()"

If you intentionally update R packages, refresh the lockfile from inside R:

    library(renv)
    renv::snapshot()

Only commit `renv.lock` when the package changes are intentional.

## Python package verification

Check that the Python package imports:

    python -c "import sc_rna_atac_regulon_benchmark; print('Package import OK')"

Check core Python dependencies:

    python -c "import pandas, numpy, matplotlib, seaborn, scipy, statsmodels; print('Python dependencies OK')"

Run tests:

    pytest

## Platform notes

| Platform | Notes |
|---|---|
| macOS Apple Silicon | Use the conda environment. Avoid mixing system R with conda R. |
| Linux | Recommended for full reproducibility and larger runs. |
| Windows | WSL2 Ubuntu is recommended. Native Windows R/conda behavior is not tested. |

## Common troubleshooting

### Conda solve is slow

Use `mamba` instead of `conda`:

    mamba env create -f environment.yml

or:

    mamba env update -f environment.yml --prune

### R package is missing

First update the conda environment:

    mamba env update -f environment.yml --prune

Then rerun:

    Rscript scripts/install_r_deps.R

### Wrong R version is being used

Check:

    which R
    R --version

The path should point inside the active conda environment.

If not, reactivate the environment:

    conda activate sc-rna-atac

### Python package is not importable

Install the local package in editable mode:

    pip install -e .

Then check:

    python -c "import sc_rna_atac_regulon_benchmark; print('OK')"

### renv appears out of sync

Check status from R:

    renv::status()

If the mismatch is expected after an intentional package update, run:

    renv::snapshot()

If the mismatch is local machine noise, avoid committing unrelated `renv.lock` changes.

## Suggested setup check before running the full pipeline

Run these commands before starting the analysis:

    conda activate sc-rna-atac
    python -c "import sc_rna_atac_regulon_benchmark; print('Python package OK')"
    pytest
    Rscript scripts/install_r_deps.R

Once these pass, continue with the workflow in the main README.