import anndata as ad
import numpy as np
import pytest

from sc_rna_atac_regulon_benchmark.preprocessing import preprocess_rna


def test_returns_anndata(tiny_adata):
    result = preprocess_rna(tiny_adata)
    assert isinstance(result, ad.AnnData)


def test_obsm_pca_present(tiny_adata):
    result = preprocess_rna(tiny_adata, n_pcs=10)
    assert "X_pca" in result.obsm


def test_pca_shape(tiny_adata):
    n_pcs = 10
    result = preprocess_rna(tiny_adata, n_pcs=n_pcs)
    assert result.obsm["X_pca"].shape == (tiny_adata.n_obs, n_pcs)


def test_obs_preserved(tiny_adata):
    original_obs = tiny_adata.obs.copy()
    result = preprocess_rna(tiny_adata)
    assert list(result.obs.columns) == list(original_obs.columns)
    assert result.n_obs == original_obs.shape[0]


def test_var_preserved(tiny_adata):
    result = preprocess_rna(tiny_adata)
    assert result.n_vars == tiny_adata.n_vars
    assert list(result.var_names) == list(tiny_adata.var_names)


def test_zero_count_cell_no_crash(tiny_adata):
    tiny_adata.X[0, :] = 0.0
    result = preprocess_rna(tiny_adata)
    assert not np.any(np.isnan(result.obsm["X_pca"]))
