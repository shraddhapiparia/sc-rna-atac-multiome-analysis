import numpy as np
import anndata as ad


def preprocess_rna(adata: ad.AnnData, n_pcs: int = 10) -> ad.AnnData:
    """Log-normalize counts and compute PCA, storing result in obsm['X_pca']."""
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = X.astype(float)

    lib_sizes = X.sum(axis=1, keepdims=True)
    lib_sizes[lib_sizes == 0] = 1.0
    X_norm = np.log1p(X / lib_sizes * 1e4)

    X_centered = X_norm - X_norm.mean(axis=0)
    n_components = min(n_pcs, adata.n_obs - 1, adata.n_vars)
    U, S, _ = np.linalg.svd(X_centered, full_matrices=False)
    adata.obsm["X_pca"] = U[:, :n_components] * S[:n_components]
    return adata
