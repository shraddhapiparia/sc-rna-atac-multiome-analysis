import numpy as np
import pandas as pd
import pytest
import anndata as ad


@pytest.fixture
def tiny_adata():
    rng = np.random.default_rng(0)
    n_obs, n_vars = 50, 100
    X = rng.poisson(5, size=(n_obs, n_vars)).astype(float)
    obs = pd.DataFrame(
        {"cell_type": ["TypeA"] * 25 + ["TypeB"] * 25},
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
    return ad.AnnData(X=X, obs=obs, var=var)
