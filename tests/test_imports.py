import importlib
import pytest

CORE_DEPS = [
    "numpy",
    "pandas",
    "matplotlib",
    "seaborn",
    "scipy",
    "statsmodels",
]


@pytest.mark.parametrize("module", CORE_DEPS)
def test_import(module):
    importlib.import_module(module)
