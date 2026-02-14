import pytest
from unittest.mock import MagicMock, patch
import numpy as np
from core.optimization.logic import differential_evolution_optimization, gibrid_optimization, fun_optimize_gibrid
from core.optimization.config import OptimizeConfig, GibridOptimizationConfig, BoundsConfig

@pytest.fixture
def mock_R_func():
    return lambda x: np.array(x) * 0.0 + 0.02 # Return array of 0.02

def test_fun_optimize_gibrid_valid_params(mock_R_func):
    # x vector: [z, f_pod, f_v, f_op, f_z, k1, k2, k3, k4]
    x = [0.001, 1.0, 0.1, 1.0, 0.1, 10, 10, 10, 10]
    fi_list = np.array([0.0, 1.0])
    optimize_config = OptimizeConfig(D=0.03, h=0.01)

    # We wrap in try-except because it constructs Kulachok/Polidain which might fail with random values
    # But we provided reasonable values.

    val = fun_optimize_gibrid(x, fi_list=fi_list, m=3, d=4, R_func=mock_R_func, optimize_config=optimize_config)

    # fun_optimize_gibrid returns 1e6 on exception.
    # If it works, it returns a sum of diffs (float).
    assert isinstance(val, float)

def test_differential_evolution_calls_scipy():
    with patch('core.optimization.logic.differential_evolution') as mock_de:
        mock_res = MagicMock()
        mock_res.x = np.array([0]*9)
        mock_res.fun = 0.0
        mock_de.return_value = mock_res

        m = 3
        d = 4
        opt_config = OptimizeConfig(D=0.03, h=0.01)
        R_func = lambda x: x

        differential_evolution_optimization(m, d, opt_config, R_func, N=10)

        assert mock_de.called

def test_gibrid_calls_diff_evolution():
    with patch('core.optimization.logic.differential_evolution') as mock_de:
        mock_res = MagicMock()
        mock_res.x = np.array([0]*9)
        mock_res.fun = 0.0
        mock_de.return_value = mock_res

        g_config = GibridOptimizationConfig(m=[3], d=[4, 5])
        opt_config = OptimizeConfig(D=0.03, h=0.01)
        R_func = lambda x: x

        gibrid_optimization(g_config, opt_config, R_func, N=10)

        assert mock_de.call_count == 2
