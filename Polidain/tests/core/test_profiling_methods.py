import pytest
import numpy as np
from pydantic import ValidationError as PydanticValidationError

from core.profiling_methods.polinmail.config import LocalPolinmailConfig, PolinmailConfig, ValidationError as PolinmailValidationError
from core.profiling_methods.polinmail.logic import PolinmailCalculator, get_matrix_coefficients, calculate_poly_coefficients

from core.profiling_methods.polidain.config import PolidainConfig, ValidationError as PolidainValidationError
from core.profiling_methods.polidain.logic import PolidainCalculator as PolidainCalcLogic, k_fun, c_fun

# --- Polinmail Tests ---

def test_local_polinmail_config_valid(local_polinmail_config):
    assert local_polinmail_config.m == 4
    assert len(local_polinmail_config.m_list) == 4

def test_local_polinmail_config_invalid():
    # m_list = [4, 5, 6, 7]
    # boundary_conditions = [-1, 0, 0, 0]
    # invalid if bc[i] == 0 and i in m_list.
    # index 1 is in m_list? No. m_list is values of powers.
    # Wait, check logic:
    # for i in range(len(bc)): if bc[i] == 0 and i in m_list: raise
    # m_list values are integers. 'i' is the index of boundary condition (derivative order).

    # If m=1, d=1. m_list=[1, 2, 3, 4] (len=4).
    # If bc=[0, 0, 0, 0].
    # i=0. 0 in [1,2,3,4]? No.
    # i=1. 1 in [1,2,3,4]? Yes. bc[1]=0. Raise.

    with pytest.raises(PolinmailValidationError):
        LocalPolinmailConfig(m=1, d=1, boundary_conditions=[0, 0, 0, 0])

def test_polinmail_config_resolve(local_polinmail_config):
    config = PolinmailConfig(config_1=local_polinmail_config)
    assert config.config_2 == local_polinmail_config
    assert config.config_3 == local_polinmail_config
    assert config.config_4 == local_polinmail_config

def test_polinmail_calculator_coefficients():
    # Simple case: y = x^2. m_list=[2]. boundary_conditions=[1] (at x=1, y=1)?
    # logic: coeffs = solve(A, B).
    # A constructed from m_list.
    # If m_list=[2]. A=[[1]]. B=[1]. coeff=[1]. y = 1*x^2.

    m_list = np.array([2])
    bc = np.array([1])
    coeffs = calculate_poly_coefficients(m_list, bc)
    assert np.allclose(coeffs, [1.0])

def test_polinmail_calculator_kinematics(polinmail_calculator):
    # Just check it runs and returns sensible types
    val = polinmail_calculator.h_phi(0.5, 1.0, 0.0, 1.0, 1)
    assert isinstance(val, float)

    val_v = polinmail_calculator.v_phi(0.5, 1.0, 0.0, 1.0, 1)
    assert isinstance(val_v, float)


# --- Polidain Tests ---

def test_polidain_config_valid():
    conf = PolidainConfig(m=3, d=4, k_1=6, k_2=6, k_3=6, k_4=6)
    assert conf.m == 3

def test_polidain_config_invalid_k():
    # k <= d
    with pytest.raises(PolidainValidationError):
        PolidainConfig(m=3, d=5, k_1=5, k_2=6, k_3=6, k_4=6) # k_1 <= d (5<=5)

def test_polidain_logic_functions():
    # Test k_fun
    # k_fun(p, m, d) -> [m, p, q, r, s]
    # q = 2p - d
    # r = 3p - 2d
    # s = 4p - 3d
    # if p=5, m=3, d=4
    # q = 10 - 4 = 6
    # r = 6 + 5 - 4 = 7
    # s = 7 + 5 - 4 = 8
    res = k_fun(5, 3, 4)
    assert res == [3, 5, 6, 7, 8]

def test_polidain_calculator_init():
    conf = PolidainConfig(m=3, d=4, k_1=6, k_2=6, k_3=6, k_4=6)
    calc = PolidainCalcLogic(conf)
    assert len(calc.c_list_1) == 5
    assert len(calc.k_list_1) == 5

    val = calc.h_phi(0.5, 1.0, 0.0, 1.0, 1)
    assert isinstance(val, float)
