import sys
import os
import pytest
import numpy as np

# Add Polidain directory to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from core.cam_geometry.config import KulachokConfig
from core.profiling_methods.polinmail.config import PolinmailConfig, LocalPolinmailConfig
from core.profiling_methods.polinmail.logic import PolinmailCalculator
from core.cam_geometry.logic import Kulachok

@pytest.fixture
def kulachok_config():
    return KulachokConfig(
        N_k=1000,
        D=30.0 * 1e-3,
        D_t=32.0 * 1e-3,
        h=9.0 * 1e-3,
        z=0.25 * 1e-3,
        f_pod=100.0 / 180 * np.pi,
        f_v=5.0 / 180 * np.pi,
        f_op=100.0 / 180 * np.pi,
        f_z=25 / 180 * np.pi,
        R_r=5 * 1e-3,
    )

@pytest.fixture
def local_polinmail_config():
    return LocalPolinmailConfig(
        m=4,
        d=1,
        boundary_conditions=[-1, 0, 0, 0]
    )

@pytest.fixture
def polinmail_config(local_polinmail_config):
    return PolinmailConfig(
        config_1=local_polinmail_config,
        config_2=local_polinmail_config,
        config_3=local_polinmail_config,
        config_4=local_polinmail_config
    )

@pytest.fixture
def polinmail_calculator(polinmail_config):
    return PolinmailCalculator(polinmail_config)

@pytest.fixture
def kulachok(kulachok_config, polinmail_calculator):
    return Kulachok(kulachok_config, polinmail_calculator)
