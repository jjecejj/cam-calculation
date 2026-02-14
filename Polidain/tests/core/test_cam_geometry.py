import pytest
import numpy as np
from pydantic import ValidationError as PydanticValidationError
from core.cam_geometry.config import KulachokConfig, ValidationError
from core.cam_geometry.logic import Kulachok, PusherDiameterError, ProfileSmoothnessError, SolvePreliminaryCalculations
from core.cam_geometry.options import CamGeometryOptions, resolve_cam_geometry_options

# --- Tests for KulachokConfig ---

def test_kulachok_config_valid(kulachok_config):
    # Should not raise
    assert kulachok_config.N_k == 1000

def test_kulachok_config_invalid_phases():
    # Sum of phases > 2*pi
    # The custom ValidationError is raised directly, not wrapped in PydanticValidationError in this setup?
    # Or maybe because it's mode='after' and raises Exception?
    # Actually, let's catch both just in case, or specifically the custom one if we saw it in logs.
    # The log showed core.cam_geometry.config.ValidationError.

    with pytest.raises(ValidationError) as excinfo:
        KulachokConfig(
            N_k=1000, D=0.03, h=0.01, z=0.001,
            f_pod=np.pi, f_v=np.pi, f_op=np.pi, f_z=0.1, # Sum > 2pi
            R_r=0.005, D_t=0.032 # Added D_t as it might be required or have default
        )
    # Check if the error message is from our validator
    assert "Ошибка в значениях фаз" in str(excinfo.value)

def test_kulachok_config_properties(kulachok_config):
    assert kulachok_config.phi_0 == 0
    assert kulachok_config.phi_1 == kulachok_config.f_z
    assert kulachok_config.omega > 0
    assert kulachok_config.T > 0
    assert kulachok_config.r0 == kulachok_config.D / 2

# --- Tests for Kulachok Logic ---

def test_kulachok_initialization(kulachok):
    assert kulachok.kulachok_data is None
    assert kulachok.tolkatel_data is None
    assert kulachok.profile_data is None
    assert kulachok.kulachok_solve_flag is False

def test_kulachok_fun_universal(kulachok):
    # Test scalar input
    val = kulachok.fun_h(0.1)
    # Check if it returns a scalar (float or numpy scalar)
    assert np.isscalar(val) or isinstance(val, (float, int))

    # Test array input
    arr = np.array([0.1, 0.2])
    val_arr = kulachok.fun_h(arr)
    assert isinstance(val_arr, np.ndarray)
    assert len(val_arr) == 2

def test_kulachok_solve_thin(kulachok):
    kulachok.solve(kulachok_type='thin', N=100)
    assert kulachok.kulachok_solve_flag
    assert kulachok.tolkatel_solve_flag
    assert kulachok.profile_solve_flag
    assert kulachok.tolkatel_solve_type == 'thin'
    assert kulachok.kulachok_data is not None
    assert kulachok.tolkatel_data is not None
    assert kulachok.profile_data is not None

def test_kulachok_solve_flat(kulachok):
    # Might raise errors depending on config
    # For default config it should pass
    kulachok.solve(kulachok_type='flat', N=100)
    assert kulachok.tolkatel_solve_type == 'flat'

def test_kulachok_solve_roller(kulachok):
    kulachok.solve(kulachok_type='roller', N=100)
    assert kulachok.tolkatel_solve_type == 'roller'

def test_profile_flat_check_raises_pusher_diameter_error(kulachok):
    # Force max velocity to be high to trigger PusherDiameterError
    # We can do this by mocking tolkatel_data or using a config that causes high velocity

    # Let's manually set data to trigger error
    kulachok.solve(kulachok_type='thin', N=100)

    # If we set D_t very small, it should fail.
    kulachok.config.D_t = 1.0 * 1e-3 # 1mm diameter

    # Re-run check
    with pytest.raises(PusherDiameterError):
        kulachok.profile_flat_check()

def test_profile_flat_check_raises_smoothness_error(kulachok):
    kulachok.solve(kulachok_type='thin', N=100)

    # H_rad + A_rad + 500*D
    # We need min_curvature <= 0
    # Let's inject negative values into curvature check

    # H_rad is property of GraphData? No, it's a field.
    # We can overwrite the array.
    kulachok.tolkatel_data.H_rad = np.array([-1000.0] * 100) # Very negative
    kulachok.tolkatel_data.A_rad = np.array([0.0] * 100)

    with pytest.raises(ProfileSmoothnessError):
        kulachok.profile_flat_check(curvature_flag=True)

def test_profile_roller_check_raises_smoothness_error(kulachok):
    kulachok.solve(kulachok_type='thin', N=100)

    # To force rho_c < 0 (and thus rho_p < 0), we can make the denominator negative.
    # den = Rc^2 + 2v^2 - Rc * a
    # If we make 'a' very large positive, Rc*a > Rc^2 + 2v^2

    kulachok.tolkatel_data.A_rad = np.array([100000.0] * 100) # Huge acceleration

    with pytest.raises(ProfileSmoothnessError):
        kulachok.profile_roller_check()

def test_preliminary_calculations_error(kulachok):
    # Reset flags
    kulachok.kulachok_solve_flag = False
    kulachok.tolkatel_solve_flag = False # Ensure both are false

    with pytest.raises(SolvePreliminaryCalculations):
        kulachok.profile_flat_check()

# --- Tests for CamGeometryOptions ---

def test_cam_geometry_options_defaults():
    opts = CamGeometryOptions()
    assert opts.calculator_type == 'polidain'
    assert opts.kulachok_type == 'thin'
    assert opts.calculator is not None # resolved by validator

def test_resolve_cam_geometry_options(kulachok_config):
    # Note: we need to handle the calculator creation which might fail if dependencies (polidain config) are default
    # but CamGeometryOptions creates default calculator if not provided.

    # Wait, CamGeometryOptions default calculator_config depends on type.
    # If type is polidain, it uses PolidainConfig.
    # If type is polinmail, it uses PolinmailConfig.

    # Let's use polinmail to avoid polidain complexity if any.
    opts = CamGeometryOptions(
        cam_config=kulachok_config,
        kulachok_type='thin',
        N=50,
        calculator_type='polinmail'
    )
    # The validator 'resolve_calculator_config' sets default PolinmailConfig
    # The validator 'resolve_calculator' creates PolinmailCalculator

    kulachok, angle = resolve_cam_geometry_options(opts)
    assert kulachok.kulachok_solve_flag
    assert kulachok.tolkatel_solve_type == 'thin'
    assert isinstance(angle, (float, int, np.floating))
