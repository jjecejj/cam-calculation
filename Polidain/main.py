import numpy as np

from core.cam_geometry import Kulachok
from core.cam_geometry.config import KulachokConfig
from core.cam_geometry.options import CamGeometryOptions
from core.profiling_methods.polinmail.config import PolinmailConfig, LocalPolinmailConfig
from core.profiling_methods.polinmail.logic import PolinmailCalculator
from options import CamSolveOptions, calculate_cam_solve
from vizualization.plotter.logic import display_graphs_kulachok, display_graphs_tolkatel, display_profile
from vizualization.plotter.options import PlotterOptions
from vizualization.rotate_animation.logic import display_animation
from vizualization.rotate_animation.options import RotateAnimationOptions

if __name__ == '__main__':
    config = KulachokConfig(
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
    local_config = LocalPolinmailConfig(m = 4, d = 1, boundary_conditions = [-1, 0, 0, 0])
    calculator_config =  PolinmailConfig(
        config_1=local_config,
    )
    calculator = PolinmailCalculator(calculator_config)

    cam_geometry_options = CamGeometryOptions(
        cam_config = config,
        calculator = calculator,
        N = 5000,
        calculate_optimal_initial_angle = True,
        kulachok_type='thin'
    )
    plotter_options = PlotterOptions(
        graphs_kulachok_flag = True,
        graphs_tolkatel_flag = True,
        graphs_profile_flag = True,
        graphs_argument_type = 't',
    )
    rotate_animation_options = RotateAnimationOptions(
        display_animation_flag = False,
        save_animation_flag = False
    )
    cso = CamSolveOptions(cam_geometry_options = cam_geometry_options,
                          plotter_options = plotter_options,
                          rotate_animation_options = rotate_animation_options)
    calculate_cam_solve(cso)
