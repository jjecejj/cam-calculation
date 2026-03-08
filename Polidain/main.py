"""
Модуль для запуска демонстрационного расчета кулачка через код.
Пример инициализации конфигурации и запуска решателя.
"""
import numpy as np

from core.cam_geometry import Kulachok
from core.cam_geometry.config import KulachokConfig
from core.cam_geometry.options import CamGeometryOptions
from core.optimization import OptimizeConfig
from core.optimization.logic import differential_evolution_optimization, get_calculator_config_from_x
from core.profiling_methods.polidain.config import PolidainConfig
from core.profiling_methods.polinmail.config import PolinmailConfig, LocalPolinmailConfig
from core.profiling_methods import PolinmailCalculator, PolidainCalculator
from core.profils import ProfileDataExtractor
from options import CamSolveOptions, calculate_cam_solve
from vizualization.plotter import display_graphs_kulachok, display_graphs_tolkatel, display_profile
from vizualization.plotter.logic import display_profile_multiplicity
from vizualization.plotter.options import PlotterOptions
from vizualization.rotate_animation import display_animation
from vizualization.rotate_animation.options import RotateAnimationOptions
from core.profils.pyzirev import data as data_pyzirev, Polidain_x
from core.profils.pyzirev import f_dif as f_dif_pyzirev
from core.profils.pyzirev import D as D_pyzirev
from core.profils.pyzirev import h as h_pyzirev
if __name__ == '__main__':
    '''opt = OptimizeConfig(
        profiling_type = 'polidain',
        D = D_pyzirev,
        h = h_pyzirev,
        R_func=ProfileDataExtractor(data_pyzirev, f_dif_pyzirev).get_H()
    )
    differential_evolution_optimization(opt, N = 2000)
    exit()'''
    config, calculator_config = get_calculator_config_from_x(Polidain_x, D_pyzirev, h_pyzirev, 'polidain')
    calculator = PolidainCalculator(calculator_config)
    kulachok = Kulachok(config = config, profile_method_calculator = calculator)
    kulachok.solve()
    display_profile_multiplicity([kulachok.profile_data, ProfileDataExtractor(data_pyzirev, f_dif_pyzirev).get_ProfileData()])

    exit()
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
    local_config = LocalPolinmailConfig(m = 4, d = 1, boundary_conditions = [1, 0, 0, 0])
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
