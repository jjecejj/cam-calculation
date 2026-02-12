import numpy as np

from core.cam_geometry import Kulachok
from core.cam_geometry.config import KulachokConfig
from core.profiling_methods.polinmail.config import default_polinmail_config, PolinmailConfig, LocalPolinmailConfig
from core.profiling_methods.polinmail.logic import PolinmailCalculator
from options import CamSolveOptions, calculate_cam_solve
from vizualization.plotter.logic import display_graphs_kulachok, display_graphs_tolkatel, display_profil
from vizualization.rotate_animation.logic import display_animation

if __name__ == '__main__':
    config = KulachokConfig(
        N_k=1000,
        D=30.0 * 1e-3,
        D_t=30.0 * 1e-3,
        h=12.0 * 1e-3,
        z=0.25 * 1e-3,
        f_pod=100.0 / 180 * np.pi,
        f_v=5.0 / 180 * np.pi,
        f_op=100.0 / 180 * np.pi,
        f_z=25 / 180 * np.pi,
        R_r=5 * 1e-3,
)
    locat_config = LocalPolinmailConfig(m = 3, d = 1, boundary_conditions = [-1, 0, 0])
    calculator_config =  PolinmailConfig(
        config_1=locat_config,
        config_2=locat_config,
        config_3=locat_config,
        config_4=locat_config,
    )
    calculator = PolinmailCalculator(calculator_config)
    kulachok = Kulachok(config, calculator)
    kulachok.solve(kulachok_type='flat', N = 200)
    type = 't'
    display_graphs_kulachok(kulachok.kulachok_data, graphs_type=type)
    display_graphs_tolkatel(kulachok.tolkatel_data, graphs_type=type)
    display_profil(kulachok.profile_data)
    display_animation(kulachok, save_flag = True)
    exit()
    cso = CamSolveOptions(kulachok_type = 'thin',
                          N = 1000,
                          graphs_kulachok_flag = True,
                          graphs_profile_flag = True)
    calculate_cam_solve(cso)
