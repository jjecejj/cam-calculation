from config import AppConfig
from core.schemas import PolidainConfig
from visualization.plotter import set_config
import numpy as np
from core.options import calculate, CamSolveOptions


if __name__ == '__main__':
    cam = PolidainConfig(
        N_k=1000,
        D=30.0 * 1e-3,
        D_t=30.0 * 1e-3,
        h=2.0 * 1e-3,
        z=0.25 * 1e-3,
        f_pod=80.0 / 180 * np.pi,
        f_v=5.0 / 180 * np.pi,
        f_op=75.0 / 180 * np.pi,
        f_z=25 / 180 * np.pi,
        m=3,
        d=4,
        k_1=6,
        k_2=6,
        k_3=6,
        k_4=6,
        R_r=5 * 1e-3,
)
    appConfig = AppConfig(cam = cam)
    set_config(appConfig.plot)
    cam_solve_options = CamSolveOptions(cam = cam,
                                        graphs_tolkatel_flag = True,
                                        graphs_kulachok_flag = False,
                                        graphs_profil_flag = True,
                                        display_animation_flag = False,
                                        save_animation_flag = True,
                                        import_dxf_flag = False,
                                        dxf_profil_name = "tolkatel",
                                        calculate_optimal_initial_angle = True,
                                        graphs_argument_type = 'degree',
                                        kulachok_type = 'roller',
                                        N = 1000,
                                        profil_and_graphs_together_flag = True,
                                        animation_profil_and_graphs_together_flag = True,
                                        animation_pause_flag = True,
    )
    calculate(cam_solve_options)
