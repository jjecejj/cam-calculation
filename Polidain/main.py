from core.cam_geometry import Kulachok_polidain, set_polidain_data
from exporters.dxf_creator import build_profile
from config import AppConfig
from core.schemas import PolidainConfig
from vizualization.plotter import set_config, display_graphs_kulachok, display_graphs_tolkatel, display_profil, display_all, display_all_comprasion, calculate_optimal_angle
from math import pi
from vizualization.rotate_animation import display_animation, set_rotate_data
from core.options import calculate, CamSolveOptions

if __name__ == '__main__':
    cam = PolidainConfig(
            N_k = 1000,
            D = 30.0 * 1e-3,
            D_t = 32.0 * 1e-3,
            h = 12.0 * 1e-3,
            z = 0.25 * 1e-3,
            f_pod = 80.0 / 180 * pi,
            f_v = 5.0 / 180 * pi,
            f_op = 75.0 / 180 * pi,
            f_z = 20 / 180 * pi,
            m = 3,
            d = 12,
            k_1 = 20,
            k_2 = 20,
            k_3 = 20,
            k_4 = 20
)
    appConfig = AppConfig(cam = cam)
    set_config(appConfig.plot)
    cam_solve_options = CamSolveOptions(cam = cam,
                                        graphs_tolkatel_flag = True,
                                        graphs_kulachok_flag = False,
                                        graphs_profil_flag = True,
                                        display_animation_flag = False,
                                        save_animation_flag = False,
                                        import_dxf_flag = False,
                                        dxf_profil_name = "kulachok"
    )
    calculate(cam_solve_options)
