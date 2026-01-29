from core.cam_geometry import Kulachok_polidain, set_polidain_data
from exporters.dxf_creator import build_profile
from config import AppConfig
from core.schemas import PolidainConfig
from vizualization.plotter import set_config, display_graphs_kulachok, display_graphs_tolkatel, display_profil, display_all, display_all_comprasion, calculate_optimal_angle
from math import pi
from vizualization.rotate_animation import display_animation, set_rotate_data
from math import sqrt

if __name__ == '__main__':
    # Исходные данные
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

    # Решение кулачка
    kulachok = Kulachok_polidain(appConfig.cam)
    kulachok.solve(kulachok_type='flat', N=100)

    # Построение графиков
    set_config(appConfig.plot)
    initial_angle = calculate_optimal_angle(kulachok)
    #display_graphs_tolkatel(kulachok.kulachok_data, initial_angle=initial_angle)
    display_profil(kulachok.profil_data, initial_angle=initial_angle)
    display_animation(set_rotate_data(kulachok, tolkatel_type=kulachok.solve_type), save_flag=True)

    # Импорт геометрии
    #build_profile(kulachok.profil_data)