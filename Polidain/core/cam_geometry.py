import numpy as np
from typing import Callable, List
from core.schemas import PolidainData, ProfilData


class CamProfileError(Exception):
    """Базовый класс для ошибок профиля кулачка"""
    pass

class PusherDiameterError(CamProfileError):
    """Ошибка: недостаточный диаметр толкателя"""
    pass

class ProfileSmoothnessError(CamProfileError):
    """Ошибка: негладкий профиль (подрез профиля)"""
    pass

class SolvePreliminaryCalculations(CamProfileError):
    """Ошибка: Не были проведены необходимые вычисления"""
    pass

def h_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    '''
    Функция радиуса кулачка от угла поворота кулачка
    :param fi: Угол поворота кулачка
    :param C_list:Массив коэффициентов при членах полиномма
    :param k_list:Массив степеней членов полинома
    :param fi_1:Конец характерного участка
    :param fi_0:Начало характерного участка
    :param h_kn_max: максимальная высота подёма кулачка на характерном участке
    '''
    temp = 1
    for i in range(0, len(C_list)):
        if k_list[i] - 0 < 0:
            continue
        temp += C_list[i] * (((fi - fi_0) / (fi_1 - fi_0)) ** k_list[i])
    return temp * h_kn_max


def v_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    '''Функция скорости от угла поворота кулачка
    :param fi: Угол поворота кулачка
    :param C_list:Массив коэффициентов при членах полиномма
    :param k_list:Массив степеней членов полинома
    :param fi_1:Конец характерного участка
    :param fi_0:Начало характерного участка
    :param h_kn_max: максимальная высота подёма кулачка на характерном участке
    '''
    temp = 0
    for i in range(0, len(C_list)):
        if k_list[i] - 1 < 0:
            continue
        temp += k_list[i] * C_list[i] * (((fi - fi_0) / (fi_1 - fi_0)) ** (k_list[i] - 1))
    # return temp * h_kn_max * (omega / (fi_1 - fi_0))
    return temp * h_kn_max


def a_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    '''Функция ускорения от угла поворота кулачка
    :param fi: Угол поворота кулачка
    :param C_list:Массив коэффициентов при членах полиномма
    :param k_list:Массив степеней членов полинома
    :param fi_1:Конец характерного участка
    :param fi_0:Начало характерного участка
    :param h_kn_max: максимальная высота подёма кулачка на характерном участке
    '''
    temp = 0
    for i in range(0, len(C_list)):
        if k_list[i] - 2 < 0:
            continue
        temp += k_list[i] * (k_list[i] - 1) * C_list[i] * (((fi - fi_0) / (fi_1 - fi_0)) ** (k_list[i] - 2))
    return temp * h_kn_max


def d_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    '''Функция рывка от угла поворота кулачка
    :param fi: Угол поворота кулачка
    :param C_list:Массив коэффициентов при членах полиномма
    :param k_list:Массив степеней членов полинома
    :param fi_1:Конец характерного участка
    :param fi_0:Начало характерного участка
    :param h_kn_max: максимальная высота подёма кулачка на характерном участке
    '''
    temp = 0
    for i in range(0, len(C_list)):
        if k_list[i] - 3 < 0:
            continue
        temp += k_list[i] * (k_list[i] - 1) * (k_list[i] - 2) * C_list[i] * (
                ((fi - fi_0) / (fi_1 - fi_0)) ** (k_list[i] - 3))
    return temp * h_kn_max


def k_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    '''Функция кракена от угла поворота кулачка
    :param fi: Угол поворота кулачка
    :param C_list:Массив коэффициентов при членах полиномма
    :param k_list:Массив степеней членов полинома
    :param fi_1:Конец характерного участка
    :param fi_0:Начало характерного участка
    :param h_kn_max: максимальная высота подёма кулачка на характерном участке
    '''
    temp = 0
    for i in range(0, len(C_list)):
        if k_list[i] - 4 < 0:
            continue
        temp += k_list[i] * (k_list[i] - 1) * (k_list[i] - 2) * (k_list[i] - 3) * C_list[i] * (
                ((fi - fi_0) / (fi_1 - fi_0)) ** (k_list[i] - 4))
    return temp * h_kn_max

def set_polidain_data(fun_list:list, N = 1000):
    fi_list = np.linspace(0, 2 * np.pi, N)
    H = fun_list[0](fi_list) * 1000
    V = fun_list[1](fi_list) * 1000
    A = fun_list[2](fi_list) * 1000
    D = fun_list[3](fi_list) * 1000
    K = fun_list[4](fi_list) * 1000
    fi_list = fi_list / np.pi * 180
    return PolidainData(H = H, V = V, A = A, D = D, K = K, fi_list = fi_list)

def set_profil_data(fun_list:list, N = 1000):
    fi_list = np.linspace(0, 2 * np.pi, N)
    X = fun_list[0](fi_list) * 1000
    Y = fun_list[1](fi_list) * 1000
    fi_list = fi_list / np.pi * 180
    return ProfilData(X = X, Y = Y, fi_list = fi_list)

def fi_list_dif(fi_array: np.ndarray, f_dif: float) -> np.ndarray:
    """
    Векторизированный сдвиг фазы
    """
    shifted = fi_array - f_dif
    return np.where(shifted <= 0, shifted + 2 * np.pi, shifted)

class Kulachok_polidain:
    def __init__(self, config):
        self.config = config
        self.kulachok_data = None
        self.tolkatel_data = None
        self.profil_data = None
        self.kulachok_solve_flag = False
        self.tolkatel_solve_flag = False
        self.profil_solve_flag = False

    def fun_universal(self, fi: float | np.ndarray, fun: Callable, sign_list: List[int],
                      const_list: List[float]) -> float | np.ndarray:
        '''
        Универсальная функция поддерживающая numpy, вычислчющая fun на всех участках
        '''

        fi_arr = np.asarray(fi)

        # Если на входе был скаляр, работаем как с 0-мерным массивом,
        is_scalar = fi_arr.ndim == 0
        if is_scalar:
            fi_arr = np.atleast_1d(fi_arr)

        # Инициализируем массив результатов нулями
        result = np.zeros_like(fi_arr, dtype=float)

        # Конфигурационные углы для краткости
        p1 = self.config.phi_1
        p2 = self.config.phi_2
        p3 = self.config.phi_3
        p4 = self.config.phi_4
        p5 = self.config.phi_5

        # --- Создаем маски для каждого участка ---
        # Участок 1: 0 <= fi < phi_1
        mask1 = (fi_arr >= 0) & (fi_arr < p1)

        # Участок 2: phi_1 <= fi < phi_2
        mask2 = (fi_arr >= p1) & (fi_arr < p2)

        # Участок 3: phi_2 <= fi < phi_3
        mask3 = (fi_arr >= p2) & (fi_arr < p3)

        # Участок 4: phi_3 <= fi < phi_4
        mask4 = (fi_arr >= p3) & (fi_arr < p4)

        # Участок 5: phi_4 <= fi < phi_5
        mask5 = (fi_arr >= p4) & (fi_arr < p5)

        # Остаток: phi_5 <= fi <= 2*pi
        mask6 = (fi_arr >= p5) & (fi_arr <= 2 * np.pi + 1e-6)  # +эпсилон для точности float

        # --- Вычисления для каждого участка ---
        # Участок 1
        if np.any(mask1):
            val = fun(fi_arr[mask1], self.config.C_list_1, self.config.k_list_1,
                      p1, fi_0=0, h_kn_max=self.config.z)
            result[mask1] = sign_list[0] * val + const_list[0]

        # Участок 2
        if np.any(mask2):
            val = fun(fi_arr[mask2], self.config.C_list_2, self.config.k_list_2,
                      p2, fi_0=p1, h_kn_max=self.config.h)
            result[mask2] = sign_list[1] * val + const_list[1]

        # Участок 3 (Выдержка - просто константа)
        if np.any(mask3):
            result[mask3] = const_list[2]

        # Участок 4 (С инверсией аргумента)
        if np.any(mask4):
            # Аргумент: phi_4 - fi + phi_3
            args = p4 - fi_arr[mask4] + p3
            val = fun(args, self.config.C_list_3, self.config.k_list_3,
                      p4, fi_0=p3, h_kn_max=self.config.h)
            result[mask4] = sign_list[2] * val + const_list[3]

        # Участок 5 (Спуск с инверсией)
        if np.any(mask5):
            args = p5 - fi_arr[mask5] + p4
            val = fun(args, self.config.C_list_4, self.config.k_list_4,
                      p5, fi_0=p4, h_kn_max=self.config.z)
            result[mask5] = sign_list[3] * val + const_list[4]

        # Остаток
        if np.any(mask6):
            result[mask6] = const_list[5]

        # Возвращаем скаляр, если на входе был скаляр, иначе массив
        if is_scalar:
            return result.item()
        return result

    def fun_h(self, fi):
        return self.fun_universal(fi, h_phi, [-1, -1, -1, -1], [self.config.r0,
                                                                self.config.r0 + self.config.h,
                                                                self.config.r0 + self.config.h,
                                                                self.config.r0 + self.config.h,
                                                                self.config.r0,
                                                                self.config.r0 - self.config.z])

    def fun_h_2(self, fi: float | np.ndarray):
        if type(fi) is np.ndarray:
            return (self.fun_universal(fi, h_phi, [-1, -1, -1, -1], [self.config.r0,
                                                                     self.config.r0 + self.config.h,
                                                                     self.config.r0 + self.config.h,
                                                                     self.config.r0 + self.config.h,
                                                                     self.config.r0,
                                                                     self.config.r0 - self.config.z]) - self.config.r0) * np.int64(
                fi >= self.config.phi_1) * np.int64(fi <= self.config.phi_4)
        return (self.fun_universal(fi, h_phi, [-1, -1, -1, -1], [self.config.r0,
                                                                 self.config.r0 + self.config.h,
                                                                 self.config.r0 + self.config.h,
                                                                 self.config.r0 + self.config.h,
                                                                 self.config.r0,
                                                                 self.config.r0 - self.config.z]) - self.config.r0) * int(
            fi >= self.config.phi_1) * int(fi <= self.config.phi_4)

    def fun_v(self, fi):
        return self.fun_universal(fi, v_phi, [-1, -1, 1, 1], [0, 0, 0, 0, 0, 0])

    def fun_v_2(self, fi):
        if type(fi) is np.ndarray:
            return self.fun_universal(fi, v_phi, [-1, -1, 1, 1], [0, 0, 0, 0, 0, 0]) * np.int64(
                fi >= self.config.phi_1) * np.int64(fi <= self.config.phi_4)
        return self.fun_universal(fi, v_phi, [-1, -1, 1, 1], [0, 0, 0, 0, 0, 0]) * int(fi >= self.config.phi_1) * int(
            fi <= self.config.phi_4)

    def fun_a(self, fi):
        return self.fun_universal(fi, a_phi, [-1, -1, -1, -1], [0, 0, 0, 0, 0, 0])

    def fun_a_2(self, fi):
        if type(fi) is np.ndarray:
            return self.fun_universal(fi, a_phi, [-1, -1, -1, -1], [0, 0, 0, 0, 0, 0]) * np.int64(
                fi >= self.config.phi_1) * np.int64(fi <= self.config.phi_4)
        return self.fun_universal(fi, v_phi, [-1, -1, 1, 1], [0, 0, 0, 0, 0, 0]) * int(fi >= self.config.phi_1) * int(
            fi <= self.config.phi_4)

    def fun_d(self, fi):
        return self.fun_universal(fi, d_phi, [-1, -1, 1, 1], [0, 0, 0, 0, 0, 0])

    def fun_d_2(self, fi):
        if type(fi) is np.ndarray:
            return self.fun_universal(fi, d_phi, [-1, -1, 1, 1], [0, 0, 0, 0, 0, 0]) * np.int64(
                fi >= self.config.phi_1) * np.int64(fi <= self.config.phi_4)
        return self.fun_universal(fi, v_phi, [-1, -1, 1, 1], [0, 0, 0, 0, 0, 0]) * int(fi >= self.config.phi_1) * int(
            fi <= self.config.phi_4)

    def fun_k(self, fi):
        return self.fun_universal(fi, k_phi, [-1, -1, -1, -1], [0, 0, 0, 0, 0, 0])

    def fun_k_2(self, fi):
        if type(fi) is np.ndarray:
            return self.fun_universal(fi, k_phi, [-1, -1, -1, -1], [0, 0, 0, 0, 0, 0]) * np.int64(
                fi >= self.config.phi_1) * np.int64(fi <= self.config.phi_4)
        return self.fun_universal(fi, v_phi, [-1, -1, 1, 1], [0, 0, 0, 0, 0, 0]) * int(fi >= self.config.phi_1) * int(
            fi <= self.config.phi_4)

    def fun_x(self, fi):
        return self.fun_h(fi) * np.cos(fi)

    def fun_y(self, fi):
        return self.fun_h(fi) * np.sin(fi)

    def set_kulachok_data(self, N = 1000):
        self.kulachok_data = set_polidain_data([self.fun_h, self.fun_v, self.fun_a, self.fun_d, self.fun_k], N = N)
        self.kulachok_solve_flag = True

    def set_tolkatel_data(self, N = 1000):
        self.tolkatel_data = set_polidain_data([self.fun_h_2, self.fun_v_2, self.fun_a_2, self.fun_d_2, self.fun_k_2], N = N)
        self.tolkatel_solve_flag = True

    def set_profil_data(self, N = 1000):
        self.profil_data = set_profil_data([self.fun_x, self.fun_y], N = N)
        self.profil_solve_flag = True

    def solve(self, N = 1000, kulachok_type = None):
        self.set_tolkatel_data(N = N)
        self.set_kulachok_data(N = N)
        if kulachok_type is None:
            self.set_profil_data(N = N)
        elif kulachok_type == 'flat':
            self.set_profil_flatpusher()
        else:
            raise ValueError('kulachok_type must be either "flat" or None')

    def profil_flatpusher_check(self, curvature_flag = None):
        if not(self.kulachok_solve_flag and self.tolkatel_solve_flag):
            raise SolvePreliminaryCalculations(f"Не были проведены предварительные вычисления закона движения толкатиля")
        max_v = np.max(self.tolkatel_data.V)
        if self.config.D_t * 1e3 <= max_v:
            raise PusherDiameterError(
                f"Недостаточный диаметр толкателя: {self.config.D_t * 1e3:.2f} <= {max_v:.2f}"
            )

        curvature_check = self.tolkatel_data.H + self.tolkatel_data.A + 500 * self.config.D
        min_curvature = np.min(curvature_check)
        if min_curvature <= 0 and curvature_flag:
            raise ProfileSmoothnessError(
                f"Негладкий профиль (min = {min_curvature:.4f}). "
                "Необходимо повысить минимальный радиус кривизны."
            )

    def set_profil_flatpusher(self, curvature_flag = False):
        self.profil_flatpusher_check(curvature_flag = curvature_flag)
        self.tolkatel_data.V = self.tolkatel_data.V
        E = self.tolkatel_data.V / self.config.omega
        fi_list = self.tolkatel_data.fi_list * np.pi / 180
        R = self.kulachok_data.H + 1e3 * self.config.D / 2
        X = R * np.cos(fi_list) - E * np.sin(fi_list)
        Y = R * np.sin(fi_list) + E * np.cos(fi_list)
        self.profil_data = ProfilData(fi_list=fi_list, X=X, Y=Y)