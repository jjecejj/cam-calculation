import numpy as np
from typing import Callable, List, Literal
from config.kulachok import KulachokConfig
from core.profiling_methods.base import BaseCalculator
from core.schemas import GraphData, ProfileData, set_graph_data, set_profile_data
from scipy.interpolate import interp1d

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

class Kulachok:
    def __init__(self, config: KulachokConfig, profile_method_calculator: BaseCalculator):
        self.config = config
        self.profile_method_calculator = profile_method_calculator

        self.kulachok_data: GraphData | None = None
        self.tolkatel_data: GraphData | None = None
        self.profile_data: ProfileData | None = None

        self.kulachok_solve_flag: bool = False
        self.tolkatel_solve_flag: bool = False
        self.profile_solve_flag: bool = False

        self.solve_type: str = None

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
            val = fun(fi_arr[mask1], p1, fi_0 = 0, h_kn_max = self.config.z, segment_number = 1)
            result[mask1] = sign_list[0] * val + const_list[0]

        # Участок 2
        if np.any(mask2):
            val = fun(fi_arr[mask2], p2, fi_0=p1, h_kn_max=self.config.h, segment_number = 2)
            result[mask2] = sign_list[1] * val + const_list[1]

        # Участок 3 (Выдержка - просто константа)
        if np.any(mask3):
            result[mask3] = const_list[2]

        # Участок 4 (С инверсией аргумента)
        if np.any(mask4):
            # Аргумент: phi_4 - fi + phi_3
            args = p4 - fi_arr[mask4] + p3
            val = fun(args, p4, fi_0=p3, h_kn_max=self.config.h, segment_number = 3)
            result[mask4] = sign_list[2] * val + const_list[3]

        # Участок 5 (Спуск с инверсией)
        if np.any(mask5):
            args = p5 - fi_arr[mask5] + p4
            val = fun(args, p5, fi_0=p4, h_kn_max=self.config.z, segment_number = 4)
            result[mask5] = sign_list[3] * val + const_list[4]

        # Остаток
        if np.any(mask6):
            result[mask6] = const_list[5]

        # Возвращаем скаляр, если на входе был скаляр, иначе массив
        if is_scalar:
            return result.item()
        return result

    def fun_h(self, fi: float | np.ndarray):
        return self.fun_universal(fi, self.profile_method_calculator.h_phi, [-1, -1, -1, -1], [self.config.r0,
                                                                self.config.r0 + self.config.h,
                                                                self.config.r0 + self.config.h,
                                                                self.config.r0 + self.config.h,
                                                                self.config.r0,
                                                                self.config.r0 - self.config.z])

    def fun_h_2(self, fi: float | np.ndarray):
        if type(fi) is np.ndarray:
            return (self.fun_universal(fi, self.profile_method_calculator.h_phi, [-1, -1, -1, -1], [self.config.r0,
                                                                     self.config.r0 + self.config.h,
                                                                     self.config.r0 + self.config.h,
                                                                     self.config.r0 + self.config.h,
                                                                     self.config.r0,
                                                                     self.config.r0 - self.config.z]) - self.config.r0) * np.int64(fi >= self.config.phi_1) * np.int64(fi <= self.config.phi_4)
        return (self.fun_universal(fi, self.profile_method_calculator.h_phi, [-1, -1, -1, -1], [self.config.r0,
                                                                 self.config.r0 + self.config.h,
                                                                 self.config.r0 + self.config.h,
                                                                 self.config.r0 + self.config.h,
                                                                 self.config.r0,
                                                                 self.config.r0 - self.config.z]) - self.config.r0) * int(fi >= self.config.phi_1) * int(fi <= self.config.phi_4)

    def fun_v(self, fi: float | np.ndarray):
        return self.fun_universal(fi, self.profile_method_calculator.v_phi, [-1, -1, 1, 1], [0, 0, 0, 0, 0, 0])

    def fun_v_2(self, fi: float | np.ndarray):
        if type(fi) is np.ndarray:
            return self.fun_universal(fi, self.profile_method_calculator.v_phi, [-1, -1, 1, 1], [0, 0, 0, 0, 0, 0]) * np.int64(
                fi >= self.config.phi_1) * np.int64(fi <= self.config.phi_4)
        return self.fun_universal(fi, self.profile_method_calculator.v_phi, [-1, -1, 1, 1], [0, 0, 0, 0, 0, 0]) * int(fi >= self.config.phi_1) * int(
            fi <= self.config.phi_4)

    def fun_a(self, fi: float | np.ndarray):
        return self.fun_universal(fi, self.profile_method_calculator.a_phi, [-1, -1, -1, -1], [0, 0, 0, 0, 0, 0])

    def fun_a_2(self, fi: float | np.ndarray):
        if type(fi) is np.ndarray:
            return self.fun_universal(fi, self.profile_method_calculator.a_phi, [-1, -1, -1, -1], [0, 0, 0, 0, 0, 0]) * np.int64(
                fi >= self.config.phi_1) * np.int64(fi <= self.config.phi_4)
        return self.fun_universal(fi, self.profile_method_calculator.a_phi, [-1, -1, 1, 1], [0, 0, 0, 0, 0, 0]) * int(fi >= self.config.phi_1) * int(
            fi <= self.config.phi_4)

    def fun_d(self, fi: float | np.ndarray):
        return self.fun_universal(fi, self.profile_method_calculator.d_phi, [-1, -1, 1, 1], [0, 0, 0, 0, 0, 0])

    def fun_d_2(self, fi: float | np.ndarray):
        if type(fi) is np.ndarray:
            return self.fun_universal(fi, self.profile_method_calculator.d_phi, [-1, -1, 1, 1], [0, 0, 0, 0, 0, 0]) * np.int64(
                fi >= self.config.phi_1) * np.int64(fi <= self.config.phi_4)
        return self.fun_universal(fi, self.profile_method_calculator.d_phi, [-1, -1, 1, 1], [0, 0, 0, 0, 0, 0]) * int(fi >= self.config.phi_1) * int(
            fi <= self.config.phi_4)

    def fun_k(self, fi: float | np.ndarray):
        return self.fun_universal(fi, self.profile_method_calculator.k_phi, [-1, -1, -1, -1], [0, 0, 0, 0, 0, 0])

    def fun_k_2(self, fi: float | np.ndarray):
        if type(fi) is np.ndarray:
            return self.fun_universal(fi, self.profile_method_calculator.k_phi, [-1, -1, -1, -1], [0, 0, 0, 0, 0, 0]) * np.int64(
                fi >= self.config.phi_1) * np.int64(fi <= self.config.phi_4)
        return self.fun_universal(fi, self.profile_method_calculator.k_phi, [-1, -1, 1, 1], [0, 0, 0, 0, 0, 0]) * int(fi >= self.config.phi_1) * int(
            fi <= self.config.phi_4)

    def fun_x(self, fi: float | np.ndarray):
        return self.fun_h(fi) * np.cos(fi)

    def fun_y(self, fi: float | np.ndarray):
        return self.fun_h(fi) * np.sin(fi)

    def set_kulachok_data(self, N: int = 1000):
        self.kulachok_data = set_graph_data([self.fun_h, self.fun_v, self.fun_a, self.fun_d, self.fun_k], self.config.omega, N = N)
        self.kulachok_solve_flag = True

    def set_tolkatel_data(self, N: int = 1000):
        self.tolkatel_data = set_graph_data([self.fun_h_2, self.fun_v_2, self.fun_a_2, self.fun_d_2, self.fun_k_2], self.config.omega, N = N)
        self.tolkatel_solve_flag = True

    def set_profile_data(self, N: int = 1000):
        self.profile_data = set_profile_data([self.fun_x, self.fun_y], N = N)
        self.profile_solve_flag = True
        self.solve_type = "thin"

    def solve(self, kulachok_type: Literal['thin', 'flat', 'roller'] = 'thin', N: int = 1000):
        self.set_tolkatel_data(N = N)
        self.set_kulachok_data(N = N)
        if kulachok_type == 'thin':
            self.set_profile_data(N = N)
        elif kulachok_type == 'flat':
            self.set_profile_flat()
        else:
            self.set_profile_roller()

    def profile_flat_check(self, curvature_flag: bool = True):
        if not(self.kulachok_solve_flag and self.tolkatel_solve_flag):
            raise SolvePreliminaryCalculations(f"Не были проведены предварительные вычисления закона движения толкатиля")

        max_v = np.max(self.tolkatel_data.V_t / self.config.omega)
        if self.config.D_t * 1e3 / 2 <= max_v:
            raise PusherDiameterError(
                f"Недостаточный диаметр толкателя: {self.config.D_t * 1e3:.2f} <= {max_v:.2f}"
            )

        curvature_check = self.tolkatel_data.H_rad + self.tolkatel_data.A_rad + 500 * self.config.D
        min_curvature = np.min(curvature_check)
        if min_curvature <= 0 and curvature_flag:
            raise ProfileSmoothnessError(
                f"Негладкий профиль (min = {min_curvature:.4f}). "
                "Необходимо повысить минимальный радиус кривизны."
            )

    def set_profile_flat(self):
        # Проверка возможности построения профиля для заданного закона движения толкателя
        self.profile_flat_check()

        # Расчёт угла откланения и эксцентроситета
        E = self.tolkatel_data.V_t / self.tolkatel_data.omega_rad
        R = np.sqrt(E**2 + (self.tolkatel_data.H_t + self.config.D * 1e3 / 2)**2)
        delta_fi = np.atan(E / (self.tolkatel_data.H_t + self.config.D * 1e3 / 2))
        fi_list_tolkatel = self.tolkatel_data.fi_list_rad
        fi_list_kulachok = (fi_list_tolkatel + delta_fi) % (2 * np.pi)

        # Выправка массивов перед интерполяцией
        index = np.argsort(fi_list_kulachok)
        R = R[index]
        fi_list_kulachok = fi_list_kulachok[index]
        unique_idx = np.unique(fi_list_kulachok, return_index=True)[1]
        fi_list_kulachok = fi_list_kulachok[unique_idx]
        R = R[unique_idx]
        fi_list_kulachok[0] = 0
        R[0] = self.config.D * 1e3 / 2
        fi_list_kulachok[-1] = 2 * np.pi
        R[-1] = self.config.D * 1e3 / 2

        # Интерполяция
        R_func = interp1d(fi_list_kulachok, R, kind="linear")
        self.profile_data = ProfileData(fi_list=fi_list_tolkatel.copy(), X=R_func(fi_list_tolkatel) * np.cos(fi_list_tolkatel), Y=R_func(fi_list_tolkatel) * np.sin(fi_list_tolkatel))
        self.profile_solve_flag = True
        self.solve_type = "flat"

    def profile_roller_check(self):
        h = self.tolkatel_data.H_rad
        v = self.tolkatel_data.V_rad
        a = self.tolkatel_data.A_rad
        Rb = self.config.D * 1e3 / 2.0
        Rr = self.config.R_r * 1e3
        Rc = Rb + Rr + h

        num = (Rc ** 2 + v ** 2) ** 1.5
        den = (Rc ** 2 + 2 * v ** 2 - Rc * a)
        den[den == 0] = 1e-9
        rho_c = num / den
        rho_p = rho_c - Rr
        if np.any(rho_p < 0):
            raise ProfileSmoothnessError(
                f"Негладкий профиль (min = {np.min(rho_p):.4f}). "
                "Необходимо повысить минимальный радиус кривизны."
            )

    def set_profile_roller(self):
        """
        Рассчитывает профиль кулачка на основе конфигурации и кинематических данных.
        ИСПРАВЛЕНО: Приведено к стандартной полярной системе (X = cos, Y = sin)
        для совместимости с анимацией.
        """

        self.profile_roller_check()

        # 1. Извлечение данных
        phi = self.tolkatel_data.fi_list_rad
        h = self.tolkatel_data.H_rad
        v = self.tolkatel_data.V_rad
        # a = self.tolkatel_data.A_rad # Не используется в координатах, только в проверке

        # Геометрические параметры
        Rb = self.config.D * 1e3 / 2.0
        Rr = self.config.R_r * 1e3

        # 2. Расчет теоретического профиля (Траектория центра ролика)
        Rc = Rb + Rr + h

        xc = Rc * np.cos(phi)
        yc = Rc * np.sin(phi)

        # 3. Расчет действительного профиля
        # Производные координат центра ролика по углу phi
        # Производная от (Rc * cos) -> Rc' * cos + Rc * (-sin) -> v * cos - Rc * sin
        dxc = v * np.cos(phi) - Rc * np.sin(phi)

        # Производная от (Rc * sin) -> Rc' * sin + Rc * cos -> v * sin + Rc * cos
        dyc = v * np.sin(phi) + Rc * np.cos(phi)

        # Модуль вектора касательной
        norm_factor = np.sqrt(dxc ** 2 + dyc ** 2)
        norm_factor[norm_factor == 0] = 1e-9

        # Координаты конструктивного профиля
        # Смещение на Rr внутрь кривой.
        # Формулы смещения эквидистанты:
        # xp = xc + Rr * dy / sqrt(...) * sign
        # yp = yc - Rr * dx / sqrt(...) * sign
        # Для внутреннего смещения (уменьшения радиуса) знаки обычно такие:
        xp = xc - Rr * (dyc / norm_factor)
        yp = yc + Rr * (dxc / norm_factor)

        # Важно: В зависимости от направления обхода (CW/CCW) знаки могут инвертироваться.
        # Если профиль "вывернется", попробуйте поменять знаки перед Rr на противоположные:
        # xp = xc - Rr * ...
        # yp = yc + Rr * ...
        # Но для стандартной математики (CCW) текущий вариант должен быть верным,
        # если мы хотим "уменьшить" профиль относительно центра ролика.

        self.profile_data = ProfileData(fi_list=self.tolkatel_data.fi_list_rad.copy(),
                                       X=xp,
                                       Y=yp)
        self.profile_solve_flag = True
        self.solve_type = "roller"