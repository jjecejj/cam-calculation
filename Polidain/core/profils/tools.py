import numpy as np
from scipy.interpolate import CubicSpline

from core.schemas import set_graph_data, set_profile_data


def fi_list_dif(fi_array: np.ndarray, f_dif: float) -> np.ndarray:
    """
    Векторизированный сдвиг фазы
    """
    shifted = fi_array - f_dif
    return np.where(shifted <= 0, shifted + 2 * np.pi, shifted)

class ProfileDataExtractor:
    def __init__(self, data: list[tuple[float, float]], f_dif: float):
        if len(data) == 0:
            raise ValueError("Data is empty")

        # Преобразование в полярные координаты
        X = []
        Y = []
        R = []
        FI = []
        for x, y in data:
            X.append(x)
            Y.append(y)
            r = np.hypot(x, y)
            if y >= 0:
                phi = np.atan2(y, x)
            else:
                phi = 2 * np.pi + np.atan2(y, x)  # угол в радианах
            R.append(r)
            FI.append(phi)
        np.append(R, R[0])  # Дублируем радиус первой точки в конец
        np.append(FI, FI[0] + 2 * np.pi)  # Дублируем угол первой точки + полный оборот

        FI = fi_list_dif(np.array(FI), f_dif)
        R = np.array(R) / 1000

        # Cортируем массивы
        sorted_indices = np.argsort(FI)
        FI = FI[sorted_indices]
        R = R[sorted_indices]
        R[-1] = R[0]

        # Интерполяция и производные с помощью CubicSpline
        cs = CubicSpline(FI, R, bc_type='periodic')

        # Создаем функции для производных
        self.H = cs  # Сама функция R(phi)
        self.V = cs.derivative(1)  # Первая производная R'(phi)
        self.A = cs.derivative(2)  # Вторая производная R''(phi)
        self.D = cs.derivative(3)  # Третья производная R'''(phi)
        self.K = cs.derivative(4)  # Четвёртая производная R''''(phi)

    def x_func(self, fi: np.ndarray | float | int):
        return self.H(fi) * np.cos(fi)

    def y_func(self, fi: np.ndarray | float | int):
        return self.H(fi) * np.sin(fi)

    def get_H(self):
        return self.H

    def get_V(self):
        return self.V

    def get_A(self):
        return self.A

    def get_D(self):
        return self.D

    def get_K(self):
        return self.K

    def get_GraphData(self, omega: float = 0.0, N: int = 1000):
        return set_graph_data([self.H, self.V, self.A, self.D, self.K], omega = omega, N = N)

    def get_ProfileData(self, N: int = 1000):
        return set_profile_data([self.x_func, self.y_func], N = N)