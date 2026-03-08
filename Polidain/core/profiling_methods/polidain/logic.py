from core.profiling_methods.polidain.config import PolidainConfig
from core.profiling_methods.base import BaseCalculator

def k_fun(p, m, d):
    """
    Вычисляет степени членов полинома.

    Args:
        p: Степень второго члена полинома.
        m: Степень первого члена полинома.
        p - d: Разность между вторым и последующими членами полинома.

    Returns:
        list: [m, p, q, r, s] - массив степеней членов полинома.
    """
    q = p + p - d
    r = q + p - d
    s = r + p - d
    return [m, p, q, r, s]

def c_fun(p, m, d):
    """
    Вычисляет коэффициенты полинома.

    Args:
        p: Степень второго члена полинома.
        m: Младшая степень полинома.
        d: Разность между степенями членов полинома.

    Returns:
        list: [C2, Cp, Cq, Cr, Cs] - массив коэффициентов при членах полинома.
    """
    temp = k_fun(p, m, d)
    q = temp[2]
    r = temp[3]
    s = temp[4]
    C2 = p * q * r * s / ((p - m) * (q - m) * (r - m) * (s - m))
    Cp = -m * q * r * s / ((p - m) * (q - p) * (r - p) * (s - p))
    Cq = m * p * s * r / ((q - m) * (r - q) * (q - p) * (s - q))
    Cr = -m * p * q * s / ((r - m) * (s - r) * (r - p) * (r - q))
    Cs = m * p * q * r / ((s - m) * (s - p) * (s - q) * (s - r))
    return [C2, Cp, Cq, Cr, Cs]

def h_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    """
    Вычисляет значение перемещения кулачка (высоту подъема).

    Args:
        fi: Текущий угол поворота кулачка.
        C_list: Массив коэффициентов при членах полинома.
        k_list: Массив степеней членов полинома.
        fi_1: Конечный угол характерного участка.
        fi_0: Начальный угол характерного участка.
        h_kn_max: Максимальная высота подъема кулачка на участке.

    Returns:
        float: Значение перемещения.
    """
    temp = 0
    for i in range(0, len(C_list)):
        if k_list[i] - 0 < 0:
            continue
        temp += C_list[i] * (((fi - fi_0) / (fi_1 - fi_0)) ** k_list[i])
    return temp * h_kn_max


def v_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    """
    Вычисляет аналог скорости (первая производная перемещения).

    Args:
        fi: Текущий угол поворота кулачка.
        C_list: Массив коэффициентов при членах полинома.
        k_list: Массив степеней членов полинома.
        fi_1: Конечный угол характерного участка.
        fi_0: Начальный угол характерного участка.
        h_kn_max: Максимальная высота подъема кулачка на участке.

    Returns:
        float: Значение аналога скорости.
    """
    temp = 0
    for i in range(0, len(C_list)):
        if k_list[i] - 1 < 0:
            continue
        temp += k_list[i] * C_list[i] * (((fi - fi_0) / (fi_1 - fi_0)) ** (k_list[i] - 1)) / (fi_1 - fi_0)
    return temp * h_kn_max


def a_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    """
    Вычисляет аналог ускорения (вторая производная перемещения).

    Args:
        fi: Текущий угол поворота кулачка.
        C_list: Массив коэффициентов при членах полинома.
        k_list: Массив степеней членов полинома.
        fi_1: Конечный угол характерного участка.
        fi_0: Начальный угол характерного участка.
        h_kn_max: Максимальная высота подъема кулачка на участке.

    Returns:
        float: Значение аналога ускорения.
    """
    temp = 0
    for i in range(0, len(C_list)):
        if k_list[i] - 2 < 0:
            continue
        temp += k_list[i] * (k_list[i] - 1) * C_list[i] * (((fi - fi_0) / (fi_1 - fi_0)) ** (k_list[i] - 2)) / ((fi_1 - fi_0)**2)
    return temp * h_kn_max

def d_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    """
    Вычисляет аналог рывка (третья производная перемещения).

    Args:
        fi: Текущий угол поворота кулачка.
        C_list: Массив коэффициентов при членах полинома.
        k_list: Массив степеней членов полинома.
        fi_1: Конечный угол характерного участка.
        fi_0: Начальный угол характерного участка.
        h_kn_max: Максимальная высота подъема кулачка на участке.

    Returns:
        float: Значение аналога рывка.
    """
    temp = 0
    for i in range(0, len(C_list)):
        if k_list[i] - 3 < 0:
            continue
        temp += k_list[i] * (k_list[i] - 1) * (k_list[i] - 2) * C_list[i] * (
                ((fi - fi_0) / (fi_1 - fi_0)) ** (k_list[i] - 3)) / ((fi_1 - fi_0)**3)
    return temp * h_kn_max

def k_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    """
    Вычисляет четвертую производную перемещения.

    Args:
        fi: Текущий угол поворота кулачка.
        C_list: Массив коэффициентов при членах полинома.
        k_list: Массив степеней членов полинома.
        fi_1: Конечный угол характерного участка.
        fi_0: Начальный угол характерного участка.
        h_kn_max: Максимальная высота подъема кулачка на участке.

    Returns:
        float: Значение четвертой производной.
    """
    temp = 0
    for i in range(0, len(C_list)):
        if k_list[i] - 4 < 0:
            continue
        temp += k_list[i] * (k_list[i] - 1) * (k_list[i] - 2) * (k_list[i] - 3) * C_list[i] * (
                ((fi - fi_0) / (fi_1 - fi_0)) ** (k_list[i] - 4)) / ((fi_1 - fi_0)**4)
    return temp * h_kn_max

class PolidainCalculator(BaseCalculator):
    """
    Калькулятор профиля методом Polidain (полиномиальная аппроксимация).
    Использует полиномы высоких степеней для обеспечения гладкости функций.
    """

    def __init__(self, config: PolidainConfig):
        """
        Инициализация калькулятора.
        Предварительно рассчитывает коэффициенты полиномов для всех участков.

        Args:
            config: Конфигурация Polidain.
        """
        BaseCalculator.__init__(self, config)

        self.c_list_1: list[float] = c_fun(self.config.k_1, self.config.m, self.config.d)
        self.c_list_2: list[float] = c_fun(self.config.k_2, self.config.m, self.config.d)
        self.c_list_3: list[float] = c_fun(self.config.k_3, self.config.m, self.config.d)
        self.c_list_4: list[float] = c_fun(self.config.k_4, self.config.m, self.config.d)

        self.k_list_1: list[int] = k_fun(self.config.k_1, self.config.m, self.config.d)
        self.k_list_2: list[int] = k_fun(self.config.k_2, self.config.m, self.config.d)
        self.k_list_3: list[int] = k_fun(self.config.k_3, self.config.m, self.config.d)
        self.k_list_4: list[int] = k_fun(self.config.k_4, self.config.m, self.config.d)

        self._segments = {
            1: (self.c_list_1, self.k_list_1),
            2: (self.c_list_2, self.k_list_2),
            3: (self.c_list_3, self.k_list_3),
            4: (self.c_list_4, self.k_list_4),
        }

    def segment_selection(self, segment_number: int):
        """
        Выбирает коэффициенты и степени полинома для заданного участка.

        Args:
            segment_number: Номер участка.

        Returns:
            Tuple[list[float], list[int]]: Кортеж (список коэффициентов, список степеней).

        Raises:
            ValueError: Если участок не настроен.
        """
        if segment_number not in self._segments:
            raise ValueError(f'Участок {segment_number} не настроен')
        return self._segments[segment_number]

    def h_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """См. BaseCalculator.h_phi"""
        c_list, k_list = self.segment_selection(segment_number)
        return h_phi(fi, c_list, k_list, fi_1, fi_0, h_kn_max)

    def v_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """См. BaseCalculator.v_phi"""
        c_list, k_list = self.segment_selection(segment_number)
        return v_phi(fi, c_list, k_list, fi_1, fi_0, h_kn_max)

    def a_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """См. BaseCalculator.a_phi"""
        c_list, k_list = self.segment_selection(segment_number)
        return a_phi(fi, c_list, k_list, fi_1, fi_0, h_kn_max)

    def d_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """См. BaseCalculator.d_phi"""
        c_list, k_list = self.segment_selection(segment_number)
        return d_phi(fi, c_list, k_list, fi_1, fi_0, h_kn_max)

    def k_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """См. BaseCalculator.k_phi"""
        c_list, k_list = self.segment_selection(segment_number)
        return k_phi(fi, c_list, k_list, fi_1, fi_0, h_kn_max)
