import numpy as np
from core.profiling_methods.polinmail.config import PolinmailConfig
from core.profiling_methods.base import BaseCalculator
from typing import List

def get_matrix_coefficients(m_list: List | np.ndarray) -> np.ndarray:
    """
    Формирует матрицу коэффициентов для системы линейных уравнений.

    Args:
        m_list: Список степеней полинома.

    Returns:
        np.ndarray: Матрица коэффициентов.
    """
    m_list = np.asarray(m_list)
    A = [np.array([1 for i in range(0, len(m_list))])]
    for i in range(0, len(m_list) - 1):
        A.append(A[-1] * (m_list - i))
    return np.array(A)

def calculate_poly_coefficients(m_list, boundary_conditions):
    """
    Вычисляет коэффициенты полинома на основе граничных условий.

    Args:
        m_list: Список степеней полинома.
        boundary_conditions: Граничные условия.

    Returns:
        np.ndarray: Коэффициенты полинома (с обратным знаком).
    """
    A = get_matrix_coefficients(m_list)
    B = boundary_conditions
    coeffs = np.linalg.solve(A, B)
    return -coeffs

def h_phi(fi: float, m_list: list | np.ndarray, coeffs: list | np.ndarray, fi_1: float, fi_0: float, h_kn_max: float):
    """
    Вычисляет значение перемещения кулачка.

    Args:
        fi: Текущий угол поворота кулачка.
        m_list: Массив степеней членов полинома.
        coeffs: Массив коэффициентов при членах полинома.
        fi_1: Конец характерного участка.
        fi_0: Начало характерного участка.
        h_kn_max: Максимальная высота подъема кулачка на характерном участке.

    Returns:
        float: Значение перемещения.
    """
    temp = 0
    for i in range(0, len(coeffs)):
        if m_list[i] - 0 < 0:
            continue
        temp += coeffs[i] * (((fi - fi_0) / (fi_1 - fi_0)) ** m_list[i])
    return temp * h_kn_max

def v_phi(fi: float, m_list: list | np.ndarray, coeffs: list | np.ndarray, fi_1: float, fi_0: float, h_kn_max: float):
    """
    Вычисляет аналог скорости (первая производная перемещения).

    Args:
        fi: Текущий угол поворота кулачка.
        m_list: Массив степеней членов полинома.
        coeffs: Массив коэффициентов при членах полинома.
        fi_1: Конец характерного участка.
        fi_0: Начало характерного участка.
        h_kn_max: Максимальная высота подъема кулачка на характерном участке.

    Returns:
        float: Значение аналога скорости.
    """
    temp = 0
    for i in range(0, len(coeffs)):
        if m_list[i] - 1 < 0:
            continue
        temp += m_list[i] * coeffs[i] * (((fi - fi_0) / (fi_1 - fi_0)) ** (m_list[i] - 1)) / (fi_1 - fi_0)
    return temp * h_kn_max


def a_phi(fi: float, m_list: list | np.ndarray, coeffs: list | np.ndarray, fi_1: float, fi_0: float, h_kn_max: float):
    """
    Вычисляет аналог ускорения (вторая производная перемещения).

    Args:
        fi: Текущий угол поворота кулачка.
        m_list: Массив степеней членов полинома.
        coeffs: Массив коэффициентов при членах полинома.
        fi_1: Конец характерного участка.
        fi_0: Начало характерного участка.
        h_kn_max: Максимальная высота подъема кулачка на характерном участке.

    Returns:
        float: Значение аналога ускорения.
    """
    temp = 0
    for i in range(0, len(coeffs)):
        if m_list[i] - 2 < 0:
            continue
        temp += m_list[i] * (m_list[i] - 1) * coeffs[i] * (((fi - fi_0) / (fi_1 - fi_0)) ** (m_list[i] - 2)) / ((fi_1 - fi_0)**2)
    return temp * h_kn_max


def d_phi(fi: float, m_list: list | np.ndarray, coeffs: list | np.ndarray, fi_1: float, fi_0: float, h_kn_max: float):
    """
    Вычисляет аналог рывка (третья производная перемещения).

    Args:
        fi: Текущий угол поворота кулачка.
        m_list: Массив степеней членов полинома.
        coeffs: Массив коэффициентов при членах полинома.
        fi_1: Конец характерного участка.
        fi_0: Начало характерного участка.
        h_kn_max: Максимальная высота подъема кулачка на характерном участке.

    Returns:
        float: Значение аналога рывка.
    """
    temp = 0
    for i in range(0, len(coeffs)):
        if m_list[i] - 3 < 0:
            continue
        temp += m_list[i] * (m_list[i] - 1) * (m_list[i] - 2) * coeffs[i] * (
                ((fi - fi_0) / (fi_1 - fi_0)) ** (m_list[i] - 3)) / ((fi_1 - fi_0)**3)
    return temp * h_kn_max


def k_phi(fi: float, m_list: list | np.ndarray, coeffs: list | np.ndarray, fi_1: float, fi_0: float, h_kn_max: float):
    """
    Вычисляет четвертую производную перемещения.

    Args:
        fi: Текущий угол поворота кулачка.
        m_list: Массив степеней членов полинома.
        coeffs: Массив коэффициентов при членах полинома.
        fi_1: Конец характерного участка.
        fi_0: Начало характерного участка.
        h_kn_max: Максимальная высота подъема кулачка на характерном участке.

    Returns:
        float: Значение четвертой производной.
    """
    temp = 0
    for i in range(0, len(coeffs)):
        if m_list[i] - 4 < 0:
            continue
        temp += m_list[i] * (m_list[i] - 1) * (m_list[i] - 2) * (m_list[i] - 3) * coeffs[i] * (
                ((fi - fi_0) / (fi_1 - fi_0)) ** (m_list[i] - 4)) / ((fi_1 - fi_0)**4)
    return temp * h_kn_max

class PolinmailCalculator(BaseCalculator):
    """
    Калькулятор профиля методом Polinmail (метод на основе решения СЛАУ для полиномов).
    Позволяет гибко настраивать степени полиномов и граничные условия.
    """

    def __init__(self, config: PolinmailConfig):
        """
        Инициализация калькулятора.
        Рассчитывает коэффициенты полиномов для всех активных конфигураций (участков).

        Args:
            config: Конфигурация Polinmail.
        """
        BaseCalculator.__init__(self, config)

        # Инициализация конфигураций
        self.config_1 = self.config.config_1
        self.config_2 = self.config.config_2
        self.config_3 = self.config.config_3
        self.config_4 = self.config.config_4

        self.m_list_1 = self.config_1.m_list
        self.m_list_2 = self.config_2.m_list
        self.m_list_3 = self.config_3.m_list
        self.m_list_4 = self.config_4.m_list

        self.coeffs_1 = calculate_poly_coefficients(self.m_list_1, self.config_1.boundary_conditions)
        self.coeffs_2 = calculate_poly_coefficients(self.m_list_2, self.config_2.boundary_conditions)
        self.coeffs_3 = calculate_poly_coefficients(self.m_list_3, self.config_3.boundary_conditions)
        self.coeffs_4 = calculate_poly_coefficients(self.m_list_4, self.config_4.boundary_conditions)

        # Словарь сегментов с данными (m_list и coeffs)
        self._segments = {
            1: (self.m_list_1, self.coeffs_1),
            2: (self.m_list_2, self.coeffs_2),
            3: (self.m_list_3, self.coeffs_3),
            4: (self.m_list_4, self.coeffs_4)
        }

    def segment_selection(self, segment_number: int):
        """
        Выбирает параметры сегмента по его номеру.

        Args:
            segment_number: Номер сегмента.

        Returns:
            Tuple[list, np.ndarray]: Кортеж (список степеней, массив коэффициентов).

        Raises:
            ValueError: Если сегмент не настроен.
        """
        if segment_number not in self._segments:
            raise ValueError(f"Сегмент {segment_number} не настроен")
        return self._segments[segment_number]

    def h_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """См. BaseCalculator.h_phi"""
        m_list, coeffs = self.segment_selection(segment_number)
        return h_phi(fi, m_list, coeffs, fi_1, fi_0, h_kn_max)

    def v_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """См. BaseCalculator.v_phi"""
        m_list, coeffs = self.segment_selection(segment_number)
        return v_phi(fi, m_list, coeffs, fi_1, fi_0, h_kn_max)

    def a_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """См. BaseCalculator.a_phi"""
        m_list, coeffs = self.segment_selection(segment_number)
        return a_phi(fi, m_list, coeffs, fi_1, fi_0, h_kn_max)

    def d_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """См. BaseCalculator.d_phi"""
        m_list, coeffs = self.segment_selection(segment_number)
        return d_phi(fi, m_list, coeffs, fi_1, fi_0, h_kn_max)

    def k_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """См. BaseCalculator.k_phi"""
        m_list, coeffs = self.segment_selection(segment_number)
        return k_phi(fi, m_list, coeffs, fi_1, fi_0, h_kn_max)
