import numpy as np

from core.profiling_methods.polinmail.config import PolinmailConfig
from core.profiling_methods.base.logic import BaseCalculator
from typing import List

def get_matrix_coefficients(m_list: List | np.ndarray) -> np.ndarray:
    m_list = np.asarray(m_list)
    A = [np.array([1 for i in range(0, len(m_list))])]
    for i in range(0, len(m_list) - 1):
        A.append(A[-1] * (m_list - i))
    return np.array(A)

def calculate_poly_coefficients(m_list, boundary_conditions):
    A = get_matrix_coefficients(m_list)
    B = boundary_conditions
    coeffs = np.linalg.solve(A, B)
    return coeffs

def h_phi(fi: float, m_list: list | np.ndarray, coeffs: list | np.ndarray, fi_1: float, fi_0: float, h_kn_max: float):
    temp = 1
    for i in range(0, len(coeffs)):
        if m_list[i] - 0 < 0:
            continue
        temp += coeffs[i] * (((fi - fi_0) / (fi_1 - fi_0)) ** m_list[i])
    return temp * h_kn_max

def v_phi(fi: float, m_list: list | np.ndarray, coeffs: list | np.ndarray, fi_1: float, fi_0: float, h_kn_max: float):
    '''Функция скорости от угла поворота кулачка
    :param fi: Угол поворота кулачка
    :param m_list: Массив степеней членов полинома
    :param coeffs: Массив коэффициентов при членах полинома
    :param fi_1: Конец характерного участка
    :param fi_0: Начало характерного участка
    :param h_kn_max: максимальная высота подёма кулачка на характерном участке
    '''
    temp = 0
    for i in range(0, len(coeffs)):
        if m_list[i] - 1 < 0:
            continue
        temp += m_list[i] * coeffs[i] * (((fi - fi_0) / (fi_1 - fi_0)) ** (m_list[i] - 1)) / (fi_1 - fi_0)
    return temp * h_kn_max


def a_phi(fi: float, m_list: list | np.ndarray, coeffs: list | np.ndarray, fi_1: float, fi_0: float, h_kn_max: float):
    '''Функция ускорения от угла поворота кулачка
    :param fi: Угол поворота кулачка
    :param m_list: Массив степеней членов полинома
    :param coeffs: Массив коэффициентов при членах полинома
    :param fi_1: Конец характерного участка
    :param fi_0: Начало характерного участка
    :param h_kn_max: максимальная высота подёма кулачка на характерном участке
    '''
    temp = 0
    for i in range(0, len(coeffs)):
        if m_list[i] - 2 < 0:
            continue
        temp += m_list[i] * (m_list[i] - 1) * coeffs[i] * (((fi - fi_0) / (fi_1 - fi_0)) ** (m_list[i] - 2)) / ((fi_1 - fi_0)**2)
    return temp * h_kn_max


def d_phi(fi: float, m_list: list | np.ndarray, coeffs: list | np.ndarray, fi_1: float, fi_0: float, h_kn_max: float):
    '''Функция рывка от угла поворота кулачка
    :param fi: Угол поворота кулачка
    :param m_list: Массив степеней членов полинома
    :param coeffs: Массив коэффициентов при членах полинома
    :param fi_1: Конец характерного участка
    :param fi_0: Начало характерного участка
    :param h_kn_max: максимальная высота подёма кулачка на характерном участке
    '''
    temp = 0
    for i in range(0, len(coeffs)):
        if m_list[i] - 3 < 0:
            continue
        temp += m_list[i] * (m_list[i] - 1) * (m_list[i] - 2) * coeffs[i] * (
                ((fi - fi_0) / (fi_1 - fi_0)) ** (m_list[i] - 3)) / ((fi_1 - fi_0)**3)
    return temp * h_kn_max


def k_phi(fi: float, m_list: list | np.ndarray, coeffs: list | np.ndarray, fi_1: float, fi_0: float, h_kn_max: float):
    '''Функция кракена от угла поворота кулачка
    :param fi: Угол поворота кулачка
    :param m_list: Массив степеней членов полинома
    :param coeffs: Массив коэффициентов при членах полинома
    :param fi_1: Конец характерного участка
    :param fi_0: Начало характерного участка
    :param h_kn_max: максимальная высота подёма кулачка на характерном участке
    '''
    temp = 0
    for i in range(0, len(coeffs)):
        if m_list[i] - 4 < 0:
            continue
        temp += m_list[i] * (m_list[i] - 1) * (m_list[i] - 2) * (m_list[i] - 3) * coeffs[i] * (
                ((fi - fi_0) / (fi_1 - fi_0)) ** (m_list[i] - 4)) / ((fi_1 - fi_0)**4)
    return temp * h_kn_max

class PolinmailCalculator(BaseCalculator):
    def __init__(self, config: PolinmailConfig):
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
        """Выбор параметров сегмента по его номеру."""
        if segment_number not in self._segments:
            raise ValueError(f"Сегмент {segment_number} не настроен")
        return self._segments[segment_number]

    def h_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """Основной расчет функции h_phi."""
        m_list, coeffs = self.segment_selection(segment_number)
        return h_phi(fi, m_list, coeffs, fi_1, fi_0, h_kn_max)

    def v_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """Основной расчет функции v_phi."""
        m_list, coeffs = self.segment_selection(segment_number)
        return v_phi(fi, m_list, coeffs, fi_1, fi_0, h_kn_max)

    def a_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """Основной расчет функции a_phi."""
        m_list, coeffs = self.segment_selection(segment_number)
        return a_phi(fi, m_list, coeffs, fi_1, fi_0, h_kn_max)

    def d_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """Основной расчет функции d_phi."""
        m_list, coeffs = self.segment_selection(segment_number)
        return d_phi(fi, m_list, coeffs, fi_1, fi_0, h_kn_max)

    def k_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """Основной расчет функции k_phi."""
        m_list, coeffs = self.segment_selection(segment_number)
        return k_phi(fi, m_list, coeffs, fi_1, fi_0, h_kn_max)