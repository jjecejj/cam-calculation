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

        # Calculate powers and coefficients once based on configuration
        self.k_list: List[int] = get_powers(self.config.n)
        self.c_list: List[float] = calculate_coefficients(self.k_list)

    def h_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        return h_phi(fi, self.m_list, self.coeffs, fi_1, fi_0, h_kn_max)

    def v_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        return v_phi(fi, self.c_list, self.k_list, fi_1, fi_0, h_kn_max)

    def a_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        return a_phi(fi, self.c_list, self.k_list, fi_1, fi_0, h_kn_max)

    def d_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        return d_phi(fi, self.c_list, self.k_list, fi_1, fi_0, h_kn_max)

    def k_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        return k_phi(fi, self.c_list, self.k_list, fi_1, fi_0, h_kn_max)