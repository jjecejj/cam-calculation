from core.profiling_methods.polinmail.polinmail_config import PolinmailConfig
from core.profiling_methods.base.base import BaseCalculator
from typing import List

def get_powers(n: int) -> List[int]:
    """
    Возвращает список степеней для полинома из n членов.
    n=4: 3-я степень (члены x^2, x^3) - аналог S=3x^2 - 2x^3 (Cubic)
    n=5: 4-я степень (члены x^3, x^4) - аналог S=4x^3 - 3x^4 (Quartic)
    n=6: 5-я степень (члены x^3, x^4, x^5) - аналог S=10x^3 - 15x^4 + 6x^5 (Quintic 3-4-5)
    """
    if n == 4:
        return [2, 3]
    elif n == 5:
        return [3, 4]
    elif n == 6:
        return [3, 4, 5]
    return []

def calculate_coefficients(powers: List[int]) -> List[float]:
    """
    Вычисляет коэффициенты C_i для полинома, удовлетворяющего условиям S(1)=1, и производным=0.
    Использует обобщенную формулу для системы Вандермонда производных.
    C_i = (-1)^(N-1) * product(p_j for j!=i) / product(p_i - p_j for j!=i)
    """
    coeffs = []
    num_terms = len(powers)
    sign_factor = (-1) ** (num_terms - 1)

    for i, p_i in enumerate(powers):
        numerator = 1.0
        denominator = 1.0
        for j, p_j in enumerate(powers):
            if i == j:
                continue
            numerator *= p_j
            denominator *= (p_i - p_j)
        coeffs.append(sign_factor * numerator / denominator)
    return coeffs

def h_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    """Функция радиуса кулачка"""
    temp = 0
    normalized_fi = (fi - fi_0) / (fi_1 - fi_0)
    for i in range(len(C_list)):
        if k_list[i] < 0: continue
        temp += C_list[i] * (normalized_fi ** k_list[i])
    return temp * h_kn_max

def v_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    """Функция скорости"""
    temp = 0
    normalized_fi = (fi - fi_0) / (fi_1 - fi_0)
    for i in range(len(C_list)):
        if k_list[i] < 1: continue
        temp += k_list[i] * C_list[i] * (normalized_fi ** (k_list[i] - 1))
    return temp * h_kn_max / (fi_1 - fi_0)

def a_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    """Функция ускорения"""
    temp = 0
    normalized_fi = (fi - fi_0) / (fi_1 - fi_0)
    for i in range(len(C_list)):
        if k_list[i] < 2: continue
        temp += k_list[i] * (k_list[i] - 1) * C_list[i] * (normalized_fi ** (k_list[i] - 2))
    return temp * h_kn_max / ((fi_1 - fi_0)**2)

def d_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    """Функция рывка (jerk)"""
    temp = 0
    normalized_fi = (fi - fi_0) / (fi_1 - fi_0)
    for i in range(len(C_list)):
        if k_list[i] < 3: continue
        temp += k_list[i] * (k_list[i] - 1) * (k_list[i] - 2) * C_list[i] * (normalized_fi ** (k_list[i] - 3))
    return temp * h_kn_max / ((fi_1 - fi_0)**3)

def k_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    """Функция пинга (ping)"""
    temp = 0
    normalized_fi = (fi - fi_0) / (fi_1 - fi_0)
    for i in range(len(C_list)):
        if k_list[i] < 4: continue
        temp += k_list[i] * (k_list[i] - 1) * (k_list[i] - 2) * (k_list[i] - 3) * C_list[i] * (normalized_fi ** (k_list[i] - 4))
    return temp * h_kn_max / ((fi_1 - fi_0)**4)


class PolinmailCalculator(BaseCalculator):
    def __init__(self, config: PolinmailConfig):
        BaseCalculator.__init__(self, config)

        # Calculate powers and coefficients once based on configuration
        self.k_list: List[int] = get_powers(self.config.n)
        self.c_list: List[float] = calculate_coefficients(self.k_list)

    def h_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        return h_phi(fi, self.c_list, self.k_list, fi_1, fi_0, h_kn_max)

    def v_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        return v_phi(fi, self.c_list, self.k_list, fi_1, fi_0, h_kn_max)

    def a_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        return a_phi(fi, self.c_list, self.k_list, fi_1, fi_0, h_kn_max)

    def d_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        return d_phi(fi, self.c_list, self.k_list, fi_1, fi_0, h_kn_max)

    def k_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        return k_phi(fi, self.c_list, self.k_list, fi_1, fi_0, h_kn_max)
