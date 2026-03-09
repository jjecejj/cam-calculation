from dataclasses import dataclass, field
from typing import Callable, Literal

import numpy as np


@dataclass
class OptimizeConfig:
    profiling_type: Literal['polidain', 'polinmail']
    R_func: Callable
    D: float  # Базовый диаметр кулака (мм)
    h: float  # Максимальное перемещение толкателя
    N_k: float = 1000 # Количество оборотов кулачка в минут
    H_func: Callable | None = None
    tolkatel_type: Literal['thin', 'flat', 'roller'] = 'thin'
    D_t: float = 0.0  # Диаметр толкателя (мм)
    R_r: float = 0.0  # Радиус ролика (мм)

@dataclass
class BoundsConfig:
    z: tuple = field(default=(0.0, 0.5))
    f_pod: tuple = field(default=(0.001, np.pi))
    f_v: tuple = field(default=(0.0, np.pi / 6))
    f_op: tuple = field(default=(0.001, np.pi))
    f_z: tuple = field(default=(0.0, np.pi / 6))

    """Опции Polidain"""
    m: tuple = field(default=(2, 100))
    d: tuple = field(default=(1, 100))
    k_1: tuple = field(default=(10, 1000))
    k_2: tuple = field(default=(10, 100))
    k_3: tuple = field(default=(10, 100))
    k_4: tuple = field(default=(10, 1000))

    """Опции Polinmail"""
    m_1: tuple = field(default=(2, 100))
    d_1: tuple = field(default=(1, 100))
    m_2: tuple = field(default=(2, 100))
    d_2: tuple = field(default=(1, 100))
    m_3: tuple = field(default=(2, 100))
    d_3: tuple = field(default=(1, 100))
    m_4: tuple = field(default=(2, 100))
    d_4: tuple = field(default=(1, 100))

    @property
    def bounds_polidain(self) -> list:
        return [self.z, self.f_pod, self.f_v, self.f_op, self.f_z, self.k_1, self.k_2, self.k_3, self.k_4, self.m, self.d]

    @property
    def integrality_polidain(self) -> list:
        return [False, False, False, False, False, True, True, True, True, True, True]

    @property
    def bounds_polinmail(self) -> list:
        return [self.z, self.f_pod, self.f_v, self.f_op, self.f_z, self.m_1, self.d_1, self.m_2, self.d_2, self.m_3, self.d_3, self.m_4, self.d_4]

    @property
    def integrality_polinmail(self) -> list:
        return [False, False, False, False, False, True, True, True, True, True, True, True, True]


@dataclass
class DifferentialEvolutionConfig:
    strategy: str = 'best1bin'
    maxiter: int = 1000
    popsize: int = 400
    workers: int = -1
    tol: float = 0.001
