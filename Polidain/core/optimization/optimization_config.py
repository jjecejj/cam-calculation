from dataclasses import dataclass, field

import numpy as np


@dataclass
class OptimizeConfig:
    D: float  # Базовый диаметр кулака (мм)
    h: float  # Максимальное перемещение толкателя
    N_k: float = 1000 # Количество оборотов кулачка в минут

@dataclass
class BoundsConfig:
    z: tuple = field(default=(0.01, 0.5))
    f_pod: tuple = field(default=(0.001, np.pi))
    f_v: tuple = field(default=(0.00001, np.pi / 6))
    f_op: tuple = field(default=(0.001, np.pi))
    f_z: tuple = field(default=(0.001, np.pi / 6))
    k_1: tuple = field(default=(10, 1000))
    k_2: tuple = field(default=(10, 100))
    k_3: tuple = field(default=(10, 100))
    k_4: tuple = field(default=(10, 1000))

    @property
    def bounds(self) -> list:
        return [self.z, self.f_pod, self.f_v, self.f_op, self.f_z, self.k_1, self.k_2, self.k_3, self.k_4]


@dataclass
class DifferentialEvolutionConfig:
    strategy: str = 'best1bin'
    maxiter: int = 1000
    popsize: int = 300
    workers: int = -1
    tol: float = 0.001
    integrality: list[bool] = field(default_factory=lambda: [False, False, False, False, False, True, True, True, True])

@dataclass
class GibridOptimizationConfig:
    m: list = field(default_factory=list)
    d: list = field(default_factory=list)
