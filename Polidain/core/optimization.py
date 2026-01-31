import numpy as np
from scipy.optimize import differential_evolution
from multiprocessing import freeze_support
from dataclasses import dataclass, field
from core.cam_geometry import Kulachok_polidain
from core.schemas import PolidainConfig
from typing import Callable


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

def fun_optimize_gibrid(x, fi_list = None, m = None, d = None, R_func = None, optimize_config = None):
    z = x[0]  # Тепловой зазор (мм)
    f_pod = x[1]  # Фаза подъёма (град)
    f_v = x[2] # Фаза выдержки (град)
    f_op = x[3]  # Фаза опускания (град)
    f_z = x[4]  # Фаза теплового зазора (град)
    k_1 = round(x[5])  # коэффициент агрессивности первого участка (выбор зазора)
    k_2 = round(x[6])  # коэффициент агрессивности второго участка (Фаза подъёма)
    k_3 = round(x[7])  # коэффициент агрессивности четвёртого участка (Фаза опускания)
    k_4 = round(x[8])  # коэффициент агрессивности пятого участка (Фаза выбора зазора)
    D = optimize_config.D + 2 * z
    h = optimize_config.h - D / 2
    try:
        config = PolidainConfig(
            N_k=optimize_config.N_k,
            D=D,
            h=h,
            z=z,
            f_pod=f_pod,
            f_v=f_v,
            f_op=f_op,
            f_z=f_z,
            m=m,
            d=d,
            k_1=k_1,
            k_2=k_2,
            k_3=k_3,
            k_4=k_4
        )
        kulachok = Kulachok_polidain(config)
        temp = np.sum(np.abs(kulachok.fun_h(fi_list) - R_func(fi_list)))
        return temp
    except:
        return 1e6

def differential_evolution_optimization(m, d, optimize_config, R_func, differential_evolution_config: DifferentialEvolutionConfig | None = None, bounds_config: BoundsConfig | None = None, N = 1000):
    if differential_evolution_config is None:
        differential_evolution_config = DifferentialEvolutionConfig()
    if bounds_config is None:
        bounds_config = BoundsConfig()
    fi_list = np.linspace(0, 2 * np.pi, N)
    result = differential_evolution(
        fun_optimize_gibrid,
        bounds_config.bounds,
        args = (fi_list, m, d, R_func, optimize_config),
        strategy=differential_evolution_config.strategy,
        maxiter=differential_evolution_config.maxiter,
        popsize=differential_evolution_config.popsize,
        tol=differential_evolution_config.tol,
        workers=differential_evolution_config.workers,
        integrality= differential_evolution_config.integrality,
        updating='deferred'
    )
    print("Best result:", result.x.tolist())
    print("Function value:", result.fun)

def gibrid_optimization(gibrid_optimization_config: GibridOptimizationConfig, optimize_config: OptimizeConfig, R_func: Callable[[np.ndarray], np.ndarray], differential_evolution_config: DifferentialEvolutionConfig | None = None, bounds_config: BoundsConfig | None = None, N: int = 1000):
    for i in gibrid_optimization_config.m:
        for j in gibrid_optimization_config.d:
            print("m:", i, "d:", j)
            differential_evolution_optimization(i, j, optimize_config, R_func, differential_evolution_config = differential_evolution_config, bounds_config = bounds_config, N = N)

if __name__ == '__main__':
    freeze_support()

    from pyzirev_profil import R_func as R_func_pyzirev
    gibrid_optimization(GibridOptimizationConfig(m = [3],
                                                 d = [4, 5, 6, 7, 8]),
                        OptimizeConfig(D = 17.8009 * 2,
                                       h = 30.8018837267781),
                        R_func_pyzirev)
