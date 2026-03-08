import numpy as np
from numpy import ndarray

from scipy.optimize import differential_evolution
from typing import Callable, Literal
from core.cam_geometry.config import KulachokConfig
from core.optimization.config import DifferentialEvolutionConfig, BoundsConfig, OptimizeConfig
from core.profiling_methods import PolinmailCalculator
from core.profiling_methods.polidain.config import PolidainConfig
from core.cam_geometry import Kulachok
from core.profiling_methods.polidain import PolidainCalculator
from core.profiling_methods.polinmail.config import PolinmailConfig, LocalPolinmailConfig
from core.profils import ProfileDataExtractor


def get_differential_evolution_function(profiling_type: str, bounds_config: BoundsConfig):
    if profiling_type == 'polidain':
        return fun_optimize_polidain, bounds_config.bounds_polidain, bounds_config.integrality_polidain
    elif profiling_type == 'polinmail':
        return fun_optimize_polinmail, bounds_config.bounds_polinmail, bounds_config.integrality_polinmail

def fun_optimize_polidain(x: ndarray, fi_list: ndarray = None, R_func: Callable | None = None, optimize_config: OptimizeConfig | None = None):
    """
    Целевая функция оптимизации (минимизируемая функция).
    Рассчитывает отклонение профиля кулачка с параметрами x от целевого профиля R_func.

    Args:
        x: Вектор варьируемых параметров [z, f_pod, f_v, f_op, f_z, k_1, k_2, k_3, k_4].
        fi_list: Массив углов для расчета отклонения.
        R_func: Целевая функция профиля (callable), возвращающая радиус-вектор от угла.
        optimize_config: Конфигурация оптимизации (базовые параметры кулачка).

    Returns:
        float: Сумма модулей разности между рассчитанным и целевым профилем.
               В случае ошибки расчета возвращает большое число (1e6).
    """
    # Основные параметры
    z = x[0]  # Тепловой зазор (мм)
    f_pod = x[1]  # Фаза подъёма (град)
    f_v = x[2] # Фаза выдержки (град)
    f_op = x[3]  # Фаза опускания (град)
    f_z = x[4]  # Фаза теплового зазора (град)

    # Параметры метода профилирования
    m = round(x[5])
    d = round(x[6])
    k_1 = round(x[7])  # коэффициент агрессивности первого участка (выбор зазора)
    k_2 = round(x[8])  # коэффициент агрессивности второго участка (Фаза подъёма)
    k_3 = round(x[9])  # коэффициент агрессивности четвёртого участка (Фаза опускания)
    k_4 = round(x[10])  # коэффициент агрессивности пятого участка (Фаза выбора зазора)

    D = optimize_config.D + 2 * z
    h = optimize_config.h
    try:
        config = KulachokConfig(
            N_k=optimize_config.N_k,
            D=D,
            h=h,
            z=z,
            f_pod=f_pod,
            f_v=f_v,
            f_op=f_op,
            f_z=f_z,
        )
        polidain_config = PolidainConfig(
            m=m,
            d=d,
            k_1=k_1,
            k_2=k_2,
            k_3=k_3,
            k_4=k_4
        )
        polidain_calculator = PolidainCalculator(polidain_config)
        kulachok = Kulachok(config, polidain_calculator)
        temp = np.sum(np.abs(kulachok.fun_h(fi_list) - R_func(fi_list)))
        print(temp)
        return temp
    except Exception as e:
        print(f"Произошла ошибка: {e}")
        return 1e6

def fun_optimize_polinmail(x: ndarray, fi_list: ndarray = None, R_func: Callable | None = None, optimize_config: OptimizeConfig | None = None):
    """
    Целевая функция оптимизации (минимизируемая функция).
    Рассчитывает отклонение профиля кулачка с параметрами x от целевого профиля R_func.

    Args:
        x: Вектор варьируемых параметров [z, f_pod, f_v, f_op, f_z, k_1, k_2, k_3, k_4].
        fi_list: Массив углов для расчета отклонения.
        R_func: Целевая функция профиля (callable), возвращающая радиус-вектор от угла.
        optimize_config: Конфигурация оптимизации (базовые параметры кулачка).

    Returns:
        float: Сумма модулей разности между рассчитанным и целевым профилем.
               В случае ошибки расчета возвращает большое число (1e6).
    """
    # Основные параметры
    z = x[0]  # Тепловой зазор (мм)
    f_pod = x[1]  # Фаза подъёма (град)
    f_v = x[2] # Фаза выдержки (град)
    f_op = x[3]  # Фаза опускания (град)
    f_z = x[4]  # Фаза теплового зазора (град)

    # Параметры метода профилирования
    m_1 = round(x[5])
    d_1 = round(x[6])
    m_2 = round(x[7])
    d_2 = round(x[8])
    m_3 = round(x[9])
    d_3 = round(x[10])
    m_4 = round(x[11])
    d_4 = round(x[12])

    D = optimize_config.D + 2 * z
    h = optimize_config.h
    try:
        config = KulachokConfig(
            N_k=optimize_config.N_k,
            D=D,
            h=h,
            z=z,
            f_pod=f_pod,
            f_v=f_v,
            f_op=f_op,
            f_z=f_z,
        )
        config_1 = LocalPolinmailConfig(
            m = m_1,
            d = d_1
        )
        config_2 = LocalPolinmailConfig(
            m=m_2,
            d=d_2
        )
        config_3 = LocalPolinmailConfig(
            m=m_3,
            d=d_3
        )
        config_4 = LocalPolinmailConfig(
            m=m_4,
            d=d_4
        )
        polinmail_config = PolinmailConfig(
            config_1=config_1,
            config_2=config_2,
            config_3=config_3,
            config_4=config_4
        )
        polinmail_calculator = PolinmailCalculator(polinmail_config)
        kulachok = Kulachok(config, polinmail_calculator)
        temp = np.sum(np.abs(kulachok.fun_h(fi_list) - R_func(fi_list)))
        print(temp)
        return temp
    except Exception as e:
        print(f"Произошла ошибка: {e}")
        return 1e6

def differential_evolution_optimization(optimize_config: OptimizeConfig,
                                        differential_evolution_config: DifferentialEvolutionConfig = DifferentialEvolutionConfig(),
                                        bounds_config: BoundsConfig = BoundsConfig(),
                                        N: int = 1000):
    """
    Запускает дифференциальную эволюцию для поиска оптимальных параметров кулачка.

    Args:
        optimize_config: Конфигурация оптимизации.
        differential_evolution_config: Настройки алгоритма дифференциальной эволюции.
        bounds_config: Границы изменения параметров.
        N: Количество точек для сравнения профилей.
    """
    fi_list = np.linspace(0, 2 * np.pi, N)
    fun, bounds, integrality, = get_differential_evolution_function(optimize_config.profiling_type, bounds_config)
    result = differential_evolution(
        fun,
        bounds,
        args = (fi_list, optimize_config.R_func, optimize_config),
        strategy=differential_evolution_config.strategy,
        maxiter=differential_evolution_config.maxiter,
        popsize=differential_evolution_config.popsize,
        tol=differential_evolution_config.tol,
        workers=differential_evolution_config.workers,
        integrality= integrality,
        updating='deferred'
    )
    print("Best result:", result.x.tolist())
    print("Function value:", result.fun)

def get_calculator_config_from_x(x, D, h, calculator_type: str):
    config = KulachokConfig(
        N_k=1000,
        D=D,
        h=h,
        z=x[0],
        f_pod=x[1],
        f_v=x[2],
        f_op=x[3],
        f_z=x[4],
    )
    if calculator_type == 'polidain':
        calculator_config = PolidainConfig(
            m=x[5],
            d=x[6],
            k_1=x[7],
            k_2=x[8],
            k_3=x[9],
            k_4=x[10]
        )
    elif calculator_type == 'polinmail':
        local_config_1 = LocalPolinmailConfig(m=x[5], d=x[6])
        local_config_2 = LocalPolinmailConfig(m=x[7], d=x[8])
        local_config_3 = LocalPolinmailConfig(m=x[9], d=x[10])
        local_config_4 = LocalPolinmailConfig(m=x[11], d=x[12])
        calculator_config = PolinmailConfig(
            config_1=local_config_1,
            config_2=local_config_2,
            config_3=local_config_3,
            config_4=local_config_4,
        )
    return config, calculator_config