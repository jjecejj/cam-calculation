import numpy as np
from Polidain.main import Kulachok_polidain, PolidainConfig
from pyzirev_profil import R_func as R_func_pyzirev
from math import pi
from multiprocessing import freeze_support

def fun_optimize(x):
    z = x[0]  # Тепловой зазор (мм)
    f_pod = x[1]  # Фаза подъёма (град)
    f_v = x[2] # Фаза выдержки (град)
    f_op = x[3]  # Фаза опускания (град)
    f_z = x[4]  # Фаза теплового зазора (град)
    m = round(x[5])  # степень при C2 только целочисленное и не меньше 3!!!
    d = round(x[6])  # разность между степенями членов полинома только целочисленное и не меньше 1!!!
    k_1 = round(x[7])  # коэффициент агрессивности первого участка (выбор зазора)
    k_2 = round(x[8])  # коэффициент агрессивности второго участка (Фаза подъёма)
    k_3 = round(x[9]) # коэффициент агрессивности четвёртого участка (Фаза опускания)
    k_4 = round(x[10])  # коэффициент агрессивности пятого участка (Фаза выбора зазора)

    # вычисляемы параметры
    D = 2 * z + 17.8009 * 2  # Базовый диаметр кулака (мм)
    h = 30.8018837267781 - D / 2  # Максимальное перемещение толкателя
    N_k = 1000  # Количество оборотов кулачка в минут

    f_dif = 1.1423973285781066
    try:
        config = PolidainConfig(
            N_k=N_k,
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
        fi_list_1 = np.linspace(0, 2 * pi, 2000)
        fi_list_2 = []
        for i in fi_list_1:
            if i - f_dif <= 0:
                fi_list_2.append(2 * pi + i - f_dif)
            else:
                fi_list_2.append(i - f_dif)
        fi_list_2 = np.array(fi_list_2)
        temp = np.sum(np.abs(np.vectorize(kulachok.fun_h)(fi_list_2) - R_func_pyzirev(fi_list_1)))
        print(temp)
        return temp
    except:
        return 1e6

from scipy.optimize import differential_evolution

# 1. Нужно задать границы для каждой переменной (минимум, максимум)
# Порядок: z, f_pod, f_v, f_op, f_z, m, d, k1, k2, k3, k4
bounds = [
    (0.01, 0.5),          # z
    (0.01, np.pi),          # f_pod
    (0.01, np.pi/6),        # f_v
    (0.01, np.pi),          # f_op
    (0.001, np.pi/2),        # f_z
    (3, 10),             # m (будет округляться внутри)
    (1, 10),              # d
    (1, 10000),            # k1
    (1, 10000),            # k2
    (1, 10000),            # k3
    (1, 10000)             # k4
]

# Ваша функция остается такой же, differential_evolution сам перебирает значения.
# strategy='best1bin' - классическая стратегия, popsize=15 - размер популяции
if __name__ == '__main__':
    freeze_support()
    result = differential_evolution(
        fun_optimize,
        bounds,
        strategy='best1bin',
        maxiter=150,
        popsize=25,
        tol=0.01,
        workers=-1  # Использовать все ядра процессора для ускорения
    )

    print("Best result:", (result.x).tolist())
    print("Function value:", result.fun)