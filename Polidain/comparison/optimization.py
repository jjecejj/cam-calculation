import os
import sys
import optuna
import numpy as np
from scipy.optimize import minimize
from pyzirev_profil import R_func as R_func_pyzirev
from math import pi
from scipy.optimize import differential_evolution
from multiprocessing import freeze_support

current_file_path = os.path.abspath(__file__)
current_dir = os.path.dirname(current_file_path)
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
from main import Kulachok_polidain, PolidainConfig

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
        fi_list_1 = np.linspace(0, 2 * pi, 1000)
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

def optuna_optimization():
    def objective(trial):
        # Optuna сама предлагает значения
        z = trial.suggest_float("z", 0.01, 0.5)
        f_pod = trial.suggest_float("f_pod", 10 / 180 * pi, 120 / 180 * pi)
        f_v = trial.suggest_float("f_v", 1 / 180 * pi, 20 / 180 * pi)
        f_op = trial.suggest_float("f_op", 10 / 180 * pi, 120 / 180 * pi)
        f_z = trial.suggest_float("f_z", 1 / 180 * pi, 60 / 180 * pi)

        # Целочисленные переменные - Optuna знает, что они целые!
        m = trial.suggest_int("m", 3, 10)
        d = trial.suggest_int("d", 1, 5)
        k_1 = trial.suggest_int("k_1", 1, 100)
        k_2 = trial.suggest_int("k_2", 1, 100)
        k_3 = trial.suggest_int("k_3", 1, 100)
        k_4 = trial.suggest_int("k_4", 1, 100)

        # Собираем x для вашей старой логики (или переписываем fun_optimize, чтобы принимать переменные напрямую)
        x = [z, f_pod, f_v, f_op, f_z, m, d, k_1, k_2, k_3, k_4]

        # Вызов вашей функции (можно немного переписать fun_optimize, чтобы не делать там round снова, но не критично)
        error = fun_optimize(x)

        return error

    # Создаем исследование
    study = optuna.create_study(direction="minimize")
    # Запускаем оптимизацию (например, на 500 попыток)
    study.optimize(objective, n_trials=5000, n_jobs=-1)

    print("Best params:", study.best_params)
    print("Best value:", study.best_value)

def differential_evolution_optimization():
    # 1. Нужно задать границы для каждой переменной (минимум, максимум)
    # Порядок: z, f_pod, f_v, f_op, f_z, m, d, k1, k2, k3, k4
    bounds = [
        (0.01, 0.5),  # z
        (0.01, np.pi),  # f_pod
        (0.01, np.pi / 6),  # f_v
        (0.01, np.pi),  # f_op
        (0.001, np.pi / 2),  # f_z
        (3, 10),  # m (будет округляться внутри)
        (1, 10),  # d
        (1, 10000),  # k1
        (1, 10000),  # k2
        (1, 10000),  # k3
        (1, 10000)  # k4
    ]
    freeze_support()
    result = differential_evolution(
        fun_optimize,
        bounds,
        strategy='best1bin',
        maxiter=1000,
        popsize=10,
        tol=0.01,
        workers=-1  # Использовать все ядра процессора для ускорения
    )
    print("Best result:", (result.x).tolist())
    print("Function value:", result.fun)

def minimize_optimization(type = "Powell"):
    x_0 = np.array([0.25, 80.0 / 180 * pi, 5.0 / 180 * pi, 75.0 / 180 * pi, 25 / 180 * pi, 3, 2, 30, 30, 30, 30])
    res = minimize(fun_optimize, x_0, method='Powell')
    print("Best result:", (res.x).tolist())
    print("Function value:", res.fun)

fi_list_1 = np.linspace(0, 2 * pi, 1000)
fi_list_2 = []
f_dif = 1.1423973285781066
for i in fi_list_1:
    if i - f_dif <= 0:
        fi_list_2.append(2 * pi + i - f_dif)
    else:
        fi_list_2.append(i - f_dif)
fi_list_2 = np.array(fi_list_2)

def fun_optimize_gibrid(x, m = None, d = None):
    z = x[0]  # Тепловой зазор (мм)
    f_pod = x[1]  # Фаза подъёма (град)
    f_v = x[2] # Фаза выдержки (град)
    f_op = x[3]  # Фаза опускания (град)
    f_z = x[4]  # Фаза теплового зазора (град)
    #m = round(x[5])  # степень при C2 только целочисленное и не меньше 3!!!
    #d = round(x[6])  # разность между степенями членов полинома только целочисленное и не меньше 1!!!
    k_1 = round(x[5])  # коэффициент агрессивности первого участка (выбор зазора)
    k_2 = round(x[6])  # коэффициент агрессивности второго участка (Фаза подъёма)
    k_3 = round(x[7])  # коэффициент агрессивности четвёртого участка (Фаза опускания)
    k_4 = round(x[8])  # коэффициент агрессивности пятого участка (Фаза выбора зазора)

    # вычисляемы параметры
    D = 2 * z + 17.8009 * 2  # Базовый диаметр кулака (мм)
    h = 30.8018837267781 - D / 2  # Максимальное перемещение толкателя
    N_k = 1000  # Количество оборотов кулачка в минут
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
        temp = np.sum(np.abs(kulachok.fun_h(fi_list_2) - R_func_pyzirev(fi_list_1)))
        return temp
    except:
        return 1e6

def differential_evolution_optimization(m, d):
    # 1. Нужно задать границы для каждой переменной (минимум, максимум)
    # Порядок: z, f_pod, f_v, f_op, f_z, m, d, k1, k2, k3, k4
    bounds = [
        (0.01, 0.5),  # z
        (0.001, np.pi),  # f_pod
        (0.00001, np.pi / 6),  # f_v
        (0.001, np.pi),  # f_op
        (0.001, np.pi / 6),  # f_z
        (1, 1000),  # k1
        (1, 100),  # k2
        (1, 100),  # k3
        (1, 1000)  # k4
    ]
    freeze_support()
    result = differential_evolution(
        fun_optimize_gibrid,
        bounds,
        args = (m, d),
        strategy='best1bin',
        maxiter=1000,
        popsize=300,
        tol=0.001,
        workers=-1  # Использовать все ядра процессора для ускорения
    )
    print("Best result:", (result.x).tolist())
    print("Function value:", result.fun)

def gibrid_optimization():
    m = [3]
    d = [2, 3, 4]
    for i in m:
        for j in d:
            print("m:", i, "d:", j)
            differential_evolution_optimization(i, j)

if __name__ == '__main__':
    gibrid_optimization()