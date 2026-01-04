import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import math
from math import pi
from pyzirev_profil import R_func as R_func_pyzirev
from pyzirev_profil import X as X_pyzirev
from pyzirev_profil import Y as Y_pyzirev

current_file_path = os.path.abspath(__file__)
current_dir = os.path.dirname(current_file_path)
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
from main import Kulachok_polidain, PolidainConfig


def create_config(x, N = 1000):
    fi_list = np.linspace(0, 2 * pi, N)
    # Варьируемые парамметры параметры
    z = x[0]  # Тепловой зазор (мм)
    f_pod = x[1]  # Фаза подъёма (град)
    f_v = x[2]  # Фаза выдержки (град)
    f_v = 0.001
    f_op = x[3]  # Фаза опускания (град)
    f_z = x[4]  # Фаза теплового зазора (град)
    m = round(x[5])  # степень при C2 только целочисленное и не меньше 3!!!
    d = round(x[6])  # разность между степенями членов полинома только целочисленное и не меньше 1!!!
    k_1 = round(x[7])  # коэффициент агрессивности первого участка (выбор зазора)
    k_2 = round(x[8])  # коэффициент агрессивности второго участка (Фаза подъёма)
    k_3 = round(x[9])  # коэффициент агрессивности четвёртого участка (Фаза опускания)
    k_4 = round(x[10])  # коэффициент агрессивности пятого участка (Фаза выбора зазора)

    # вычисляемы параметры
    D = 2 * z + 17.8009 * 2  # Базовый диаметр кулака (мм)
    h = np.max(R_func_pyzirev(fi_list)) - D / 2  # Максимальное перемещение толкателя
    N_k = 1000  # Количество оборотов кулачка в минут

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
    return config

def graf_profil(x, f_dif = 1.1423973285781066, N = 1000):
    kulachok = Kulachok_polidain(create_config(x))
    fi_list_1 = np.linspace(0, 2 * pi, N)
    fi_list_2 = []
    for i in fi_list_1:
        if i - f_dif <= 0:
            fi_list_2.append(2 * pi + i - f_dif)
        else:
            fi_list_2.append(i - f_dif)
    fi_list_2 = np.array(fi_list_2)
    X = []
    Y = []
    for i in range(0, len(fi_list_1)):
        if fi_list_1[i] <= np.pi:
            X.append(kulachok.fun_h(fi_list_2[i]) * np.cos(fi_list_1[i]))
            Y.append(kulachok.fun_h(fi_list_2[i]) * np.sin(fi_list_1[i]))
        else:
            X.append(kulachok.fun_h(fi_list_2[i]) * np.cos(2 * np.pi - fi_list_1[i]))
            Y.append(-kulachok.fun_h(fi_list_2[i]) * np.sin(2 * np.pi - fi_list_1[i]))

    plt.figure(figsize=(6, 6))
    plt.plot(X, Y)
    plt.plot(X_pyzirev, Y_pyzirev)
    plt.scatter([0], [0])
    plt.xlabel('X, м')
    plt.ylabel('Y, м')
    #plt.xlim(np.min(X_pyzirev) - 0.01, np.max(X_pyzirev) + 0.01)
    #plt.ylim(np.min(Y) - 0.01, np.max(Y) + 0.01)
    plt.grid(True)
    plt.title("Профиль кулачка")
    plt.show()

def graf_h(x, f_dif = 1.1423973285781066, N = 1000):
    kulachok = Kulachok_polidain(create_config(x))
    fi_list_1 = np.linspace(0, 2 * pi, N)
    fi_list_2 = []
    for i in fi_list_1:
        if i - f_dif <= 0:
            fi_list_2.append(2 * pi + i - f_dif)
        else:
            fi_list_2.append(i - f_dif)
    fi_list_2 = np.array(fi_list_2)

    plt.plot(fi_list_1, kulachok.fun_h(fi_list_2), c='r')
    plt.plot(fi_list_1, R_func_pyzirev(fi_list_1), c='b')
    plt.grid(True)
    plt.show()


def print_parameters(x):
    # Распаковка параметров (как в твоем коде)
    z = x[0]
    f_pod = x[1]
    f_v = x[2]
    f_op = x[3]
    f_z = x[4]
    m = round(x[5])
    d = round(x[6])
    k_1 = round(x[7])
    k_2 = round(x[8])
    k_3 = round(x[9])
    k_4 = round(x[10])

    # Вывод в консоль с форматированием
    print("=" * 75)
    print(f"{'Параметр':<8} | {'Значение (Raw)':<14} | {'Доп. инфо':<15} | {'Описание'}")
    print("=" * 75)

    # Линейные параметры
    print(f"{'z':<8} | {z:<14.5f} | {'-':<15} | Тепловой зазор (мм)")

    # Угловые параметры (добавляем перевод в градусы)
    print(f"{'f_pod':<8} | {f_pod:<14.5f} | {math.degrees(f_pod):6.2f}°        | Фаза подъёма")
    print(f"{'f_v':<8} | {f_v:<14.5f} | {math.degrees(f_v):6.2f}°        | Фаза выдержки")
    print(f"{'f_op':<8} | {f_op:<14.5f} | {math.degrees(f_op):6.2f}°        | Фаза опускания")
    print(f"{'f_z':<8} | {f_z:<14.5f} | {math.degrees(f_z):6.2f}°        | Фаза теплового зазора")

    # Целочисленные параметры (структурные)
    print("-" * 75)
    print(f"{'m':<8} | {m:<14} | {'int':<15} | Степень полинома при C2")
    print(f"{'d':<8} | {d:<14} | {'int':<15} | Разность степеней")

    # Коэффициенты агрессивности
    print("-" * 75)
    print(f"{'k_1':<8} | {k_1:<14} | {'int':<15} | Агрессивность (выбор зазора)")
    print(f"{'k_2':<8} | {k_2:<14} | {'int':<15} | Агрессивность (подъём)")
    print(f"{'k_3':<8} | {k_3:<14} | {'int':<15} | Агрессивность (опускание)")
    print(f"{'k_4':<8} | {k_4:<14} | {'int':<15} | Агрессивность (выбор зазора 2)")
    print("=" * 75)

if __name__ == '__main__':
    import matplotlib_settings.profile_1

    x = [0.2158609144273846, 1.4752676726321512, 0.3659194493259062, 1.4734837012705113, 0.34118203634644295, 3.1956975577607416, 5.213885042319334, 6475.819912613473, 8222.264738671533, 5111.141550658319, 5282.439747507636]
    x = [0.2605034729471177, 1.4625876049540705, 0.07314946299807315, 1.7692621299837572, 0.35350901825132053, 3.084325638718698, 9.736116587649615, 6169.004306957386, 6564.0020771292075, 13.259848424733718, 8079.336211740891]
    x = [0.27281413896589235, 1.6266298614374053, 1e-05, 1.6263067704268104, 0.37264120880113805, 3, 2, 982.0397766773281, 17.651854389975107, 17.52268688338083, 979.5145395794359]
    print_parameters(x)
    graf_profil(x)
    graf_h(x)



