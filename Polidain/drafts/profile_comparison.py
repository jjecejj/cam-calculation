import numpy as np
from attr import dataclass
from scipy import optimize
from scipy.optimize import minimize

from main import Kulachok_polidain, PolidainConfig
from pyzirev_profil import R_func as R_func_pyzirev
from math import pi

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


if __name__ == '__main__':
    '''
    x_0 = np.array([0.25, 80.0 / 180 * pi, 5.0 / 180 * pi, 75.0 / 180 * pi, 25 / 180 * pi, 3, 2, 30, 30, 30, 30])
    res = minimize(fun_optimize, x_0, method='Powell')
    print(res.x)
    '''
    #x = res.x
    x = [2.62528480e-01, 1.57484468e+00, 1.41363763e-01, 1.55687019e+00, 3.66956177e-01,
         2.91138287e+00, 1.25323327e+01, 9.33149535e+13,
 2.87852183e+01, 3.57689408e+01, 2.24511846e+14]
    fi_list = np.linspace(0, 2 * pi, 1000)
    # Варьируемые парамметры параметры
    z = x[0]  # Тепловой зазор (мм)
    f_pod = x[1]  # Фаза подъёма (град)
    f_v = x[2]  # Фаза выдержки (град)
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

    # разность начальной фазы
    f_dif = 1.1423973285781066

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

    import matplotlib.pyplot as plt
    from pyzirev_profil import X as X_pyzirev
    from pyzirev_profil import Y as Y_pyzirev

    plt.plot(fi_list_1, np.vectorize(kulachok.fun_h)(fi_list_2), c='r')
    plt.plot(fi_list_1, R_func_pyzirev(fi_list_1), c='b')
    plt.show()

    X = np.vectorize(kulachok.fun_h)(fi_list_2) * np.sin(fi_list_2 - (2 * pi - f_dif))
    Y = np.vectorize(kulachok.fun_h)(fi_list_2) * np.cos(fi_list_2 - (2 * pi - f_dif))

    plt.figure(figsize=(6, 6))
    plt.plot(X, Y)
    plt.plot(X_pyzirev, Y_pyzirev)
    plt.scatter([0], [0])
    plt.xlabel('X, м')
    plt.ylabel('Y, м')
    plt.xlim(np.min(X) - 0.01, np.max(X) + 0.01)
    plt.ylim(np.min(Y) - 0.01, np.max(Y) + 0.01)
    plt.grid(True)
    plt.title("Профиль кулачка")
    plt.show()



