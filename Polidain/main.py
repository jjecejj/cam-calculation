from math import sin, cos, pi
import numpy as np
from pydantic import BaseModel, Field, model_validator
from typing import Callable, List
import matplotlib.pyplot as plt

class ValidationError(Exception):
    pass

def k_fun(p, m, d):
    """
        :param m: Младшая степень полинома
        :param p: Степень второго члена полинома
        :param d: Разность между степенями членов полинома
        :return: [m, p, q, r, s] - массив степеней членов полинома
        """
    q = p + p - d
    r = q + p - d
    s = r + p - d
    return [m, p, q, r, s]

def C_fun(p, m, d):
    """
    :param m: Младшая степень полинома
    :param p: Степень второго члена полинома
    :param d: Разность между степенями членов полинома
    :return: [C2, Cp, Cq, Cr, Cs] - массив коэффициентов при членах полиномма
    """
    temp = k_fun(p, m, d)
    q = temp[2]
    r = temp[3]
    s = temp[4]
    C2 = -p * q * r * s / ((p - m) * (q - m) * (r - m) * (s - m))
    Cp = m * q * r * s / ((p - m) * (q - p) * (r - p) * (s - p))
    Cq = -m * p * s * r / ((q - m) * (r - q) * (q - p) * (s - q))
    Cr = m * p * q * s / ((r - m) * (s - r) * (r - p) * (r - q))
    Cs = -m * p * q * r / ((s - m) * (s - p) * (s - q) * (s - r))
    return [C2, Cp, Cq, Cr, Cs]

def h_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    '''
    Функция перемещения от угла поворота кулачка
    :param fi: Угол поворота кулачка
    :param C_list:Массив коэффициентов при членах полиномма
    :param k_list:Массив степеней членов полинома
    :param fi_1:Конец характерного участка
    :param fi_0:Начало характерного участка
    :param h_kn_max: максимальная высота подёма кулачка на характерном участке
    '''
    temp = 1
    for i in range(0, len(C_list)):
        if k_list[i] - 0 < 0:
            continue
        temp += C_list[i] * (((fi - fi_0) / (fi_1 - fi_0)) ** k_list[i])
    return temp * h_kn_max

def v_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    '''Функция скорости от угла поворота кулачка
    :param fi: Угол поворота кулачка
    :param C_list:Массив коэффициентов при членах полиномма
    :param k_list:Массив степеней членов полинома
    :param fi_1:Конец характерного участка
    :param fi_0:Начало характерного участка
    :param h_kn_max: максимальная высота подёма кулачка на характерном участке
    '''
    temp = 0
    for i in range(0, len(C_list)):
        if k_list[i] - 1 < 0:
            continue
        temp += k_list[i] * C_list[i] * (((fi - fi_0) / (fi_1 - fi_0)) ** (k_list[i] - 1))
    # return temp * h_kn_max * (omega / (fi_1 - fi_0))
    return temp * h_kn_max

def a_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    '''Функция ускорения от угла поворота кулачка
    :param fi: Угол поворота кулачка
    :param C_list:Массив коэффициентов при членах полиномма
    :param k_list:Массив степеней членов полинома
    :param fi_1:Конец характерного участка
    :param fi_0:Начало характерного участка
    :param h_kn_max: максимальная высота подёма кулачка на характерном участке
    '''
    temp = 0
    for i in range(0, len(C_list)):
        if k_list[i] - 2 < 0:
            continue
        temp += k_list[i] * (k_list[i] - 1) * C_list[i] * (((fi - fi_0) / (fi_1 - fi_0)) ** (k_list[i] - 2))
    return temp * h_kn_max

def d_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    '''Функция рывка от угла поворота кулачка
    :param fi: Угол поворота кулачка
    :param C_list:Массив коэффициентов при членах полиномма
    :param k_list:Массив степеней членов полинома
    :param fi_1:Конец характерного участка
    :param fi_0:Начало характерного участка
    :param h_kn_max: максимальная высота подёма кулачка на характерном участке
    '''
    temp = 0
    for i in range(0, len(C_list)):
        if k_list[i] - 3 < 0:
            continue
        temp += k_list[i] * (k_list[i] - 1) * (k_list[i] - 2) * C_list[i] * (
                    ((fi - fi_0) / (fi_1 - fi_0)) ** (k_list[i] - 3))
    return temp * h_kn_max

def k_phi(fi, C_list, k_list, fi_1, fi_0, h_kn_max):
    '''Функция кракена от угла поворота кулачка
    :param fi: Угол поворота кулачка
    :param C_list:Массив коэффициентов при членах полиномма
    :param k_list:Массив степеней членов полинома
    :param fi_1:Конец характерного участка
    :param fi_0:Начало характерного участка
    :param h_kn_max: максимальная высота подёма кулачка на характерном участке
    '''
    temp = 0
    for i in range(0, len(C_list)):
        if k_list[i] - 4 < 0:
            continue
        temp += k_list[i] * (k_list[i] - 1) * (k_list[i] - 2) * (k_list[i] - 3) * C_list[i] * (
                    ((fi - fi_0) / (fi_1 - fi_0)) ** (k_list[i] - 4))
    return temp * h_kn_max

class PolidainConfig(BaseModel):
    """
    Класс-валидатор входных данных
    Проверяет типы и ограничения
    Передаваемы данные:
        'N_k':Количество оборотов кулачка в минут\n
        'D':Базовый диаметр кулака (м)\n
        'h':Максимальное перемещение толкателя (м)\n
        'z':Тепловой зазор (м)\n
        'f_pod':Фаза подъёма (рад)\n
        'f_v':Фаза выдержки (рад)\n
        'f_op':Фаза опускания (рад)\n
        'f_z':Фаза теплового зазора (рад)
        'm':Степень при C2 только целочисленное и не меньше 3\n
        'd':Разность между степенями членов полинома только целочисленное и не меньше 1\n
        'k_1':Коэффициент агрессивности первого участка (выбор зазора)\n
        'k_2':Коэффициент агрессивности второго участка (Фаза подъёма)\n
        'k_3':Коэффициент агрессивности четвёртого участка (Фаза опускания)\n
        'k_4':Коэффициент агрессивности пятого участка (Фаза выбора зазора)
    """
    # Общие параметры
    N_k: float = Field(..., gt=0, description="Обороты в минуту")
    D: float = Field(..., gt=0, description="Базовый диаметр")
    h: float = Field(..., gt=0, description="Максимальное перемещение")
    z: float = Field(..., ge=0, description="Тепловой зазор")

    # Фазовые углы
    f_pod: float = Field(..., gt=0, description="Фаза подъёма (рад)")
    f_v: float = Field(..., gt=0, description=" Фаза выдержки (рад)")
    f_op: float = Field(..., gt=0, description="Фаза опускания (рад)")
    f_z: float = Field(..., gt=0, description="Фаза теплового зазора (рад)")

    # Параметры полинома
    m: int = Field(..., ge=3, description="Степень при C2, не меньше 3")
    d: int = Field(..., ge=1, description="Разность степеней, не меньше 1")

    # Коэффициенты агрессивности
    k_1: int = Field(..., gt=0, description="Коэффициент агрессивности первого участка (выбор зазора)")
    k_2: int = Field(..., gt=0, description="Коэффициент агрессивности второго участка (Фаза подъёма)")
    k_3: int = Field(..., gt=0, description="Коэффициент агрессивности четвёртого участка (Фаза опускания)")
    k_4: int = Field(..., gt=0, description="Коэффициент агрессивности пятого участка (Фаза выбора зазора)")

    @model_validator(mode='after')
    def check_k_less_than_d(self):
        """Проверяет, что все коэффициенты k строго меньше d"""
        # Создаем словарь для удобной проверки в цикле
        k_values = {
            'k_1': self.k_1,
            'k_2': self.k_2,
            'k_3': self.k_3,
            'k_4': self.k_4
        }

        # Проходим по всем k
        for name, value in k_values.items():
            # Если условие нарушено (k <= d)
            if value <= self.d:
                raise ValueError(
                    f"Ошибка в параметре {name}: значение ({value}) "
                    f"должно быть строго больше d ({self.d})"
                )

        return self

    # --- Вычисляемые свойства (Автоматический расчет) ---
    @property
    def phi_0(self) -> float:
        return 0

    @property
    def phi_1(self) -> float:
        return self.phi_0 + self.f_z

    @property
    def phi_2(self) -> float:
        return self.phi_1 + self.f_pod

    @property
    def phi_3(self) -> float:
        return self.phi_2 + self.f_v

    @property
    def phi_4(self) -> float:
        return self.phi_3 + self.f_op

    @property
    def phi_5(self) -> float:
        return self.phi_4 + self.f_z

    @property
    def omega(self) -> float:
        """Угловая скорость (рад/с)"""
        return self.N_k * 2 * pi / 60

    @property
    def T(self) -> float:
        """Период оборота (с)"""
        return 60 / self.N_k

    @property
    def r0(self) -> float:
        """Базовый радиус (м)"""
        return self.D / 2

    @property
    def C_list_1(self) -> list:
        return C_fun(self.k_1, self.m, self.d)

    @property
    def k_list_1(self) -> list:
        return k_fun(self.k_1, self.m, self.d)

    @property
    def C_list_2(self) -> list:
        return C_fun(self.k_2, self.m, self.d)

    @property
    def k_list_2(self) -> list:
        return k_fun(self.k_2, self.m, self.d)

    @property
    def C_list_3(self) -> list:
        return C_fun(self.k_3, self.m, self.d)

    @property
    def k_list_3(self) -> list:
        return k_fun(self.k_3, self.m, self.d)

    @property
    def C_list_4(self) -> list:
        return C_fun(self.k_4, self.m, self.d)

    @property
    def k_list_4(self) -> list:
        return k_fun(self.k_4, self.m, self.d)


class Kulachok_polidain:
    def __init__(self, config:PolidainConfig):
        self.config = config

    def fun_universal(self, fi: float, fun: Callable, sign_list: List[int], const_list: List[float]) -> float:
        """
        Универсальная функция для расчета кинематики.
        :param fi: Текущий угол (рад)
        :param fun: Ссылка на мат. функцию (h_phi, v_phi...)
        :param sign_list: Список знаков для участков [уч1, уч2, уч4, уч5]
        :param const_list: Список констант (смещений) для каждого участка + остаток
        """
        # Валидация входного угла
        if fi < 0 or fi > 2 * pi:
            raise ValidationError(f"Недопустимое значение угла: {fi:.4f} рад")

        # Участок 1: Выбор зазора
        elif fi < self.config.phi_1:
            val = fun(fi, self.config.C_list_1, self.config.k_list_1,
                      self.config.phi_1, fi_0=0, h_kn_max=self.config.z)
            return sign_list[0] * val + const_list[0]

        # Участок 2: Подъем
        elif fi < self.config.phi_2:
            val = fun(fi, self.config.C_list_2, self.config.k_list_2,
                      self.config.phi_2, fi_0=self.config.phi_1, h_kn_max=self.config.h)
            return sign_list[1] * val + const_list[1]

        # Участок 3: Выдержка
        elif fi < self.config.phi_3:
            return const_list[2]

        # Участок 4: Опускание (с инверсией времени)
        elif fi < self.config.phi_4:
            # Аргумент времени: phi_4 - fi + phi_3
            args = self.config.phi_4 - fi + self.config.phi_3
            val = fun(args, self.config.C_list_3, self.config.k_list_3,
                      self.config.phi_4, fi_0=self.config.phi_3, h_kn_max=self.config.h)
            return sign_list[2] * val + const_list[3]

        # Участок 5: Выбор зазора (спуск)
        elif fi < self.config.phi_5:
            args = self.config.phi_5 - fi + self.config.phi_4
            val = fun(args, self.config.C_list_4, self.config.k_list_4,
                      self.config.phi_5, fi_0=self.config.phi_4, h_kn_max=self.config.z)
            return sign_list[3] * val + const_list[4]

        # Остаток цикла
        elif fi <= 2 * pi:
            return const_list[5]

    def fun_h(self, fi):
        return self.fun_universal(fi, h_phi, [-1, -1, -1, -1], [self.config.r0,
                                                                     self.config.r0 + self.config.h,
                                                                     self.config.r0 + self.config.h,
                                                                     self.config.r0 + self.config.h,
                                                                     self.config.r0,
                                                                     self.config.r0])
    def fun_v(self, fi):
        return self.fun_universal(fi, v_phi, [1, 1, 1, 1], [0, 0, 0, 0, 0, 0])

    def fun_a(self, fi):
        return self.fun_universal(fi, a_phi, [-1, -1, -1, -1], [0, 0, 0, 0, 0, 0])

    def fun_d(self, fi):
        return self.fun_universal(fi, d_phi, [1, 1, 1, 1], [0, 0, 0, 0, 0, 0])

    def fun_k(self, fi):
        return self.fun_universal(fi, k_phi, [-1, -1, -1, -1], [0, 0, 0, 0, 0, 0])

    def display_graphs(self, N = 1000):
        fi_list = np.linspace(0, 2 * pi, N)
        H = np.vectorize(self.fun_h)(fi_list)
        V = np.vectorize(self.fun_v)(fi_list)
        A = np.vectorize(self.fun_a)(fi_list)
        D = np.vectorize(self.fun_d)(fi_list)
        K = np.vectorize(self.fun_k)(fi_list)
        fi_list = fi_list / pi * 180

        fig, axs = plt.subplots(5, 1, figsize=(8, 20))

        # --- 1. Координата ---
        axs[0].plot(fi_list, H)
        axs[0].scatter(fi_list[H.argmax()], H.max(), color='r')
        axs[0].set_xlabel(r'$/phi$, град')
        axs[0].set_ylabel('Координата толкателя, мм')
        axs[0].grid(True)

        # --- 2. Скорость ---
        axs[1].plot(fi_list, V)
        axs[1].scatter(fi_list[V.argmax()], V.max(), color='r')
        axs[1].scatter(fi_list[V.argmin()], V.min(), color='r')
        axs[1].set_xlabel(r'$/phi$, град')
        axs[1].set_ylabel('Скорость, $мм/с$')
        axs[1].grid(True)

        # --- 3. Ускорение ---
        axs[2].plot(fi_list, A)
        axs[2].scatter(fi_list[A.argmax()], A.max(), color='r')
        axs[2].scatter(fi_list[A.argmin()], A.min(), color='r')
        axs[2].set_xlabel(r'$/phi$, град')
        axs[2].set_ylabel('Ускорение, $мм/с^2$')
        axs[2].grid(True)

        # --- 4. Рывок ---
        axs[3].plot(fi_list, D)
        axs[3].scatter(fi_list[D.argmax()], D.max(), color='r')
        axs[3].scatter(fi_list[D.argmin()], D.min(), color='r')
        axs[3].set_xlabel(r'$/phi$, град')
        axs[3].set_ylabel('Рывок, $мм/с^3$')
        axs[3].grid(True)

        # --- 5. "Кракен" ---
        axs[4].plot(fi_list, K)
        axs[4].scatter(fi_list[K.argmax()], K.max(), color='r')
        axs[4].scatter(fi_list[K.argmin()], K.min(), color='r')
        axs[4].set_xlabel(r'$/phi$, град')
        axs[4].set_ylabel('Кракен, $мм/с^4$')
        axs[4].grid(True)

        plt.tight_layout()
        plt.show()

'''
    def fun_h(self, fi):
        if fi < 0 or fi > 2 * pi:
            raise ValidationError("Недопустимое значение угла поворота кулачка")
        elif fi < self.config.phi_1:
            return -h_phi(fi, self.config.C_list_1, self.config.k_list_1, self.config.phi_1, fi_0=0,h_kn_max=self.config.z) + self.config.r0
        elif fi < self.config.phi_2:
            return -h_phi(fi, self.config.C_list_2, self.config.k_list_2, self.config.phi_2, fi_0=self.config.phi_1, h_kn_max=self.config.h) + self.config.r0 + self.config.h
        elif fi < self.config.phi_3:
            return self.config.r0 + self.config.h
        elif fi < self.config.phi_4:
            return -h_phi(self.config.phi_4 - fi + self.config.phi_3, self.config.C_list_3, self.config.k_list_3, self.config.phi_4, fi_0=self.config.phi_3, h_kn_max=self.config.h) + self.config.r0 + self.config.h
        elif fi < self.config.phi_5:
            return -h_phi(self.config.phi_5 - fi + self.config.phi_4, self.config.C_list_4, self.config.k_list_4, self.config.phi_5, fi_0=self.config.phi_4, h_kn_max = self.config.z) + self.config.r0
        elif fi < pi * 2:
            return self.config.r0

    def fun_v(self, fi):
        if fi < 0 or fi > 2 * pi:
            raise ValidationError("Недопустимое значение угла поворота кулачка")
        elif fi < self.config.phi_1:
            return v_phi(fi, self.config.C_list_1, self.config.k_list_1, self.config.phi_1, fi_0=0,h_kn_max=self.config.z)
        elif fi < self.config.phi_2:
            return v_phi(fi, self.config.C_list_2, self.config.k_list_2, self.config.phi_2, fi_0=self.config.phi_1, h_kn_max=self.config.h)
        elif fi < self.config.phi_3:
            return 0
        elif fi < self.config.phi_4:
            return v_phi(self.config.phi_4 - fi + self.config.phi_3, self.config.C_list_3, self.config.k_list_3, self.config.phi_4, fi_0=self.config.phi_3, h_kn_max=self.config.h)
        elif fi < self.config.phi_5:
            return v_phi(self.config.phi_5 - fi + self.config.phi_4, self.config.C_list_4, self.config.k_list_4, self.config.phi_5, fi_0=self.config.phi_4, h_kn_max = self.config.z)
        elif fi < pi * 2:
            return 0

    def fun_a(self, fi):
        if fi < 0 or fi > 2 * pi:
            raise ValidationError("Недопустимое значение угла поворота кулачка")
        elif fi < self.config.phi_1:
            return -a_phi(fi, self.config.C_list_1, self.config.k_list_1, self.config.phi_1, fi_0=0,h_kn_max=self.config.z)
        elif fi < self.config.phi_2:
            return -a_phi(fi, self.config.C_list_2, self.config.k_list_2, self.config.phi_2, fi_0=self.config.phi_1, h_kn_max=self.config.h)
        elif fi < self.config.phi_3:
            return 0
        elif fi < self.config.phi_4:
            return -a_phi(self.config.phi_4 - fi + self.config.phi_3, self.config.C_list_3, self.config.k_list_3, self.config.phi_4, fi_0=self.config.phi_3, h_kn_max=self.config.h)
        elif fi < self.config.phi_5:
            return -a_phi(self.config.phi_5 - fi + self.config.phi_4, self.config.C_list_4, self.config.k_list_4, self.config.phi_5, fi_0=self.config.phi_4, h_kn_max = self.config.z)
        elif fi < pi * 2:
            return 0

    def fun_d(self, fi):
        if fi < 0 or fi > 2 * pi:
            raise ValidationError("Недопустимое значение угла поворота кулачка")
        elif fi < self.config.phi_1:
            return d_phi(fi, self.config.C_list_1, self.config.k_list_1, self.config.phi_1, fi_0=0,h_kn_max=self.config.z)
        elif fi < self.config.phi_2:
            return d_phi(fi, self.config.C_list_2, self.config.k_list_2, self.config.phi_2, fi_0=self.config.phi_1, h_kn_max=self.config.h)
        elif fi < self.config.phi_3:
            return 0
        elif fi < self.config.phi_4:
            return d_phi(self.config.phi_4 - fi + self.config.phi_3, self.config.C_list_3, self.config.k_list_3, self.config.phi_4, fi_0=self.config.phi_3, h_kn_max=self.config.h)
        elif fi < self.config.phi_5:
            return d_phi(self.config.phi_5 - fi + self.config.phi_4, self.config.C_list_4, self.config.k_list_4, self.config.phi_5, fi_0=self.config.phi_4, h_kn_max = self.config.z)
        elif fi < pi * 2:
            return 0

    def fun_k(self, fi):
        if fi < 0 or fi > 2 * pi:
            raise ValidationError("Недопустимое значение угла поворота кулачка")
        elif fi < self.config.phi_1:
            return -k_phi(fi, self.config.C_list_1, self.config.k_list_1, self.config.phi_1, fi_0=0,h_kn_max=self.config.z)
        elif fi < self.config.phi_2:
            return -k_phi(fi, self.config.C_list_2, self.config.k_list_2, self.config.phi_2, fi_0=self.config.phi_1, h_kn_max=self.config.h)
        elif fi < self.config.phi_3:
            return 0
        elif fi < self.config.phi_4:
            return -k_phi(self.config.phi_4 - fi + self.config.phi_3, self.config.C_list_3, self.config.k_list_3, self.config.phi_4, fi_0=self.config.phi_3, h_kn_max=self.config.h)
        elif fi < self.config.phi_5:
            return -k_phi(self.config.phi_5 - fi + self.config.phi_4, self.config.C_list_4, self.config.k_list_4, self.config.phi_5, fi_0=self.config.phi_4, h_kn_max = self.config.z)
        elif fi < pi * 2:
            return 0
'''
if __name__ == '__main__':
    # Исходные данные
    N_k = 1000
    D = 30.0 * 1e-3
    h = 12.0 * 1e-3
    z = 0.25 * 1e-3
    f_pod = 80.0 / 180 * pi
    f_v = 5.0 / 180 * pi
    f_op = 75.0 / 180 * pi
    f_z = 25 / 180 * pi

    m = 3
    d = 2
    k_1 = 30
    k_2 = 30
    k_3 = 30
    k_4 = 30
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
        print("Конфигурация успешно создана!")

    except ValidationError as e:
        print("Ошибка в данных:", e)

    kulachok = Kulachok_polidain(config)
    kulachok.display_graphs()

