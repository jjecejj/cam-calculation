from math import sin, cos, pi
import numpy as np
from pydantic import BaseModel, Field, model_validator, ConfigDict


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

class PolidainConfig(BaseModel):
    """
    Класс-валидатор входных данных
    Проверяет типы и ограничения
    Передаваемы данные:
        'N_k':Количество оборотов кулачка в минут\n
        'D':Базовый диаметр кулака (м)\n
        'D_t:Диаметр толкателя (м)\n'
        'D_к:Диаметр ролика (м)\n'
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

    # Опциональные параметры
    D_t: float = Field(default=0.0, gt=0,  description="Диаметр толкателя")
    R_r: float = Field(default=0.0, gt=0, description="Радиус ролика")

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


class PolidainData(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    fi_list: np.ndarray[float | np.ndarray] | list[float | np.ndarray]
    H: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10)
    V: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10)
    A: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10)
    D: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10)
    K: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10)

class ProfilData(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    fi_list: np.ndarray[float | np.ndarray] | list[float | np.ndarray]
    X: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10)
    Y: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10)