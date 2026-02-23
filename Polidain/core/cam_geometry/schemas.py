import numpy as np
from typing import Callable
from pydantic import BaseModel, Field, model_validator, ConfigDict

def set_graph_data(fun_list: list[Callable], omega: float, N: int = 1000):
    """
    Рассчитывает данные графиков кинематических параметров.

    Args:
        fun_list: Список функций [H, V, A, D, K].
        omega: Угловая скорость (рад/с).
        N: Количество точек расчета.

    Returns:
        GraphData: Объект с рассчитанными данными (в радианах и метрах/с).
    """
    fi_list = np.linspace(0, 2 * np.pi, N)
    H = fun_list[0](fi_list) * 1000
    V = fun_list[1](fi_list) * 1000
    A = fun_list[2](fi_list) * 1000
    D = fun_list[3](fi_list) * 1000
    K = fun_list[4](fi_list) * 1000
    return GraphData(H_rad = H, V_rad = V, A_rad = A, D_rad = D, K_rad = K, fi_list_rad = fi_list, omega_rad = omega)

def set_profile_data(fun_list: list[Callable], N: int = 1000):
    """
    Рассчитывает координаты профиля кулачка.

    Args:
        fun_list: Список функций [X, Y].
        N: Количество точек расчета.

    Returns:
        ProfileData: Объект с координатами профиля.
    """
    fi_list = np.linspace(0, 2 * np.pi, N)
    X = fun_list[0](fi_list) * 1000
    Y = fun_list[1](fi_list) * 1000
    fi_list = fi_list / np.pi * 180
    return ProfileData(X = X, Y = Y, fi_list = fi_list)

class GraphData(BaseModel):
    """
    Модель данных для хранения кинематических графиков (S, V, A, J и др.).
    Хранит базовые значения в радианах и предоставляет свойства для получения
    значений в градусах или во времени.
    """
    model_config = ConfigDict(arbitrary_types_allowed=True)

    fi_list_rad: np.ndarray[float] | list[float] = Field(..., description="Массив углов поворота (рад)")
    H_rad: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10, description="Перемещение (аналог S)")
    V_rad: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10, description="Аналог скорости")
    A_rad: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10, description="Аналог ускорения")
    D_rad: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10, description="Аналог рывка")
    K_rad: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10, description="Четвертая производная")
    omega_rad: float  = Field(ge = 0.0, description="Угловая скорость (рад/с)")

    @model_validator(mode='after')
    def check_lists(self):
        """Преобразует списки в numpy массивы после инициализации."""
        self.fi_list_rad = np.asarray(self.fi_list_rad)
        self.H_rad = np.asarray(self.H_rad)
        self.V_rad = np.asarray(self.V_rad)
        self.A_rad = np.asarray(self.A_rad)
        self.D_rad = np.asarray(self.D_rad)
        self.K_rad = np.asarray(self.K_rad)
        return self

    @property
    def fi_list_degree(self) -> np.ndarray[float]:
        """Массив углов в градусах."""
        return  np.degrees(self.fi_list_rad)
    @property
    def omega_degree(self) -> float:
        """Угловая скорость в градусах/с."""
        return self.omega_rad / np.pi * 180

    @property
    def H_degree(self) -> np.ndarray[float]:
        """Перемещение (аргумент - градусы)."""
        return self.H_rad

    @property
    def V_degree(self) -> np.ndarray[float]:
        """Скорость (аргумент - градусы)."""
        return self.V_rad * (np.pi / 180)

    @property
    def A_degree(self) -> np.ndarray[float]:
        """Ускорение (аргумент - градусы)."""
        return self.A_rad * (np.pi / 180) ** 2

    @property
    def D_degree(self) -> np.ndarray[float]:
        """Рывок (аргумент - градусы)."""
        return self.D_rad * (np.pi / 180) ** 3

    @property
    def K_degree(self) -> np.ndarray[float]:
        """Четвертая производная (аргумент - градусы)."""
        return self.K_rad * (np.pi / 180) ** 4

    @property
    def t_list(self) -> np.ndarray[float]:
        """Массив времени (с)."""
        return self.fi_list_rad / self.omega_rad

    @property
    def H_t(self) -> np.ndarray[float]:
        """Перемещение (м)."""
        return self.H_rad

    @property
    def V_t(self) -> np.ndarray[float]:
        """Скорость (м/с)."""
        return self.V_rad * self.omega_rad

    @property
    def A_t(self) -> np.ndarray[float]:
        """Ускорение (м/с^2)."""
        return self.A_rad * self.omega_rad**2

    @property
    def D_t(self) -> np.ndarray[float]:
        """Рывок (м/с^3)."""
        return self.D_rad * self.omega_rad**3

    @property
    def K_t(self) -> np.ndarray[float]:
        """Четвертая производная (м/с^4)."""
        return self.K_rad * self.omega_rad**4

class ProfileData(BaseModel):
    """
    Модель данных для хранения координат профиля кулачка.
    """
    model_config = ConfigDict(arbitrary_types_allowed=True)
    fi_list: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., description="Массив углов")
    X: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10, description="Координаты X профиля")
    Y: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10, description="Координаты Y профиля")
