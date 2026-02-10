import numpy as np
from jedi.inference.gradual.typing import Callable
from pydantic import BaseModel, Field, model_validator, ConfigDict

def set_graph_data(fun_list: list[Callable], omega: float, N: int = 1000):
    fi_list = np.linspace(0, 2 * np.pi, N)
    H = fun_list[0](fi_list) * 1000
    V = fun_list[1](fi_list) * 1000
    A = fun_list[2](fi_list) * 1000
    D = fun_list[3](fi_list) * 1000
    K = fun_list[4](fi_list) * 1000
    return GraphData(H_rad = H, V_rad = V, A_rad = A, D_rad = D, K_rad = K, fi_list_rad = fi_list, omega_rad = omega)

def set_profile_data(fun_list: list[Callable], N: int = 1000):
    fi_list = np.linspace(0, 2 * np.pi, N)
    X = fun_list[0](fi_list) * 1000
    Y = fun_list[1](fi_list) * 1000
    fi_list = fi_list / np.pi * 180
    return ProfileData(X = X, Y = Y, fi_list = fi_list)

class GraphData(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    fi_list_rad: np.ndarray[float] | list[float]
    H_rad: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10)
    V_rad: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10)
    A_rad: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10)
    D_rad: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10)
    K_rad: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10)
    omega_rad: float  = Field(ge = 0.0)

    @model_validator(mode='after')
    def check_lists(self):
        self.fi_list_rad = np.asarray(self.fi_list_rad)
        self.H_rad = np.asarray(self.H_rad)
        self.V_rad = np.asarray(self.V_rad)
        self.A_rad = np.asarray(self.A_rad)
        self.D_rad = np.asarray(self.D_rad)
        self.K_rad = np.asarray(self.K_rad)
        return self

    @property
    def fi_list_degree(self) -> np.ndarray[float]:
        return  np.degrees(self.fi_list_rad)
    @property
    def omega_degree(self) -> float:
        return self.omega_rad / np.pi * 180

    @property
    def H_degree(self) -> np.ndarray[float]:
        return self.H_rad

    @property
    def V_degree(self) -> np.ndarray[float]:
        return self.V_rad * (np.pi / 180)

    @property
    def A_degree(self) -> np.ndarray[float]:
        return self.A_rad * (np.pi / 180) ** 2

    @property
    def D_degree(self) -> np.ndarray[float]:
        return self.D_rad * (np.pi / 180) ** 3

    @property
    def K_degree(self) -> np.ndarray[float]:
        return self.K_rad * (np.pi / 180) ** 4

    @property
    def t_list(self) -> np.ndarray[float]:
        return self.fi_list_rad / self.omega_rad

    @property
    def H_t(self) -> np.ndarray[float]:
        return self.H_rad

    @property
    def V_t(self) -> np.ndarray[float]:
        return self.V_rad * self.omega_rad

    @property
    def A_t(self) -> np.ndarray[float]:
        return self.A_rad * self.omega_rad**2

    @property
    def D_t(self) -> np.ndarray[float]:
        return self.D_rad * self.omega_rad**3

    @property
    def K_t(self) -> np.ndarray[float]:
        return self.K_rad * self.omega_rad**4

class ProfileData(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    fi_list: np.ndarray[float | np.ndarray] | list[float | np.ndarray]
    X: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10)
    Y: np.ndarray[float | np.ndarray] | list[float | np.ndarray] = Field(..., min_length=10)