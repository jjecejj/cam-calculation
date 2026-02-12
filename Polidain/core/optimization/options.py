from typing import Callable
from pydantic import BaseModel, ConfigDict, Field

from core.optimization import GibridOptimizationConfig, OptimizeConfig, BoundsConfig, DifferentialEvolutionConfig, gibrid_optimization


class CamOptimizationOptions(BaseModel):
    """
    Параметры для решения задачи оптимизации параметров кулачка.
    """
    model_config = ConfigDict(arbitrary_types_allowed=True, extra='forbid')

    optimiz_config: OptimizeConfig
    gibrid_optimization_config: GibridOptimizationConfig
    bounds_config: BoundsConfig = Field(default_factory=BoundsConfig, description="Границы переменных оптимизации")
    differential_evolution_config: DifferentialEvolutionConfig = Field(default_factory=DifferentialEvolutionConfig, description="Настройки алгоритма дифференциальной эволюции")
    R_func: Callable = Field(description="Функция расчета радиус-вектора")
    display_compression_profile_flag: bool = Field(default=True, description="Сравнить полученный профиль с эталонным")
    N: int = Field(default=1000, ge=10, description="Разрешение профиля при оптимизации")

def calculate_cam_optimization(cam_optimization_options: CamOptimizationOptions):
    """
    Запускает процесс гибридной оптимизации кулачкового механизма.
    """
    gibrid_optimization(
        cam_optimization_options.gibrid_optimization_config,
        cam_optimization_options.optimiz_config,
        cam_optimization_options.R_func,
        differential_evolution_config=cam_optimization_options.differential_evolution_config,
        bounds_config=cam_optimization_options.bounds_config,
        N=cam_optimization_options.N
    )