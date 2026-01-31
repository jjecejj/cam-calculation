from pydantic import BaseModel, Field
from typing import Literal, Callable, Union
from core.cam_geometry import Kulachok_polidain
from core.schemas import PolidainConfig
from core.optimization import (
    OptimizeConfig, BoundsConfig, DifferentialEvolutionConfig,
    GibridOptimizationConfig, gibrid_optimization
)
from vizualization.plotter import (
    display_graphs_kulachok, display_graphs_tolkatel,
    display_profil, calculate_optimal_angle
)
from vizualization.rotate_animation import display_animation, set_rotate_data
from exporters.dxf_creator import build_profile
from core.pyzirev_profil import R_func as R_func_pyzirev


class CamSolveOptions(BaseModel):
    """
    Параметры для решения прямой задачи расчета геометрии кулачка.
    """
    cam: PolidainConfig
    kulachok_type: Literal['thin', 'flat', 'roller'] = Field(
        default='thin', description="Тип толкателя: остроконечный, плоский или роликовый"
    )
    N: int = Field(default=1000, ge=10, description="Количество точек расчета")
    initial_angle: float = Field(default=0.0, ge=0.0, le=360.0, description="Начальный угол поворота (град)")
    calculate_optimal_initial_angle: bool = Field(default=True, description="Авторасчет оптимального начального угла")

    # Флаги визуализации
    graphs_tolkatel_flag: bool = Field(default=False, description="Показать графики движения толкателя")
    graphs_kulachok_flag: bool = Field(default=False, description="Показать графики характеристик кулачка")
    graphs_profil_flag: bool = Field(default=False, description="Показать профиль кулачка")

    # Настройки анимации
    display_animation_flag: bool = Field(default=False, description="Запустить анимацию работы механизма")
    save_animation_flag: bool = Field(default=False, description="Сохранить анимацию в файл")
    animation_intarval: int = Field(default=50, description="Интервал между кадрами (мс)")
    animation_name_file: str = Field(default="animation", description="Имя файла для сохранения анимации")

    # Экспорт
    import_dxf_flag: bool = Field(default=False, description="Экспортировать профиль в формат DXF")
    dxf_profil_name: str = Field(default="kulachok_1", description="Имя файла DXF")
    dxf_line_type: Literal["spline", "line"] = Field(default="spline", description="Тип геометрии в DXF")


class CamOptimizationOptions(BaseModel):
    """
    Параметры для решения задачи оптимизации параметров кулачка.
    """
    optimiz_config: OptimizeConfig
    gibrid_optimization_config: GibridOptimizationConfig
    bounds_config: BoundsConfig = Field(default_factory=BoundsConfig, description="Границы переменных оптимизации")
    differential_evolution_config: DifferentialEvolutionConfig = Field(
        default_factory=DifferentialEvolutionConfig, description="Настройки алгоритма дифференциальной эволюции"
    )
    R_func: Callable = Field(default_factory=lambda: R_func_pyzirev, description="Функция расчета радиус-вектора")
    display_comprasion_profil_flag: bool = Field(default=True, description="Сравнить полученный профиль с эталонным")
    N: int = Field(default=1000, ge=10, description="Разрешение профиля при оптимизации")


def calculate_cam_solve(cam_solve_options: CamSolveOptions):
    """
    Выполняет расчет геометрии кулачка и запускает выбранные методы вывода (графики, анимация, DXF).
    """
    # 1. Инициализация и расчет
    kulachok = Kulachok_polidain(cam_solve_options.cam)
    kulachok.solve(kulachok_type=cam_solve_options.kulachok_type, N=cam_solve_options.N)

    # 2. Определение начального угла
    if cam_solve_options.calculate_optimal_initial_angle:
        cam_solve_options.initial_angle = calculate_optimal_angle(kulachok)

    # 3. Отрисовка графиков
    if cam_solve_options.graphs_tolkatel_flag:
        display_graphs_tolkatel(kulachok.kulachok_data, initial_angle=cam_solve_options.initial_angle)

    if cam_solve_options.graphs_kulachok_flag:
        display_graphs_kulachok(kulachok.kulachok_data, initial_angle=cam_solve_options.initial_angle)

    if cam_solve_options.graphs_profil_flag:
        display_profil(kulachok.profil_data, initial_angle=cam_solve_options.initial_angle)

    # 4. Анимация
    if cam_solve_options.display_animation_flag:
        rotate_data = set_rotate_data(kulachok, tolkatel_type=kulachok.solve_type)
        display_animation(
            rotate_data,
            interval=cam_solve_options.animation_intarval,
            save_flag=cam_solve_options.save_animation_flag,
            name_file=cam_solve_options.animation_name_file
        )

    # 5. Экспорт в DXF
    if cam_solve_options.import_dxf_flag:
        build_profile(
            kulachok.profil_data,
            profil_name=cam_solve_options.dxf_profil_name,
            line_type=cam_solve_options.dxf_line_type
        )


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


def calculate(cam_options: Union[CamSolveOptions, CamOptimizationOptions]):
    """
    Единая точка входа для запуска расчетов.
    Определяет тип задачи (расчет или оптимизация) и вызывает соответствующий обработчик.

    Args:
        cam_options: Объект настроек (CamSolveOptions или CamOptimizationOptions).

    Raises:
        ValueError: Если передан объект неподдерживаемого типа.
    """
    if isinstance(cam_options, CamSolveOptions):
        return calculate_cam_solve(cam_options)
    elif isinstance(cam_options, CamOptimizationOptions):
        return calculate_cam_optimization(cam_options)
    else:
        raise ValueError(
            f'cam_options должен быть объектом одного из классов: '
            f'[CamSolveOptions, CamOptimizationOptions]. Получен: {type(cam_options)}'
        )