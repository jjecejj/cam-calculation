from pydantic import BaseModel, Field, ConfigDict
from typing import Union

from core.cam_geometry.options import CamGeometryOptions
from core.optimization.options import CamOptimizationOptions, calculate_cam_optimization
from core.cam_geometry import Kulachok
from vizualization.plotter.options import PlotterOptions
from vizualization.rotate_animation.logic import display_animation, set_rotate_data, display_dashboard_animation
from vizualization.rotate_animation.options import RotateAnimationOptions
from exporters.dxf_creator.logic import build_profile
from exporters.dxf_creator.options import DxfCreatorOptions


class CamSolveOptions(BaseModel):
    """
    Параметры для решения прямой задачи расчета геометрии кулачка.
    """
    model_config = ConfigDict(arbitrary_types_allowed=True, extra='forbid')

    cam_geometry_options: CamGeometryOptions = Field(default_factory=CamGeometryOptions)
    plotter_options: PlotterOptions = Field(default_factory=PlotterOptions)
    rotate_animation_options: RotateAnimationOptions = Field(default_factory=RotateAnimationOptions)
    dxf_creator_options: DxfCreatorOptions = Field(default_factory=DxfCreatorOptions)

def calculate_cam_solve(cam_solve_options: CamSolveOptions):
    """
    Выполняет расчет геометрии кулачка и запускает выбранные методы вывода (графики, анимация, DXF).
    """
    # 1. Инициализация и расчет
    kulachok = Kulachok(cam_solve_options.cam_config, cam_solve_options.calculator)
    kulachok.solve(kulachok_type=cam_solve_options.kulachok_type, N=cam_solve_options.N)

    # 2. Определение начального угла
    if cam_solve_options.calculate_optimal_initial_angle:
        initial_angle = calculate_optimal_angle(kulachok)
    else:
        initial_angle = cam_solve_options.initial_angle

    # 3. Отрисовка графиков
    if cam_solve_options.profile_and_graphs_together_flag:
        if cam_solve_options.graphs_kulachok_flag:
            display_dashboard(kulachok, initial_angle=initial_angle, graphs_type= cam_solve_options.graphs_argument_type, target='kulachok')

        if cam_solve_options.graphs_tolkatel_flag:
            display_dashboard(kulachok, initial_angle=initial_angle, graphs_type=cam_solve_options.graphs_argument_type,  target='tolkatel')

    else:
        if cam_solve_options.graphs_tolkatel_flag:
            display_graphs_tolkatel(kulachok.kulachok_data, initial_angle=initial_angle, graphs_type= cam_solve_options.graphs_argument_type)

        if cam_solve_options.graphs_kulachok_flag:
            display_graphs_kulachok(kulachok.kulachok_data, initial_angle=initial_angle, graphs_type= cam_solve_options.graphs_argument_type)

        if cam_solve_options.graphs_profile_flag:
            display_profil(kulachok.profile_data, initial_angle=initial_angle)

    # 4. Анимация
    if cam_solve_options.display_animation_flag:
        rotate_data = set_rotate_data(kulachok, tolkatel_type=kulachok.tolkatel_solve_type)
        display_animation(
            rotate_data,
            interval=cam_solve_options.animation_intarval,
            save_flag=cam_solve_options.save_animation_flag,
            name_file=cam_solve_options.profile_animation_name_file,
            pause_flag=cam_solve_options.animation_pause_flag,
        )
    if cam_solve_options.animation_profile_and_graphs_together_flag:
        display_dashboard_animation(kulachok, tolkatel_type=kulachok.tolkatel_solve_type,
                                    interval=cam_solve_options.animation_intarval,
                                    save_flag=cam_solve_options.save_animation_flag,
                                    name_file=cam_solve_options.dashboard_animation_name_file,
                                    graphs_type=cam_solve_options.animation_graphs_argument_type,
                                    pause_flag=cam_solve_options.animation_pause_flag, )

    # 5. Экспорт в DXF
    if cam_solve_options.import_dxf_flag:
        build_profile(
            kulachok.profile_data,
            profile_name=cam_solve_options.dxf_profile_name,
            line_type=cam_solve_options.dxf_line_type
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