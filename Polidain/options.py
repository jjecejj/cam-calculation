from pydantic import BaseModel, Field, ConfigDict
from typing import Union

from core.cam_geometry.options import CamGeometryOptions, resolve_cam_geometry_options
from core.optimization.options import CamOptimizationOptions, calculate_cam_optimization
from vizualization.plotter.options import PlotterOptions, resolve_plotter_options
from vizualization.rotate_animation.options import RotateAnimationOptions, resolve_rotate_animation_options
from exporters.dxf_creator.options import DxfCreatorOptions, resolve_dxf_creator_options


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
    kulachok, initial_angle = resolve_cam_geometry_options(cam_solve_options.cam_geometry_options)
    resolve_plotter_options(cam_solve_options.plotter_options, kulachok, initial_angle)
    resolve_rotate_animation_options(cam_solve_options.rotate_animation_options, kulachok)
    resolve_dxf_creator_options(cam_solve_options.dxf_creator_options, kulachok)

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