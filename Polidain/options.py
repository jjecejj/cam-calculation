from pydantic import BaseModel, Field, ConfigDict
from typing import Union

from core.cam_geometry.options import CamGeometryOptions, resolve_cam_geometry_options
from core.optimization.options import CamOptimizationOptions
from vizualization.plotter.options import PlotterOptions, resolve_plotter_options
from vizualization.rotate_animation.options import RotateAnimationOptions, resolve_rotate_animation_options
from exporters.dxf_creator.options import DxfCreatorOptions, resolve_dxf_creator_options


class CamSolveOptions(BaseModel):
    """
    Параметры для решения прямой задачи расчета геометрии кулачка.
    Агрегирует настройки геометрии, графиков, анимации и экспорта.
    """
    model_config = ConfigDict(arbitrary_types_allowed=True, extra='forbid')

    cam_geometry_options: CamGeometryOptions = Field(default_factory=CamGeometryOptions, description="Настройки геометрии и метода расчета")
    plotter_options: PlotterOptions = Field(default_factory=PlotterOptions, description="Настройки отображения графиков")
    rotate_animation_options: RotateAnimationOptions = Field(default_factory=RotateAnimationOptions, description="Настройки анимации")
    dxf_creator_options: DxfCreatorOptions = Field(default_factory=DxfCreatorOptions, description="Настройки экспорта в DXF")

def calculate_cam_solve(cam_solve_options: CamSolveOptions):
    """
    Выполняет полный цикл расчета геометрии кулачка.
    1. Рассчитывает профиль (resolve_cam_geometry_options).
    2. Строит графики (resolve_plotter_options).
    3. Запускает анимацию (resolve_rotate_animation_options).
    4. Экспортирует DXF (resolve_dxf_creator_options).

    Args:
        cam_solve_options: Объект с полными настройками расчета.
    """
    kulachok, initial_angle = resolve_cam_geometry_options(cam_solve_options.cam_geometry_options)
    resolve_plotter_options(cam_solve_options.plotter_options, kulachok, initial_angle)
    resolve_rotate_animation_options(cam_solve_options.rotate_animation_options, kulachok)
    resolve_dxf_creator_options(cam_solve_options.dxf_creator_options, kulachok)

def calculate(cam_options: Union[CamSolveOptions, CamOptimizationOptions]):
    """
    Единая точка входа для запуска расчетов.
    Определяет тип задачи (прямой расчет или оптимизация) и вызывает соответствующий обработчик.

    Args:
        cam_options: Объект настроек (CamSolveOptions или CamOptimizationOptions).

    Raises:
        ValueError: Если передан объект неподдерживаемого типа.
    """
    if isinstance(cam_options, CamSolveOptions):
        return calculate_cam_solve(cam_options)
    elif isinstance(cam_options, CamOptimizationOptions):
        return None
        #return calculate_cam_optimization(cam_options)
    else:
        raise ValueError(
            f'cam_options должен быть объектом одного из классов: '
            f'[CamSolveOptions, CamOptimizationOptions]. Получен: {type(cam_options)}'
        )
