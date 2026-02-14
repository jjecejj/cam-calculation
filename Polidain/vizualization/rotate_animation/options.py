from pydantic import BaseModel, Field, ConfigDict
from typing import Literal

from core.cam_geometry import Kulachok
from vizualization.rotate_animation.logic import set_rotate_data, display_animation, display_dashboard_animation


class RotateAnimationOptions(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True, extra='forbid')

    # Настройки анимации
    display_animation_flag: bool = Field(default=False, description="Запустить анимацию работы механизма")
    animation_profile_and_graphs_together_flag: bool = Field(default=False, description="Запустить анимацию работы механизма c графиками")
    save_animation_flag: bool = Field(default=False, description="Сохранить анимацию в файл")
    animation_intarval: int = Field(default=50, description="Интервал между кадрами (мс)")
    profile_animation_name_file: str = Field(default="animation_profile", description="Имя файла для сохранения анимации профиля")
    dashboard_animation_name_file: str = Field(default="animation_dashboard", description="Имя файла для сохранения анимации работы механизма c графиками")
    animation_graphs_argument_type: Literal['degree', 'rad', 't'] = Field(default='degree', description="Аргумент графиков анимации")
    animation_pause_flag: bool = Field(default=False, description="Поддержка паузы во время анимации")

def resolve_rotate_animation_options(config: RotateAnimationOptions, kulachok: Kulachok):
    if config.display_animation_flag:
        display_animation(
            kulachok,
            interval=config.animation_intarval,
            save_flag=config.save_animation_flag,
            name_file=config.profile_animation_name_file,
            pause_flag=config.animation_pause_flag,
        )
    if config.animation_profile_and_graphs_together_flag:
        display_dashboard_animation(kulachok,
                                    interval=config.animation_intarval,
                                    save_flag=config.save_animation_flag,
                                    name_file=config.dashboard_animation_name_file,
                                    graphs_type=config.animation_graphs_argument_type,
                                    pause_flag=config.animation_pause_flag, )
    return None