from pydantic import BaseModel, Field
from typing import Literal


class RotateAnimationOptions(BaseModel):
    # Настройки анимации
    display_animation_flag: bool = Field(default=False, description="Запустить анимацию работы механизма")
    animation_profile_and_graphs_together_flag: bool = Field(default=False, description="Запустить анимацию работы механизма c графиками")
    save_animation_flag: bool = Field(default=False, description="Сохранить анимацию в файл")
    animation_intarval: int = Field(default=50, description="Интервал между кадрами (мс)")
    profile_animation_name_file: str = Field(default="animation_profile", description="Имя файла для сохранения анимации профиля")
    dashboard_animation_name_file: str = Field(default="animation_dashboard", description="Имя файла для сохранения анимации работы механизма c графиками")
    animation_graphs_argument_type: Literal['degree', 'rad', 't'] = Field(default='degree', description="Аргумент графиков анимации")
    animation_pause_flag: bool = Field(default=False, description="Поддержка паузы во время анимации")
