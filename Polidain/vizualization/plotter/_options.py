from pydantic import BaseModel, Field, model_validator
from typing import Literal

from .logic import (
    display_graphs_kulachok, display_graphs_tolkatel,
    display_profil, calculate_optimal_angle, display_dashboard
)

class PlotterOptions(BaseModel):
    graphs_tolkatel_flag: bool = Field(default=False, description="Показать графики движения толкателя")
    graphs_kulachok_flag: bool = Field(default=False, description="Показать графики характеристик кулачка")
    graphs_argument_type: Literal['degree', 'rad', 't'] = Field(default='degree', description="Аргумент графиков")
    graphs_profile_flag: bool = Field(default=False, description="Показать профиль кулачка")
    profile_and_graphs_together_flag: bool = Field(default=False, description="Показать профиль кулачка и графики в одном окне")

    @model_validator(mode='after')
    def resolve_conflicting_flags(self):
        """
        Автоматически разрешает конфликты флагов визуализации.
        Если выбран режим 'together', он имеет приоритет.
        """
        if self.profile_and_graphs_together_flag:
            self.graphs_profile_flag = True

            if not self.graphs_tolkatel_flag and not self.graphs_kulachok_flag:
                self.graphs_kulachok_flag = True

        return self
