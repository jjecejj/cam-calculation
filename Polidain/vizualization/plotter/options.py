from pydantic import BaseModel, Field, model_validator, ConfigDict
from typing import Literal

from core.cam_geometry import Kulachok
from vizualization.plotter import display_all, display_graphs_tolkatel, display_graphs_kulachok, display_profile


class PlotterOptions(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True, extra='forbid')

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

def resolve_plotter_options(config: PlotterOptions, kulachok: Kulachok, initial_angle: float):
    if config.profile_and_graphs_together_flag:
        if config.graphs_kulachok_flag:
            display_all(kulachok, initial_angle=initial_angle, graphs_type=config.graphs_argument_type,
                              target='kulachok')

        if config.graphs_tolkatel_flag:
            display_all(kulachok, initial_angle=initial_angle, graphs_type=config.graphs_argument_type,
                              target='tolkatel')

    else:
        if config.graphs_tolkatel_flag:
            display_graphs_tolkatel(kulachok.tolkatel_data, initial_angle=initial_angle,
                                    graphs_type=config.graphs_argument_type)

        if config.graphs_kulachok_flag:
            display_graphs_kulachok(kulachok.kulachok_data, initial_angle=initial_angle,
                                    graphs_type=config.graphs_argument_type)

        if config.graphs_profile_flag:
            display_profile(kulachok.profile_data, initial_angle=initial_angle)
