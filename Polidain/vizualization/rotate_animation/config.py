import numpy as np
from pydantic import BaseModel, ConfigDict, model_validator


class RotateProfileData(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    fi_list: np.ndarray[float | np.ndarray] | list[float | np.ndarray]
    movement_data: list[tuple[np.ndarray | list, np.ndarray | list]]
    tolkatel_data: np.ndarray[float | np.ndarray] | list[float | np.ndarray]
    tolkatel_type: str
    tolkatel_D_t: float | None = None
    tolkatel_R_r: float | None = None

    @model_validator(mode='after')
    def tolkatel_check(self):
        if self.tolkatel_type == "flat":
            if self.tolkatel_D_t is None or self.tolkatel_D_t <= 0:
                raise ValueError("tolkatel_D_t must be greater than zero")
        elif self.tolkatel_type == "roller":
            if self.tolkatel_R_r is None or self.tolkatel_R_r <= 0:
                raise ValueError("tolkatel_R_r must be greater than zero")
        elif self.tolkatel_type == "thin":
            pass
        else:
            raise ValueError("tolkatel_type must be either flat or roller or thin")
        return self

    @model_validator(mode='after')
    def data_check(self):
        if len(self.movement_data) != len(self.fi_list) or len(self.tolkatel_data) != len(self.fi_list):
            raise ValueError("Массивы с данными о перемещении элементов должны иметь одинаковую размерность!")
        return self