from typing import List
from pydantic import Field, model_validator, BaseModel

from core.profiling_methods.base.config import MethodConfig


class ValidationError(Exception):
    pass

class LocalPolinmailConfig(BaseModel):
    m: int = Field(..., ge=1, description="Степень при первом члене")
    d: int = Field(..., ge=1, description="Разность степеней, не меньше 1")
    boundary_conditions: List[float] = Field(default=[-1, 0, 0, 0, 0], description="Граничные условия")

    @property
    def m_list(self) -> List[int]:
        '''Список степеней членов полинома'''
        return [self.m + self.d * i for i in range(0, len(self.boundary_conditions))]

    @model_validator(mode='after')
    def check_boundary_conditions(self):
        m_list = self.m_list
        for i in range(len(self.boundary_conditions)):
            if (self.boundary_conditions[i] == 0) and (i in m_list):
                raise ValidationError("Значение m не соответствует граничным условиям")
        return self

class PolinmailConfig(MethodConfig):
    config_1: LocalPolinmailConfig = Field(default_factory=LocalPolinmailConfig)
    config_2: LocalPolinmailConfig | None = None
    config_3: LocalPolinmailConfig | None = None
    config_4: LocalPolinmailConfig | None = None

    @model_validator(mode='after')
    def config_resolve(self):
        if self.config_2 is None:
            self.config_2 = self.config_1
        if self.config_3 is None:
            self.config_3 = self.config_1
        if self.config_4 is None:
            self.config_4 = self.config_1
        return self

default_local_polinmail_config = LocalPolinmailConfig(
    m = 5,
    d = 1,
    boundary_conditions = [1, 0, 0, 0, 0]
)

default_polinmail_config = PolinmailConfig(
    config_1 = default_local_polinmail_config,
    config_2 = default_local_polinmail_config,
    config_3 = default_local_polinmail_config,
    config_4 = default_local_polinmail_config,
)