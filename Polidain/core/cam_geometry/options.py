from typing import Literal

from pydantic import BaseModel, Field, model_validator

from core.cam_geometry.config import KulachokConfig, default_kulachok_config
from core.profiling_methods.base import BaseCalculator
from core.profiling_methods.base.config import MethodConfig
from core.profiling_methods.polidain import PolidainCalculator
from core.profiling_methods.polidain.config import default_polidain_config, PolidainConfig
from core.profiling_methods.polinmail.config import default_polinmail_config, PolinmailConfig
from core.profiling_methods.polinmail.logic import PolinmailCalculator


class CamGeometryOptions(BaseModel):
    cam_config: KulachokConfig = Field(default=default_kulachok_config, description="Геометрические и кинематические параметры кулачка")
    calculator_config: MethodConfig | None = Field(default=None, description="Параметры профилирования кулачка")
    calculator: BaseCalculator | None = Field(default=None, description="Решатель для профилирования кулачка")
    calculator_type: Literal['polidain', 'polinmail'] = Field(default='polidain', description="Тип метода профилирования кулачка")
    kulachok_type: Literal['thin', 'flat', 'roller'] = Field(default='thin', description="Тип толкателя: остроконечный, плоский или роликовый")
    N: int = Field(default=1000, ge=10, description="Количество точек расчета")
    initial_angle: float = Field(default=0.0, ge=0.0, le=360.0, description="Начальный угол поворота (град)")
    calculate_optimal_initial_angle: bool = Field(default=True, description="Авторассчет оптимального начального угла")

    @model_validator(mode='after')
    def resolve_calculator_config(self):
        if self.calculator_config is None:
            if self.calculator_type == 'polidain':
                self.calculator_config = default_polidain_config
            elif self.calculator_type == 'polinmail':
                self.calculator_config = default_polinmail_config
        return self

    @model_validator(mode='after')
    def resolve_calculator(self):
        if self.calculator is None:
            if type(self.calculator_config) is PolidainConfig:
                self.calculator = PolidainCalculator(self.calculator_config)
            elif type(self.calculator_config) is PolinmailConfig:
                self.calculator = PolinmailCalculator(self.calculator_config)
            else:
                raise ValueError("calculator_config неправильного типа")
        return self