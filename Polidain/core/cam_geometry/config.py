import numpy as np
from pydantic import BaseModel, Field, model_validator

class ValidationError(Exception):
    pass

class KulachokConfig(BaseModel):
    """
    Класс-валидатор входных данных
    Проверяет типы и ограничения
    Передаваемы данные:
        'N_k':Количество оборотов кулачка в минут\n
        'D':Базовый диаметр кулака (м)\n
        'D_t:Диаметр толкателя (м)\n'
        'D_r:Диаметр ролика (м)\n'
        'h':Максимальное перемещение толкателя (м)\n
        'z':Тепловой зазор (м)\n
        'f_pod':Фаза подъёма (рад)\n
        'f_v':Фаза выдержки (рад)\n
        'f_op':Фаза опускания (рад)\n
        'f_z':Фаза теплового зазора (рад)
    """
    # Общие параметры
    N_k: float = Field(..., gt=0, description="Обороты в минуту")
    D: float = Field(..., gt=0, description="Базовый диаметр")
    h: float = Field(..., gt=0, description="Максимальное перемещение")
    z: float = Field(..., ge=0, description="Тепловой зазор")

    # Фазовые углы
    f_pod: float = Field(..., gt=0, description="Фаза подъёма (рад)")
    f_v: float = Field(..., gt=0, description=" Фаза выдержки (рад)")
    f_op: float = Field(..., gt=0, description="Фаза опускания (рад)")
    f_z: float = Field(..., gt=0, description="Фаза теплового зазора (рад)")

    # Опциональные параметры
    D_t: float = Field(default=0.0, ge=0,  description="Диаметр толкателя")
    R_r: float = Field(default=0.0, ge=0, description="Радиус ролика")

    @model_validator(mode='after')
    def check_f(self):
        if self.f_pod + self.f_v + self.f_op + self.f_z > 2 * np.pi:
            raise ValidationError(f"Ошибка в значениях фаз, их сумма больше 2 * pi")
        return self

    @model_validator(mode='after')
    def check_z(self):
        if self.z >= self.D / 2:
            raise ValidationError(f"Ошибка в значении теплового зазора z, его значение больше радиуса кулачка R")
        return self
    # --- Вычисляемые свойства (Автоматический расчет) ---
    @property
    def phi_0(self) -> float:
        return 0

    @property
    def phi_1(self) -> float:
        return self.phi_0 + self.f_z

    @property
    def phi_2(self) -> float:
        return self.phi_1 + self.f_pod

    @property
    def phi_3(self) -> float:
        return self.phi_2 + self.f_v

    @property
    def phi_4(self) -> float:
        return self.phi_3 + self.f_op

    @property
    def phi_5(self) -> float:
        return self.phi_4 + self.f_z

    @property
    def omega(self) -> float:
        """Угловая скорость (рад/с)"""
        return self.N_k * 2 * np.pi / 60

    @property
    def T(self) -> float:
        """Период оборота (с)"""
        return 60 / self.N_k

    @property
    def r0(self) -> float:
        """Базовый радиус (м)"""
        return self.D / 2

default_kulachok_config = KulachokConfig(
        N_k=1000,
        D=30.0 * 1e-3,
        D_t=30.0 * 1e-3,
        h=2.0 * 1e-3,
        z=0.25 * 1e-3,
        f_pod=80.0 / 180 * np.pi,
        f_v=5.0 / 180 * np.pi,
        f_op=75.0 / 180 * np.pi,
        f_z=25 / 180 * np.pi,
        R_r=5 * 1e-3,
)