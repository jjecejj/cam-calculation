from pydantic import Field
from core.profiling_methods.base.base_config import MethodConfig

class PolinmailConfig(MethodConfig):
    """
    Класс-валидатор входных данных для полинома 4-6 членов
    Передаваемые данные:
        'n': Количество членов полинома (4, 5 или 6)
        'k_1': Коэффициент агрессивности (для совместимости, не используется)
        'k_2': ...
        'k_3': ...
        'k_4': ...
    """
    n: int = Field(..., ge=4, le=6, description="Количество членов полинома (4, 5 или 6)")

    k_1: int = Field(6, gt=0, description="Коэффициент агрессивности первого участка")
    k_2: int = Field(6, gt=0, description="Коэффициент агрессивности второго участка")
    k_3: int = Field(6, gt=0, description="Коэффициент агрессивности третьего участка")
    k_4: int = Field(6, gt=0, description="Коэффициент агрессивности четвертого участка")

default_polinmail_config = PolinmailConfig(
    n=6,
    k_1=6,
    k_2=6,
    k_3=6,
    k_4=6
)