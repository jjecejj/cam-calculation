from pydantic import Field, model_validator
from core.profiling_methods.base.config import MethodConfig

class ValidationError(Exception):
    pass

class PolidainConfig(MethodConfig):
    """
    Класс-валидатор входных данных
    Проверяет типы и ограничения
    Передаваемы данные:
        'm':Степень при C2 только целочисленное и не меньше 3\n
        'd':Разность между степенями членов полинома только целочисленное и не меньше 1\n
        'k_1':Коэффициент агрессивности первого участка (выбор зазора)\n
        'k_2':Коэффициент агрессивности второго участка (Фаза подъёма)\n
        'k_3':Коэффициент агрессивности четвёртого участка (Фаза опускания)\n
        'k_4':Коэффициент агрессивности пятого участка (Фаза выбора зазора)
    """
    # Параметры полинома
    m: int = Field(..., ge=2, description="Степень первого члена полинома")
    d: int = Field(..., ge=1, description="k - d: Разность между вторым и последующими членами полинома")

    # Коэффициенты агрессивности
    k_1: int = Field(..., gt=0, description="Степень второго члена полинома первого участка (выбор зазора)")
    k_2: int = Field(..., gt=0, description="Степень второго члена полинома второго участка (Фаза подъёма)")
    k_3: int = Field(..., gt=0, description="Степень второго члена полинома четвёртого участка (Фаза опускания)")
    k_4: int = Field(..., gt=0, description="Степень второго члена полинома пятого участка (Фаза выбора зазора)")

    @model_validator(mode='after')
    def check_k(self):
        """Проверяет, что все коэффициенты k строго больше m"""
        # Создаем словарь для удобной проверки в цикле
        k_values = {
            'k_1': self.k_1,
            'k_2': self.k_2,
            'k_3': self.k_3,
            'k_4': self.k_4
        }

        # Проходим по всем k
        for name, value in k_values.items():
            # Если условие нарушено (k <= d)
            if value <= self.m:
                raise ValidationError(
                    f"Ошибка в параметре {name}: значение ({value}) "
                    f"должно быть строго больше m ({self.m})"
                )
            if value <= self.d:
                raise ValidationError(
                    f"Ошибка в параметре {name}: значение ({value}) "
                    f"должно быть строго больше d ({self.m})"
                )
        return self

default_polidain_config = PolidainConfig(
        m=3,
        d=4,
        k_1=6,
        k_2=6,
        k_3=6,
        k_4=6
)