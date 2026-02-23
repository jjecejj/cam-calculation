from typing import Literal
from pydantic import BaseModel, model_validator
import ezdxf

from core.cam_geometry.schemas import ProfileData


class ProfileExportData(BaseModel):
    """
    Модель данных для конфигурации и валидации профиля кулачка перед экспортом.
    Используется для хранения координат, проверки замкнутости контура
    и подготовки данных в формате, удобном для библиотеки ezdxf.
    """
    X: list
    Y: list

    @model_validator(mode='after')
    def check_XY(self):
        """
        Валидатор замкнутости контура.
        Проверяет, совпадает ли первая точка с последней. Если контур не замкнут,
        принудительно изменяет координаты ПЕРВОЙ точки, делая их равными
        координатам ПОСЛЕДНЕЙ точки.
        """
        if self.X[0] != self.X[-1] or self.Y[0] != self.Y[-1]:
            self.X[0] = self.X[-1]
            self.Y[0] = self.Y[-1]
        return self

    @property
    def N(self):
        """Количество точек в профиле."""
        return len(self.X)

    @property
    def fit_points(self):
        """Список кортежей (x, y), подготовленный для построения сплайна."""
        return [(self.X[i], self.Y[i]) for i in range(0, self.N)]


def build_profile(profile_data: ProfileData, profile_name: str = "kulachok", line_type: str = "spline"):
    """
    Фасадная функция (обертка) для запуска процесса генерации DXF.
    Извлекает сырые данные из объекта профиля, упаковывает их в валидируемую
    модель ProfileExportData и вызывает функцию создания файла.

    Args:
        profile_data: Объект ProfileData, содержащий атрибуты X и Y (массивы координат).
        profile_name: Имя выходного файла (без расширения).
        line_type: Тип линии ("spline" или "line").
    """
    config = ProfileExportData(X=profile_data.X, Y=profile_data.Y)
    create_profile(config, profil_name=profile_name, line_type=line_type)


def create_profile(profile_export_data: ProfileExportData, profil_name: str = "kulachok", line_type: Literal['line', "spline"] = "spline"):
    """
    Генерирует и сохраняет DXF-файл с профилем кулачка.
    Использует библиотеку ezdxf для создания чертежа версии R2018.
    Сохраняет файл по пути: 'data/output/DXF_profils/<profil_name>.dxf'.

    Args:
        profile_export_data: Валидированный объект конфигурации с координатами.
        profil_name: Имя файла.
        line_type: Метод отрисовки профиля.
            "line" - соединяет точки прямыми отрезками.
            "spline" - строит сглаженный сплайн по точкам (fit points).
    """
    # Создаем новый DXF-документ.
    doc = ezdxf.new(dxfversion='R2018', setup=True)

    # Получаем доступ к пространству модели (Modelspace)
    msp = doc.modelspace()

    if line_type == "line":
        for i in range(1, profile_export_data.N):
            msp.add_line(
                start=(profile_export_data.X[i - 1], profile_export_data.Y[i - 1]),
                end=(profile_export_data.X[i], profile_export_data.Y[i])
            )
    elif line_type == "spline":
        msp.add_spline(fit_points=profile_export_data.fit_points)

    # Сохранение файла
    doc.saveas('data\\output\\DXF_profils\\' + profil_name + ".dxf")
