from pydantic import BaseModel, model_validator
import ezdxf


class ProfilConfig(BaseModel):
    """
    Модель данных для конфигурации и валидации профиля кулачка перед экспортом.

    Используется для хранения координат, проверки замкнутости контура
    и подготовки данных в формате, удобном для библиотеки ezdxf.

    Attributes:
        X (list): Список координат X точек профиля (в мм).
        Y (list): Список координат Y точек профиля (в мм).
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

        Returns:
            self: Экземпляр модели с гарантированно замкнутым контуром.
        """
        if self.X[0] != self.X[-1] or self.Y[0] != self.Y[-1]:
            self.X[0] = self.X[-1]
            self.Y[0] = self.Y[-1]
        return self

    @property
    def N(self):
        """
        int: Количество точек в профиле.
        """
        return len(self.X)

    @property
    def fit_points(self):
        """
        list[tuple]: Список кортежей (x, y), подготовленный для построения сплайна.
        """
        return [(self.X[i], self.Y[i]) for i in range(0, self.N)]


def build_profile(profil_data, profil_name: str = "kulachok", line_type: str = "spline"):
    """
    Фасадная функция (обертка) для запуска процесса генерации DXF.

    Извлекает сырые данные из объекта профиля, упаковывает их в валидируемую
    конфигурацию ProfilConfig и вызывает функцию создания файла.

    Args:
        profil_data: Объект, содержащий атрибуты X и Y (массивы координат).
        profil_name (str, optional): Имя выходного файла (без расширения).
            По умолчанию "kulachok".
        line_type (str, optional): Тип линии ("spline" или "line").
            По умолчанию "spline".
    """
    config = ProfilConfig(X=profil_data.X, Y=profil_data.Y)
    create_profil(config, profil_name=profil_name, line_type=line_type)


def create_profil(config: ProfilConfig, profil_name: str = "kulachok", line_type: str = "spline"):
    """
    Генерирует и сохраняет DXF-файл с профилем кулачка.

    Использует библиотеку ezdxf для создания чертежа версии R2018.
    Сохраняет файл по пути: 'data/output/DXF_profils/<profil_name>.dxf'.

    Args:
        config (ProfilConfig): Валидированный объект конфигурации с координатами.
        profil_name (str, optional): Имя файла. По умолчанию "kulachok".
        line_type (str, optional): Метод отрисовки профиля.
            - "line": Соединяет точки прямыми отрезками.
              ВНИМАНИЕ: В этом режиме координаты делятся на 1000 (перевод мм -> м).
            - "spline": Строит сглаженный сплайн по точкам (fit points).
              В этом режиме координаты используются как есть (без масштабирования).
    """
    # Создаем новый DXF-документ.
    doc = ezdxf.new(dxfversion='R2018', setup=True)

    # Получаем доступ к пространству модели (Modelspace)
    msp = doc.modelspace()

    if line_type == "line":
        for i in range(1, config.N):
            # print(config.X[i - 1], config.Y[i - 1]) # Отладка
            msp.add_line(
                start=(config.X[i - 1], config.Y[i - 1]),
                end=(config.X[i], config.Y[i])
            )
    elif line_type == "spline":
        msp.add_spline(fit_points=config.fit_points)

    # Сохранение файла
    # Лучше использовать os.path.join, но оставляем как в оригинале для совместимости
    doc.saveas('data\\output\\DXF_profils\\' + profil_name + ".dxf")