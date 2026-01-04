import os
import sys
import sys
import ezdxf
import numpy as np
from pydantic import BaseModel, Field, model_validator
from math import pi

current_file_path = os.path.abspath(__file__)
current_dir = os.path.dirname(current_file_path)
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
from main import Kulachok_polidain, PolidainConfig

class ValidationError(Exception):
    pass

class ProfilConfig(BaseModel):
    X:list
    Y:list

    @model_validator(mode='after')
    def check_XY(self):
        if self.X[0] != self.X[-1] or self.Y[0] != self.Y[-1]:
            raise ValidationError(f"Ошибка в массивах X и Y:"
                                  f"Массивы должны быть замкнуты!")
        return self

    @property
    def N(self):
        return len(self.X)

    @property
    def fit_points(self):
        return [(X[i], Y[i]) for i in range(0, self.N)]

def create_profil(config: ProfilConfig, profil_name = "kulachok", line_type="spline"):
    # 1. Создаем новый DXF-документ.
    # setup=True настраивает стандартные стили линий и текста.
    doc = ezdxf.new(dxfversion='R2018', setup=True)

    # 2. Получаем доступ к пространству модели (Modelspace)
    msp = doc.modelspace()

    if line_type == "line":
        for i in range(1, config.N):
            print(config.X[i - 1], config.Y[i - 1])
            msp.add_line(start=(config.X[i - 1] / 1000, config.Y[i - 1] / 1000), end=(config.X[i] / 1000, config.Y[i] / 1000))
    elif line_type == "spline":
        msp.add_spline(fit_points=config.fit_points)
    doc.saveas('profils\\' + profil_name + ".dxf")

if __name__ == '__main__':
    # Исходные данные
    N_k = 1000  # Количество оборотов кулачка в минут
    z = 0.27281 * 1e-3  # Тепловой зазор (мм)
    f_pod = 1.62663  # Фаза подъёма (град)
    f_v = 0.00001  # Фаза выдержки (град)
    f_op = 1.62631  # Фаза опускания (град)
    f_z = 0.37264  # Фаза теплового зазора (град)

    D = 2 * z + 17.8009 * 2  # Базовый диаметр кулака (мм)
    h = 30.8018837267781 - D / 2  # Максимальное перемещение толкателя

    m = 3  # степень при C2 только целочисленное и не меньше 3!!!
    d = 2  # разность между степенями членов полинома только целочисленное и не меньше 1!!!
    k_1 = 982  # коэффициент агрессивности первого участка (выбор зазора)
    k_2 = 18  # коэффициент агрессивности второго участка (Фаза подъёма)
    k_3 = 18  # коэффициент агрессивности четвёртого участка (Фаза опускания)
    k_4 = 980  # коэффициент агрессивности пятого участка (Фаза выбора зазора)

    # Настройки генерации
    N = 360  # Точность профиля (количество точек)
    profil_name = "kulachok_360"

    try:
        config = PolidainConfig(
            N_k=N_k,
            D=D,
            h=h,
            z=z,
            f_pod=f_pod,
            f_v=f_v,
            f_op=f_op,
            f_z=f_z,
            m=m,
            d=d,
            k_1=k_1,
            k_2=k_2,
            k_3=k_3,
            k_4=k_4
        )
        print("Конфигурация успешно создана!")
    except ValidationError as e:
        print("Ошибка в данных:", e)
        sys.exit(1)

    kulachok = Kulachok_polidain(config)
    fi_list = np.linspace(0, 2 * pi, N)
    X = kulachok.fun_x(fi_list).tolist() + [kulachok.fun_x(fi_list)[0]]
    Y = kulachok.fun_y(fi_list).tolist() + [kulachok.fun_y(fi_list)[0]]
    profil_config = ProfilConfig(X=X, Y=Y)
    create_profil(profil_config, profil_name=profil_name, line_type="spline")