import numpy as np
from numpy import ndarray
from pydantic import BaseModel, Field, ConfigDict, model_validator
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle

class RotateProfileData(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    fi_list: np.ndarray[float | np.ndarray] | list[float | np.ndarray]
    movement_data: list[tuple[ndarray | list, ndarray | list]]
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

def rotate_profile_data(X, Y, angle):
    angle = np.radians(-angle) + np.pi / 2
    x_new = X * np.cos(angle) - Y * np.sin(angle)
    y_new = X * np.sin(angle) + Y * np.cos(angle)
    return (x_new, y_new)

def set_rotate_data(kulachok, tolkatel_type):
    if not(kulachok.kulachok_solve_flag and kulachok.tolkatel_solve_flag and kulachok.profil_solve_flag):
        raise ValueError(f"Не были проведены предварительные вычисления кулачка")

    fi_list = kulachok.kulachok_data.fi_list.copy()
    movement_data = []
    for i in fi_list:
        movement_data.append(rotate_profile_data(kulachok.profil_data.X, kulachok.profil_data.Y, i))
    tolkatel_data = kulachok.tolkatel_data.H.copy() + kulachok.config.D * 1e3 / 2
    tolkatel_D_t = kulachok.config.D_t * 1e3
    tolkatel_R_r = kulachok.config.R_r * 1e3

    return RotateProfileData(fi_list = fi_list,
                             movement_data = movement_data,
                             tolkatel_data = tolkatel_data,
                             tolkatel_type = tolkatel_type,
                             tolkatel_D_t = tolkatel_D_t,
                             tolkatel_R_r = tolkatel_R_r)


def display_animation(data: RotateProfileData, interval: int = 50, save_flag: bool = False, name_file: str | None = None):
    """
    Анимирует вращение кулачка и движение толкателя.

    :param data: Объект RotateProfileData с рассчитанными данными
    :param interval: Интервал между кадрами в миллисекундах
    """

    # 1. Подготовка фигуры
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect('equal')
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.set_title("Анимация кулачкового механизма")
    ax.set_xlabel("X [мм]")
    ax.set_ylabel("Y [мм]")

    # 2. Инициализация графических элементов
    # Линия профиля кулачка
    cam_line, = ax.plot([], [], 'b-', lw=2, label='Кулачок')
    # Точка центра вращения
    ax.plot(0, 0, 'k+', markersize=10, markeredgewidth=2)

    # Элементы толкателя (зависят от типа)
    follower_element = None

    # Для тарельчатого толкателя (линия)
    if data.tolkatel_type == "flat":
        # Создаем линию тарелки
        follower_element, = ax.plot([], [], 'r-', lw=3, label='Толкатель')

    # Для роликового толкателя (круг)
    elif data.tolkatel_type == "roller":
        # Создаем круг (Patch), начальная позиция (0,0) будет обновлена
        follower_element = Circle((0, 0), radius=data.tolkatel_R_r, color='r', fill=True, alpha=0.5, label='Ролик')
        ax.add_patch(follower_element)

    # Для острого толкателя (маркер)
    else:  # thin
        follower_element, = ax.plot([], [], 'rv', markersize=10, label='Толкатель')

    # 3. Вычисление границ графика (чтобы оси не "скакали")
    # Собираем все координаты X и Y кулачка для определения максимума
    all_x = []
    all_y = []
    # Берем каждый 10-й кадр для ускорения расчета границ (можно и все)
    for coords in data.movement_data:
        all_x.extend(coords[0])
        all_y.extend(coords[1])

    max_val = max(np.max(np.abs(all_x)), np.max(np.abs(all_y)))

    # Учитываем вылет толкателя для верхней границы
    max_h = np.max(data.tolkatel_data)
    if data.tolkatel_type == "flat":
        width = data.tolkatel_D_t / 2
        limit = max(max_val, max_h + 10, width + 10)
    elif data.tolkatel_type == "roller":
        limit = max(max_val, max_h + data.tolkatel_R_r + 10)
    else:
        limit = max(max_val, max_h + 10)

    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.legend(loc='upper right')

    # 4. Функция обновления кадра
    def update(frame_idx):
        # -- Обновление кулачка --
        # movement_data[frame_idx] это кортеж (X_array, Y_array)
        x_cam, y_cam = data.movement_data[frame_idx]
        cam_line.set_data(x_cam, y_cam)

        # -- Обновление толкателя --
        current_h = data.tolkatel_data[frame_idx]

        if data.tolkatel_type == "flat":
            # Рисуем горизонтальную линию длиной D_t на высоте current_h
            half_w = data.tolkatel_D_t / 2
            follower_element.set_data([-half_w, half_w], [current_h, current_h])

        elif data.tolkatel_type == "roller":
            # Перемещаем круг. Центр круга по X=0 (обычно), по Y = current_h
            # Важно: current_h обычно координаты центра ролика.
            # Если это координата касания, нужно добавить радиус: (0, current_h + R_r)
            follower_element.set_center((0, current_h))

        else:  # thin
            # Просто точка/маркер
            follower_element.set_data([0], [current_h])

        return cam_line, follower_element

    # 5. Запуск анимации
    ani = FuncAnimation(
        fig,
        update,
        frames=len(data.fi_list),
        interval=interval,
        blit=False
    )

    if save_flag:
        if name_file is None:
            name_file = 'cam_animation.gif'
        else:
            if not name_file.lower().endswith('.gif'):
                name_file += '.gif'
        ani.save('data\\output\\cam_animations\\' + name_file, writer='pillow', fps=round(1000 / interval))

    plt.show()
