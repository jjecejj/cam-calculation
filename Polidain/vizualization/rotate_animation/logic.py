from numpy import ndarray
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle
import matplotlib.gridspec as gridspec
from typing import Literal

from core.cam_geometry import Kulachok
from vizualization.rotate_animation.config import RotateProfileData


def rotate_profile_data(X: float | ndarray, Y: float | ndarray, angle: float):
    angle = -angle - np.pi / 2
    x_new = X * np.cos(angle) - Y * np.sin(angle)
    y_new = X * np.sin(angle) + Y * np.cos(angle)
    return (x_new, y_new)


def set_rotate_data(kulachok: Kulachok):
    if not (kulachok.kulachok_solve_flag and kulachok.tolkatel_solve_flag and kulachok.profile_solve_flag):
        raise ValueError(f"Не были проведены предварительные вычисления кулачка")

    fi_list = kulachok.kulachok_data.fi_list_rad.copy()
    movement_data = []
    for i in fi_list:
        movement_data.append(rotate_profile_data(kulachok.profile_data.X, kulachok.profile_data.Y, i))

    tolkatel_data = kulachok.tolkatel_data.H_rad.copy() + kulachok.config.D * 1e3 / 2
    tolkatel_D_t = kulachok.config.D_t * 1e3
    tolkatel_R_r = kulachok.config.R_r * 1e3

    return RotateProfileData(fi_list=fi_list,
                             movement_data=movement_data,
                             tolkatel_data=tolkatel_data,
                             tolkatel_type=kulachok.tolkatel_solve_type,
                             tolkatel_D_t=tolkatel_D_t,
                             tolkatel_R_r=tolkatel_R_r)


def display_animation(kulachok: Kulachok, interval: int = 50, save_flag: bool = False,
                      name_file: str | None = None, pause_flag=False):
    """
    Анимирует вращение кулачка и движение толкателя.
    """

    data = set_rotate_data(kulachok)

    # 1. Подготовка фигуры
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect('equal')
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.set_title("Анимация кулачкового механизма")
    ax.set_xlabel("X [мм]")
    ax.set_ylabel("Y [мм]")

    # 2. Инициализация графических элементов
    cam_line, = ax.plot([], [], 'b-', lw=2, label='Кулачок')
    ax.plot(0, 0, 'k+', markersize=10, markeredgewidth=2)

    follower_element = None

    if data.tolkatel_type == "flat":
        follower_element, = ax.plot([], [], 'r-', lw=3, label='Толкатель')

    elif data.tolkatel_type == "roller":
        follower_element = Circle((0, 0), radius=data.tolkatel_R_r, color='r', fill=True, alpha=0.5, label='Ролик')
        ax.add_patch(follower_element)

    else:  # thin
        follower_element, = ax.plot([], [], 'rv', markersize=10, label='Толкатель')

    # 3. Вычисление границ
    all_x = []
    all_y = []
    for coords in data.movement_data[::10]:  # Сэмплинг для скорости
        all_x.extend(coords[0])
        all_y.extend(coords[1])

    max_val = max(np.max(np.abs(all_x)), np.max(np.abs(all_y)))
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
        x_cam, y_cam = data.movement_data[frame_idx]
        cam_line.set_data(x_cam, y_cam)

        # current_h - это расстояние от центра до поверхности толкателя (положительное число)
        # Так как толкатель снизу, координата Y будет отрицательной: -current_h
        current_h = data.tolkatel_data[frame_idx]

        if data.tolkatel_type == "flat":
            half_w = data.tolkatel_D_t / 2
            # !!! ИЗМЕНЕНИЕ: ставим -current_h
            follower_element.set_data([-half_w, half_w], [-current_h, -current_h])

        elif data.tolkatel_type == "roller":
            # !!! ИЗМЕНЕНИЕ: Центр круга уходит вниз: -(H + R_r)
            follower_element.set_center((0, -(current_h + data.tolkatel_R_r)))

        else:  # thin
            # !!! ИЗМЕНЕНИЕ: Y = -current_h
            follower_element.set_data([0], [-current_h])

        return cam_line, follower_element

    if pause_flag:
        anim_state = {'paused': False}

        def toggle_pause(event):
            if event.inaxes is not None or event.canvas is not None:
                if anim_state['paused']:
                    ani.resume()
                    anim_state['paused'] = False
                    fig.suptitle("Анимация ЗАПУЩЕНА (клик для паузы)", fontsize=10, color='green')
                else:
                    ani.pause()
                    anim_state['paused'] = True
                    fig.suptitle("Анимация ПРИОСТАНОВЛЕНА (клик для запуска)", fontsize=10, color='red')
                plt.draw()

        fig.canvas.mpl_connect('button_press_event', toggle_pause)

    ani = FuncAnimation(
        fig, update, frames=len(data.fi_list), interval=interval, blit=False
    )

    if save_flag:
        if name_file is None:
            name_file = 'cam_animation.gif'
        elif not name_file.lower().endswith('.gif'):
            name_file += '.gif'
        ani.save('data\\output\\cam_animations\\' + name_file, writer='pillow', fps=round(1000 / interval))

    plt.show()


def display_dashboard_animation(kulachok: Kulachok,
                                interval: int = 50,
                                save_flag: bool = False,
                                pause_flag: bool = False,
                                name_file: str | None = None,
                                graphs_type: Literal['degree', 'rad', 't'] = 'degree'):
    """
    Анимирует приборную панель: слева бегущий курсор по графикам, справа вращение механизма.
    """
    # 1. Подготовка данных
    t_type = kulachok.tolkatel_solve_type
    rotate_data = set_rotate_data(kulachok, tolkatel_type=t_type)
    data_kin = kulachok.tolkatel_data

    if graphs_type == 'degree':
        x_axis = data_kin.fi_list_degree
        xlabel = r'$\phi$, град'
        y_datasets = [data_kin.H_rad, data_kin.V_degree, data_kin.A_degree, data_kin.D_degree, data_kin.K_degree]
        ylabels = ['S, мм', 'V, мм/град', 'A, мм/град²', 'J, мм/град³', 'K, мм/град⁴']
    elif graphs_type == 'rad':
        x_axis = data_kin.fi_list_rad
        xlabel = r'$\phi$, рад'
        y_datasets = [data_kin.H_rad, data_kin.V_rad, data_kin.A_rad, data_kin.D_rad, data_kin.K_rad]
        ylabels = ['S, мм', 'V, мм/рад', 'A, мм/рад²', 'J, мм/рад³', 'K, мм/рад⁴']
    else:  # 't'
        x_axis = data_kin.t_list
        xlabel = "t, c"
        y_datasets = [data_kin.H_t, data_kin.V_t, data_kin.A_t, data_kin.D_t, data_kin.K_t]
        ylabels = ['S, мм', 'V, мм/с', 'A, мм/с²', 'J, мм/с³', 'K, мм/с⁴']

    # 2. Настройка фигуры
    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(5, 3, figure=fig)

    # --- ЛЕВАЯ ЧАСТЬ: ГРАФИКИ ---
    ax_graphs = []
    cursors = []
    dots = []

    for i in range(5):
        ax = fig.add_subplot(gs[i, 0:2])
        ax.plot(x_axis, y_datasets[i], 'b-', alpha=0.6)
        ax.set_ylabel(ylabels[i], fontsize=9)
        ax.grid(True, linestyle=':')

        line = ax.axvline(x_axis[0], color='r', alpha=0.8, lw=1)
        dot, = ax.plot([x_axis[0]], [y_datasets[i][0]], 'ro', markersize=4)

        cursors.append(line)
        dots.append(dot)
        ax_graphs.append(ax)

        if i < 4:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel(xlabel)

    ax_graphs[0].set_title("Кинематика толкателя", fontsize=14)

    # --- ПРАВАЯ ЧАСТЬ: МЕХАНИЗМ ---
    ax_mech = fig.add_subplot(gs[:, 2])
    ax_mech.set_aspect('equal')
    ax_mech.grid(True, linestyle='--', alpha=0.6)
    # !!! ИЗМЕНЕНИЕ: Новое название
    ax_mech.set_title("Кулачок")

    cam_line, = ax_mech.plot([], [], 'b-', lw=2, label='Кулачок')
    ax_mech.plot(0, 0, 'k+', markersize=10)

    follower_element = None
    if rotate_data.tolkatel_type == "flat":
        follower_element, = ax_mech.plot([], [], 'r-', lw=3, label='Толкатель')
    elif rotate_data.tolkatel_type == "roller":
        follower_element = Circle((0, 0), radius=rotate_data.tolkatel_R_r, color='r', fill=True, alpha=0.5,
                                  label='Ролик')
        ax_mech.add_patch(follower_element)
    else:  # thin
        follower_element, = ax_mech.plot([], [], 'rv', markersize=10, label='Толкатель')

    # Границы
    all_x_cam = []
    all_y_cam = []
    for coords in rotate_data.movement_data[::10]:
        all_x_cam.extend(coords[0])
        all_y_cam.extend(coords[1])

    max_val = max(np.max(np.abs(all_x_cam)), np.max(np.abs(all_y_cam)))
    max_h = np.max(rotate_data.tolkatel_data)

    if rotate_data.tolkatel_type == "roller":
        limit = max(max_val, max_h + rotate_data.tolkatel_R_r + 10)
    else:
        limit = max(max_val, max_h + 10)

    ax_mech.set_xlim(-limit, limit)
    ax_mech.set_ylim(-limit, limit)

    # 3. Функция обновления
    def update(frame_idx):
        updated_artists = []

        # Графики
        current_x = x_axis[frame_idx]
        for i in range(5):
            current_y = y_datasets[i][frame_idx]
            cursors[i].set_xdata([current_x, current_x])
            updated_artists.append(cursors[i])
            dots[i].set_data([current_x], [current_y])
            updated_artists.append(dots[i])

        # Механизм
        x_c, y_c = rotate_data.movement_data[frame_idx]
        cam_line.set_data(x_c, y_c)
        updated_artists.append(cam_line)

        current_h = rotate_data.tolkatel_data[frame_idx]

        if rotate_data.tolkatel_type == "flat":
            half_w = rotate_data.tolkatel_D_t / 2
            # !!! ИЗМЕНЕНИЕ: -current_h
            follower_element.set_data([-half_w, half_w], [-current_h, -current_h])
            updated_artists.append(follower_element)

        elif rotate_data.tolkatel_type == "roller":
            # !!! ИЗМЕНЕНИЕ: -(H + R_r)
            follower_element.set_center((0, -(current_h + rotate_data.tolkatel_R_r)))
            updated_artists.append(follower_element)

        else:
            # !!! ИЗМЕНЕНИЕ: -current_h
            follower_element.set_data([0], [-current_h])
            updated_artists.append(follower_element)

        return updated_artists

    # Логика паузы
    anim_state = {'paused': False}

    def toggle_pause(event):
        if event.inaxes is not None or event.canvas is not None:
            if anim_state['paused']:
                ani.resume()
                anim_state['paused'] = False
                fig.suptitle("Анимация ЗАПУЩЕНА (клик для паузы)", fontsize=10, color='green')
            else:
                ani.pause()
                anim_state['paused'] = True
                fig.suptitle("Анимация ПРИОСТАНОВЛЕНА (клик для запуска)", fontsize=10, color='red')
            plt.draw()

    if pause_flag:
        fig.canvas.mpl_connect('button_press_event', toggle_pause)

    ani = FuncAnimation(
        fig, update, frames=len(x_axis), interval=interval, blit=False
    )

    if save_flag:
        if name_file is None:
            name_file = 'dashboard_animation.gif'
        elif not name_file.lower().endswith('.gif'):
            name_file += '.gif'
        ani.save('data\\output\\cam_animations\\' + name_file, writer='pillow', fps=round(1000 / interval))
        print(f"Анимация сохранена в {name_file}")

    plt.tight_layout()
    plt.show()