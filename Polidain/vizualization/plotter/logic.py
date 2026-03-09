import numpy as np
import matplotlib.pyplot as plt
from typing import Literal, Union, List

from matplotlib import gridspec
from numpy import ndarray

from vizualization.plotter.config import PlotConfig
from core.cam_geometry import Kulachok
from core.cam_geometry.schemas import ProfileData, GraphData

fig_size: tuple = (12, 12)

def set_config(config: PlotConfig):
    """
    Применяет глобальные настройки стилей matplotlib из конфигурации.

    Args:
        config: Объект конфигурации PlotConfig с параметрами шрифтов и стилей.
    """
    rc = config.rc
    plt.rcParams.update(rc)
    plt.rcParams["font.serif"] = config.font_serif
    plt.rcParams["font.family"] = config.font_family
    plt.rcParams['font.size'] = config.font_size
    plt.rcParams['mathtext.fontset'] = config.mathtext_fontset
    plt.rcParams['mathtext.rm'] = config.mathtext_rm
    plt.rcParams['mathtext.it'] = config.mathtext_it


def _plot_component(ax, x: list | ndarray, y: list | ndarray, ylabel: str, xlabel: str):
    """
    Вспомогательная функция для построения одного графика с отметкой минимума и максимума.

    Args:
        ax: Объект осей matplotlib.
        x: Данные оси X.
        y: Данные оси Y.
        ylabel: Подпись оси Y.
        xlabel: Подпись оси X.
    """
    x_arr = np.asarray(x)
    y_arr = np.asarray(y)

    ax.plot(x_arr, y_arr)

    idx_max = np.argmax(y_arr)
    idx_min = np.argmin(y_arr)

    ax.scatter(x_arr[idx_max], y_arr[idx_max], color='r')
    ax.scatter(x_arr[idx_min], y_arr[idx_min], color='r')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True)


'''
def _plot_vertical_lines(ax, phi_0, phi_1, phi_2, phi_3, phi_4, phi_5, initial_angle: float = 0.0):
    ax.axvline(x=phi_0 - initial_angle, color='red', linestyle='--')
    ax.axvline(x=phi_1 - initial_angle, color='red', linestyle='--')
    ax.axvline(x=phi_2 - initial_angle, color='red', linestyle='--')
    ax.axvline(x=phi_3 - initial_angle, color='red', linestyle='--')
    ax.axvline(x=phi_4 - initial_angle, color='red', linestyle='--')
    ax.axvline(x=phi_5 - initial_angle, color='red', linestyle='--')
'''


def _plot_vertical_lines(ax, phi_0, phi_1, phi_2, phi_3, phi_4, phi_5, initial_angle: float = 0.0):
    """
    Рисует вертикальные линии фаз и подписывает их.
    """
    phis = [phi_0, phi_1, phi_2, phi_3, phi_4, phi_5]
    labels = [r'$\phi_0$', r'$\phi_1$', r'$\phi_2$', r'$\phi_3$', r'$\phi_4$', r'$\phi_5$']

    # transform=ax.get_xaxis_transform() позволяет задавать X в координатах данных графика,
    # а Y в относительных координатах осей (где 0 - самый низ, 1 - самый верх)
    trans = ax.get_xaxis_transform()

    for phi, label in zip(phis, labels):
        x_val = phi - initial_angle
        ax.axvline(x=x_val, color='red', linestyle='--', alpha=0.6)
        ax.text(x_val, 0.95, f' {label} ', color='darkred', transform=trans,
                ha='left', va='top', fontsize=11,
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=1))

def set_initial_angle_degree(_fi_list: list | np.ndarray, initial_angle: float | int = 0.0):
    """
    Вычисляет индексы сортировки для сдвига массива углов (в градусах).

    Args:
        _fi_list: Массив углов.
        initial_angle: Начальный угол сдвига.

    Returns:
        np.ndarray: Индексы сортировки.
    """
    __fi_list = np.asarray(_fi_list)
    arr = np.mod(__fi_list - initial_angle, 360)
    i_list = np.argsort(arr)
    return i_list

def set_initial_angle_rad(_fi_list: list | np.ndarray, initial_angle: float | int = 0.0):
    """
    Вычисляет индексы сортировки для сдвига массива углов (в радианах).

    Args:
        _fi_list: Массив углов.
        initial_angle: Начальный угол сдвига.

    Returns:
        np.ndarray: Индексы сортировки.
    """
    __fi_list = np.asarray(_fi_list)
    arr = np.mod(__fi_list - initial_angle, 2 * np.pi)
    i_list = np.argsort(arr)
    return i_list

def set_initial_angle_t(_fi_list: list | np.ndarray, initial_angle: float | int = 0.0):
    """
    Вычисляет индексы сортировки для сдвига массива времени.

    Args:
        _fi_list: Массив времени.
        initial_angle: Начальное время сдвига.

    Returns:
        np.ndarray: Индексы сортировки.
    """
    __fi_list = np.asarray(_fi_list)
    arr = np.mod(__fi_list - initial_angle, np.max(_fi_list))
    i_list = np.argsort(arr)
    return i_list

def set_initial_angle(_fi_list: list | np.ndarray, initial_angle: float | int = 0.0, angle_type = Literal['t', 'rad', 'degree']):
    """
    Универсальная функция для сдвига массива аргументов (фазировка графика).

    Args:
        _fi_list: Массив аргументов (углы или время).
        initial_angle: Значение сдвига.
        angle_type: Тип аргумента ('t', 'rad', 'degree').

    Returns:
        np.ndarray: Индексы сортировки.
    """
    if angle_type == 'rad':
        return set_initial_angle_rad(_fi_list, initial_angle = initial_angle)
    elif angle_type == 'degree':
        return set_initial_angle_degree(_fi_list, initial_angle = initial_angle)
    else:
        return set_initial_angle_t(_fi_list, initial_angle = initial_angle)

def calculate_optimal_angle(kulachok: Kulachok):
    """
    Рассчитывает оптимальный начальный угол для отображения графиков (центрирование рабочего хода).

    Args:
        kulachok: Объект кулачка.

    Returns:
        float: Оптимальный угол в градусах.
    """
    return np.degrees(-(2 * np.pi - kulachok.config.phi_5) / 2)


def display_graphs_kulachok(data: GraphData, initial_angle: float | int = 0,
                            graphs_type: Literal['t', 'rad', 'degree'] = 'degree'):
    """
    Отображает графики кинематики для кулачка (H, V, A, J, K).

    Args:
        data: Данные графиков (GraphData).
        initial_angle: Начальный угол сдвига.
        graphs_type: Тип аргумента ('t', 'rad', 'degree').
    """
    fig, axs = plt.subplots(5, 1, figsize=fig_size)
    fig.suptitle("Графики для кулачка", fontsize=16)

    if graphs_type == 'degree':
        x_data = data.fi_list_degree
        xl = r'$\phi$, град'
        i_list = set_initial_angle(x_data, initial_angle=initial_angle, angle_type = graphs_type)
        _plot_component(axs[0], x_data, data.H_rad[i_list], 'Координата точки контакта, мм', xl)
        _plot_component(axs[1], x_data, data.V_degree[i_list], 'Скорость, $мм/град$', xl)
        _plot_component(axs[2], x_data, data.A_degree[i_list], 'Ускорение, $мм/град^2$', xl)
        _plot_component(axs[3], x_data, data.D_degree[i_list], 'Рывок, $мм/град^3$', xl)
        _plot_component(axs[4], x_data, data.K_degree[i_list], 'Кракен, $мм/град^4$', xl)
        for ax in axs:
            _plot_vertical_lines(ax, data.phi_0_degree, data.phi_1_degree, data.phi_2_degree, data.phi_3_degree, data.phi_4_degree, data.phi_5_degree, initial_angle=initial_angle)

    elif graphs_type == 'rad':
        x_data = data.fi_list_rad
        xl = r'$\phi$, рад'
        initial_angle_rad = initial_angle * np.pi / 180
        i_list = set_initial_angle(x_data, initial_angle=initial_angle_rad, angle_type = graphs_type)
        _plot_component(axs[0], x_data, data.H_rad[i_list], 'Координата точки контакта, мм', xl)
        _plot_component(axs[1], x_data, data.V_rad[i_list], 'Скорость, $мм/рад$', xl)
        _plot_component(axs[2], x_data, data.A_rad[i_list], 'Ускорение, $мм/рад^2$', xl)
        _plot_component(axs[3], x_data, data.D_rad[i_list], 'Рывок, $мм/рад^3$', xl)
        _plot_component(axs[4], x_data, data.K_rad[i_list], 'Кракен, $мм/рад^4$', xl)
        for ax in axs:
            _plot_vertical_lines(ax, data.phi_0_rad, data.phi_1_rad, data.phi_2_rad, data.phi_3_rad, data.phi_4_rad, data.phi_5_rad, initial_angle=initial_angle_rad)

    elif graphs_type == 't':
        x_data = data.t_list
        xl = "t, c"
        initial_time = initial_angle * np.pi / 180 / data.omega_rad if initial_angle != 0 else 0
        i_list = set_initial_angle(x_data, initial_angle=initial_time, angle_type = graphs_type)
        _plot_component(axs[0], x_data, data.H_t[i_list], 'Координата точки контакта, мм', xl)
        _plot_component(axs[1], x_data, data.V_t[i_list], 'Скорость, $мм/с$', xl)
        _plot_component(axs[2], x_data, data.A_t[i_list], 'Ускорение, $мм/с^2$', xl)
        _plot_component(axs[3], x_data, data.D_t[i_list], 'Рывок, $мм/с^3$', xl)
        _plot_component(axs[4], x_data, data.K_t[i_list], 'Кракен, $мм/с^4$', xl)
        for ax in axs:
            _plot_vertical_lines(ax, data.phi_0_t, data.phi_1_t, data.phi_2_t, data.phi_3_t, data.phi_4_t, data.phi_5_t, initial_angle=initial_time)


    plt.tight_layout()
    plt.show()


def display_graphs_tolkatel(data: GraphData, initial_angle: float | int = 0,
                            graphs_type: Literal['t', 'rad', 'degree'] = 'degree'):
    """
    Отображает графики кинематики для толкателя.

    Args:
        data: Данные графиков (GraphData).
        initial_angle: Начальный угол сдвига.
        graphs_type: Тип аргумента ('t', 'rad', 'degree').
    """
    fig, axs = plt.subplots(5, 1, figsize=fig_size)
    fig.suptitle("Графики для толкателя", fontsize=16)

    if graphs_type == 'degree':
        x_data = data.fi_list_degree
        xl = r'$\phi$, град'
        i_list = set_initial_angle(x_data, initial_angle=initial_angle, angle_type=graphs_type)
        _plot_component(axs[0], x_data, data.H_rad[i_list], 'Перемещение толкателя, мм', xl)
        _plot_component(axs[1], x_data, data.V_t[i_list], 'Скорость, $мм/град$', xl)
        _plot_component(axs[2], x_data, data.A_t[i_list], 'Ускорение, $мм/град^2$', xl)
        _plot_component(axs[3], x_data, data.D_t[i_list], 'Рывок, $мм/град^3$', xl)
        _plot_component(axs[4], x_data, data.K_t[i_list], 'Кракен, $мм/град^4$', xl)

        for ax in axs:
            _plot_vertical_lines(ax, data.phi_0_degree, data.phi_1_degree, data.phi_2_degree, data.phi_3_degree,
                                 data.phi_4_degree, data.phi_5_degree, initial_angle=initial_angle)

    elif graphs_type == 'rad':
        x_data = data.fi_list_rad
        xl = r'$\phi$, рад'
        initial_angle_rad = initial_angle * np.pi / 180
        i_list = set_initial_angle(x_data, initial_angle=initial_angle_rad, angle_type=graphs_type)
        _plot_component(axs[0], x_data, data.H_rad[i_list], 'Перемещение толкателя, мм', xl)
        _plot_component(axs[1], x_data, data.V_t[i_list], 'Скорость, $мм/рад$', xl)
        _plot_component(axs[2], x_data, data.A_t[i_list], 'Ускорение, $мм/рад^2$', xl)
        _plot_component(axs[3], x_data, data.D_t[i_list], 'Рывок, $мм/рад^3$', xl)
        _plot_component(axs[4], x_data, data.K_t[i_list], 'Кракен, $мм/рад^4$', xl)

        for ax in axs:
            _plot_vertical_lines(ax, data.phi_0_rad, data.phi_1_rad, data.phi_2_rad, data.phi_3_rad, data.phi_4_rad,
                                 data.phi_5_rad, initial_angle=initial_angle_rad)

    elif graphs_type == 't':
        x_data = data.t_list
        xl = "t, c"
        initial_time = initial_angle * np.pi / 180 / data.omega_rad if initial_angle != 0 else 0
        i_list = set_initial_angle(x_data, initial_angle=initial_time, angle_type=graphs_type)
        _plot_component(axs[0], x_data, data.H_t[i_list], 'Перемещение толкателя, мм', xl)
        _plot_component(axs[1], x_data, data.V_t[i_list], 'Скорость, $мм/с$', xl)
        _plot_component(axs[2], x_data, data.A_t[i_list], 'Ускорение, $мм/с^2$', xl)
        _plot_component(axs[3], x_data, data.D_t[i_list], 'Рывок, $мм/с^3$', xl)
        _plot_component(axs[4], x_data, data.K_t[i_list], 'Кракен, $мм/с^4$', xl)

        for ax in axs:
            _plot_vertical_lines(ax, data.phi_0_t, data.phi_1_t, data.phi_2_t, data.phi_3_t, data.phi_4_t, data.phi_5_t,
                                 initial_angle=initial_time)

    plt.tight_layout()
    plt.show()

def display_profile(data: ProfileData, initial_angle: float | int = 0):
    """
    Отображает профиль кулачка.

    Args:
        data: Данные профиля (ProfileData).
        initial_angle: Начальный угол сдвига.
    """
    plt.figure(figsize=(6, 6))
    i_list = set_initial_angle(data.fi_list, initial_angle=initial_angle)
    X = np.asarray(data.X)
    Y = np.asarray(data.Y)
    X_sorted = X[i_list]
    Y_sorted = Y[i_list]
    X_plot = np.append(X_sorted, X_sorted[0])
    Y_plot = np.append(Y_sorted, Y_sorted[0])

    plt.plot(X_plot, Y_plot)
    plt.scatter([0], [0], color='black', marker='x', label='Центр вращения')
    plt.xlabel('X, мм')
    plt.ylabel('Y, мм')
    margin = 0.01
    plt.xlim(np.min(X) - margin, np.max(X) + margin)
    plt.ylim(np.min(Y) - margin, np.max(Y) + margin)
    plt.grid(True)
    plt.title("Профиль кулачка")
    plt.axis('equal')
    plt.legend()
    plt.show()


def display_profile(data: ProfileData, initial_angle: float | int = 0):
    """
    Отображает профиль кулачка с пунктирными линиями фаз.

    Args:
        data: Данные профиля (ProfileData).
        initial_angle: Начальный угол сдвига.
    """
    plt.figure(figsize=(6, 6))

    # Сортировка точек для отрисовки профиля
    i_list = set_initial_angle(data.fi_list, initial_angle=initial_angle)
    X = np.asarray(data.X)
    Y = np.asarray(data.Y)
    X_sorted = X[i_list]
    Y_sorted = Y[i_list]
    X_plot = np.append(X_sorted, X_sorted[0])
    Y_plot = np.append(Y_sorted, Y_sorted[0])

    plt.plot(X_plot, Y_plot, linewidth=2, label='Профиль')
    plt.scatter([0], [0], color='black', marker='x', s=100, label='Центр вращения')

    # --- Отрисовка линий фаз ---
    phis_rad = [data.phi_0_rad, data.phi_1_rad, data.phi_2_rad, data.phi_3_rad, data.phi_4_rad, data.phi_5_rad]
    labels = [r'$\phi_0$', r'$\phi_1$', r'$\phi_2$', r'$\phi_3$', r'$\phi_4$', r'$\phi_5$']

    # Массив углов профиля (в градусах)
    fi_list_deg = np.asarray(data.fi_list)

    for phi_rad, label in zip(phis_rad, labels):
        phi_deg = np.degrees(phi_rad)

        # Находим ближайший к фазе индекс в массиве углов (с учетом перехода через 360)
        diff = np.abs(np.mod(fi_list_deg - phi_deg + 180, 360) - 180)
        idx = np.argmin(diff)

        # Координаты точки на самом профиле
        x_end = X[idx]
        y_end = Y[idx]

        # Рисуем пунктирную линию от центра до края профиля
        plt.plot([0, x_end], [0, y_end], color='gray', linestyle='--', alpha=0.7)

        # Добавляем текстовую подпись чуть дальше от края (выносим на 15%)
        length = np.hypot(x_end, y_end)
        if length > 0:
            text_x = x_end * 1.15
            text_y = y_end * 1.15
            plt.text(text_x, text_y, label, color='darkred', ha='center', va='center',
                     fontsize=12, bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', pad=1))
    # ---------------------------

    plt.xlabel('X, мм')
    plt.ylabel('Y, мм')

    # Увеличиваем отступы, чтобы поместились текстовые подписи фаз
    margin_x = (np.max(X) - np.min(X)) * 0.15
    margin_y = (np.max(Y) - np.min(Y)) * 0.15
    plt.xlim(np.min(X) - margin_x, np.max(X) + margin_x)
    plt.ylim(np.min(Y) - margin_y, np.max(Y) + margin_y)

    plt.grid(True)
    plt.title("Профиль кулачка\n")
    plt.axis('equal')
    plt.legend(loc='upper right')
    plt.show()


def display_profile_multiplicity(data: List[ProfileData], initial_angle: float | int = 0):
    plt.figure(figsize=(6, 6))

    all_X = []
    all_Y = []

    for profile in data:
        i_list = set_initial_angle(profile.fi_list, initial_angle=initial_angle)
        X = np.asarray(profile.X)
        Y = np.asarray(profile.Y)
        X_sorted = X[i_list]
        Y_sorted = Y[i_list]
        X_plot = np.append(X_sorted, X_sorted[0])
        Y_plot = np.append(Y_sorted, Y_sorted[0])

        plt.plot(X_plot, Y_plot)
        plt.scatter([0], [0], color='black', marker='x')

        all_X.extend(X)
        all_Y.extend(Y)

        # --- Отрисовка линий фаз для каждого профиля ---
        phis_rad = [profile.phi_0_rad, profile.phi_1_rad, profile.phi_2_rad,
                    profile.phi_3_rad, profile.phi_4_rad, profile.phi_5_rad]
        labels = [r'$\phi_0$', r'$\phi_1$', r'$\phi_2$', r'$\phi_3$', r'$\phi_4$', r'$\phi_5$']
        fi_list_deg = np.asarray(profile.fi_list)

        for phi_rad, label in zip(phis_rad, labels):
            phi_deg = np.degrees(phi_rad)
            diff = np.abs(np.mod(fi_list_deg - phi_deg + 180, 360) - 180)
            idx = np.argmin(diff)
            x_end = X[idx]
            y_end = Y[idx]

            plt.plot([0, x_end], [0, y_end], color='gray', linestyle='--', alpha=0.5)

            length = np.hypot(x_end, y_end)
            if length > 0:
                plt.text(x_end * 1.15, y_end * 1.15, label, color='darkred', ha='center', va='center',
                         fontsize=12, bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', pad=1))
        # -----------------------------------------------

    plt.xlabel('X, мм')
    plt.ylabel('Y, мм')

    # Динамические отступы с учетом всех профилей
    margin_x = (np.max(all_X) - np.min(all_X)) * 0.15
    margin_y = (np.max(all_Y) - np.min(all_Y)) * 0.15
    plt.xlim(np.min(all_X) - margin_x, np.max(all_X) + margin_x)
    plt.ylim(np.min(all_Y) - margin_y, np.max(all_Y) + margin_y)

    plt.grid(True)
    plt.title("Профили кулачка")
    plt.axis('equal')
    plt.show()


def display_all(kulachok: Kulachok, initial_angle: Union[float, int, str] = 0,
                graphs_type: Literal['t', 'rad', 'degree'] = 'degree'):
    """
    Отображает все графики: кинематику кулачка, кинематику толкателя и профиль.

    Args:
        kulachok: Объект кулачка.
        initial_angle: Начальный угол сдвига (или 'auto').
        graphs_type: Тип аргумента ('t', 'rad', 'degree').
    """
    if initial_angle == "auto":
        target_angle = calculate_optimal_angle(kulachok)
    elif isinstance(initial_angle, (int, float)):
        target_angle = initial_angle
    else:
        raise ValueError("initial_angle должен быть числом или 'auto'")

    display_graphs_kulachok(kulachok.kulachok_data, initial_angle=target_angle, graphs_type=graphs_type)
    display_graphs_tolkatel(kulachok.tolkatel_data, initial_angle=target_angle, graphs_type=graphs_type)
    display_profile(kulachok.profile_data, initial_angle=target_angle)


def display_dashboard(kulachok: Kulachok,
                      initial_angle: Union[float, int, str] = 0,
                      graphs_type: Literal['t', 'rad', 'degree'] = 'degree',
                      target: Literal['tolkatel', 'kulachok'] = 'tolkatel'):
    """
    Отображает единую панель (Dashboard): слева 5 графиков кинематики, справа профиль кулачка.

    Args:
        kulachok: Объект механизма.
        initial_angle: Начальный угол (число или 'auto').
        graphs_type: Тип оси X ('degree', 'rad', 't').
        target: Чьи графики кинематики строить: 'tolkatel' (толкатель) или 'kulachok' (кулачок).
    """
    # 1. Определяем угол сдвига
    if initial_angle == "auto":
        target_angle = calculate_optimal_angle(kulachok)
    elif isinstance(initial_angle, (int, float)):
        target_angle = initial_angle
    else:
        raise ValueError("initial_angle должен быть числом или 'auto'")

    # 2. Выбираем данные для кинематики (Толкатель или Кулачок)
    if target == 'tolkatel':
        data = kulachok.tolkatel_data
        base_title = "Кинематика толкателя"
        ylabels = ['Перемещение, мм', 'Скорость', 'Ускорение', 'Рывок', 'Кракен']
    else:
        data = kulachok.kulachok_data
        base_title = "Кинематика кулачка"
        ylabels = ['Координата, мм', 'Скорость', 'Ускорение', 'Рывок', 'Кракен']

    # 3. Подготовка данных в зависимости от типа оси X (t, rad, degree)
    if graphs_type == 'degree':
        x_data = data.fi_list_degree
        xlabel = r'$\phi$, град'
        units = ['', r'$мм/град$', r'$мм/град^2$', r'$мм/град^3$', r'$мм/град^4$']
        i_list = set_initial_angle(x_data, initial_angle=target_angle, angle_type=graphs_type)
        y_datasets = [data.H_degree, data.V_degree, data.A_degree, data.D_degree, data.K_degree]

        # Данные для вертикальных линий
        phis = (data.phi_0_degree, data.phi_1_degree, data.phi_2_degree, data.phi_3_degree, data.phi_4_degree,
                data.phi_5_degree)
        shift_val = target_angle

    elif graphs_type == 'rad':
        x_data = data.fi_list_rad
        xlabel = r'$\phi$, рад'
        units = ['', r'$мм/рад$', r'$мм/рад^2$', r'$мм/рад^3$', r'$мм/рад^4$']
        initial_angle_rad = target_angle * np.pi / 180
        i_list = set_initial_angle(x_data, initial_angle=initial_angle_rad, angle_type=graphs_type)
        y_datasets = [data.H_rad, data.V_rad, data.A_rad, data.D_rad, data.K_rad]

        phis = (data.phi_0_rad, data.phi_1_rad, data.phi_2_rad, data.phi_3_rad, data.phi_4_rad, data.phi_5_rad)
        shift_val = initial_angle_rad

    elif graphs_type == 't':
        x_data = data.t_list
        xlabel = "t, c"
        units = ['', r'$мм/с$', r'$мм/с^2$', r'$мм/с^3$', r'$мм/с^4$']
        initial_time = target_angle * np.pi / 180 / data.omega_rad if target_angle != 0 else 0
        i_list = set_initial_angle(x_data, initial_angle=initial_time, angle_type=graphs_type)
        y_datasets = [data.H_t, data.V_t, data.A_t, data.D_t, data.K_t]

        phis = (data.phi_0_t, data.phi_1_t, data.phi_2_t, data.phi_3_t, data.phi_4_t, data.phi_5_t)
        shift_val = initial_time

    # --- СОЗДАНИЕ ГРАФИЧЕСКОГО ОКНА ---
    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(5, 3, figure=fig)

    # === ЛЕВАЯ ЧАСТЬ: Кинематика ===
    ax_kinematics = []
    for i in range(5):
        ax = fig.add_subplot(gs[i, 0:2])
        full_ylabel = f"{ylabels[i]} {units[i]}" if i > 0 else ylabels[i]
        _plot_component(ax, x_data, y_datasets[i][i_list], full_ylabel, xlabel)

        if i < 4:
            ax.set_xlabel("")
            ax.tick_params(labelbottom=False)
        ax_kinematics.append(ax)

    ax_kinematics[0].set_title(base_title, fontsize=14)

    # Добавляем вертикальные линии на все графики кинематики
    for ax in ax_kinematics:
        _plot_vertical_lines(ax, *phis, initial_angle=shift_val)

    # === ПРАВАЯ ЧАСТЬ: Профиль ===
    ax_profile = fig.add_subplot(gs[:, 2])

    profile_data = kulachok.profile_data

    i_list_prof = set_initial_angle(profile_data.fi_list, initial_angle=target_angle)
    X = np.asarray(profile_data.X)
    Y = np.asarray(profile_data.Y)
    X_sorted = X[i_list_prof]
    Y_sorted = Y[i_list_prof]
    X_plot = np.append(X_sorted, X_sorted[0])
    Y_plot = np.append(Y_sorted, Y_sorted[0])

    ax_profile.plot(X_plot, Y_plot, linewidth=2, label='Профиль')
    ax_profile.scatter([0], [0], color='black', marker='x', s=100, label='Центр вращения')

    # --- Отрисовка линий фаз на дашборде ---
    phis_rad_prof = [profile_data.phi_0_rad, profile_data.phi_1_rad, profile_data.phi_2_rad,
                     profile_data.phi_3_rad, profile_data.phi_4_rad, profile_data.phi_5_rad]
    labels_prof = [r'$\phi_0$', r'$\phi_1$', r'$\phi_2$', r'$\phi_3$', r'$\phi_4$', r'$\phi_5$']
    fi_list_deg_prof = np.asarray(profile_data.fi_list)

    for phi_rad, label in zip(phis_rad_prof, labels_prof):
        phi_deg = np.degrees(phi_rad)
        diff = np.abs(np.mod(fi_list_deg_prof - phi_deg + 180, 360) - 180)
        idx = np.argmin(diff)
        x_end = X[idx]
        y_end = Y[idx]

        ax_profile.plot([0, x_end], [0, y_end], color='gray', linestyle='--', alpha=0.7)

        length = np.hypot(x_end, y_end)
        if length > 0:
            ax_profile.text(x_end * 1.15, y_end * 1.15, label, color='darkred', ha='center', va='center',
                            fontsize=12, bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', pad=1))
    # ---------------------------------------

    ax_profile.set_xlabel('X, мм')
    ax_profile.set_ylabel('Y, мм')
    ax_profile.set_title("Профиль кулачка", fontsize=14)
    ax_profile.axis('equal')
    ax_profile.grid(True)
    ax_profile.legend(loc='upper right')

    # Динамические отступы
    margin_x = (np.max(X) - np.min(X)) * 0.15
    margin_y = (np.max(Y) - np.min(Y)) * 0.15
    ax_profile.set_xlim(np.min(X) - margin_x, np.max(X) + margin_x)
    ax_profile.set_ylim(np.min(Y) - margin_y, np.max(Y) + margin_y)

    plt.tight_layout()
    plt.show()