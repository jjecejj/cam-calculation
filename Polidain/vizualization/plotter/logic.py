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
        i_list = set_initial_angle(x_data, initial_angle=initial_angle, angle_type = graphs_type)
        _plot_component(axs[0], x_data, data.H_rad[i_list], 'Перемещение толкателя, мм', xl)
        _plot_component(axs[1], x_data, data.V_t[i_list], 'Скорость, $мм/град$', xl)
        _plot_component(axs[2], x_data, data.A_t[i_list], 'Ускорение, $мм/град^2$', xl)
        _plot_component(axs[3], x_data, data.D_t[i_list], 'Рывок, $мм/град^3$', xl)
        _plot_component(axs[4], x_data, data.K_t[i_list], 'Кракен, $мм/град^4$', xl)

    elif graphs_type == 'rad':
        x_data = data.fi_list_rad
        xl = r'$\phi$, рад'
        initial_angle_rad = initial_angle * np.pi / 180
        i_list = set_initial_angle(x_data, initial_angle=initial_angle_rad, angle_type = graphs_type)
        _plot_component(axs[0], x_data, data.H_rad[i_list], 'Перемещение толкателя, мм', xl)
        _plot_component(axs[1], x_data, data.V_t[i_list], 'Скорость, $мм/рад$', xl)
        _plot_component(axs[2], x_data, data.A_t[i_list], 'Ускорение, $мм/рад^2$', xl)
        _plot_component(axs[3], x_data, data.D_t[i_list], 'Рывок, $мм/рад^3$', xl)
        _plot_component(axs[4], x_data, data.K_t[i_list], 'Кракен, $мм/рад^4$', xl)

    elif graphs_type == 't':
        x_data = data.t_list
        xl = "t, c"
        initial_time = initial_angle * np.pi / 180 / data.omega_rad if initial_angle != 0 else 0
        i_list = set_initial_angle(x_data, initial_angle=initial_time, angle_type = graphs_type)
        _plot_component(axs[0], x_data, data.H_t[i_list], 'Перемещение толкателя, мм', xl)
        _plot_component(axs[1], x_data, data.V_t[i_list], 'Скорость, $мм/с$', xl)
        _plot_component(axs[2], x_data, data.A_t[i_list], 'Ускорение, $мм/с^2$', xl)
        _plot_component(axs[3], x_data, data.D_t[i_list], 'Рывок, $мм/с^3$', xl)
        _plot_component(axs[4], x_data, data.K_t[i_list], 'Кракен, $мм/с^4$', xl)

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

def display_profile_multiplicity(data: List[ProfileData], initial_angle: float | int = 0):
    plt.figure(figsize=(6, 6))
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

    plt.xlabel('X, мм')
    plt.ylabel('Y, мм')
    margin = 0.01
    plt.xlim(np.min(X) - margin, np.max(X) + margin)
    plt.ylim(np.min(Y) - margin, np.max(Y) + margin)
    plt.grid(True)
    plt.title("Профили кулачка")
    plt.axis('equal')
    plt.show()

def display_all(kulachok: Kulachok, initial_angle: Union[float, int, str] = 0, graphs_type: Literal['t', 'rad', 'degree'] = 'degree'):
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
        # Заголовки для Y осей
        ylabels = ['Перемещение, мм', 'Скорость', 'Ускорение', 'Рывок', 'Кракен']
    else:
        data = kulachok.kulachok_data
        base_title = "Кинематика кулачка"
        ylabels = ['Координата, мм', 'Скорость', 'Ускорение', 'Рывок', 'Кракен']

    # 3. Подготовка данных в зависимости от типа оси X (t, rad, degree)
    if graphs_type == 'degree':
        x_data = data.fi_list_degree
        xlabel = r'$\phi$, град'
        # Единицы измерения для подписей
        units = ['', r'$мм/град$', r'$мм/град^2$', r'$мм/град^3$', r'$мм/град^4$']
        i_list = set_initial_angle(x_data, initial_angle=target_angle, angle_type=graphs_type)
        # Данные Y
        y_datasets = [data.H_degree, data.V_degree, data.A_degree, data.D_degree, data.K_degree]  # H_rad обычно в мм, ок

    elif graphs_type == 'rad':
        x_data = data.fi_list_rad
        xlabel = r'$\phi$, рад'
        units = ['', r'$мм/рад$', r'$мм/рад^2$', r'$мм/рад^3$', r'$мм/рад^4$']
        initial_angle_rad = target_angle * np.pi / 180
        i_list = set_initial_angle(x_data, initial_angle=initial_angle_rad, angle_type=graphs_type)
        y_datasets = [data.H_rad, data.V_rad, data.A_rad, data.D_rad, data.K_rad]

    elif graphs_type == 't':
        x_data = data.t_list
        xlabel = "t, c"
        units = ['', r'$мм/с$', r'$мм/с^2$', r'$мм/с^3$', r'$мм/с^4$']
        initial_time = target_angle * np.pi / 180 / data.omega_rad if target_angle != 0 else 0
        i_list = set_initial_angle(x_data, initial_angle=initial_time, angle_type=graphs_type)
        y_datasets = [data.H_t, data.V_t, data.A_t, data.D_t, data.K_t]

    # --- СОЗДАНИЕ ГРАФИЧЕСКОГО ОКНА ---
    fig = plt.figure(figsize=(16, 10))
    # Сетка: 5 строк, 3 колонки.
    # Графики занимают первые 2 колонки (слева), Профиль занимает 3-ю колонку (справа)
    gs = gridspec.GridSpec(5, 3, figure=fig)

    # === ЛЕВАЯ ЧАСТЬ: Кинематика ===
    ax_kinematics = []
    for i in range(5):
        ax = fig.add_subplot(gs[i, 0:2])  # Занимаем 2 колонки из 3
        # Формируем полный заголовок оси Y с единицами измерения
        full_ylabel = f"{ylabels[i]} {units[i]}" if i > 0 else ylabels[i]
        _plot_component(ax, x_data, y_datasets[i][i_list], full_ylabel, xlabel)

        # Убираем подписи оси X у всех графиков кроме нижнего для чистоты
        if i < 4:
            ax.set_xlabel("")
            ax.tick_params(labelbottom=False)
        ax_kinematics.append(ax)

    ax_kinematics[0].set_title(base_title, fontsize=14)

    # === ПРАВАЯ ЧАСТЬ: Профиль ===
    ax_profile = fig.add_subplot(gs[:, 2])  # Занимаем последнюю колонку, все строки

    profile_data = kulachok.profile_data

    # Сортировка для профиля всегда через set_initial_angle (она внутри понимает градусы/радианы по умолчанию)
    i_list_prof = set_initial_angle(profile_data.fi_list, initial_angle=target_angle)
    X = np.asarray(profile_data.X)
    Y = np.asarray(profile_data.Y)
    X_sorted = X[i_list_prof]
    Y_sorted = Y[i_list_prof]
    X_plot = np.append(X_sorted, X_sorted[0])
    Y_plot = np.append(Y_sorted, Y_sorted[0])

    ax_profile.plot(X_plot, Y_plot, linewidth=2)
    ax_profile.scatter([0], [0], color='black', marker='x', s=100, label='Центр вращения')

    # Оформление профиля
    ax_profile.set_xlabel('X, мм')
    ax_profile.set_ylabel('Y, мм')
    ax_profile.set_title("Профиль кулачка", fontsize=14)
    ax_profile.axis('equal')  # Важно: сохраняем пропорции
    ax_profile.grid(True)
    ax_profile.legend()

    # Добавляем отступы
    margin = 2  # мм
    ax_profile.set_xlim(np.min(X) - margin, np.max(X) + margin)
    ax_profile.set_ylim(np.min(Y) - margin, np.max(Y) + margin)

    plt.tight_layout()
    plt.show()