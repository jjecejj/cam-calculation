import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from typing import Literal, Union

from numpy import ndarray

from vizualization.plotter.config import PlotConfig
from core.cam_geometry import Kulachok
from core.cam_geometry.schemas import ProfileData, GraphData

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
    fig, axs = plt.subplots(5, 1, figsize=(8, 20))
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
    fig, axs = plt.subplots(5, 1, figsize=(8, 20))
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
