import numpy as np
import matplotlib.pyplot as plt
from typing import Union
from numpy import ndarray
from core.cam_geometry import set_polidain_data, set_profil_data


def set_config(config):
    """
    Применяет настройки конфигурации к глобальным параметрам matplotlib.

    Обновляет шрифты, размеры текста и настройки отображения математических формул
    в соответствии с переданным объектом конфигурации.

    Args:
        config: Объект конфигурации, содержащий атрибуты стиля (rc, font_serif и др.).
    """
    rc = config.rc
    plt.rcParams.update(rc)
    plt.rcParams["font.serif"] = config.font_serif
    plt.rcParams["font.family"] = config.font_family
    plt.rcParams['font.size'] = config.font_size
    plt.rcParams['mathtext.fontset'] = config.mathtext_fontset
    plt.rcParams['mathtext.rm'] = config.mathtext_rm
    plt.rcParams['mathtext.it'] = config.mathtext_it


def _plot_component(ax, x, y, ylabel):
    """
    Универсальная вспомогательная функция для построения графика с экстремумами.

    Рисует график, находит глобальный минимум и максимум, отмечает их точками.

    Args:
        ax (matplotlib.axes.Axes): Ось (subplot), на которой нужно построить график.
        x (list | ndarray): Данные для оси абсцисс (угол поворота).
        y (list | ndarray): Данные для оси ординат (значение кинематического параметра).
        ylabel (str): Подпись для оси ординат (например, 'Скорость, мм/с').
    """
    # Гарантируем, что работаем с numpy массивами
    x_arr = np.asarray(x)
    y_arr = np.asarray(y)

    # Основной график
    ax.plot(x_arr, y_arr)

    # Поиск индексов максимума и минимума
    idx_max = np.argmax(y_arr)
    idx_min = np.argmin(y_arr)

    # Точки экстремумов
    ax.scatter(x_arr[idx_max], y_arr[idx_max], color='r')
    ax.scatter(x_arr[idx_min], y_arr[idx_min], color='r')

    # Оформление
    ax.set_xlabel(r'$\phi$, град')
    ax.set_ylabel(ylabel)
    ax.grid(True)


def set_initial_angle(_fi_list: list | ndarray, initial_angle: float | int = 0.0):
    """
    Вычисляет индексы сортировки для фазового сдвига массива углов.

    Позволяет переупорядочить данные так, чтобы график начинался с заданного
    начального угла, с учетом цикличности (360 градусов).

    Args:
        _fi_list (list | ndarray): Исходный массив углов.
        initial_angle (float | int, optional): Угол начала отсчета. По умолчанию 0.0.

    Returns:
        ndarray: Массив индексов, который нужно использовать для сортировки данных.
    """
    __fi_list = np.asarray(_fi_list)
    arr = np.mod(__fi_list - initial_angle, 360)
    i_list = np.argsort(arr)
    return i_list


def calculate_optimal_angle(kulachok):
    """
    Внутренняя функция для автоматического определения оптимального угла сдвига.

    Рассчитывает угол симметрии на основе параметра phi_5 из конфигурации.

    Args:
        kulachok: Объект механизма.

    Returns:
        float: Рассчитанный угол в градусах.
    """
    return np.degrees(-(2 * np.pi - kulachok.config.phi_5) / 2)


def display_graphs_kulachok(data, initial_angle: float | int = 0):
    """
    Строит и отображает 5 кинематических графиков для кулачка.

    Визуализирует радиус-вектор, скорость, ускорение, рывок и "кракен".

    Args:
        data: Объект данных кулачка (атрибуты fi_list, H, V, A, D, K).
        initial_angle (float | int, optional): Угол сдвига фазы. По умолчанию 0.
    """
    fig, axs = plt.subplots(5, 1, figsize=(8, 20))
    fig.suptitle("Графики для кулачка", fontsize=16)

    i_list = set_initial_angle(data.fi_list, initial_angle=initial_angle)

    _plot_component(axs[0], data.fi_list, data.H[i_list], 'Координата точки контакта, мм')
    _plot_component(axs[1], data.fi_list, data.V[i_list], 'Скорость, $мм/с$')
    _plot_component(axs[2], data.fi_list, data.A[i_list], 'Ускорение, $мм/с^2$')
    _plot_component(axs[3], data.fi_list, data.D[i_list], 'Рывок, $мм/с^3$')
    _plot_component(axs[4], data.fi_list, data.K[i_list], 'Кракен, $мм/с^4$')

    plt.tight_layout()
    plt.show()


def display_graphs_tolkatel(data, initial_angle: float | int = 0):
    """
    Строит и отображает 5 кинематических графиков для толкателя.

    Визуализирует перемещение, скорость, ускорение, рывок и "кракен".

    Args:
        data: Объект данных толкателя (атрибуты fi_list, H, V, A, D, K).
        initial_angle (float | int, optional): Угол сдвига фазы. По умолчанию 0.
    """
    fig, axs = plt.subplots(5, 1, figsize=(8, 20))
    fig.suptitle("Графики для толкателя", fontsize=16)

    i_list = set_initial_angle(data.fi_list, initial_angle=initial_angle)

    _plot_component(axs[0], data.fi_list, data.H[i_list], 'Перемещение толкателя, мм')
    _plot_component(axs[1], data.fi_list, data.V[i_list], 'Скорость, $мм/с$')
    _plot_component(axs[2], data.fi_list, data.A[i_list], 'Ускорение, $мм/с^2$')
    _plot_component(axs[3], data.fi_list, data.D[i_list], 'Рывок, $мм/с^3$')
    _plot_component(axs[4], data.fi_list, data.K[i_list], 'Кракен, $мм/с^4$')

    plt.tight_layout()
    plt.show()


def display_profil(data, initial_angle: float | int = 0):
    """
    Отображает геометрический профиль кулачка в координатах X-Y.

    Строит замкнутый контур профиля, отмечает центр вращения.

    Args:
        data: Объект с координатами профиля (атрибуты X, Y и fi_list).
        initial_angle (float | int, optional): Угол сдвига фазы для отрисовки.
            По умолчанию 0.
    """
    plt.figure(figsize=(6, 6))

    i_list = set_initial_angle(data.fi_list, initial_angle=initial_angle)
    X = np.asarray(data.X)
    Y = np.asarray(data.Y)

    # Используем i_list для упорядочивания точек профиля при отрисовке линии
    plt.plot(X[i_list], Y[i_list])
    plt.scatter([0], [0], color='black', marker='x', label='Центр вращения')

    plt.xlabel('X, мм')
    plt.ylabel('Y, мм')

    # Отступы для красивого отображения
    margin = 0.01
    plt.xlim(np.min(X) - margin, np.max(X) + margin)
    plt.ylim(np.min(Y) - margin, np.max(Y) + margin)

    plt.grid(True)
    plt.title("Профиль кулачка")
    plt.axis('equal')
    plt.legend()
    plt.show()


def display_all(kulachok, initial_angle: Union[float, int, str] = 0):
    """
    Комплексная функция для отображения всех графиков механизма.

    Последовательно строит графики кинематики кулачка, толкателя и профиля.
    Поддерживает автоматический расчет угла сдвига фазы.

    Args:
        kulachok: Главный объект механизма, содержащий данные кулачка, толкателя и профиля.
        initial_angle (float | int | str, optional): Угол начала отсчета.
            - Если число: используется как фиксированный угол сдвига.
            - Если "auto": угол рассчитывается автоматически через __calculate_optimal_angle.
            - По умолчанию 0.

    Raises:
        ValueError: Если initial_angle имеет недопустимый тип.
    """
    if initial_angle == "auto":
        target_angle = calculate_optimal_angle(kulachok)
    elif isinstance(initial_angle, (int, float)):
        target_angle = initial_angle
    else:
        raise ValueError("initial_angle должен быть числом или 'auto'")

    display_graphs_kulachok(kulachok.kulachok_data, initial_angle=target_angle)
    display_graphs_tolkatel(kulachok.tolkatel_data, initial_angle=target_angle)
    display_profil(kulachok.profil_data, initial_angle=target_angle)

def display_graphs_comprasion(data_1, data_2, initial_angle: float | int = 0):
    if data_1.fi_list.shape[0] != data_2.fi_list.shape[0]:
        raise ValueError("Оба кулачка должны иметь одинаковую размерность массивов")

    fig, axs = plt.subplots(5, 1, figsize=(8, 20))
    fig.suptitle("Графики для толкателя", fontsize=16)

    i_list = set_initial_angle(data_1.fi_list, initial_angle=initial_angle)

    _plot_component(axs[0], data_1.fi_list, data_1.H[i_list], 'Радиус кулачка, мм')
    _plot_component(axs[0], data_2.fi_list, data_2.H[i_list], 'Радиус кулачка, мм')

    _plot_component(axs[1], data_1.fi_list, data_1.V[i_list], 'Скорость, $мм/с$')
    _plot_component(axs[1], data_2.fi_list, data_2.V[i_list], 'Скорость, $мм/с$')

    _plot_component(axs[2], data_1.fi_list, data_1.A[i_list], 'Ускорение, $мм/с^2$')
    _plot_component(axs[2], data_2.fi_list, data_2.A[i_list], 'Ускорение, $мм/с^2$')

    _plot_component(axs[3], data_1.fi_list, data_1.D[i_list], 'Рывок, $мм/с^3$')
    _plot_component(axs[3], data_2.fi_list, data_2.D[i_list], 'Рывок, $мм/с^3$')

    _plot_component(axs[4], data_1.fi_list, data_1.K[i_list], 'Кракен, $мм/с^4$')
    _plot_component(axs[4], data_2.fi_list, data_2.K[i_list], 'Кракен, $мм/с^4$')

    plt.tight_layout()
    plt.show()

def display_profil_comprasion(data_1, data_2, initial_angle: float | int = 0):
    plt.figure(figsize=(6, 6))

    i_list = set_initial_angle(data_1.fi_list, initial_angle=initial_angle)
    X_1 = np.asarray(data_1.X)
    Y_1 = np.asarray(data_1.Y)
    X_2 = np.asarray(data_2.X)
    Y_2 = np.asarray(data_2.Y)

    # Используем i_list для упорядочивания точек профиля при отрисовке линии
    plt.plot(X_1[i_list], Y_1[i_list])
    plt.plot(X_2[i_list], Y_2[i_list])
    plt.scatter([0], [0], color='black', marker='x', label='Центр вращения')

    plt.xlabel('X, мм')
    plt.ylabel('Y, мм')

    # Отступы для красивого отображения
    margin = 0.01
    plt.xlim(np.min(X_1) - margin, np.max(X_1) + margin)
    plt.ylim(np.min(Y_1) - margin, np.max(Y_1) + margin)
    plt.xlim(np.min(X_2) - margin, np.max(X_2) + margin)
    plt.ylim(np.min(Y_2) - margin, np.max(Y_2) + margin)

    plt.grid(True)
    plt.title("Профиль кулачка")
    plt.axis('equal')
    plt.legend()
    plt.show()

def display_all_comprasion(kulachok, fun_list: list, initial_angle: Union[float, int, str] = 0):
    if initial_angle == "auto":
        target_angle = calculate_optimal_angle(kulachok)
    elif isinstance(initial_angle, (int, float)):
        target_angle = initial_angle
    else:
        raise ValueError("initial_angle должен быть числом или 'auto'")

    display_graphs_comprasion(kulachok.kulachok_data, set_polidain_data(fun_list[0:5], N = kulachok.kulachok_data.fi_list.shape[0]),  initial_angle=target_angle)
    display_profil_comprasion(kulachok.profil_data, set_profil_data(fun_list[5:7], N = kulachok.profil_data.fi_list.shape[0]), initial_angle=target_angle)