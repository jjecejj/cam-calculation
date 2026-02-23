from core.profiling_methods.base.config import MethodConfig
from abc import ABC, abstractmethod


class BaseCalculator(ABC):
    """
    Абстрактный базовый класс для реализации методов расчета профилей кулачков.

    Определяет интерфейс для вычисления кинематических характеристик (перемещение, скорость,
    ускорение и т.д.) в зависимости от угла поворота.
    """

    def __init__(self, config: MethodConfig):
        """
        Инициализирует калькулятор с заданной конфигурацией метода.

        Args:
            config: Объект конфигурации метода расчета (наследник MethodConfig).
        """
        self.config: MethodConfig = config

    @abstractmethod
    def h_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """
        Вычисляет безразмерное перемещение (аналог S) для заданного угла.

        Args:
            fi: Текущий угол поворота (рад).
            fi_1: Конечный угол участка (рад).
            fi_0: Начальный угол участка (рад).
            h_kn_max: Максимальное перемещение на данном участке.
            segment_number: Номер текущего участка профиля.

        Returns:
            Значение перемещения.
        """
        pass

    @abstractmethod
    def v_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """
        Вычисляет аналог скорости (первая производная перемещения по углу) для заданного угла.

        Args:
            fi: Текущий угол поворота (рад).
            fi_1: Конечный угол участка (рад).
            fi_0: Начальный угол участка (рад).
            h_kn_max: Максимальное перемещение на данном участке.
            segment_number: Номер текущего участка профиля.

        Returns:
            Значение аналога скорости.
        """
        pass

    @abstractmethod
    def a_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """
        Вычисляет аналог ускорения (вторая производная перемещения по углу) для заданного угла.

        Args:
            fi: Текущий угол поворота (рад).
            fi_1: Конечный угол участка (рад).
            fi_0: Начальный угол участка (рад).
            h_kn_max: Максимальное перемещение на данном участке.
            segment_number: Номер текущего участка профиля.

        Returns:
            Значение аналога ускорения.
        """
        pass

    @abstractmethod
    def d_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """
        Вычисляет аналог рывка (третья производная перемещения по углу) для заданного угла.

        Args:
            fi: Текущий угол поворота (рад).
            fi_1: Конечный угол участка (рад).
            fi_0: Начальный угол участка (рад).
            h_kn_max: Максимальное перемещение на данном участке.
            segment_number: Номер текущего участка профиля.

        Returns:
            Значение аналога рывка.
        """
        pass

    @abstractmethod
    def k_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        """
        Вычисляет аналог производной рывка (четвертая производная перемещения по углу) для заданного угла.

        Args:
            fi: Текущий угол поворота (рад).
            fi_1: Конечный угол участка (рад).
            fi_0: Начальный угол участка (рад).
            h_kn_max: Максимальное перемещение на данном участке.
            segment_number: Номер текущего участка профиля.

        Returns:
            Значение четвертой производной.
        """
        pass
