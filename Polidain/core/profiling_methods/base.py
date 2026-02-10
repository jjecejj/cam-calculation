from config.profiling_methods.base import MethodConfig
from abc import ABC, abstractmethod

class  BaseCalculator(ABC):
    def __init__(self, config: MethodConfig):
        self.config: MethodConfig = config

    @abstractmethod
    def h_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        pass

    @abstractmethod
    def v_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        pass

    @abstractmethod
    def a_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        pass

    @abstractmethod
    def d_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        pass

    @abstractmethod
    def k_phi(self, fi: float, fi_1: float, fi_0: float, h_kn_max: float, segment_number: int):
        pass