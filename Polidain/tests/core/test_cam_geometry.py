from collections import namedtuple
from unittest.mock import patch, MagicMock
import numpy as np
import pytest

from core.cam_geometry.config import default_kulachok_config
from core.cam_geometry.logic import Kulachok
from core.profiling_methods.polidain import PolidainCalculator
from core.profiling_methods.polidain.config import default_polidain_config


CamParams = namedtuple('CamParams', ['phi_0', 'phi_1', 'phi_2', 'phi_3', 'phi_4', 'phi_5', 'h', 'z', 'r0'])

@pytest.fixture
def params(kulachok):
    """Фикстура, извлекающая все ключевые параметры кулачка в один объект."""
    c = kulachok.config
    return CamParams(c.phi_0,
                     c.phi_1,
                     c.phi_2,
                     c.phi_3,
                     c.phi_4,
                     c.phi_5,
                     c.h,
                     c.z,
                     c.r0)

@pytest.fixture
def kulachok():
    """Фикстура, которая создает экземпляр кулачка перед каждым тестом."""
    config = default_kulachok_config
    profile_method_calculator = PolidainCalculator(default_polidain_config)
    return Kulachok(config, profile_method_calculator)

def test_fun_h(kulachok, params):
    assert kulachok.fun_h(params.phi_0) == pytest.approx(params.r0 - params.z, abs=1e-9)
    assert kulachok.fun_h(params.phi_1) == pytest.approx(params.r0, abs=1e-9)
    assert kulachok.fun_h(params.phi_2) == pytest.approx(params.r0 + params.h, abs=1e-9)
    assert kulachok.fun_h(params.phi_3) == pytest.approx(params.r0 + params.h, abs=1e-9)
    assert kulachok.fun_h(params.phi_4) == pytest.approx(params.r0, abs=1e-9)
    assert kulachok.fun_h(params.phi_5) == pytest.approx(params.r0 - params.z, abs=1e-9)

def test_fun_h_tolkatel(kulachok, params):
    assert kulachok.fun_h_tolkatel(params.phi_0) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_h_tolkatel(params.phi_1) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_h_tolkatel(params.phi_2) == pytest.approx(params.h, abs=1e-9)
    assert kulachok.fun_h_tolkatel(params.phi_3) == pytest.approx(params.h, abs=1e-9)
    assert kulachok.fun_h_tolkatel(params.phi_4) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_h_tolkatel(params.phi_5) == pytest.approx(0, abs=1e-9)

def test_fun_v(kulachok, params):
    """Тестирование скорости толкателя (первая производная перемещения)."""
    # На границах всех кинематических участков скорость должна падать до нуля,
    # чтобы обеспечить безударный переход между фазами (выстой -> движение -> выстой).
    assert kulachok.fun_v(params.phi_0) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_v(params.phi_1) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_v(params.phi_2) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_v(params.phi_3) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_v(params.phi_4) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_v(params.phi_5) == pytest.approx(0, abs=1e-9)

    # Дополнительная проверка: в середине фазы подъема (между phi_1 и phi_2)
    # скорость должна быть больше нуля
    phi_mid_rise = (params.phi_1 + params.phi_2) / 2
    assert kulachok.fun_v(phi_mid_rise) > 0

    # В середине фазы опускания (между phi_3 и phi_4) скорость должна быть меньше нуля
    phi_mid_fall = (params.phi_3 + params.phi_4) / 2
    assert kulachok.fun_v(phi_mid_fall) < 0

def test_fun_v_tolkatel(kulachok, params):
    """Тестирование скорости толкателя (первая производная перемещения)."""
    # На границах всех кинематических участков скорость должна падать до нуля,
    # чтобы обеспечить безударный переход между фазами (выстой -> движение -> выстой).
    assert kulachok.fun_v_tolkatel(params.phi_0) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_v_tolkatel(params.phi_1) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_v_tolkatel(params.phi_2) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_v_tolkatel(params.phi_3) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_v_tolkatel(params.phi_4) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_v_tolkatel(params.phi_5) == pytest.approx(0, abs=1e-9)

    # Дополнительная проверка: в середине фазы подъема (между phi_1 и phi_2)
    # скорость должна быть больше нуля
    phi_mid_rise = (params.phi_1 + params.phi_2) / 2
    assert kulachok.fun_v_tolkatel(phi_mid_rise) > 0

    # В середине фазы опускания (между phi_3 и phi_4) скорость должна быть меньше нуля
    phi_mid_fall = (params.phi_3 + params.phi_4) / 2
    assert kulachok.fun_v_tolkatel(phi_mid_fall) < 0

def test_fun_a(kulachok, params):
    """Тестирование ускорения толкателя (вторая производная перемещения)."""
    assert kulachok.fun_a(params.phi_0) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_a(params.phi_1) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_a(params.phi_2) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_a(params.phi_3) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_a(params.phi_4) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_a(params.phi_5) == pytest.approx(0, abs=1e-9)

def test_fun_a_tolkatel(kulachok, params):
    """Тестирование ускорения толкателя (вторая производная перемещения)."""
    assert kulachok.fun_a_tolkatel(params.phi_0) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_a_tolkatel(params.phi_1) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_a_tolkatel(params.phi_2) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_a_tolkatel(params.phi_3) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_a_tolkatel(params.phi_4) == pytest.approx(0, abs=1e-9)
    assert kulachok.fun_a_tolkatel(params.phi_5) == pytest.approx(0, abs=1e-9)

# ==========================================
# Тесты для установки данных и флагов
# ==========================================

def test_set_data_flags(kulachok):
    """Проверяем, что методы генерации данных корректно устанавливают флаги."""
    kulachok.set_kulachok_data(N=10)
    assert kulachok.kulachok_solve_flag is True
    assert hasattr(kulachok, 'kulachok_data')

    kulachok.set_tolkatel_data(N=10)
    assert kulachok.tolkatel_solve_flag is True
    assert hasattr(kulachok, 'tolkatel_data')

    kulachok.set_profile_data(N=10)
    assert kulachok.profile_solve_flag is True
    assert kulachok.tolkatel_solve_type == "thin"
    assert hasattr(kulachok, 'profile_data')


# ==========================================
# Тесты для роутинга метода solve()
# ==========================================

def test_solve_routing(kulachok):
    """Проверяем, что функция solve правильно распределяет вызовы методов."""
    # С помощью patch.object мы "перехватываем" вызовы методов, чтобы они
    # не выполняли реальные расчеты, а просто фиксировали факт своего вызова.
    with patch.object(kulachok, 'set_profile_data') as mock_thin, \
         patch.object(kulachok, 'set_profile_flat') as mock_flat, \
         patch.object(kulachok, 'set_profile_roller') as mock_roller, \
         patch.object(kulachok, 'set_tolkatel_data'), \
         patch.object(kulachok, 'set_kulachok_data'):

        # Проверяем ветку 'thin'
        kulachok.solve(kulachok_type='thin', N=10)
        mock_thin.assert_called_once_with(N=10)
        mock_flat.assert_not_called()
        mock_roller.assert_not_called()

        mock_thin.reset_mock()

        # Проверяем ветку 'flat'
        kulachok.solve(kulachok_type='flat', N=10)
        mock_flat.assert_called_once()
        mock_thin.assert_not_called()

        mock_flat.reset_mock()

        # Проверяем ветку 'roller'
        kulachok.solve(kulachok_type='roller', N=10)
        mock_roller.assert_called_once()


# ==========================================
# Тесты для проверок (Checkers) и исключений
# ==========================================

def test_profile_flat_check_preliminary_error(kulachok):
    """Проверка ошибки: предварительные вычисления не выполнены."""
    kulachok.kulachok_solve_flag = False
    kulachok.tolkatel_solve_flag = False

    with pytest.raises(Exception, match="Не были проведены предварительные вычисления"):
        kulachok.profile_flat_check()

def test_profile_flat_check_diameter_error(kulachok):
    """Проверка валидации диаметра плоского толкателя (PusherDiameterError)."""
    # Имитируем, что предварительные расчеты выполнены
    kulachok.kulachok_solve_flag = True
    kulachok.tolkatel_solve_flag = True

    # Создаем фиктивные данные, чтобы спровоцировать математическую ошибку
    mock_data = MagicMock()
    kulachok.config.D_t = 0.01  # Диаметр 10 мм
    mock_data.V_t = np.array([1000.0])  # Очень высокая скорость (вызовет заклинивание)
    kulachok.tolkatel_data = mock_data

    # Замени Exception на PusherDiameterError
    with pytest.raises(Exception, match="Недостаточный диаметр толкателя"):
        kulachok.profile_flat_check()

def test_profile_roller_check_smoothness_error(kulachok):
    """Проверка гладкости роликового профиля (ProfileSmoothnessError)."""
    # Искусственно задаем массивы с огромным отрицательным ускорением,
    # что гарантированно приведет к пересечению профиля (отрицательному радиусу кривизны)
    mock_data = MagicMock()
    mock_data.H_rad = np.array([0.0])
    mock_data.V_rad = np.array([0.0])
    mock_data.A_rad = np.array([1000.0])
    kulachok.tolkatel_data = mock_data

    kulachok.config.D = 0.01
    kulachok.config.R_r = 0.05

    # Замени Exception на ProfileSmoothnessError
    with pytest.raises(Exception, match="Негладкий профиль"):
        kulachok.profile_roller_check()