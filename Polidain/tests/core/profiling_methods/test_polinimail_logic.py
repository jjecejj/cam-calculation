import pytest

from core.profiling_methods.polinmail.config import default_polinmail_config
from core.profiling_methods.polinmail.logic import calculate_poly_coefficients, PolinmailCalculator


@pytest.fixture
def calculator():
    """Фикстура, которая создает экземпляр калькулятора перед каждым тестом."""
    config = default_polinmail_config
    return PolinmailCalculator(config)

def test_segment_selection_success(calculator):
    """Проверяем, что существующий участок возвращает кортеж с коэффициентами."""
    m_list, coeffs = calculator.segment_selection(1)

    assert len(m_list) > 0
    assert len(coeffs) > 0


def test_segment_selection_error(calculator):
    """Проверяем, что при запросе несуществующего участка выбрасывается ValueError."""
    with pytest.raises(ValueError) as exc_info:
        calculator.segment_selection(5)


@pytest.mark.parametrize("segment", [1, 2, 3, 4])
def test_h_phi(calculator, segment):
    fi_0 = 0
    fi_1 = 1
    h_kn_max = 10

    # Проверка в начале участка (fi = 0)
    assert calculator.h_phi(fi_0, fi_1, fi_0, h_kn_max, segment) == pytest.approx(0)

    # Проверка в середине участка (fi = 0.5)
    h_mid = calculator.h_phi(0.5, fi_1, fi_0, h_kn_max, segment)
    assert 0 <= h_mid <= h_kn_max

    # Проверка в конце участка (fi = 1)
    assert calculator.h_phi(fi_1, fi_1, fi_0, h_kn_max, segment) == pytest.approx(h_kn_max)


@pytest.mark.parametrize("segment", [1, 2, 3, 4])
def test_v_phi(calculator, segment):
    fi_0 = 0
    fi_1 = 1
    h_kn_max = 10

    # В начале участка (fi = 0) скорость должна быть равна 0
    v_start = calculator.v_phi(fi_0, fi_1, fi_0, h_kn_max, segment)
    assert v_start == pytest.approx(0, abs=1e-9), f"Скорость в начале {segment} участка не равна 0"

    # В середине участка (fi = 0.5) кулачок движется, скорость не должна быть нулевой
    v_mid = calculator.v_phi(0.5, fi_1, fi_0, h_kn_max, segment)
    assert isinstance(v_mid, float)
    assert v_mid != 0, "В середине подъема скорость не должна быть равна нулю"

    # В конце участка (fi = 1) скорость снова падает до 0
    v_end = calculator.v_phi(fi_1, fi_1, fi_0, h_kn_max, segment)
    assert v_end == pytest.approx(0, abs=1e-9), f"Скорость в конце {segment} участка не равна 0"

@pytest.mark.parametrize("segment", [1, 2, 3, 4])
def test_a_phi(calculator, segment):
    fi_0 = 0
    fi_1 = 1
    h_kn_max = 10

    # Проверяем ускорение в начале участка
    a_start = calculator.a_phi(fi_0, fi_1, fi_0, h_kn_max, segment)
    assert a_start == pytest.approx(0, abs=1e-9)

    # Проверяем ускорение в середине участка
    a_mid = calculator.a_phi(0.5, fi_1, fi_0, h_kn_max, segment)
    assert isinstance(a_mid, float)

    # Проверяем ускорение в конце участка
    a_end = calculator.a_phi(fi_1, fi_1, fi_0, h_kn_max, segment)
    assert a_end == pytest.approx(0, abs=1e-9)