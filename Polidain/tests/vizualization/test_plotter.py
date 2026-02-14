import pytest
import numpy as np
from unittest.mock import MagicMock, patch

from vizualization.plotter.logic import (
    display_graphs_kulachok,
    display_graphs_tolkatel,
    display_profile,
    display_all,
    display_dashboard
)
from vizualization.rotate_animation.logic import display_animation, rotate_profile_data

@pytest.fixture
def mock_plt():
    with patch('vizualization.plotter.logic.plt') as mock:
        yield mock

@pytest.fixture
def mock_plt_anim():
    with patch('vizualization.rotate_animation.logic.plt') as mock:
        yield mock

@pytest.fixture
def mock_animation():
    with patch('vizualization.rotate_animation.logic.FuncAnimation') as mock:
        yield mock

def test_rotate_profile_data():
    # Simple check
    # if angle = -pi/2.
    # angle_in_func = -(-pi/2) - pi/2 = pi/2 - pi/2 = 0.
    # x_new = X*cos(0) - Y*sin(0) = X
    # y_new = X*sin(0) + Y*cos(0) = Y
    x, y = rotate_profile_data(10.0, 5.0, -np.pi/2)
    assert np.isclose(x, 10.0)
    assert np.isclose(y, 5.0)

def test_display_graphs_kulachok(kulachok, mock_plt):
    # Setup mock
    fig_mock = MagicMock()
    axs_mock = MagicMock()
    # Allow indexing on axs
    axs_mock.__getitem__.return_value = MagicMock()
    mock_plt.subplots.return_value = (fig_mock, axs_mock)

    # Setup data
    kulachok.solve(kulachok_type='thin', N=100)

    # Run function
    display_graphs_kulachok(kulachok.kulachok_data, initial_angle=0)

    # Assertions
    assert mock_plt.subplots.called
    assert mock_plt.show.called

def test_display_graphs_tolkatel(kulachok, mock_plt):
    # Setup mock
    fig_mock = MagicMock()
    axs_mock = MagicMock()
    axs_mock.__getitem__.return_value = MagicMock()
    mock_plt.subplots.return_value = (fig_mock, axs_mock)

    kulachok.solve(kulachok_type='thin', N=100)
    display_graphs_tolkatel(kulachok.tolkatel_data, initial_angle=0)
    assert mock_plt.subplots.called
    assert mock_plt.show.called

def test_display_profile(kulachok, mock_plt):
    kulachok.solve(kulachok_type='thin', N=100)
    display_profile(kulachok.profile_data, initial_angle=0)
    assert mock_plt.figure.called or mock_plt.subplots.called # display_profile uses plt.figure()
    assert mock_plt.show.called

def test_display_all(kulachok, mock_plt):
    # Setup mock
    fig_mock = MagicMock()
    axs_mock = MagicMock()
    axs_mock.__getitem__.return_value = MagicMock()
    mock_plt.subplots.return_value = (fig_mock, axs_mock)

    kulachok.solve(kulachok_type='thin', N=100)
    display_all(kulachok, initial_angle=0)
    # This calls other display functions which call subplots/figure/show.
    assert mock_plt.subplots.called or mock_plt.figure.called
    assert mock_plt.show.called

def test_display_dashboard(kulachok, mock_plt):
    kulachok.solve(kulachok_type='thin', N=100)
    display_dashboard(kulachok, initial_angle=0)
    assert mock_plt.figure.called
    assert mock_plt.show.called

def test_display_animation(kulachok, mock_plt_anim, mock_animation):
    kulachok.solve(kulachok_type='thin', N=100)

    # Mock ax.plot to return list of lines
    fig_mock = MagicMock()
    ax_mock = MagicMock()
    mock_plt_anim.subplots.return_value = (fig_mock, ax_mock)

    line_mock = MagicMock()
    ax_mock.plot.return_value = [line_mock]

    display_animation(kulachok, interval=50, save_flag=False)

    assert mock_plt_anim.subplots.called
    assert mock_animation.called
    assert mock_plt_anim.show.called

def test_display_animation_save(kulachok, mock_plt_anim, mock_animation):
    kulachok.solve(kulachok_type='thin', N=100)

    fig_mock = MagicMock()
    ax_mock = MagicMock()
    mock_plt_anim.subplots.return_value = (fig_mock, ax_mock)
    line_mock = MagicMock()
    ax_mock.plot.return_value = [line_mock]

    anim_mock = MagicMock()
    mock_animation.return_value = anim_mock

    display_animation(kulachok, interval=50, save_flag=True, name_file="test")

    assert anim_mock.save.called
