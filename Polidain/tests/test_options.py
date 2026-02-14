import pytest
from unittest.mock import MagicMock, patch
from core.optimization.options import CamOptimizationOptions, calculate_cam_optimization
from core.optimization.config import OptimizeConfig, GibridOptimizationConfig
from vizualization.plotter.options import PlotterOptions, resolve_plotter_options
from vizualization.rotate_animation.options import RotateAnimationOptions, resolve_rotate_animation_options

def test_cam_optimization_options():
    opt_config = OptimizeConfig(D=0.03, h=0.01)
    gibrid_config = GibridOptimizationConfig(m=[3], d=[4])
    options = CamOptimizationOptions(
        optimiz_config=opt_config,
        gibrid_optimization_config=gibrid_config,
        R_func=lambda x: x
    )
    assert options.N == 1000

@patch('core.optimization.options.gibrid_optimization')
def test_calculate_cam_optimization(mock_gibrid):
    opt_config = OptimizeConfig(D=0.03, h=0.01)
    gibrid_config = GibridOptimizationConfig(m=[3], d=[4])
    options = CamOptimizationOptions(
        optimiz_config=opt_config,
        gibrid_optimization_config=gibrid_config,
        R_func=lambda x: x
    )
    calculate_cam_optimization(options)
    assert mock_gibrid.called

def test_plotter_options_validation():
    opts = PlotterOptions(profile_and_graphs_together_flag=True)
    # Validator should set graphs_profile_flag and graphs_kulachok_flag
    assert opts.graphs_profile_flag
    assert opts.graphs_kulachok_flag

@patch('vizualization.plotter.options.display_dashboard')
@patch('vizualization.plotter.options.display_graphs_tolkatel')
@patch('vizualization.plotter.options.display_graphs_kulachok')
@patch('vizualization.plotter.options.display_profile')
def test_resolve_plotter_options(mock_profile, mock_kulachok, mock_tolkatel, mock_dashboard, kulachok):
    kulachok.solve(kulachok_type='thin', N=100)

    # Test separate flags
    opts = PlotterOptions(graphs_tolkatel_flag=True, graphs_profile_flag=True)
    resolve_plotter_options(opts, kulachok, initial_angle=0)
    assert mock_tolkatel.called
    assert mock_profile.called
    assert not mock_dashboard.called

    # Reset mocks
    mock_tolkatel.reset_mock()
    mock_profile.reset_mock()

    # Test together flag
    opts_together = PlotterOptions(profile_and_graphs_together_flag=True)
    resolve_plotter_options(opts_together, kulachok, initial_angle=0)
    assert mock_dashboard.called

def test_rotate_animation_options():
    opts = RotateAnimationOptions(display_animation_flag=True)
    assert opts.animation_intarval == 50

@patch('vizualization.rotate_animation.options.display_animation')
@patch('vizualization.rotate_animation.options.display_dashboard_animation')
def test_resolve_rotate_animation_options(mock_dash, mock_anim, kulachok):
    kulachok.solve(kulachok_type='thin', N=100)

    opts = RotateAnimationOptions(display_animation_flag=True, animation_profile_and_graphs_together_flag=True)
    resolve_rotate_animation_options(opts, kulachok)

    assert mock_anim.called
    assert mock_dash.called
