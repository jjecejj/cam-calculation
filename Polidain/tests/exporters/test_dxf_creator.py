import pytest
from unittest.mock import MagicMock, patch
from exporters.dxf_creator.logic import ProfileExportData, create_profil, build_profile
from core.cam_geometry.schemas import ProfileData
import numpy as np

def test_profile_export_data_validation():
    # Test valid closed loop
    X = [0, 1, 1, 0, 0]
    Y = [0, 0, 1, 1, 0]
    data = ProfileExportData(X=X, Y=Y)
    assert data.X[0] == data.X[-1]

    # Test open loop -> should be closed automatically
    X_open = [0, 1, 1, 0]
    Y_open = [0, 0, 1, 1]
    data_open = ProfileExportData(X=X_open, Y=Y_open)
    # The validation logic says:
    # if self.X[0] != self.X[-1] or self.Y[0] != self.Y[-1]:
    #     self.X[0] = self.X[-1]
    #     self.Y[0] = self.Y[-1]

    # In X_open: X[0]=0, X[-1]=0. Match.
    # In Y_open: Y[0]=0, Y[-1]=1. Mismatch.
    # So logic executes.
    # X[0] becomes X[-1] (0).
    # Y[0] becomes Y[-1] (1).

    assert data_open.X[0] == 0
    assert data_open.Y[0] == 1

@patch('exporters.dxf_creator.logic.ezdxf.new')
def test_create_profil_spline(mock_new):
    mock_doc = MagicMock()
    mock_msp = MagicMock()
    mock_new.return_value = mock_doc
    mock_doc.modelspace.return_value = mock_msp

    X = [0, 1, 1, 0, 0]
    Y = [0, 0, 1, 1, 0]
    data = ProfileExportData(X=X, Y=Y)

    create_profil(data, profil_name="test", line_type="spline")

    assert mock_msp.add_spline.called
    assert mock_doc.saveas.called

@patch('exporters.dxf_creator.logic.ezdxf.new')
def test_create_profil_line(mock_new):
    mock_doc = MagicMock()
    mock_msp = MagicMock()
    mock_new.return_value = mock_doc
    mock_doc.modelspace.return_value = mock_msp

    X = [0, 1000, 1000, 0, 0]
    Y = [0, 0, 1000, 1000, 0]
    data = ProfileExportData(X=X, Y=Y)

    create_profil(data, profil_name="test", line_type="line")

    assert mock_msp.add_line.called
    # It should call add_line N-1 times (4 times)
    assert mock_msp.add_line.call_count == 4

@patch('exporters.dxf_creator.logic.create_profil')
def test_build_profile(mock_create_profil):
    # ProfileData requires min_length=10
    X = [0] * 10
    Y = [1] * 10
    fi = [0] * 10
    profile_data = ProfileData(X=X, Y=Y, fi_list=fi)
    build_profile(profile_data, profile_name="test")
    assert mock_create_profil.called
