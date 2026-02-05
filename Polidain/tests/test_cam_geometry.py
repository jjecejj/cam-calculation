import unittest
import numpy as np
import sys
import os

# Ensure the parent directory is in the path to import core modules
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from core.schemas import PolidainConfig, PolidainData, ProfilData
from core.cam_geometry import Kulachok_polidain

class TestCamGeometry(unittest.TestCase):
    def setUp(self):
        # Configuration similar to polidainConfig_default_thin but with adjusted parameters
        # to ensure validity for flat and roller followers during testing
        # Use parameters from polidainConfig_default_flat which are likely optimized
        self.config = PolidainConfig(
            N_k=1000,
            D=30.0 * 1e-3,
            D_t=50.0 * 1e-3, # Slightly increased D_t
            h=2.0 * 1e-3,
            z=0.25 * 1e-3,
            f_pod=80.0 / 180 * np.pi,
            f_v=5.0 / 180 * np.pi,
            f_op=75.0 / 180 * np.pi,
            f_z=25 / 180 * np.pi,
            m=3,
            d=4,
            k_1=6,
            k_2=6,
            k_3=6,
            k_4=6,
            R_r=5 * 1e-3
        )
        self.cam = Kulachok_polidain(self.config)

    def test_solve_thin(self):
        self.cam.solve(N=100, kulachok_type='thin')
        self.assertIsInstance(self.cam.kulachok_data, PolidainData)
        self.assertIsInstance(self.cam.tolkatel_data, PolidainData)
        self.assertIsInstance(self.cam.profil_data, ProfilData)
        self.assertEqual(len(self.cam.kulachok_data.H_rad), 100)
        self.assertEqual(self.cam.solve_type, 'thin')

        # Check basic properties
        self.assertEqual(self.cam.kulachok_data.omega_rad, self.config.omega)

    def test_solve_flat(self):
        self.cam.solve(N=100, kulachok_type='flat')
        self.assertIsInstance(self.cam.profil_data, ProfilData)
        self.assertEqual(self.cam.solve_type, 'flat')
        self.assertEqual(len(self.cam.profil_data.X), 100)

    def test_solve_roller(self):
        self.cam.solve(N=100, kulachok_type='roller')
        self.assertIsInstance(self.cam.profil_data, ProfilData)
        self.assertEqual(self.cam.solve_type, 'roller')
        self.assertEqual(len(self.cam.profil_data.X), 100)

    def test_invalid_type(self):
        with self.assertRaises(ValueError):
            self.cam.solve(N=100, kulachok_type='invalid')

if __name__ == '__main__':
    unittest.main()
