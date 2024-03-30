import unittest
import numpy as np
from goph420_lab02.functions import root_newton_raphson


class TestNewtonRaphson(unittest.TestCase):

    def setUp(self):
        radius = 1.75 #[m]
        solid_density = 925 #[km/m^3]
        water_density = 1000 #[kg.m^3]
        self.h = 0.6
        self.f = lambda h: ((h ** 2) / (3 * radius)) - ((4 / 3) * ((radius ** 2) / h) * ((solid_density / water_density) - 1)) - h
        self.dfdx = lambda h : ((2 * h) / (3 * radius)) + (((4 / 3) * ( radius ** 2) / (h ** 2)) * ((solid_density / water_density) - 1)) - 1
    
    # Quiz 2 comparing values
    def test_value(self):
        excel_root = 0.5872119182
        excel_itr = 4
        excel_error = np.array([2.206525816E-02, 2.813744239E-4, 4.758997602E-08])
        root, itr, error = root_newton_raphson(self.h, self.f, self.dfdx)
        self.assertAlmostEqual(root, excel_root)
        self.assertAlmostEqual(itr, excel_itr)
        self.assertTrue(np.allclose(error, excel_error, atol = 1e-8))