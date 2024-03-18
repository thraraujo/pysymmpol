import os 
import sys
sys.path.insert(0, os.path.abspath('..'))

import unittest
import sympy as sp
import pysymmpol as sy
import pysymmpol.utils.tools as ut

class TestHallLittlewood(unittest.TestCase):


    def test_hl_2_2coord(self) -> None:
        part = (2,)
        yg = sy.YoungDiagram(part)

        Q = sp.Symbol('Q')

        x = ut.create_x_coord(2)
        t = ut.tx_power_sum(yg.boxes, len(x))

        mp = sy.MonomialPolynomial(yg).explicit(x) 
        sc = sy.SchurPolynomial(yg).explicit(t)
        hl0 = sy.HallLittlewoodPolynomial(yg).explicit(x, 0)
        hl1 = sy.HallLittlewoodPolynomial(yg).explicit(x, 1)

        self.assertEqual((hl0 - sc).simplify(), 0)
        self.assertEqual((hl1 - mp).simplify(), 0)


    def test_hl_3_2coord(self) -> None:
        part = (3,)
        yg = sy.YoungDiagram(part)

        Q = sp.Symbol('Q')

        x = ut.create_x_coord(2)
        t = ut.tx_power_sum(yg.boxes, len(x))

        mp = sy.MonomialPolynomial(yg).explicit(x) 
        sc = sy.SchurPolynomial(yg).explicit(t)
        hl0 = sy.HallLittlewoodPolynomial(yg).explicit(x, 0)
        hl1 = sy.HallLittlewoodPolynomial(yg).explicit(x, 1)

        self.assertEqual((hl0 - sc).simplify(), 0)
        self.assertEqual((hl1 - mp).simplify(), 0)


    def test_hl_31_2coord(self) -> None:
        part = (3,1)
        yg = sy.YoungDiagram(part)

        Q = sp.Symbol('Q')

        x = ut.create_x_coord(2)
        t = ut.tx_power_sum(yg.boxes, len(x))

        mp = sy.MonomialPolynomial(yg).explicit(x) 
        sc = sy.SchurPolynomial(yg).explicit(t)
        hl0 = sy.HallLittlewoodPolynomial(yg).explicit(x, 0)
        hl1 = sy.HallLittlewoodPolynomial(yg).explicit(x, 1)

        self.assertEqual((hl0 - sc).simplify(), 0)
        self.assertEqual((hl1 - mp).simplify(), 0)


    def test_hl_2_3coord(self) -> None:
        part = (2,)
        yg = sy.YoungDiagram(part)

        Q = sp.Symbol('Q')

        x = ut.create_x_coord(3)
        t = ut.tx_power_sum(yg.boxes, len(x))

        mp = sy.MonomialPolynomial(yg).explicit(x) 
        sc = sy.SchurPolynomial(yg).explicit(t)
        hl0 = sy.HallLittlewoodPolynomial(yg).explicit(x, 0)
        hl1 = sy.HallLittlewoodPolynomial(yg).explicit(x, 1)

        self.assertEqual((hl0 - sc).simplify(), 0)
        self.assertEqual((hl1 - mp).simplify(), 0)


    def test_hl_21_3coord(self) -> None:
        part = (2,1)
        yg = sy.YoungDiagram(part)

        Q = sp.Symbol('Q')

        x = ut.create_x_coord(3)
        t = ut.tx_power_sum(yg.boxes, len(x))

        mp = sy.MonomialPolynomial(yg).explicit(x) 
        sc = sy.SchurPolynomial(yg).explicit(t)
        hl0 = sy.HallLittlewoodPolynomial(yg).explicit(x, 0)
        hl1 = sy.HallLittlewoodPolynomial(yg).explicit(x, 1)

        self.assertEqual((hl0 - sc).simplify(), 0)
        self.assertEqual((hl1 - mp).simplify(), 0)


    def test_hl_22_3coord(self) -> None:
        part = (2,2)
        yg = sy.YoungDiagram(part)

        Q = sp.Symbol('Q')

        x = ut.create_x_coord(3)
        t = ut.tx_power_sum(yg.boxes, len(x))

        mp = sy.MonomialPolynomial(yg).explicit(x) 
        sc = sy.SchurPolynomial(yg).explicit(t)
        hl0 = sy.HallLittlewoodPolynomial(yg).explicit(x, 0)
        hl1 = sy.HallLittlewoodPolynomial(yg).explicit(x, 1)

        self.assertEqual((hl0 - sc).simplify(), 0)
        self.assertEqual((hl1 - mp).simplify(), 0)


    def test_hl_42_2coord(self) -> None:
        part = (4,2)
        yg = sy.YoungDiagram(part)

        Q = sp.Symbol('Q')

        x = ut.create_x_coord(2)

        hl = sy.HallLittlewoodPolynomial(yg).explicit(x, Q)

        self.assertEqual((hl -  ( x[0]**4 * x[1]**2 + x[0]**2 * x[1]**4 + (1 - Q) * x[0]**3 * x[1]**3 )).simplify(), 0)


if __name__ == '__main__':
    unittest.main()
        
