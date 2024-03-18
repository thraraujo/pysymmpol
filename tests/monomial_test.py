import os 
import sys
sys.path.insert(0, os.path.abspath('..'))

import unittest
import sympy as sp
import pysymmpol as sy
import pysymmpol.utils.tools as ut

class TestMonomial(unittest.TestCase):


    def test_partition21_2coord(self) -> None:
        x = ut.create_x_coord(2)
        xx = dict(enumerate(x,1))
        part = (2,1)
        yg = sy.YoungDiagram(part)
        mon_symm_pol = sy.MonomialPolynomial(yg) 
        self.assertEqual(mon_symm_pol.explicit(x).expand(), (xx[1]**2 * xx[2] + xx[1]* xx[2]**2).expand())


    def test_partition21_1coord(self) -> None:
        x = ut.create_x_coord(1)
        xx = dict(enumerate(x,1))
        part = (2,1)
        yg = sy.YoungDiagram(part)
        mon_symm_pol = sy.MonomialPolynomial(yg) 
        self.assertEqual(mon_symm_pol.explicit(x), 0)


    def test_partition21_3coord(self) -> None:
        x = ut.create_x_coord(3)
        xx = dict(enumerate(x,1))
        part = (2,1)
        yg = sy.YoungDiagram(part)
        mon_symm_pol = sy.MonomialPolynomial(yg) 
        self.assertEqual(mon_symm_pol.explicit(x).expand(), (xx[1]**2 * xx[2] + xx[1]**2 * xx[3] + xx[1]* xx[2]**2 + xx[3]* xx[2]**2 +  xx[3]** 2 * xx[1] +  xx[3]** 2 * xx[2] ).expand())
        

    def test_partition3311_3coord(self) -> None:
        x = ut.create_x_coord(3)
        xx = dict(enumerate(x,1))
        part = (3,3,1,1)
        yg = sy.YoungDiagram(part)
        mon_symm_pol = sy.MonomialPolynomial(yg) 
        self.assertEqual(mon_symm_pol.explicit(x), 0)


    def test_partition3311_4coord(self) -> None:
        x = ut.create_x_coord(4)
        xx = dict(enumerate(x,1))
        part = (3,3,1,1)
        yg = sy.YoungDiagram(part)
        mon_symm_pol = sy.MonomialPolynomial(yg) 
        self.assertEqual((mon_symm_pol.explicit(x) -
                          ( xx[1]**3* xx[2]**3 * xx[3] * xx[4] + xx[1]**3 * xx[2]* xx[3] ** 3 * xx[4] + xx[1]**3 * xx[2] * xx[3] * xx[4]**3
                            + xx[1] * xx[2]**3 * xx[3]**3 * xx[4] + xx[1] * xx[2]**3 * xx[3] * xx[4]**3 + xx[1] * xx[2] * xx[3]**3 * xx[4]**3)).simplify(), 0)


if __name__ == '__main__':
    unittest.main()
