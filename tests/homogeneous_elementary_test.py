import os 
import sys
sys.path.insert(0, os.path.abspath('..'))

import unittest
import sympy as sp
import pysymmpol as sy
import pysymmpol.utils.tools as ut

class TestHomogeneousElementary(unittest.TestCase):


    def test_initialization(self) -> None:
        t = sp.symbols('t1:2')
        with self.assertRaises(TypeError):
            h3 = sy.HomogeneousPolynomial(3)
            h3.explicit(t)
        with self.assertRaises(TypeError):
            e4 = sy.ElementaryPolynomial(4)
            e4.explicit(t)


    def test_homogeneous_hp6(self) -> None:
        t = sp.symbols('t1:10')
        hp6 = t[0]**6 / 720 + t[0]**4 * t[1]/ 24 + t[0]**3 * t[2 ]/ 6 + t[0]**2 * t[1]**2/ 4 + t[0]**2 * t[3] /2 +t[0]*t[1]*t[2] + t[0]*t[4] + t[1]**3/ 6 + t[2]**2/ 2 + t[1]*t[3] + t[5]
        hp = sy.HomogeneousPolynomial(6).explicit(t)
        s = (hp6 - hp).simplify()
        self.assertEqual(s, 0)


    def test_elementary_ep5(self) -> None:
        t = sp.symbols('t1:10')
        ep5 = t[0]**5 / 120 - t[0]**3 * t[1] / 6 + t[0]**2 * t[2]/2 + t[0]*t[1]**2 / 2 - t[0] * t[3] - t[1]* t[2] + t[4]
        ep = sy.ElementaryPolynomial(5).explicit(t)
        s = (ep5 - ep).simplify()
        self.assertEqual(s, 0)


    def test_homogeneous_coordinates(self) -> None:
        n = 15
        homogeneous = [sy.HomogeneousPolynomial(i) for i in range(n)]
        t = sp.symbols(f't1:{n}')
        t_dict = ut.create_miwa(n)
        for a in homogeneous:
            self.assertEqual(a.explicit(t) - a.explicit(t_dict), 0)


    def test_elementary_coordinates(self) -> None:
        n = 15
        elementary = [sy.ElementaryPolynomial(i) for i in range(n)]
        t = sp.symbols(f't1:{n}')
        t_dict = ut.create_miwa(n)
        for a in elementary:
            self.assertEqual(a.explicit(t) - a.explicit(t_dict), 0)


if __name__ == '__main__':
    unittest.main()
        
