import os 
import sys
sys.path.insert(0, os.path.abspath('..'))

import unittest
import sympy as sp
import pysymmpol as sy
import pysymmpol.utils.tools as ut

class TestSchurPolynomial(unittest.TestCase):
    '''
    We have two different implementations of the Schur polynomials.
    The first, faster, uses the determinant of the homogeneous polynomials.
    The second is defined via de characters of the symmetric group.
    In this test we calculate the equality of these polynomials for some partitions,
    including two edge cases: partitions = (7, 5, 3) and (7, 5, 3, 1)
    '''

    def test_schur_3(self) -> None:
        partition  = (3,)
        yg = sy.YoungDiagram(partition)
        t = ut.create_miwa(yg.boxes)
        sch = sy.SchurPolynomial(yg)
        self.assertEqual(sch.explicit(t), sch._schur_characters(t))


    def test_schur_2111(self) -> None:
        partition  = (2, 1, 1, 1)
        yg = sy.YoungDiagram(partition)
        t = ut.create_miwa(yg.boxes)
        sch = sy.SchurPolynomial(yg)
        self.assertEqual(sch.explicit(t), sch._schur_characters(t))


    def test_schur_421(self) -> None:
        partition  = (4, 2, 1)
        yg = sy.YoungDiagram(partition)
        t = ut.create_miwa(yg.boxes)
        sch = sy.SchurPolynomial(yg)
        self.assertEqual(sch.explicit(t), sch._schur_characters(t))

    def test_schur_43(self) -> None:
        partition  = (4, 3)
        yg = sy.YoungDiagram(partition)
        t = ut.create_miwa(yg.boxes)
        sch = sy.SchurPolynomial(yg)
        self.assertEqual(sch.explicit(t), sch._schur_characters(t))


    def test_schur_753(self) -> None:
        partition  = (7, 5, 3)
        yg = sy.YoungDiagram(partition)
        t = ut.create_miwa(yg.boxes)
        sch = sy.SchurPolynomial(yg)
        self.assertEqual(sch.explicit(t) - sch._schur_characters(t), 0)


    def test_schur_753321(self) -> None:
        partition  = (7, 5, 3, 1)
        yg = sy.YoungDiagram(partition)
        t = ut.create_miwa(yg.boxes)
        sch = sy.SchurPolynomial(yg)
        self.assertEqual(sch.explicit(t) - sch._schur_characters(t), 0)



if __name__ == '__main__':
    unittest.main()
        
