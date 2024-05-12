import os 
import sys
sys.path.insert(0, os.path.abspath('..'))

import unittest
import numpy as np
import sympy as sp
import pysymmpol as sy
from pysymmpol.utils.inner import _conjugate


class TestYoung(unittest.TestCase):

    def test_initialization_young_list(self) -> None:
        with self.assertRaises(TypeError):
            sy.YoungDiagram([4,3,2])
            sy.YoungDiagram([5,3,3,1])
            sy.YoungDiagram([8])
        with self.assertRaises(ValueError):
            sy.YoungDiagram((2,3,4))
            sy.YoungDiagram((4,3,3,2,1,6))
            sy.YoungDiagram((4,3,3,5,1,1))


    def test_properties_young(self) -> None:
        young1 = sy.YoungDiagram((4,3,2))
        self.assertEqual(young1.partition, (4,3,2))
        self.assertEqual(young1.boxes, 9)
        self.assertEqual(young1.rows, 3)
        self.assertEqual(young1.columns, 4)


    def test_methods_young(self) -> None:
        young1 = sy.YoungDiagram((4,3,2))
        self.assertEqual(young1.transpose().partition, (3,3,2,1))
        self.assertEqual(young1.conjugacy_partition(), {1: 0, 2: 1, 3: 1, 4: 1})
        self.assertEqual(young1.count_diagonal(), 2)


    def test_contains_interlacing_young(self) -> None:
        mu = sy.YoungDiagram((7,5,3,1))
        nu1 = sy.YoungDiagram((6,))
        nu2 = sy.YoungDiagram((6,3,2))
        nu3 = sy.YoungDiagram((5,4,2,1))
        nu4 = sy.YoungDiagram((5,2,2,1))
        nu5 = sy.YoungDiagram((5,2,2,1,1))
        self.assertEqual(mu.contains(nu1), True)
        self.assertEqual(mu.contains(nu2), True)
        self.assertEqual(mu.contains(nu3), True)
        self.assertEqual(mu.contains(nu4), True)
        self.assertEqual(mu.contains(nu5), False)
        self.assertEqual(mu.interlaces(nu1), False)
        self.assertEqual(mu.interlaces(nu2), True)
        self.assertEqual(mu.interlaces(nu3), True)
        self.assertEqual(mu.interlaces(nu4), False)
        self.assertEqual(mu.interlaces(nu5), False)


    def test_diagonal_transpose_young(self) -> None:
        mu1 = sy.YoungDiagram((7,5,3,1))
        mu2 = sy.YoungDiagram((6,))
        mu3 = sy.YoungDiagram((6,3,2))
        mu4 = sy.YoungDiagram((5,4,2,1))
        mu5 = sy.YoungDiagram((7,7,6,5,5,1,1))

        self.assertEqual(mu1.count_diagonal(), 3)
        self.assertEqual(mu2.count_diagonal(), 1)
        self.assertEqual(mu3.count_diagonal(), 2)
        self.assertEqual(mu4.count_diagonal(), 2)
        self.assertEqual(mu5.count_diagonal(), 5)

        self.assertEqual(mu1.count_diagonal(), 3)
        self.assertEqual(mu2.count_diagonal(), 1)
        self.assertEqual(mu3.count_diagonal(), 2)
        self.assertEqual(mu4.count_diagonal(), 2)
        self.assertEqual(mu5.count_diagonal(), 5)

        self.assertEqual(mu1.transpose().partition, (4,3,3,2,2,1,1))
        self.assertEqual(mu2.transpose().partition, (1,1,1,1,1,1))
        self.assertEqual(mu3.transpose().partition, (3,3,2,1,1,1))
        self.assertEqual(mu4.transpose().partition, (4,3,2,2,1))
        self.assertEqual(mu5.transpose().partition, (7,5,5,5,5,3,2))

        self.assertEqual(mu1.transpose().partition, _conjugate((7,5,3,1)))
        self.assertEqual(mu2.transpose().partition, _conjugate((6,)))
        self.assertEqual(mu3.transpose().partition, _conjugate((6,3,2)))
        self.assertEqual(mu4.transpose().partition, _conjugate((5,4,2,1)))
        self.assertEqual(mu5.transpose().partition, _conjugate((7,7,6,5,5,1,1)))


    def test_diagonal_hook_length(self) -> None:
        mu = [
            sy.YoungDiagram(np.array((10, 9, 9, 8, 8, 8, 7, 7, 7, 7, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4,
                                      3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))),
            sy.YoungDiagram(np.array((54,53,52,51,50,49,48,47,46,45,44,43,42,41,40))), 
            sy.YoungDiagram(np.array((33,30,28,26,26,22,20,18,18,17,15,12,11))), 
            sy.YoungDiagram(np.array((12,12,12,12,12,12,12,12,12,12,12,12))), 
            sy.YoungDiagram(np.array((10,))), 
            sy.YoungDiagram(np.array((1,1,1,1,1,1,1,1,1,1,1))), 
            ]

        for m in mu:
            self.assertEqual(m.count_diagonal(), m.boxes - np.sum(m.frobenius_coordinates(False)))


    def test_properties_conjugacy(self) -> None:
        conjugacy1 = sy.ConjugacyClass({1: 0, 2: 1, 3: 1, 4: 1})
        self.assertEqual(conjugacy1.conjugacy, (0, 1, 1, 1))
        self.assertEqual(conjugacy1.boxes, 9)
        self.assertEqual(conjugacy1.rows, 3)
        self.assertEqual(conjugacy1.columns, 4)


    def test_methods_conjugacy(self) -> None:
        conjugacy1 = sy.ConjugacyClass({1: 0, 2: 1, 3: 1, 4: 1})
        self.assertEqual(conjugacy1.conjugacy_partition(), (4,3,2))


    def test_equavalence_objects(self) -> None:
        young2 = sy.YoungDiagram((10, 9, 9, 8, 8, 8, 7, 7, 7, 7, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4,
                                  3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
        conjugacy2 = sy.ConjugacyClass({1: 10, 2: 9, 3: 8, 4: 7, 5: 6, 6: 5, 7: 4, 8: 3, 9: 2, 10: 1})
        self.assertEqual(sy.ConjugacyClass(young2.conjugacy_partition()), conjugacy2)
        self.assertEqual(sy.YoungDiagram(conjugacy2.conjugacy_partition()), young2)
        self.assertEqual(young2.boxes, conjugacy2.boxes)
        self.assertEqual(young2.rows, conjugacy2.rows)
        self.assertEqual(young2.columns, conjugacy2.columns)


if __name__ == '__main__':
    unittest.main()
        
