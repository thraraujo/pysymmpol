import sys
sys.path.append('/home/thiago/Documents/projects/work/pySymmetricPolynomials/src/')

import unittest
import partition as pt

class TestYoung(unittest.TestCase):

    def test_properties_young(self) -> None:
        young1 = pt.YoungDiagram((4,3,2))
        self.assertEqual(young1.partition, (4,3,2))
        self.assertEqual(young1.boxes, 9)
        self.assertEqual(young1.rows, 3)
        self.assertEqual(young1.columns, 4)


    def test_methods_young(self) -> None:
        young1 = pt.YoungDiagram((4,3,2))
        self.assertEqual(young1.transpose().partition, (3,3,2,1))
        self.assertEqual(young1.conjugacy_partition(), {1: 0, 2: 1, 3: 1, 4: 1})
        self.assertEqual(young1.count_diagonal(), 2)


    def test_properties_conjugacy(self) -> None:
        conjugacy1 = pt.ConjugacyClass({1: 0, 2: 1, 3: 1, 4: 1})
        self.assertEqual(conjugacy1.conjugacy, (0, 1, 1, 1))
        self.assertEqual(conjugacy1.boxes, 9)
        self.assertEqual(conjugacy1.rows, 3)
        self.assertEqual(conjugacy1.columns, 4)


    def test_methods_conjugacy(self) -> None:
        conjugacy1 = pt.ConjugacyClass({1: 0, 2: 1, 3: 1, 4: 1})
        self.assertEqual(conjugacy1.conjugacy_partition(), (4,3,2))


    def test_equavalence_objects(self) -> None:
        young2 = pt.YoungDiagram((10, 9, 9, 8, 8, 8, 7, 7, 7, 7, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5))
        conjugacy2 = pt.ConjugacyClass({1: 0, 2: 0, 3: 0, 4: 0, 5: 6, 6: 5, 7: 4, 8: 3, 9: 2, 10: 1})
        self.assertEqual(pt.ConjugacyClass(young2.conjugacy_partition()), conjugacy2)
        self.assertEqual(pt.YoungDiagram(conjugacy2.conjugacy_partition()), young2)
        self.assertEqual(young2.boxes, conjugacy2.boxes)
        self.assertEqual(young2.rows, conjugacy2.rows)
        self.assertEqual(young2.columns, conjugacy2.columns)



if __name__ == '__main__':
    unittest.main()
        
