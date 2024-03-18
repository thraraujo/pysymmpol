import os 
import sys
sys.path.insert(0, os.path.abspath('..'))

import unittest
import sympy as sp
import pysymmpol as sy

class TestYoung(unittest.TestCase):

    def test_initialization_young_list(self) -> None:
        with self.assertRaises(ValueError):
            sy.State(-1)
        with self.assertRaises(TypeError):
            sy.State('the_dude')
            sy.State([3])



    def test_conjugacy_partition_states(self) -> None:
        state1 = sy.State(50)
        state2 = sy.State(60)
        self.assertEqual(len(state1.conjugacy_states()), 204226)
        self.assertEqual(len(state1.conjugacy_states()), len(state1.partition_states()))
        self.assertEqual(len(state2.conjugacy_states()), len(state2.partition_states()))
        self.assertEqual(len(state2.conjugacy_states()), 966467)


if __name__ == '__main__':
    unittest.main()
        
