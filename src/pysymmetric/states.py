import sympy as sp
import numpy as np

class State:
    '''
    For a given level n, I want to find the
    bosonic and fermionic states associated to
    YoungDiagrams and conjugacy classes.
    '''

    def __init__(self, level: int) -> list:
        self._level = level


    @property
    def level(self):
        return self._level

    def _accel_asc(self):
        '''
        This is a fast algorithm to generate
        integer partitions.  The author of this beauty is Jerome
        Kelleher. He argues that it is the fasted algorithm in the
        market today, so we use it. See more here:
        https://jeromekelleher.net/category/combinatorics.html
        '''
        n = self.level

        a = [0 for i in range(n + 1)]
        k = 1
        y = n - 1
        while k != 0:
            x = a[k - 1] + 1
            k -= 1
            while 2 * x <= y:
                a[k] = x
                y -= x
                k += 1
            l = k + 1
            while x <= y:
                a[k] = x
                a[l] = y
                yield a[:k + 2]
                x += 1
                y -= 1
            a[k] = x + y
            y = x + y - 1
            yield a[:k + 1]


    @property
    def conjugacy_states(self) -> list:
        '''
        For a level n, this function gives all vector k
        (the conjugacy class states) that belong to this subspace.
        Remember that these states are built with the operators J_{-m} 
        of the Heisenberg algebra. In other words, these are
        bosonic states. 
        '''
        lev = self.level

        vectors_k = [] # Here I create a list to save the vectors

        for a in self._accel_asc():
            vec = [0]*lev
            for i in range(len(a)):
                vec[a[i]-1] += 1
            vectors_k.append(vec)

        return vectors_k


    @property
    def partition_states(self) -> list:
        '''
        For the level n, this function gives the partitions that
        belong to this subspace. These states are built from the
        the fermionic operators. 
        '''
        lev = self._level

        partitions = [] # Here I create a list to save the tuples

        for a in self._accel_asc():
            vec = [0]*lev
            for i in range(len(a)):
                vec[i] += a[i]
            vec.sort(reverse=True)
            partitions.append(vec)

        return partitions
