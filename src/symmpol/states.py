from dataclasses import dataclass
import sympy as sp
import numpy as np


@dataclass
class State:
    '''
    For a given level n, I want to find the
    bosonic and fermionic states associated to
    YoungDiagrams and conjugacy classes.
    '''
    _level: int


    def __post_init__(self) -> None:
        if self._level < 0:
            raise ValueError("Level must be non-negative")


    @property
    def level(self) -> int:
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


    def _conjugacy_states(self) -> tuple:
        '''
        For a level n, this function gives all vector k
        (the conjugacy class states) that belong to this subspace.
        Remember that these states are built with the operators J_{-m} 
        of the Heisenberg algebra. In other words, these are
        bosonic states. This nonpublic method will make our life easier.
        '''
        lev = self.level

        if lev == 0:
            return tuple()
        else:
            vectors_k = [] # Here I create a list to save the vectors
            vectors_k_tuple = [] # Here I create a list to save the tuples after all manipulations

            for a in self._accel_asc():
                vec = [0]*lev
                for i in range(len(a)):
                    vec[a[i]-1] += 1
                vectors_k.append(vec)

            for a in vectors_k:
                vectors_k_tuple.append(tuple(a))

            return tuple(vectors_k_tuple)


    def conjugacy_states(self) -> tuple:
        '''
        Here I return the previous method in
        the dictionary form.
        '''
        if self.level == 0:
            return tuple()
        else:
            vectors_k = self._conjugacy_states()
            states = [] # Here is a list to collect all states

            for a in vectors_k:
                states.append(dict(enumerate(a,1)))

            return tuple(states)


    def partition_states(self) -> tuple:
        '''
        For the level n, this function gives the partitions that
        belong to this subspace. These states are built from the
        the fermionic operators. 
        '''
        lev = self._level

        if lev == 0:
            return tuple( )
        else:
            partitions = [] # Here I create a list to save the states
            partitions_tuple = [] # Here I create a list to save the tuples after all manipulations

            for a in self._accel_asc():
                vec = [0]*lev
                for i in range(len(a)):
                    vec[i] += a[i]
                vec.sort(reverse=True)
                partitions.append(vec)

            for a in partitions:
                partitions_tuple.append(tuple(a))

            return tuple(partitions_tuple)
