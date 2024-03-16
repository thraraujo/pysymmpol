from dataclasses import dataclass
import sympy as sp
import numpy as np

from .utils import _accel_asc


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
            raise ValueError("Level must be a non-negative integer.")
        if not isinstance(self._level, int):
            raise TypeError("Level must be a non-negative integer.")


    @property
    def level(self) -> int:
        return self._level


    def _conjugacy_states(self) -> tuple:
        '''
        For a level n, this function gives all vector k
        (the conjugacy class states) that belong to this subspace.
        Remember that these states are built with the operators J_{-m} 
        of the Heisenberg algebra. In other words, these are
        bosonic states. This nonpublic method can make our life easier.
        '''
        lev = self.level

        if lev == 0:
            return tuple()
        else:
            vectors_k = [] # Here I create a list to save the vectors
            vectors_k_tuple = [] # Here I create a list to save the tuples after all manipulations

            for a in _accel_asc(lev):
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

            for a in _accel_asc(lev):
                vec = [0]*lev
                for i in range(len(a)):
                    vec[i] += a[i]
                vec.sort(reverse=True)
                partitions.append(vec)

            for a in partitions:
                partitions_tuple.append(tuple(a))

            return tuple(partitions_tuple)
