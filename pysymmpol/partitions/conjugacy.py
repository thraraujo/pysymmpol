from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from .young import YoungDiagram


@dataclass(frozen=True)
class ConjugacyClass: 
    '''
    Represents the conjugacy class of partitions. For example,
    the partition L = (5, 3, 2, 2, 2, 1, 1) can represented by
    L = (1: 2, 2: 3, 3: 1, 4: 0, 5: 1) that means 2 rows of
    length 1, 3 rows of length 2, 2 rows of length 3, 0 rows of
    length 4 and 1 row of length 5.

    To avoid cluttering, we could simply write the vector
    k = (2, 3, 1, 0, 1) to describe the same partition.
    We do it in an internal method below, but we will return
    dictionaries to avoid confusion with the partitions. 

    The dictionary must be in the form {1: k1, 2: k2, ..., n: kn}

    The vectors k = (k1,k2, ..., kn) represents one conjugacy class of
    the symmetric group S_n, and for this reason, we call it conjugacy
    class vectors.

    In terms of the Heisenberg algebra, the components of these vectors
    denote the power of the Heisenberg operator J_{-n} when acting on the
    vacuum state |0>
    '''
    _conjugacy_vector: dict


    def __post_init__(self) -> None:
        '''
        Validates the form of the dictionary: {1: k1, 2: k2, ..., n: kn}
        '''

        keys = tuple(self._conjugacy_vector.keys())
        A = all([keys[i] - keys[i-1] == 1 for i in range(1, len(keys))])

        if keys[0] != 1 or not A:
           raise TypeError("Argument must be a dictionary in the form {1: k1, 2: k2, ..., n: kn}")


    @property
    def conjugacy(self) -> tuple:
        '''
        Getter for the conjugacy class vector
        '''
        return tuple(self._conjugacy_vector.values())


    @property
    def rows(self) -> int:
        '''
        Gives the number of rows in this diagram.
        This is the sum over the conjugacy. 
        '''
        return sum(self.conjugacy)


    @property
    def columns(self) -> int:
        '''
        Gives the number of columns in this diagram.
        for the vector (1: k1, ..., n: kn), the number of
        columns is n: the last entry in the dictionary.  
        (max iterates over the keys). 
        '''
        return max(self._conjugacy_vector)


    @property
    def boxes(self) -> int:
        '''
        Gives the number of boxes in the diagram.
        '''
        box = 0 
        for x, y in self._conjugacy_vector.items():
            box += x*y
        return box


    def draw_diagram(self, n=0) -> None:
        '''
        Here we have a pictorial representation of the Young diagram
        associated to the conjugacy class in French notation. 
        Here I just call the function for the partition states. 
        '''
        young = YoungDiagram(self.conjugacy_partition())

        young.draw_diagram(n)


    def conjugacy_partition(self) -> tuple:
        '''
        Converts the conjugacy class vector to the partition notation.
        Example:
        [1,2,0,4] is equivalent to the Young
        diagram [4,4,4,4,2,2,1].

        These vectors also denote the bosonic partition state that we can 
        build with the bosonic states
        power of the operators J_{-n} that acts on the
        vaccum state |0>.

        The algorithm works as follows:
        We first create an empty array partition = []
        and we define the range for the loop to be the length of
        the vector k: In our case [1,2,0,3], i = 0,1,2,3.

        
        - loop 01: i = 0 k[0] = 1: 
                row = [1]
                partition = [] U [1] = [1]
        - loop 02: i = 1 k[1] = 2
                row = [2,2]
                partition = [1] U [2,2] = [1,2,2]
        - loop 03: i = 2 k[2] = 0
                row = []
                partition = [1,2,2] U [] = [1,2,2]
        - loop 04: i = 3 k[3] = 4
                row = [4,4,4,4]
                partition = [1,2,2] U [4,4,4,4] = [1,2,2,4,4,4,4]

        Then when we put it in decreasing order we find [4,4,4,4,2,2,1]
        '''

        partition = np.array([])

        for i in range(len(self.conjugacy)):
            row = np.repeat(i+1, self.conjugacy[i])
            partition = np.concatenate((partition, row), axis=0)

        partition = -np.sort(-partition, axis=0) # Sorts the numpy array in decreasing order.
        partition = partition.astype(int) # Converts entries into integers

        return tuple(partition)
