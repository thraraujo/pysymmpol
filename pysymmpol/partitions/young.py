from __future__ import annotations

from dataclasses import dataclass
from itertools import pairwise
from numbers import Number

import numpy as np
import sympy as sp

'''
In this module we define the class of Young (or Ferrers) Diagrams
using the usual partition notation.
'''

@dataclass(frozen=True)
class YoungDiagram:
    '''
    Represents Young diagrams in the standard
    partition notation. It is a monotonic decreasing sequence
    L = (L1, L2, L3, ...,Ln) with L1 >= L2 >=L3 >= ... >= Ln.

    Example: (3,2,2,1).
    '''
    _partition: tuple


    def __post_init__(self) -> None:

        '''
        Validade Young diagram:
        1) The argument is a tuple or a numpy array.
        2) The argument must be monotonic decreasing.
        '''

        # Validade its type.
        val_A = isinstance(self._partition, (tuple, np.ndarray))
        val_B = all(isinstance(m, Number) for m in self._partition )
        if val_A and val_B:
            par = self._partition
        else:
            raise TypeError(f"Argument must be a tuple or a numpy array with numeric entries.")

        # Validade if it is a monotonic decreasing sequence. 
        val_C = all(x >= y for x,y in pairwise(par))

        if not val_C:
            raise ValueError("Argument must be a monotonic decreasing sequence.")
        

    @property
    def partition(self) -> tuple:
        '''
        Gives the partition notation for the Young diagram.
        '''
        return self._partition


    @property
    def rows(self) -> int:
        '''
        This gives the number of rows in the diagram.
        '''
        return len(self.partition)


    @property
    def columns(self) -> int:
        '''
        This gives the number of columns in the diagram.
        '''
        return self.partition[0]


    @property
    def boxes(self) -> int: 
        '''
        This gives the number of boxes in the diagram.
        '''
        return sum(self.partition)


    def count_diagonal(self) -> int: 
        '''
        Gives the number of boxes in the diagonal.
        '''

        return len(self.frobenius_coordinates())


    def draw_diagram(self, n: int=0) -> None:
        '''
        Pictorial representation of the Young/Ferrers diagram in
        French notation. 
        '''
        emo = ['■', '•', '🎲', '🎯', '#']
        for i in range(self.rows):
            if self.partition[-i-1] > 0:
                print(f"{emo[n]} " * self.partition[-i-1])


    def conjugacy_partition(self) -> dict:
        '''
        Converts the partition notation to 
        the conjugacy class notation
        Example:
        [4,4,4,2,2,1] to {1: 1, 2: 2, 3: 0, 4: 3}.

        The algorithm works as follows:

        It creates a list of zeros with length equal
        to the largest row. In our case, 4
        conjugacy = [0, 0, 0, 0].

        It iterates over the partition, adding 1s to the
        corresponding slot in the conjugacy vector, for example
        partition (4,4,4,2,2,1) gives
        - loop 01: i = 4
            conjugacy[3] = 1
            conjugacy = [0,0,0,1]
        - loop 02: i = 4
            conjugacy[3] = 2
            conjugacy = [0,0,0,2]
        - loop 03: i = 4
            conjugacy[3] = 3
            conjugacy = [0,0,0,3]
        - loop 04: i = 2
            conjugacy[1] = 1
            conjugacy = [0,1,0,3]
        - loop 05: i = 2
            conjugacy[1] = 2
            conjugacy = [0,2,0,3]
        - loop 06: i = 1
            conjugacy[0] = 1
            conjugacy = [1,2,0,3]

        At the end we convert to a dictionary, and
        conjugacy = {1: 1, 2: 2, 3: 0, 4: 3}
        '''

        length = self.partition[0]
        conjugacy = [0]*length

        for i in self.partition:
            conjugacy[i-1] += 1

        # We finally create a dictionary for the conjugacy class
        conjugacy_class = dict(enumerate(conjugacy,1)) 

        return conjugacy_class


    def transpose(self) -> YoungDiagram:
        '''
        Gives the transposed (or conjugate) Young diagram. 
        For example, the conjugate of [3,2] is [2,2,1].

        Here I follow an interesting property that Knuth mentions in
        TAOCP, volume 4A, equation (11) of section 7.2.1.4 - Other
        representations of partitions. He claims that for any partiton
        L = (L1, L2, ...), its coeficcients satisfy

                            Li - L(i+1) = LTci

        where LTc is the conjugate of L in the conjugacy class notation. 
        '''

        m = tuple(self.partition) + (0,) # pad a zero at the end of the partition tuple

        m_transpose_un = np.array([])

        for j in range(len(m)-1):
            rows = m[j] - m[j+1] # This gives the number of rows of length j+1
            m_transpose_un = np.append(m_transpose_un, np.array(rows * [j+1]))


        m_transpose = -np.sort(-m_transpose_un, axis=0) # Sorts the numpy array in decreasing order.
        m_transpose = m_transpose.astype(int) # Converts entries into integers

        return YoungDiagram(tuple(m_transpose))


    def frobenius_coordinates(self, fermionic: bool=True) -> list:
        '''
        Frobenius coordinates for the diagrams are determined as follows:

        Given a partition L = (L1, L2, ...), the Frobenius coordinates
        are defined by (a_n , b_n), where a_n = L_n -n and b_n - L'_n - n,
        (the prime denotes the conjugate diagram). Note that
        we subtract 1 because Python lists start at 0. For fermionic
        representations (which are the default), we need to add 1/2 because
        indices are half-integers. Therefore, an overall offset of -1/2 is
        required for fermionic representation, while -1 suffices for the
        standard representation.

        Additionally, I prefer to consider the representation where all
        negative sites are occupied. Thus, in the fermionic
        representation, I adjust the notation to -b for clarity.
        '''

        if fermionic:
            sign = -1
            offset = - sp.Rational(1,2)
        else: 
            sign = 1
            offset = - 1

        partition = self.partition
        conjugate = self.transpose().partition
        FrobCoor = []

        short = min([partition, conjugate], key=len)

        for m in range(len(short)):

            alpha = partition[m] - m + offset
            beta  = conjugate[m] - m + offset

            if alpha < 0 or beta < 0: 
                break

            FrobCoor.append((alpha, sign * beta))

        return FrobCoor


    def contains(self, other_young: YoungDiagram) -> bool:
        '''
        Checks if the original partition contains the new one.
        '''

        # Create some lists to concatenate these objects. 
        partition = self.partition
        other_partition = other_young.partition

        # Pad some zeros to write the two partitions in the same form (length).
        if len(other_partition) < len(partition):
            other_partition += (0,)*abs(len(partition) - len(other_partition))
        elif len(other_partition) > len(partition):
            partition += (0,)*abs(len(partition) - len(other_partition))

        for x,y in zip(partition, other_partition):
            if x < y:
                return False

        return True


    def interlaces(self, other_young: YoungDiagram) -> bool:
        '''
        Checks if the original partition interlaces the new one.
        '''

        # Create some lists because I need to concatenate these objects. 
        partition = self.partition
        other_partition = other_young.partition

        # Here I pad some zeros to write the two partitions in the same form (length).
        if len(other_partition) < len(partition):
            other_partition += (0,)*abs(len(partition) - len(other_partition))
        elif len(other_partition) > len(partition):
            partition += (0,)*abs(len(partition) - len(other_partition))

        if not self.contains(other_young):
            return False
        else:
            for x,y in zip(other_partition, partition[1:]):
                if x < y:
                    return False

        return True
            

