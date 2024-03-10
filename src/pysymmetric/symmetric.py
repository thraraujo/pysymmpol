import sympy as sp
import numpy as np
from .young import ConjugacyClass
from .states import State

'''
Here we have a class to the calculation of the
complete homogeneous polynomials. We Have two implementations:
The first uses the Miwa coordinates and the second uses the
coordinates X
'''

class HomogeneousPolynomial:
    '''
    Here we have a class that defines the 
    complete homogeneous polynomials.
    '''

    def __init__(self, level: int) -> None:
        self._level = level


    def _states(self) -> list:
        '''
        For this level, we find the conjugacy class
        states.
        '''
        states = State(self._level)
        return states.conjugacy_states
        

    def expand(self, t: list): # t are the Miwa coordinates
        '''
        This function gives the expansion of the complete
        symmetric polynomials.
        '''

        if len(t) < self._level:
            raise TypeError('''The list t must have, at least, as many coordinates
                as the level of the conjugacy class''')

        else:
            if self._level < 0:
                return 0

            elif self._level == 0:
                return 1

            else:
                # first of all, given the level above, we need to find the vectors k

                vectors_k = self._states() 
                polynomial = 0

                for vector in vectors_k:
                    A = _Monomial(ConjugacyClass(vector))
                    polynomial += A._monomial(t)

                return polynomial



class ElementaryPolynomial:
    '''
    Here we have a class that defines the 
    elementary symmetric polynomials. Their
    implementation are not completelly independent,
    that is the reason why we include them in this file.
    '''


    def __init__(self, level: int) -> None:
        self._level = level


    def expand(self, t: list):
        '''
        This function gives the expansion of the elementary
        symmetric polynomials.
        '''

        if len(t) < self._level:
            raise TypeError('''The list t must have, at least, as many coordinates
                as the level of the conjugacy class''')
        else:
            if self._level < 0:
                return 0

            elif self._level == 0:
                return 1

            else:
                n = self._level
                H = HomogeneousPolynomial(n)

                # We use that elementary and homogeneous polynomials are related via: E_k(t) = (-1) H_k(-t) 
                t = -np.array(t)

                polynomial = ((-1)**n) * (H.expand(t))

                return polynomial
            


class _Monomial:
    '''
    Here we have a class to the calculation of the
    monomials necessary for the calculation of the
    Homogeneous Symmetric Polynomials.
    '''

    def __init__(self, conjugacy_class: ConjugacyClass) -> None:
        self._conjugacy_class = conjugacy_class
        self._vector_k = conjugacy_class.conjugacy 


    @property
    def level(self) -> int:
        '''
        This function gives the level of the conjugacy class vector k, 
        that is, the number given by sum_i i k_i for a given bosonic 
        state k = (k_1, k_2, ...). This corresponds to the number of
        boxes in the partition described by this conjugacy class vector. 
        '''

        length = len(self._vector_k)
        lev = 0
        for i in range(length):
            lev += (i+1)*self._vector_k[i]
        # Test if the level is the number of boxes equals the level
        assert lev == self._conjugacy_class.boxes 
        return lev


    def _monomial(self, t: list): # t are the Miwa coordinates
        '''
        This function gives the monomial in the 
        definition of the complete symmetric polynomials
        associated with a given conjugacy class vector k.
        '''

        if len(t) < self.level:
            raise TypeError('''The list t must have, at least, as many coordinates
                as the level of the conjugacy class''')
        else:
            vector = self._vector_k
            product = 1

            for n in range(len(vector)):
                product *= sp.Rational(1, sp.factorial(vector[n])) * t[n]**(vector[n])

            return product

