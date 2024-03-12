import sympy as sp
from .conjugacy import ConjugacyClass
from .states import State
from ._monomial import _Monomial

class HomogeneousPolynomial:
    '''
    Here we have a class that defines the 
    complete homogeneous polynomials.
    '''


    def __init__(self, level: int) -> None:
        self._level = level

    @property
    def level(self) -> int:
        return self._level

    def _states(self) -> tuple:
        '''
        For this level, we find the conjugacy class
        states.
        '''
        states = State(self._level)
        return states.conjugacy_states()
        

    def explicit(self, t: tuple, pol: bool=False): # t are the Miwa coordinates
        '''
        This function gives the expansion of the complete
        symmetric polynomials.
        '''

        if isinstance(t, dict):
            t = tuple(t.values())

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

                if pol:
                    return sp.poly(polynomial, domain='QQ')
                else:
                    return polynomial
