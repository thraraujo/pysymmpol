import sympy as sp
from ..partitions.conjugacy import ConjugacyClass
from ..partitions.states import State

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


class _Monomial:
    '''
    Here we have a class to the calculation of the
    monomials necessary for the calculation of the
    Homogeneous Symmetric Polynomials.
    This is not related to the Monomial Symmetric Polynomials. 
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


    def _monomial(self, t: tuple): # t are the Miwa coordinates
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
