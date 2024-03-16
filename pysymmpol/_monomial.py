import sympy as sp
from .conjugacy import ConjugacyClass


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
