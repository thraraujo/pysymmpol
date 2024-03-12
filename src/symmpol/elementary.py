import numpy as np
from .homogeneous import HomogeneousPolynomial


class ElementaryPolynomial:
    '''
    Here we have a class that defines the 
    elementary symmetric polynomials. Their
    implementation is built from the homogeneous polynimials.
    In particular, we use the relations: e_k(t) = (-1) h_k(-t) 
    '''


    def __init__(self, level: int) -> None:
        self._level = level


    @property
    def level(self) -> int:
        return self._level


    def explicit(self, t, pol=False):
        '''
        This function gives the expansion of the elementary
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
                n = self._level
                H = HomogeneousPolynomial(n)

                tl = tuple(-np.array(t))

                polynomial = ((-1)**n) * (H.explicit(tl, pol))

                return polynomial
