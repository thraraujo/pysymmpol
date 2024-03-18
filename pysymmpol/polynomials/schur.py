import sympy as sp
import numpy as np
from itertools import product

from .homogeneous import HomogeneousPolynomial, _Monomial
from .elementary import ElementaryPolynomial

from ..partitions.states import State
from ..partitions.young import YoungDiagram
from ..partitions.conjugacy import ConjugacyClass
from ..utils.tools import character


class SchurPolynomial:
    '''
    Implementations of the Schur polynomials.
    1. We calculate these polynomials using the determinant of
    the Homogeneous polynomials.
    2. We also calculate them using the characters. We use this second
    implementation to test our results.
    '''


    def __init__(self, young: YoungDiagram) -> None:
        '''
        Initialization of the Schur polynomia. It depends on the
        Young Diagram.
        '''
        self._young = young
        self._partition = young.partition


    def explicit(self, t: tuple, pol: bool=False, other_young: YoungDiagram=YoungDiagram((0,))):
        '''
        Calculates the Schur polynomial in terms of Miwa coordinates using
        the determinant formula.

        First argument is the set of Miwa coordinates. The second argument is a
        boolean to define the sympy polynomial. The third argument gives the 
        skew Schur polynomials.

        There is a method below, skew_schur, to make the calculation of
        skew-Schur polynomials more explicit.
        '''

        if isinstance(t, dict):
            t = tuple(t.values())

        if not self._young.contains(other_young):
            return 0

        a = self._partition
        b = other_young._partition

        # Here I write the two partitions in the same form
        if a < b:
            while len(a) < len(b):
                a += (0,)
        else: 
            while len(a) > len(b):
                b += (0,)

        if self._partition[0] == 0:
            return 1
        else:
            l = self._young.rows - a.count(0)
            level = self._young.boxes

            H = sp.zeros(l, l)

            for (n,m) in product(tuple(range(l)), repeat=2):
                h = HomogeneousPolynomial(a[n] - b[m] - n + m)
                H[n,m] = h.explicit(t)

            if pol:
                return sp.Poly(H.det(), domain='QQ')
            else:
                return H.det()


    def skew_schur(self, t: tuple, other_young: YoungDiagram, pol: bool=False):
        '''
        This method returns the Skew Schur polynomials. 
        It is a wrap of the explicit method. 
        '''
        return self.explicit(t, pol, other_young)


    def _schur_characters(self, t: tuple, pol: bool=False):
        '''
        This method returts the Schur polynomial in 
        terms of Miwa coordinates using the characters expansion.
        This method is slower, but it is used to test the implementation. 
        It adds another safety layer to this code. 
        '''

        if isinstance(t, dict):
            t = tuple(t.values())

        level = self._young.boxes

        if level == 0:
            return 1
        else:
            schur = 0
            vectors = State(level).conjugacy_states() 

            for vector in vectors:
                vector = ConjugacyClass(vector)
                mono = _Monomial(vector)
                schur += mono._monomial(t) * character(self._young ,vector)

            if pol:
                return sp.poly(schur, domain='QQ')
            else:
                schur_pol = schur
                return schur_pol

