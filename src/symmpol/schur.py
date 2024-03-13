import sympy as sp
import numpy as np
from itertools import product

from .young import YoungDiagram
from .conjugacy import ConjugacyClass
from .homogeneous import HomogeneousPolynomial
from .elementary import ElementaryPolynomial
from .states import State
from ._monomial import _Monomial
from .utils import create_x_coord

class SchurPolynomial:
    '''
    Here we have two implementations of the Schur polynomials.
    1. We calculate these polynomials using the determinant of
    the Homogeneous polynomials.
    2. We also calculate them using the characters. We use this second
    implementation to test our results.
    '''


    def __init__(self, young: YoungDiagram) -> None:
        self._young = young
        self._partition = young.partition


    def schur(self, t: tuple, pol: bool=False, other_young: YoungDiagram = YoungDiagram((0,))):
        '''
        Here we calculate the Schur polynomial in 
        terms of Miwa coordinates using the determinant formula.
        If we Enter an second argument, it calculates the
        Skew-Schur Polynomials.
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


    def skew_schur(self, t, other_young: YoungDiagram, pol: bool=False):
        '''
        Here we calculate the Skew Schur polynomials.
        It is a wrap of the schur method. 
        '''
        return self.schur(t, pol, other_young)



    ####### CHARACTER IMPLEMENTATION ####### 
    # I use this implementation for tests. But it also have
    # the calculation of characters, and it might be useful later. 

    def _vandermonde(self, x: tuple): # Vandermonde determinant
        '''
        Here I want to calculate the Vandermonde polynomial. 
        '''
        m = len(x)

        vandermonde = 1

        for i in range(m):
            for j in range(m):
                if i < j:
                    vandermonde *= (x[i] - x[j])
        return vandermonde


    def _power_sum(self, x: tuple, j: int): #Power Sum
        '''
        In order to calculate the Newton polynomial,
        we need the power sums P_j(\vec{x}). These objects
        are equivalent to 'j*t_j' where t_j are the Miwa coordinates.
        The partition gives the number of rows.
        '''

        y = np.array(x)
        yj = y ** j

        return sum(yj)


    def _newton_polynomial(self, x: tuple, vector: ConjugacyClass):
        '''
        Calculation of the Newton polynomials themselves.
        '''

        k = vector.conjugacy
        
        r = len(k)
        newton = 1

        for j in range(r):
            newton *= self._power_sum(x, j+1) ** (k[j])
        return newton


    def _character(self, vector: ConjugacyClass):
        '''
        Here I calculate the characters using the Frobenius Character Formula.
        '''

        m = len(self._partition)

        x = create_x_coord(m)

        l = []    # These are the powers (l1, l2, ..., lm) defined in the formula.
        for i in range(m):
            l.append(self._partition[i] + m - i -1) # the minus comes from the fact that python lists start at 0

        power = 1 # This is the coefficient in x I need to extract
        for j in range(len(l)):
            power *= x[j]**l[j]

        polynomial = sp.poly( self._vandermonde(x) * self._newton_polynomial(x, vector) )
        coeff = polynomial.coeff_monomial(power)

        return coeff


    def _schur_characters(self, t: tuple, pol: bool=False):
        '''
        Here we calculate the Schur polynomial in 
        terms of Miwa coordinates using the characters.
        This expression will be used to test our code. 
        '''

        if isinstance(t, dict):
            t = tuple(t.values())

        level = self._young.boxes

        if level == 0:
            return 1
        else:
            schur = 0
            vectors = State(level).conjugacy_states() # This method gives the stattes as a tuple
            #t = sp.symbols('t1:{}'.format(level + 1))

            for vector in vectors:
                vector = ConjugacyClass(vector)
                mono = _Monomial(vector)
                schur += mono._monomial(t) * self._character(vector)

            if pol:
                return sp.poly(schur, domain='QQ')
            else:
                return schur.expand().simplify()

