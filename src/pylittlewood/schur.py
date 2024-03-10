import sympy as sp
import numpy as np
from itertools import product
from .young import YoungDiagram, ConjugacyClass
from .symmetric import HomogeneousPolynomial, _Monomial
from .states import State

class SchurPolynomial:
    '''
    Here is an implementation of the Schur polynomials.
    We calculate these polynomials using the determinant of
    the Homogeneous polynomials. But we assert that they are
    correct using the calculation of the same objects
    from the characters. 
    '''

    def __init__(self, partition: YoungDiagram) -> None:
        self._young = partition
        self._partition = partition.partition


    def _vandermonde(self):
        '''
        Here I want to calculate the Vandermonde polynomial. 
        '''
        m = len(self._partition) # Number of rows of the partition
        vandermonde = 1
        x = np.array(sp.symbols('x1:{}'.format(m+1)))

        for i in range(m):
            for j in range(m):
                if i < j:
                    vandermonde *= (x[i] - x[j])
        return vandermonde


    def _power_sum(self, j: int): #NewtonFactor
        '''
        In order to calculate the Newton polynomial,
        we need the power sums P_j(\vec{x}). These objects
        are equivalent to j t_j where t_j are the Miwa coordinates.
        The partition gives the number of rows.
        '''
        m = len(self._partition) # number of rows

        x = np.array(sp.symbols('x1:{}'.format(m+1)))
        xj = x**j

        power = 0
        for j in range(m):
            power += xj[j]
        return power


    def _newton_polynomial(self, vector: ConjugacyClass):
        '''
        Calculation of the Newton polynomials themselves.
        '''

        k = vector.conjugacy
        
        m = len(self._partition) # number of rows
        r = len(k)
        newton = 1
        for j in range(r):
            newton *= self._power_sum(j+1)**(k[j])
        return newton


    def _characters(self, vector: ConjugacyClass):
        '''
        Here I calculate the characters using the Frobenius Character Formula.
        '''

        m = len(self._partition)

        x = np.array(sp.symbols('x1:{}'.format(m+1)))

        l = []    # These are the powers (l1, l2, ..., lm) defined in the formula.
        for i in range(m):
            l.append(self._partition[i] + m - i -1) # the minus comes from the fact that python lists start at 0

        power = 1 # This is the coefficient in x I need to extract
        for j in range(len(l)):
            power *= x[j]**l[j]

        polynomial = sp.poly( self._vandermonde() * self._newton_polynomial(vector) )
        coeff = polynomial.coeff_monomial(power)

        return coeff


    def _schur_characters(self):
        '''
        Here we calculate the Schur polynomial in 
        terms of Miwa coordinates using the characters.
        This expression will be used to test our code. 
        '''

        level = self._young.boxes
        if level == 0:
            return 1
        else:
            schur = 0
            vectors = State(level).conjugacy_states
            t = sp.symbols('t1:{}'.format(level + 1))

            for vector in vectors:
                vector = ConjugacyClass(vector)
                mono = _Monomial(vector)
                schur += mono._monomial(t) * self._characters(vector)
            return schur


    def schur(self, other_young: YoungDiagram = YoungDiagram([0])):
        '''
        Here we calculate the Schur polynomial in 
        terms of Miwa coordinates using the determinant formula.
        If we Enter an second argument, it calculates the
        Skew-Schur Polynomials.
        '''

        # Here I write the two partitions in the same form
        while len(other_young._partition) < len(self._partition):
            other_young._partition.append(0)

        if self._partition[0] == 0:
            return 1
        else:
            l = len(self._partition) # length = number of rows in the partition
            level = self._young.boxes
            t = sp.symbols('t1:{}'.format(level + 1))

            H = sp.zeros(l, l)

            for (n,m) in product(tuple(range(l)), repeat = 2):
                h = HomogeneousPolynomial(self._partition[n] - other_young._partition[m] - n + m)
                H[n,m] = h.expand(t)

            # Here I test to see if the Schur polynomials are
            # correct. Probably I can extend this test to the skew Schur
            # polynomials as well.
            if other_young._partition == [0]*len(self._partition):
                assert H.det() == self._schur_characters()


            return H.det()

    def _xschur(self, N=1):
        '''
        Here I calculate the schur polynomials in the x coordinates
        x = (x1, ..., xn). I am gonna use to assert the results of the
        Hall-Littlewood polynomials. 
        '''
        l = len(self._partition) # length = number of rows in the partition
        x = np.array(sp.symbols('x1:{}'.format(N+1)))

        level = self._young.boxes
        t = [sum(x**(j+1))/(j+1) for j in range(level+1)]

        H = sp.zeros(l, l)

        for (n,m) in product(tuple(range(l)), repeat = 2):
            h = HomogeneousPolynomial(self._partition[n] - n + m)
            H[n,m] = h.expand(t)

        return H.det().simplify()

    def skew_schur(self, other_young: YoungDiagram):
        '''
        Here we calculate the Skew Schur polynomials.
        It is a wrap of the schur method. 
        '''
        return self.schur(other_young)
