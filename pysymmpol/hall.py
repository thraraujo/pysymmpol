import sympy as sp
import numpy as np
from .young import YoungDiagram
from .schur import SchurPolynomial
from itertools import permutations, product
from .utils import create_miwa, tx_power_sum


class HallLittlewoodPolynomial:
    '''
    Here is an implementation of the Hall-Littlewood polynomials.
    '''

    def __init__(self, young: YoungDiagram, Q: object) -> None:
        '''
        Initialization of the Hall-Littlewood polynomials.
        It depends on a partition and on the coordinates
        x = [x1, x2, ..., xn]. 

        Hall-Littlewood polynomials are given by (LaTeX)

        P_{lambda}(x_1, dots, x_n; Q) =
                    prod_{i geq 0} prod_{j=1}^{p(i)} frac{1-Q}{1-Q^j} times
                    sum_{omega in mathfrak{S}_n} omega left( x_1^{lambda_1} cdots x_n^{lambda_n}
        prod_{i<j} frac{x_i - Q x_j}{x_i - x_j} right)

        For practical reasons, we need to calculate these terms independently. 
        '''

        self._young = young
        self._partition = young.partition
        self._Q = Q


    def _factor(self, partition):
        '''
        Let us first calculate the prefactor
        prod_{i >= 0} prod_{j=1}^{p(i)} frac{(1- Q)}{(1- Q^j)}
        This term does not depend  on the coordinates (x).
        '''
        prod = 1

        for i in range(max(partition)+1):
            if self._Q  == 1:
                for j in range(1, partition.count(i)+1):
                    prod *= sp.Rational(1, j)
            elif self._Q  == 0:
                return 1
            else:
                for j in range(1, partition.count(i)+1):
                    prod *= (1 - self._Q ) / (1 - self._Q  ** j)

        if type(prod) is float and abs(prod) > 1:
            return int(prod)
        else:
            return prod
        

    def _quotient(self, x, i, j): 
        '''
        Next we need to consider the product that is inside the sum.
        prod_{i < j} frac{xi - Q xj}{xi - xj}. Let us first calculate
        the quotient. Observe that the denominator is the Vandermonde determinant. 
        '''
        return (x[i] - self._Q  * x[j]) / (x[i] - x[j])


    def _xproducts(self, x, partition):
        '''
        We also need to calculate the products of the coordinates
        power the partition arms (or legs, if you prefer).
        '''

        prod1 = 1

        for i in range(1, len(partition) +1):
            prod1 *= x[i]**partition[i]

        return prod1
    

    def explicit(self, x: tuple, pol: bool=False): 
        '''
        Finally, we need put all these definitions together
        to calculate the Hall-Littlewood themselves. 
        '''

        # Here I write partitions and coordinates in the same length
        n = len(x)
        m = len(self._partition)

        diff = n - m 

        if diff < 0:
            # Because we have less coordinates than the partitions.
            # So we pad some zeros to the coordinates. 
            return 0 
        elif diff > 0:
            # You passed more coordinates than the length of the partition.
            # Padding zeros to the partition...
            _x = dict(enumerate(x,1)) 
            partition = self._partition + (0,)*abs(len(self._partition) - n)
        else:
            _x = dict(enumerate(x,1)) 
            partition = self._partition


        # Write the partition as a dictionary. This makes the permutations easier to handle. 
        L = dict(enumerate(partition, 1)) 

        # Here we create a dictionary of permutations of the n indices - starting with 1
        sigma = dict(enumerate(permutations(range(1,len(_x)+1)),1)) 

        sum = 0

        for perm in tuple(sigma.values()): # This is the sum over all permutations

            prod1 = 1
            prod2 = 1

            xx = dict(zip(perm, _x.values())) # Dictionary for the x-coordinates with keys given by the permutations

            prod1 = self._xproducts(xx, L)

            for (i, j) in product(range(1, len(xx)+1), repeat=2):
                if i < j:
                    prod2 *= self._quotient(xx, i, j)

            sum += prod1 * prod2

        hl = (self._factor(partition) * sum).simplify()

        if pol:
            return sp.Poly(hl, domain='QQ')
        else:
            return hl


class MonomialPolynomial(HallLittlewoodPolynomial):

    def __init__(self, young: YoungDiagram) -> None:

        '''
        Initialization of the Monomial Symmetric Polynomials. 

        Since they are Hall-Littlewood polynomials with Q=0, we 
        Defined them through an inheritance of the
        Halllittlewoodpolynomial class. 

        They depends on a partition and on the coordinates
        x = [x1, x2, ..., xn]. 
        '''
        self._young = young
        self._partition = young.partition
        self._Q = 1


    def explicit(self, x: tuple, pol: bool=False): 
        '''
        Here I want to particularize to the case Q = 1. 
        '''
        return super().explicit(x, pol)

