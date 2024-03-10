import sympy as sp
import numpy as np
from .young import YoungDiagram
from .schur import SchurPolynomial
from itertools import permutations, product
#from .symmetric import HomogeneousPolynomial, _Monomial
#from .states import State

class HallLittlewood:
    '''
    Here is an implementation of the Hall-Littlewood polynomials.
    '''

    def __init__(self, young: YoungDiagram, n: int = 2, Q=0) -> None:
        '''
        Initialization of the Hall-Littlewood polynomials.
        It depends on a partition and on the coordinates
        x = [x1, x2, ..., xn]. 
        '''

        self._Q = Q
        self._young = young

        # Here I write partitions and coordinates in the same length
        if len(young.partition) < n:
            self._x = dict(enumerate(sp.symbols('x1:{}'.format(n+1)),1)) 
            self._partition = young.partition + [0]*abs(len(young.partition) - n)
        else:
            self._x = dict(enumerate(sp.symbols('x1:{}'.format(len(young.partition)+1)),1))
            self._partition = young.partition


    def _factor(self):
        '''
        The Hall-Littlewood polynomials are built from some differents terms.
        The first is the product
        prod_{i >= 0} prod_{j=1}^{p(i)} frac{(1- Q)}{(1- Q^j)}
        Observe that this factor depends only on the partition. 
        '''
        prod = 1

        for i in range(max(self._partition)+1):
            if self._Q == 1:
                for j in range(1, self._partition.count(i)+1):
                    prod *= sp.Rational(1, j)
            else:
                for j in range(1, self._partition.count(i)+1):
                    prod *= (1 - self._Q) / (1 - self._Q**j)

        if type(prod) is float and abs(prod) > 1:
            return int(prod)
        else:
            return prod
        

    def _quotient(self, x, i, j): 
        '''
        Next we need to consider the product that is inside the sum.
        prod_{i < j} frac{xi - Q xj}{ xi - xj}. Let us first calculate
        the quotient. 
        '''
        return (x[i] - self._Q * x[j]) / (x[i] - x[j])


    def _xproducts(self, x, partition):
        '''
        We also need to calculate the products of the coordinates
        power the partition legs.
        '''

        prod1 = 1

        for i in range(1, len(partition) +1):
            prod1 *= x[i]**partition[i]

        return prod1
    

    def hall_littlewood(self): 
        '''
        Finally, we need to calculate the Hall-Littlewood polynomials themselves.
        '''

        L = dict(enumerate(self._partition, 1)) # The partition

        # Here we create a dictionary of permutations of the n indices - starting with 1
        sigma = dict(enumerate(permutations(range(1,len(self._x)+1)),1)) 

        sum = 0

        if self._Q == 0:
            '''
            If we want to use the _test below, I need to comment this line and calculate
            the Hall-Littlewood polynomial using the else condition below. Here I
            particularize the case Q=0 because the final answer has the denominator is simplified. 
            '''
            hl = SchurPolynomial(self._young)._xschur(len(self._x)).expand().simplify()

        else:
            for perm in tuple(sigma.values()): # This is the sum over all permutations

                prod1 = 1
                prod2 = 1

                x = dict(zip(perm, self._x.values())) # With this line, I build a dictionary with keys given by the permutations

                prod1 = self._xproducts(x, L)

                for (i, j) in product(range(1, len(x)+1), repeat=2):
                    if i < j:
                        prod2 *= self._quotient(x, i,j)

                sum += prod1 * prod2

            hl = (self._factor() * sum).simplify()

        return hl

    def _test(self):
        '''
        Here I want to test if the definitions are correct. In particular,
        I want to see if at Q = 0, the Hall-Littlewood polynomials
        give the Schur polynomials, as we expect. 
        '''
        
        test = HallLittlewood(self._young, len(self._x), 0).hall_littlewood()
        schur = SchurPolynomial(self._young)._xschur(len(self._x)).simplify()

        assert (test - schur).simplify() == 0
        return 0
