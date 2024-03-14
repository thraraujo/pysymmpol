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

    def __init__(self, young: YoungDiagram) -> None:
        '''
        Initialization of the Hall-Littlewood polynomials.
        It depends on a partition and on the coordinates
        x = [x1, x2, ..., xn]. 
        '''

        self._young = young
        self._partition = young.partition


    def _factor(self, partition, Q: object=0):
        '''
        The Hall-Littlewood polynomials are built from some differents terms.
        The first is the product
        prod_{i >= 0} prod_{j=1}^{p(i)} frac{(1- Q)}{(1- Q^j)}
        Observe that this factor depends only on the partition. 
        '''
        prod = 1

        for i in range(max(partition)+1):
            if Q == 1:
                for j in range(1, partition.count(i)+1):
                    prod *= sp.Rational(1, j)
            else:
                for j in range(1, partition.count(i)+1):
                    prod *= (1 - Q) / (1 - Q**j)

        if type(prod) is float and abs(prod) > 1:
            return int(prod)
        else:
            return prod
        

    def _quotient(self, x, i, j, Q: object=0): 
        '''
        Next we need to consider the product that is inside the sum.
        prod_{i < j} frac{xi - Q xj}{ xi - xj}. Let us first calculate
        the quotient. 
        '''
        return (x[i] - Q * x[j]) / (x[i] - x[j])


    def _xproducts(self, x, partition):
        '''
        We also need to calculate the products of the coordinates
        power the partition legs.
        '''

        prod1 = 1

        for i in range(1, len(partition) +1):
            prod1 *= x[i]**partition[i]

        return prod1
    

    def explicit(self, x: tuple, Q: object=0, pol: bool=False): 
        '''
        Finally, we need to calculate the Hall-Littlewood polynomials themselves.
        '''

        # Here I write partitions and coordinates in the same length
        # Problems with this conditional
        n = len(x)
        if len(self._partition) < n:
            _x = dict(enumerate(sp.symbols('x1:{}'.format(n+1)),1)) 
            partition = self._partition + (0,)*abs(len(self._partition) - n)
        elif len(self._partition) > n:
            _x = dict(enumerate(sp.symbols('x1:{}'.format(len(self._partition)+1)),1))
            partition = self._partition
        else:
            _x = dict(enumerate(x,1)) 
            partition = self._partition


        L = dict(enumerate(partition, 1)) # The partition lambda as a dictionary

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
                    prod2 *= self._quotient(xx, i, j, Q)

            sum += prod1 * prod2

        hl = (self._factor(partition, Q) * sum).simplify()

        #if self._Q == 0:
            #'''
            #If we want to use the _test below, I need to comment this line and calculate
            #the Hall-Littlewood polynomial using the else condition below. Here I
            #particularize the case Q=0 because the final answer has the denominator is simplified. 
            #'''
            #xx = tx_power_sum(self._young.boxes, len(x))
            #assert (hl - SchurPolynomial(self._young).explicit(xx)).simplify() == 0
            #hl = SchurPolynomial(self._young).explicit(xx, pol)

        if pol:
            return sp.Poly(hl, domain='QQ')
        else:
            return hl


    def _test(self):
        
        '''
        MEANINGLESS
        Here I want to test if the definitions are correct. In particular,
        I want to see if at Q = 0, the Hall-Littlewood polynomials
        give the Schur polynomials, as we expect. 
        '''
        
        test = HallLittlewood(self._young, len(self._x), 0).hall_littlewood()
        schur = SchurPolynomial(self._young)._xschur(len(self._x)).simplify()

        assert (test - schur).simplify() == 0
        return 0
