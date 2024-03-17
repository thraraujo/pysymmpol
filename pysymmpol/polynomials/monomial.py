from .hall import HallLittlewoodPolynomial
from ..partitions.young import YoungDiagram

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
