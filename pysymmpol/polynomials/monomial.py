from .hall import HallLittlewoodPolynomial
from ..partitions.young import YoungDiagram

class MonomialPolynomial:

    def __init__(self, young: YoungDiagram) -> None:

        '''
        Initialization of the Monomial Symmetric Polynomials. 
        These polynomial are defined from the Hall-Littlewood
        polynomials at the point Q = 1.

        They depends on a partition and on the coordinates
        x = [x1, x2, ..., xn]. 
        '''
        self._young = young
        self._partition = young.partition
        self._hall = HallLittlewoodPolynomial(self._young)


    @property
    def partition(self):
        '''
        Getter for the partition associated to this monomial.
        '''
        return self._partition


    def explicit(self, x: tuple, pol: bool=False): 
        '''
        Determine the explicit expression for the Monomial
        Symmetric Polynomials. 
        '''
        return self._hall.explicit(x, 1, pol)
