import sympy as sp
import numpy as np

from ..partitions.young import YoungDiagram
from ..partitions.conjugacy import ConjugacyClass


def vandermonde(x: tuple): # Vandermonde determinant
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


def newton_polynomial(x: tuple, vector: ConjugacyClass):
    '''
    Calculation of the Newton polynomials themselves.
    '''

    k = vector.conjugacy

    r = len(k)
    newton = 1

    for j in range(r):
        newton *= _power_sum(x, j+1) ** (k[j])
    return newton


def character(young_diagram: YoungDiagram, vector: ConjugacyClass):
    '''
    Here I calculate the characters using the Frobenius Character Formula.
    '''

    young = young_diagram.partition
    m = len(young)

    x = create_x_coord(m)

    l = []    # These are the powers (l1, l2, ..., lm) defined in the formula.
    for i in range(m):
        l.append(young[i] + m - i -1) # the minus comes from the fact that python lists start at 0

    power = 1 # This is the coefficient in x I need to extract
    for j in range(len(l)):
        power *= x[j]**l[j]

    polynomial = sp.poly( vandermonde(x) * newton_polynomial(x, vector) )
    coeff = polynomial.coeff_monomial(power)

    return coeff


def partitions_of(n: int) -> int:
    '''
    Small function that returns the
    partition of the integer n.
    '''
    return sum(1 for x in _accel_asc(n))


def create_miwa(n: int) -> dict:
    '''
    This function creates the appropriate Miwa coordinates as a dictionary.
    '''
    return dict(enumerate(sp.symbols(f't1:{n+1}'), 1))


def create_x_coord(m) -> tuple:
    '''
    This function creates the x coordinates.
    '''
    x = np.array(sp.symbols(f'x1:{m+1}'))

    return tuple(x)


def tx_power_sum(n: int, m: int=1) -> tuple:
    '''
    This function creates the x coordinates.
    The integer n is the length of the Miwa coordinates.
    The integer m is the length of the x coordinates. 
    '''
    x = np.array(sp.symbols(f'x1:{m+1}'))
    t = [sum(x ** j) / j  for j in range(1, n+1)] 

    return tuple(t)
