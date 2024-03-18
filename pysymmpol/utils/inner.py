import sympy as sp
import numpy as np

from ..partitions.young import YoungDiagram
from ..partitions.conjugacy import ConjugacyClass


def _accel_asc(n):
    '''
    This fast algorithm to generate integer
    partitions is the heart of this project.

    The author of this beauty is Jerome Kelleher, and 
    He argues that it is the fasted algorithm available nowadays.
    See more here: https://jeromekelleher.net/category/combinatorics.html
    '''
    a = [0 for i in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]


def _power_sum(x: tuple, j: int): #Power Sum
    '''
    This is the calculation of power sums. 
    In order to calculate the Newton polynomial P_j(\vec{x}).
    These objects are equivalent to 'j*t_j' where t_j
    are the Miwa coordinates.
    '''
    y = np.array(x)
    yj = y ** j
    return sum(yj)
