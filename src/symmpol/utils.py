import sympy as sp
import numpy as np


def _accel_asc(n):
    '''
    This fast algorithm to generate integer
    partitions is the heart of this project.

    The author of this beauty is Jerome Kelleher, and 
    He argues that it is the fasted algorithm in the
    market today, so we use it. See more here:
    https://jeromekelleher.net/category/combinatorics.html
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
