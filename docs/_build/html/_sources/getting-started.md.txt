# Getting Started


## Dependencies

This package was tested on:
- Python 3.12
- Sympy 1.12
- Numpy 1.26.4

## Installation

The package can be installed with pip:

```bash
$ pip install pysymmpol
```

## Basic Usage

`PySymmPol` has seven main classes for manipulation of different 
symmetric polynomials.

### YoungDiagram and ConjugacyClass

For the construction and manipulation of Young diagrams, we need to import 
the YoungDiagram and the ConjugacyClass classes. 
```python
from pysymmpol import YoungDiagram, ConjugacyClass
```
The difference between these two classes are the different representations of 
the diagrams. YoungDiagram represents them using a monotonic decreasing 
sequence; and ConjugacyClass as a sequence representing the cycle of the symmetric 
group S_n. Let us create the partition (3,2,1) represented as a
tuple in the YoungDiagram class, and as a dictionary in the ConjugacyClass {1: 1, 2: 1, 3: 1},
```python
young = YoungDiagram((3,2,1))
conjugacy = ConjugacyClass({1: 1, 2: 1, 3: 1})
```
Both objects describe the same mathematical entity, the partition 6=3+2+1. In fact, 
we have the usual pictorial representation 
```python
young.draw_diagram(4)

conjugacy.draw_diagram(4)
```
that gives the output (the argument 4 means that we draw the octothorpe, 
there are 4 other symbols).

```python
#
# #
# # #

#
# #
# # #
```
Description of the other functionalities can be seen in the tutorial. 

### Homogeneous and Elementary Polynomials

These classes can be initialized as 
```python
from pysymmpol import HomogeneousPolynomial, ElementaryPolynomial
from pysymmpol.utils import create_miwa
```
We also imported the function `create_miwa` in the `utils` module, for convenience. 
Let us create the polynomials at the level n=3. We can instanciate them, and find 
the explicit expression using the method `explicit(t)`, there `t` are the 
Miwa coordinates.
```python
n = 3
t = create_miwa(n)

homogeneous = HomogeneousPolynomial(n)
elementary = ElementaryPolynomial(n)
print(f"homogeneous: {homogeneous.explicit(t)}")
print(f"elementary: {elementary.explicit(t)})
```
that gives the output 
```
homogeneous: t1**3/6 + t1*t2 + t3
elementary: t1**3/6 - t1*t2 + t3
```

### Schur Polynomials

For Schur polynomials, we need to instanciate a partition before the polynomial itself. 
Let us use the Young diagram we considered a few lines above,
```python
from pysymmpol import YoungDiagram
from pysymmpol import SchurPolynomial
from pysymmpol.utils import create_miwa
```
The class Young diagram has a getter for the number of boxes in the diagram. 
We use it to built the Miwa coordinates. The class SchurPolynomial is 
instanciated using the young diagram.
```python
young = YoungDiagram((3,2,1))
t = create_miwa(young.boxes)

schur = SchurPolynomial(young)

print(f"schur: {schur.explicit(t)}")
```
that gives the output 
```
schur: t1**6/45 - t1**3*t3/3 + t1*t5 - t3**2
```
In the documentation and tutorial, you can find some examples to 
find skew-Schur polynomials.


### Monomial Symmetric Polynomials

For Monomial symmetric polynomials, we have a similar structure. 
```python
from pysymmpol import YoungDiagram
from pysymmpol import MonomialPolynomial
from pysymmpol.utils import create_x_coord
```
The only difference is that we import the function `create_x_coord` from the `utils` module.
```python
young = YoungDiagram((3,2,1))

n = 3
x = create_x_coord(n)

monomial = MonomialPolynomial(young)

print(f"monomial: {monomial.explicit(x)}")
```
that gives the output 
```
monomial: x1*x2*x3*(x1**2*x2 + x1**2*x3 + x1*x2**2 + x1*x3**2 + x2**2*x3 + x2*x3**2)
```

### Hall-Littlewood Polynomials

Finally, for the Hall-Littlewood polynomials, besides the partitions, we also 
need the deformation parameter Q (because t has been used to denote the Miwa coordinates). 
```python
from sympy import Symbol
from pysymmpol import YoungDiagram
from pysymmpol import HallLittlewoodPolynomial
from pysymmpol.utils import create_x_coord
```
The method `explicit(x, Q)` needs another argument. 
```python
Q = Symbol('Q')
young = YoungDiagram((3,2,1))

n = 3
x = create_x_coord(n)

hall_littlewood = HallLittlewoodPolynomial(young)

print(f"hall-littlewood: {hall_littlewood.explicit(x, Q)}")
```
that gives the output 
```
hall-littlewood: x1*x2*x3*(-Q**2*x1*x2*x3 - Q*x1*x2*x3 + x1**2*x2 + x1**2*x3 + x1*x2**2 + 2*x1*x2*x3 + x1*x3**2 + x2**2*x3 + x2*x3**2)
```
