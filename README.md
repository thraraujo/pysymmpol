# pySymmPol: Symmetric Polynomials

![Static Badge](https://img.shields.io/badge/3.12-green?style=plastic&logo=python&logoColor=yellow&label=python)
![Static Badge](https://img.shields.io/badge/Lab-blue?style=plastic&logo=jupyter&logoColor=yellow&label=Jupyter)
![Static Badge](https://img.shields.io/badge/1.26-orange?style=plastic&logo=numpy&logoColor=green&label=Numpy)
![Static Badge](https://img.shields.io/badge/1.12-blue?style=plastic&logo=sympy&logoColor=green&label=Sympy)
![Static Badge](https://img.shields.io/badge/os-Linux?style=plastic&logo=Linux&logoColor=white&label=GNU%2FLinux)
[![License: GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This Python package is designed for the manipulation of various
symmetric polynomials or functions. It includes the following types:

    1. Complete homogeneous Symmetric Functions
    2. Elementary symmetric polynomials
    3. Monomial symmetric polynomials 
    4. Schur polynomials
    5. Hall-Littlewood polynomials

Additionally, the package contains a module with basic functionalities
for manipulating integer partitions and Young diagrams.

Read our Statement of need [here](STATEMENT-OF-NEED.md).

Tutorials can be found in the [documentation page](https://thraraujo.github.io/pysymmpol/index.html).

## Dependencies

This package has been tested with the following versions:
- Python >= 3.11
- SymPy >= 1.11
- NumPy >= 1.26.2

## Installation

The package can be installed with pip:
```bash
$ pip install pysymmpol
```

## Basic Usage

The `PySymmPol` package has seven main classes for manipulating various symmetric polynomials.

### YoungDiagram and ConjugacyClass

For the construction and manipulation of Young diagrams, we need to import 
the YoungDiagram and the ConjugacyClass classes. 
```python
from pysymmpol import YoungDiagram, ConjugacyClass
```
The distinction between these two classes lies in the representations of the diagrams 
they handle. The YoungDiagram class represents diagrams using a monotonic decreasing sequence, 
while the ConjugacyClass class represents them as a sequence representing the cycle 
of the symmetric group SnSn​. For example, let's consider the partition 
(3,2,1), which is represented as a tuple in the YoungDiagram class and 
as a dictionary in the ConjugacyClass class {1: 1, 2: 1, 3: 1}, respectively.
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
that give the same output (the argument 4 means that we draw the octothorpe, and
there are 4 other symbols available).

```python
#
# #
# # #

#
# #
# # #
```
Further details on the other functionalities can be found in the tutorial.

### Homogeneous and Elementary Polynomials

These classes can be initialized as 
```python
from pysymmpol import HomogeneousPolynomial, ElementaryPolynomial
from pysymmpol.utils import create_miwa
```
We also imported the function `create_miwa` from the `utils` module for convenience. 
Now, let's create the polynomials at level n=3. We can instantiate 
them and find their explicit expressions using the `explicit(t)` method, where `t` 
represents the Miwa coordinates. The block
```python
n = 3
t = create_miwa(n)

homogeneous = HomogeneousPolynomial(n)
elementary = ElementaryPolynomial(n)
print(f"homogeneous: {homogeneous.explicit(t)}")
print(f"elementary: {elementary.explicit(t)}")
```
gives the output 
```
homogeneous: t1**3/6 + t1*t2 + t3
elementary: t1**3/6 - t1*t2 + t3
```

### Schur Polynomials

To create Schur polynomials, we first need to instantiate a partition 
before defining the polynomial itself. Let's use the Young diagram we 
considered a few lines above, then
```python
from pysymmpol import YoungDiagram
from pysymmpol import SchurPolynomial
from pysymmpol.utils import create_miwa
```
The YoungDiagram class includes a getter method for retrieving the number 
of boxes in the diagram, which we utilize to construct the Miwa coordinates. 
Subsequently, the SchurPolynomial class is instantiated using the Young diagram. Then
```python
young = YoungDiagram((3,2,1))
t = create_miwa(young.boxes)

schur = SchurPolynomial(young)

print(f"schur: {schur.explicit(t)}")
```
gives
```
schur: t1**6/45 - t1**3*t3/3 + t1*t5 - t3**2
```
The documentation and tutorial contain examples demonstrating how to find 
skew-Schur polynomials.

### Monomial Symmetric Polynomials

For Monomial symmetric polynomials, we have a similar structure. 
```python
from pysymmpol import YoungDiagram
from pysymmpol import MonomialPolynomial
from pysymmpol.utils import create_x_coord
```
The only difference is the function `create_x_coord` from the `utils` module. Therefore,
```python
young = YoungDiagram((3,2,1))

n = 3
x = create_x_coord(n)

monomial = MonomialPolynomial(young)

print(f"monomial: {monomial.explicit(x)}")
```
gives the output 
```
monomial: x1*x2*x3*(x1**2*x2 + x1**2*x3 + x1*x2**2 + x1*x3**2 + x2**2*x3 + x2*x3**2)
```

### Hall-Littlewood Polynomials

In addition to partitions, for the Hall-Littlewood polynomials, 
we also require the deformation parameter Q (as t has been used to 
denote the Miwa coordinates).
```python
from sympy import Symbol
from pysymmpol import YoungDiagram
from pysymmpol import HallLittlewoodPolynomial
from pysymmpol.utils import create_x_coord
```
The method `explicit(x, Q)` needs another argument. Finally, the code
```python
Q = Symbol('Q')
young = YoungDiagram((3,2,1))

n = 3
x = create_x_coord(n)

hall_littlewood = HallLittlewoodPolynomial(young)

print(f"hall-littlewood: {hall_littlewood.explicit(x, Q)}")
```
gives
```
hall-littlewood: x1*x2*x3*(-Q**2*x1*x2*x3 - Q*x1*x2*x3 + x1**2*x2 + x1**2*x3 + x1*x2**2 + 2*x1*x2*x3 + x1*x3**2 + x2**2*x3 + x2*x3**2)
```

# References

Here are some recommended resources covering symmetric polynomials,
combinatorics, and their significance in theoretical physics:

- [Symmetric Functions and Hall Polynomials - Ian Macdonald](https://books.google.com.br/books/about/Symmetric_Functions_and_Hall_Polynomials.html?id=srv90XiUbZoC&redir_esc=y)
- [An Introduction to Symmetric Functions and Their Combinatorics - Eric S. Egge](https://bookstore.ams.org/stml-91)
- [Counting with Symmetric Functions - Mendes and Remmel](https://link.springer.com/book/10.1007/978-3-319-23618-6)
- [Solitons: Differential Equations, Symmetries and Infinite Dimensional Algebras - Miwa, Jimbo and Date](https://books.google.com.br/books/about/Solitons.html?id=kQDw1ZcqLjUC&redir_esc=y)
- [The uses of random partitions - Andrei Okounkov](https://arxiv.org/abs/math-ph/0309015)

# Citation & Contributing

If you found this package useful in your research, please consider citing the companion 
paper available here: [arxiv.org/abs/2403.13580](https://arxiv.org/abs/2403.13580). 
```
@article{Araujo:2024piv,
    author = "Araujo, Thiago",
    title = "{PySymmPol: Symmetric Polynomials in Python}",
    eprint = "2403.13580",
    archivePrefix = "arXiv",
    primaryClass = "math.CO",
    month = "3",
    year = "2024"
}
```
Feeling like contributing? Fork the project and open a pull request with your modifications. 
Found a bug? Just open a GitHub issue.
