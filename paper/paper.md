---
title: 'PySymmPol - Symmetric Polynomials'
tags:
  - Python
  - Combinatorics
  - Physics
authors:
  - name: Thiago Araujo
    orcid: 0000-0001-5792-2530
    affiliation: "1, 2"
affiliations:
 - name: Institute for Theoretical Physics, UNESP, R. Dr. Bento Teobaldo Ferraz, 271, Bloco II,  Barra-Funda, CEP 01140-070, São Paulo/SP, Brazil.
   index: 1
 - name: Institute of Physics, São Paulo University, Rua do Matão 1371 - CEP 05508-090 Cidade Universitária, São Paulo/SP, Brazil
   index: 2

date: 21 March 2024
bibliography: paper.bib
---

# Summary

**PySymmPol** is a Python package designed for efficient manipulation of
symmetric polynomials. It provides functionalities for working with
various types of symmetric polynomials, including elementary,
homogeneous, monomial symmetric, (skew-) Schur, and Hall-Littlewood
polynomials. In addition to polynomial operations, **PySymmPol** offers
tools to explore key properties of integer partitions and Young
diagrams, such as transposition, Frobenius coordinates, characters of
symmetric groups and others. 

This package originated from research conducted in the field of
integrable systems applied to string theory, and the AdS/CFT
correspondence. **PySymmPol** aims to facilitate computational tasks
related to symmetric polynomials and their applications in diverse
fields.

# Statement of need 

Symmetric polynomials play a crucial role across various domains of
mathematics and theoretical physics due to their rich structure and
broad applications. They arise naturally in combinatorics, algebraic
geometry, representation theory, and mathematical
physics [@Macdonald:1998], [@Fulton:2004], [@Babelon:2003],
and [@Wheeler:2010]. These polynomials encode
essential information about symmetries and patterns, making them
indispensable in the study of symmetric functions and their
connections to diverse mathematical structures. Moreover, symmetric
polynomials find extensive applications in theoretical physics,
particularly in quantum mechanics, statistical mechanics, and quantum
field theory. Their utility extends to areas such as algebraic
combinatorics, where they serve as powerful tools for solving
combinatorial problems and understanding intricate relationships
between different mathematical objects. Thus, tools and libraries like
PySymmPol provide researchers and practitioners with efficient means
to explore and manipulate symmetric polynomials, facilitating
advancements in both theoretical studies and practical applications.

# Core features and functionalities 

**PySymmPol** is composed of several modules in two main packages. 

1. Partitions
2. Polynomials

Our keen interest in exploring these representations stems
from their relevance in the realm of two-dimensional Conformal Field
Theories (CFTs). In this context, bosonic states are labeled by
conjugacy class vectors, highlighting the importance of
understanding and manipulating such representations. Conversely,
fermions are represented using the standard representation,
emphasizing the need to bridge the conceptual gap between these
distinct frameworks. Investigating these representations not only
sheds light on the fundamental properties of fermion states within
CFTs but also provides valuable insights for theoretical developments
in quantum field theory and related fields.

Further details on the other functionalities can be found in the tutorial.
The *ACCELASC* algorithm [@Kelleher:2009] greatly improved the
speed of the methods associated to these calculations. 

One of the main goals of the package is to provide the definitions of the 
polynomials in terms of the Miwa coordinates, or power sums, 
$$ t_j = \frac{1}{j} \sum_{i=1}^N x_i^j $$ where $\vec{x} = 
(x_1, \dots, x_N)$ and $\mathbf{t} = (t_1, t_2, \dots)$.

# Acknowledgements

I would like to thank Fapesp, São Paulo Research Foundation, for
finantial support, grant \textbf{2022/06599-0}. I would also like to
thank the Free and Open Source Software community.

# References
