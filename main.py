from pysymmpol import YoungDiagram, ConjugacyClass
from pysymmpol import SchurPolynomial
from pysymmpol import MonomialPolynomial
from pysymmpol import HallLittlewoodPolynomial
from pysymmpol import HomogeneousPolynomial, ElementaryPolynomial
from pysymmpol.utils import create_miwa, create_x_coord
from sympy import Symbol

def main():
    # YOUNG  DIAGRAMS
    young = YoungDiagram((3,2,1))
    conjugacy = ConjugacyClass({1: 1, 2: 1, 3: 1})

    young.draw_diagram(4)
    conjugacy.draw_diagram(4)


    # HOMOGENEOUS AND ELEMENTARY POLYNOMIALS
    n = 3
    t = create_miwa(n)

    homogeneous = HomogeneousPolynomial(n)
    elementary = ElementaryPolynomial(n)
    print(f"=> homogeneous: {homogeneous.explicit(t)}")
    print(f"=> elementary: {elementary.explicit(t)}")

    # SCHUR POLYNOMIALS
    t = create_miwa(young.boxes)

    schur = SchurPolynomial(young)

    print(f"=> schur: {schur.explicit(t)}")

    n = 3
    x = create_x_coord(n)

    monomial = MonomialPolynomial(young)

    print(f"=> monomial: {monomial.explicit(x)}")

    Q = Symbol('Q')
    young = YoungDiagram((3,2,1))

    n = 3
    x = create_x_coord(n)

    hall_littlewood = HallLittlewoodPolynomial(young)

    print(f"=> hall-littlewood: {hall_littlewood.explicit(x, Q)}")






if __name__=="__main__": 
    main() 
