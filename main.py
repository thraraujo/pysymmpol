import pysymmpol as sy
import pysymmpol.utils as ut

def main():

    a = sy.YoungDiagram((1,1,1))
    a.draw_diagram()
    print(a.boxes)

    b = sy.State(4)
    print(f'conjugacy states: {b.conjugacy_states()}')
    print(f'parition states: {b.partition_states()}')


    t = ut.create_miwa(a.boxes)
    sc = sy.SchurPolynomial(a)
    print(sc.explicit(t))

if __name__ == "__main__":main()
