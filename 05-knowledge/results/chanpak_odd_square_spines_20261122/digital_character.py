"""Exact mod-four bias and its six-step radical on the two E61 charts."""
from itertools import product

from rational_carries import chart, address_height, TARGET
from fermat_receiver import p as read


def advance(state):
    a, b, c, d = state
    return (a+c, a-c, b+d, d-b)


def polar_nullity(n):
    # Rank of the binary adjacency matrix of the square of the n-vertex path.
    basis = {}
    for j in range(n):
        row = sum(1<<i for i in (j-2,j-1,j+1,j+2) if 0 <= i < n)
        while row:
            bit = row.bit_length()-1
            if bit not in basis:
                basis[bit] = row
                break
            row ^= basis[bit]
    return n-len(basis)


def bias(n):
    table = (1,0,-2,0,0,-8,-8,0,-16,-32,0,0)
    return (-64)**(n//12)*table[n%12]


def check():
    for j in range(4):
        state = tuple(int(i==j) for i in range(4))
        base = state
        for _ in range(12):
            state = advance(state)
        assert state == tuple(-64*s for s in base)
    for n in range(1,121):
        assert polar_nullity(n) == (0,1,0,1,2,1)[n%6]
    for n in range(17):
        assert bias(n) == sum((-1)**((read(x)//2+x)%2) for x in range(1<<n))
    for spine in (4,8):
        data = chart(spine)
        sign = (-1)**(((TARGET-data[-1])//2)%2)
        for m,n in product(range(7),repeat=2):
            expected = ((1<<(m+n))+sign*bias(m)*bias(n))//2
            actual = sum((address_height(data,x,y)-TARGET)%4==0
                         for x,y in product(range(1<<m),range(1<<n)))
            assert actual == expected
        m,n=61,61
        print("spine",spine,"full-width addresses 0..2^61-1 in each coordinate:",
              ((1<<(m+n))+sign*bias(m)*bias(n))//2,"target residues mod 4")
    print("Matrix identity, nullities, one-word characters and two-word heights passed.")


if __name__ == "__main__":
    check()
