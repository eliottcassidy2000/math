#!/usr/bin/env python3
"""Direct tiny tournament controls, plus exact polynomial lift/one-read inverse."""
from fractions import Fraction
from itertools import permutations, product
from math import comb, factorial


def require(test, label):
    if not test:
        raise RuntimeError(label)


def port_pair(g, active):
    n = len(g)
    x, y = n, n+1
    return (tuple(a | (1 << y) for a in g)
            + (((1 << n)-1) | ((not active) << y), active << x))


def h(g):
    n = len(g)
    T = [[0]*n for _ in range(1 << n)]
    for i in range(n):
        T[1 << i][i] = 1
    for S in range(1, 1 << n):
        for i, v in enumerate(T[S]):
            if not v:
                continue
            J = g[i] & ~S
            while J:
                j = J & -J
                T[S | j][j.bit_length()-1] += v
                J -= j
    return sum(T[-1])


def op(q, active):
    z = [0, 0]+q
    if active:
        for d, a in enumerate(q):
            for j in range(1, d+1, 2):
                z[d-j+1] += 2*a*comb(d,j)
    return z


def inverse_op(z, active):
    c = z[:]
    q = [0]*(len(z)-2)
    for j in range(len(q)-1, -1, -1):
        q[j] = c[j+2]
        if active:
            for k in range(1, j+1, 2):
                c[j-k+1] -= 2*q[j]*comb(j,k)
    require(op(q,active)==z,"inverse operator")
    return q


def wrapper_weights(R, t, k):
    v = [0, 1] + [0]*(R-1)
    for _ in range(k):
        w = [0]*(R+1)
        for r in range(1, R+1):
            for s, c in ((r,r*r),(r+1,r),(r-1,-r)):
                if 1 <= s <= R:
                    w[s] += v[r]*c*pow(s,2*t)
        v = w
    return v[1:]


def bin_ceil(a):
    return (a-1).bit_length()


def read(q, r):
    return sum(a*pow(r,i) for i,a in enumerate(q))


def nearest(x):
    return (2*x.numerator+x.denominator)//(2*x.denominator)


def genetic_check():
    checked = 0
    for width in range(1, 4):
        for word in product((0,1),repeat=width):
            q = [0,0,1]
            g = (2,0)
            for b in word:
                q = op(q,b)
                g = port_pair(g,b)
            # Read one active wrapper with a short literal padding. This
            # confirms operator composition on graphs, without checking
            # huge DP state spaces or treating this as the robust decoder.
            gp = port_pair(port_pair(g,0),1)
            require(h(g)==read(q,1),"literal seed")
            require(h(gp)==read(q,1)+4*read(q,2),"literal wrapper")
            checked += 1
    return checked


def decode_check():
    checked = 0
    for width in range(1, 9):
        n = 2+2*width
        R = (n+1)//2
        k = R-1
        C = pow(3,k)*pow(R,2*k)
        E = factorial(n)*pow(R,n)
        t = R*bin_ceil(8*C*E)
        ws = wrapper_weights(R,t,k)
        out = set()
        for word in product((0,1),repeat=width):
            q = [0,0,1]
            for b in word:
                q = op(q,b)
            ev = [read(q,x) for x in range(1,R+1)]
            N = sum(c*x for c,x in zip(ws,ev))
            out.add(N)
            for j in range(R-1,-1,-1):
                x = nearest(Fraction(N,ws[j]))
                require(x==ev[j],"single read decoder")
                N -= ws[j]*x
            require(N==0,"decoder remainder")
            for b in reversed(word):
                bb = q[2] != 0
                require(bb==bool(b),"low-degree digit")
                q = inverse_op(q,bb)
            require(q==[0,0,1],"recovered seed")
            checked += 1
        require(len(out)==1<<width,"actual count collision")
    return checked


def fixed_odd(g, sigma):
    c = next(i for i,j in enumerate(sigma) if i==j)
    pairs = [(i,j) for i,j in enumerate(sigma) if i<j]
    return sum(all(g[x]>>y&1 for x,y in zip(a,a[1:]+(c,)))
               for representatives in product(*pairs)
               for a in permutations(representatives))


def odd_controls():
    checked = 0
    for base, seed in [((6,4,0),1),((2,4,1),3)]:
        for word in product((0,1),repeat=3):
            g, sigma = base, (2,1,0)
            for b in word:
                x=len(g)
                g=port_pair(g,b)
                sigma+=(x+1,x)
            require(fixed_odd(g,sigma)==1,"fixed-channel preservation")
            require(h(g)%4==seed,"residue-sheet preservation")
            checked += 1
    return checked


if __name__ == "__main__":
    print("literal tournament configurations:",genetic_check())
    print("weighted readout and branch inverses:",decode_check())
    print("odd fixed-channel controls:",odd_controls())
