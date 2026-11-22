#!/usr/bin/env python3
"""Arithmetic identities and short literal checks for the two-port inverse."""
from itertools import permutations, product
from fractions import Fraction

A = ((7,1,2,1),(4,2,0,1),(4,2,1,2),(5,1,1,1))
B = ((6,3,1,4),(3,0,1,2),(4,0,1,2),(5,2,1,2))
I = tuple(tuple(int(i == j) for j in range(4)) for i in range(4))
M = (A, B)
C = (1,0,1,0)
V = (-1,5,6,-6)
SEED = (2,0,1,0)


def dot(x,y):
    return sum(a*b for a,b in zip(x,y))


def row(x,a):
    return tuple(dot(x,y) for y in zip(*a))


def col(a,y):
    return tuple(dot(x,y) for x in a)


def mat(a,b):
    return tuple(row(x,b) for x in a)


def determinant(a):
    total = 0
    for p in permutations(range(len(a))):
        sign = (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
        z = sign
        for i,j in enumerate(p):
            z *= a[i][j]
        total += z
    return total


def mpow(a,n):
    out = I
    while n:
        if n % 2:
            out = mat(out,a)
        a = mat(a,a)
        n //= 2
    return out


def inv(a):
    d = determinant(a)
    def minor(i,j):
        return tuple(tuple(x for l,x in enumerate(w) if l != j)
                     for k,w in enumerate(a) if k != i)
    return tuple(tuple(Fraction((-1)**(i+j)*determinant(minor(j,i)),d)
                       for j in range(4)) for i in range(4))


def f(z):
    return z**4-11*z**3+13*z*z-15*z+4


def g(z):
    return (z**3-2*z*z-3*z+2, 4*z*z+6*z-6,
            z**3-6*z*z+19*z-8, 6*z*z-6*z+4)


def need(condition):
    if not condition:
        raise ArithmeticError("two-port certificate mismatch")


def poly_mod3_remainder(p,q):
    p = list(p)
    for k in range(len(p)-len(q),-1,-1):
        lead = p[k+len(q)-1] % 3
        for j,x in enumerate(q):
            p[k+j] = (p[k+j]-lead*x) % 3
    return p[:len(q)-1]


def literal(word):
    # Published root-child controller: entry 1, return words 041/403.
    tape, flip = [], 0
    for macro in ((1,),) + tuple(((0,4,1),(4,0,3))[b] for b in word):
        for a in macro:
            root, drop = divmod(a,3)
            tape.append(3*(root^flip)+(drop^flip if drop < 2 else 2))
        if macro[-1] == 3:
            flip ^= 1
    T = (2,4,1)
    ports = [0,1]
    for a in tape[:-1]:
        root, drop = divmod(a,3)
        v, n = ports[root], len(T)
        T = tuple(x | ((1<<n) if i == v else 0) for i,x in enumerate(T))
        T += (((1<<n)-1) ^ (1<<v),)
        if drop < 2:
            ports[drop] = n
    return T


def hamilton_count(T):
    n = len(T)
    dp = [[0]*n for _ in range(1<<n)]
    for i in range(n):
        dp[1<<i][i] = 1
    for S in range(1,1<<n):
        for i, value in enumerate(dp[S]):
            if value:
                flags = T[i] & ~S
                while flags:
                    b = flags & -flags
                    flags ^= b
                    dp[S|b][b.bit_length()-1] += value
    return sum(dp[-1])


def counts_check(m):
    suffix = col(mpow(A,100+16*m),C)
    low, high = min(suffix), 3*14**m*max(suffix)
    Ai, Bi = inv(A), inv(B)
    seen = set()
    for w in product((0,1),repeat=m):
        x = SEED
        for b in w:
            x = row(x,M[b])
        h = dot(x,suffix)
        need(h not in seen and low <= h <= high)
        seen.add(h)
        if m <= 2:
            need(hamilton_count(literal(w)) == dot(x,C))
            need(hamilton_count(literal(w+(0,))) == dot(row(x,A),C))
        for b in reversed(w):
            actual = int(dot(x,V)<0)
            need(actual == b)
            x = row(x, (Ai,Bi)[actual])
        need(x == SEED)
    return high


def main():
    need(determinant(A)==4 and determinant(B)==-4)
    need(col(A,V)==(4,0,0,0) and col(B,V)==(-9,-9,-10,-1))
    need(determinant(tuple(zip(*(col(mpow(A,i),C) for i in range(4)))))==496)
    # Both sides are polynomial of degree at most four.
    for z in range(5):
        test = tuple(tuple(z*I[i][j]-A[i][j] for j in range(4)) for i in range(4))
        need(determinant(test)==f(z))
        need(col(test,g(z))==tuple(f(z)*x for x in C))
    need(all(f(z)%3 for z in range(3)))
    need(all(any(poly_mod3_remainder((1,0,1,1,1),q))
             for q in ((1,0,1),(2,1,1),(2,2,1))))
    last = 0
    for m in (1,2,4,8):
        need(min(col(mpow(A,100+16*m),C))>last)
        last = counts_check(m)
    print("matrix constants; 12 short literal count checks; "
          "distinct/inverted states and counts at widths 1,2,4,8")


if __name__ == "__main__":
    main()
