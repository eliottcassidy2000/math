#!/usr/bin/env python3
"""Check exact moments against exhaustive directed signed half-path counts for m <= 4."""
from fractions import Fraction
from itertools import combinations, permutations, product
from math import factorial


def formula(m):
    total = Fraction(0)
    for r in range(m // 2 + 1):
        k = m - 2*r
        p = [Fraction(1)]+[Fraction(0)]*r
        for _ in range(k+1):
            p = [p[j]+sum(2*p[j-i]/4**i for i in range(1,j+1))
                 for j in range(r+1)]
        total += factorial(k)*factorial(m)*p[r]
    return total


def enumerate_signs(m):
    # Fix the terminal signs to 0: independent gauge flips preserve the
    # distribution of all cross-arc patterns for this endpoint condition.
    orders = list(permutations(range(m)))
    pairs = list(combinations(range(m), 2))
    s = s2 = n = 0
    for cross in product(range(4), repeat=len(pairs)):
        A = [[0] * (2*m) for _ in range(2*m)]
        for (i,j), code in zip(pairs,cross):
            a,b = code & 1, code >> 1
            for x,y,w in ((0,0,a),(0,1,b),(1,0,1-b),(1,1,1-a)):
                A[2*i+x][2*j+y] = w
                A[2*j+y][2*i+x] = 1-w
        count = 0
        for order in orders:
            for signs in product(range(2),repeat=m-1):
                path = [2*v+x for v,x in zip(order,signs+(0,))]
                count += all(A[a][b] for a,b in zip(path,path[1:]))
        n += 1
        s += count
        s2 += count*count
    return Fraction(s,n),Fraction(s2,n)


if __name__ == '__main__':
    for m in range(1,5):
        mean, second = enumerate_signs(m)
        if mean != factorial(m) or second != formula(m):
            raise RuntimeError((m,mean,second))
        print(m, mean, second)
