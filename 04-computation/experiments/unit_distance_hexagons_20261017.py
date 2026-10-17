#!/usr/bin/env python3
"""Hexagons (6-cycles) in the finite unit-distance graphs G_p = Cay(F_(p^2), mu_(p+1)), p prime, by direct count (S39).

For each prime 5 <= p <= PMAX: the number c_0 of 6-cycles through the vertex 0 by depth-first search (vertex-transitivity
gives C_6 = p^2 c_0 / 6), the closed 6-walk count W_6 from the exact Kloosterman moments, and the number of non-simple
closed 6-walks per vertex D = W_6 - 2 c_0, compared with the derived count

    D(p) = (p+1)(11p^2 - 2p + 1)    (p = 1 mod 3; G_p triangle-free)
    D(p) = (p+1)(11p^2 + 10p - 15)  (p = 2 mod 3; every vertex in p+1 triangles)

(tree-like walks d(5d^2 - 6d + 2), 4-cycles with one inserted backtrack (p^2 - 1)(5d - 4), lollipops (p+1)(p-1)^2,
and for p = 2 mod 3 the triangle shapes 4p(p+1) + 8(p+1)(p-1) minus the 8(p+1) walks counted twice; d = p + 1),
and the closed form C_6 = p^2 (W_6 - D) / 12 with W_6 = p^4 + 11p^3 + 10p^2 + p + 1 - (p+1) b_p.
Usage: python unit_distance_hexagons_20261017.py moments.json PMAX
Session opus-2026-10-17-S39.
"""
import sys, json
import numpy as np
from numba import njit
from sympy import jacobi_symbol, primerange


def kr(d, n):
    return int(jacobi_symbol(d % n, n)) if d % n else 0


@njit(cache=True)
def count_c6_through0(nbr, adj0):
    deg = nbr.shape[1]
    cnt = 0
    for i1 in range(deg):
        a = nbr[0, i1]
        for i2 in range(deg):
            b = nbr[a, i2]
            if b == 0:
                continue
            for i3 in range(deg):
                c = nbr[b, i3]
                if c == 0 or c == a:
                    continue
                for i4 in range(deg):
                    d = nbr[c, i4]
                    if d == 0 or d == a or d == b:
                        continue
                    for i5 in range(deg):
                        e = nbr[d, i5]
                        if e == 0 or e == a or e == b or e == c:
                            continue
                        if adj0[e]:
                            cnt += 1
    return cnt // 2


def graph(p):
    n = next(x for x in range(2, p) if kr(x, p) == -1)
    steps = [(x, y) for x in range(p) for y in range(p) if (x * x - n * y * y) % p == 1]
    assert len(steps) == p + 1
    nbr = np.zeros((p * p, p + 1), dtype=np.int64)
    for x in range(p):
        for y in range(p):
            for j, (dx, dy) in enumerate(steps):
                nbr[x * p + y, j] = ((x + dx) % p) * p + (y + dy) % p
    return nbr


def eta46(N):
    s = np.zeros(N + 1, dtype=object); s[0] = 1
    for m in (1, 2, 3, 6):
        for n in range(1, N // m + 1):
            for _ in range(2):
                s[m * n:] = s[m * n:] - s[:N + 1 - m * n].copy()
    return s


if __name__ == '__main__':
    M = {int(k): v for k, v in json.load(open(sys.argv[1]))['moments'].items()}
    PMAX = int(sys.argv[2])
    P46 = eta46(PMAX + 1)
    print(' p  (p/3)          C_6        W_6          D   D derived   C_6 closed form')
    for p in primerange(5, PMAX + 1):
        nbr = graph(p)
        adj0 = np.zeros(p * p, dtype=np.bool_); adj0[nbr[0]] = True
        c0 = count_c6_through0(nbr, adj0)
        assert (p * p * c0) % 6 == 0
        C6 = p * p * c0 // 6
        W6 = ((p + 1) ** 6 + (p + 1) * M[p][5]) // (p * p)
        D = W6 - 2 * c0
        Dd = (p + 1) * (11 * p * p - 2 * p + 1) if p % 3 == 1 else (p + 1) * (11 * p * p + 10 * p - 15)
        b = int(P46[p - 1])
        W6c = p ** 4 + 11 * p ** 3 + 10 * p * p + p + 1 - (p + 1) * b
        C6c = p * p * (W6c - Dd) // 12
        print('%3d  %2d  %14d %10d %10d %10d %16d' % (p, kr(p, 3), C6, W6, D, Dd, C6c), flush=True)
        assert W6 == W6c and D == Dd and C6 == C6c, p
    print('ALL CHECKS PASSED')
