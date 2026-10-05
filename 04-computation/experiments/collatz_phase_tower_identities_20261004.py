#!/usr/bin/env python3
"""
Exact identities of the window phases of the q-adic Syracuse frequency recursion (opus, 2026-10-04; Fourier note 4j).
theta_{n,d} := (u 2^-d mod q^n)/q^n in Q/Z (the phase of the recursion at level n, depth d, unit u).
 (i)   q-th root:  q theta_{n+1,d} = theta_{n,d}            (consecutive levels)
 (ii)  odometer:   theta_{n,d-1}   = 2 theta_{n,d}           (consecutive depths; the inverse doubling map)
 (iii) Pascal tower (q = 3 = 2 + 1):  theta_{N-m,d} = sum_i C(m,i) theta_{N,d-i}   -- every level is a binomial
       transform of the deepest level's phases; for q = 2^j + 1 the same with step j: sum_i C(m,i) theta_{N,d-ji}.
All checked exactly in rational arithmetic (u = 1, 5, 7, 11; n <= 12; d <= 24; m <= 6; and q = 5 with step 2).
"""
from fractions import Fraction as F
from math import comb

def theta(u, n, d, q=3):
    mod = q ** n
    return F((u * pow(2, -d, mod)) % mod, mod)

def fr(x):
    return x - (x.numerator // x.denominator)

if __name__ == "__main__":
    ok = True
    for u in (1, 5, 7, 11):
        for n in range(2, 12):
            for d in range(1, 25):
                if fr(3 * theta(u, n + 1, d)) != theta(u, n, d): ok = False
                if d >= 2 and fr(2 * theta(u, n, d)) != theta(u, n, d - 1): ok = False
        N = 14
        for m in range(1, 7):
            for d in range(m + 1, 25):
                if theta(u, N - m, d) != fr(sum(comb(m, i) * theta(u, N, d - i) for i in range(m + 1))): ok = False
    print("q = 3: (i) cube root, (ii) odometer, (iii) Pascal tower:", "ALL OK" if ok else "FAIL")
    ok5 = True
    for u in (1, 3):
        N = 10
        for m in range(1, 5):
            for d in range(2 * m + 1, 25):
                if theta(u, N - m, d, 5) != fr(sum(comb(m, i) * theta(u, N, d - 2 * i, 5) for i in range(m + 1))): ok5 = False
    print("q = 5: stretched Pascal tower (step 2):", "ALL OK" if ok5 else "FAIL")
