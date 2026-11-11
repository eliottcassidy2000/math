#!/usr/bin/env python3
"""exploration (not a deliverable): exact characteristic polynomials of T_k through integer companion blocks."""
import sys

import mpmath as mp
from flint import fmpz_mat, fmpz_poly
from sympy import cyclotomic_poly, Symbol, Poly, totient

mp.mp.dps = 60


def companion(m):
    z = Symbol("z")
    phi = Poly(cyclotomic_poly(m, z), z)
    d = phi.degree()
    co = [int(c) for c in phi.all_coeffs()][::-1]      # co[i] z^i, monic
    C = [[0] * d for _ in range(d)]
    for i in range(1, d):
        C[i][i - 1] = 1
    for i in range(d):
        C[i][d - 1] = -co[i]
    return C, d


def matmul(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def big_matrix(m, k):
    C, d = companion(m)
    pw = [[[int(i == j) for j in range(d)] for i in range(d)]]
    for _ in range(m - 1):
        pw.append(matmul(pw[-1], C))
    inv2 = pow(2, -1, m)
    n = 2 * m
    B = [[0] * (n * d) for _ in range(n * d)]
    for eta in range(m):
        for beta in range(2):
            i = 2 * eta + beta
            j0 = 2 * ((2 * eta) % m)
            j1 = 2 * ((2 * eta + 1) % m) + 1
            e = (k * ((eta + beta) * inv2 % m)) % m
            for a in range(d):
                B[i * d + a][j0 * d + a] += 1
                for b in range(d):
                    B[i * d + a][j1 * d + b] += pw[e][a][b]
    return B


def main():
    for m in [int(a) for a in sys.argv[1:]] or [3, 5, 7]:
        B = big_matrix(m, 1)
        M = fmpz_mat(B)
        cp = M.charpoly()
        print(f"m = {m}: integer block matrix {len(B)} x {len(B)}; norm char poly degree {cp.degree()}")
        c, facs = cp.factor()
        best = None
        for f, mult in facs:
            co = [int(x) for x in f.coeffs()]
            roots = mp.polyroots(co[::-1], maxsteps=2000, extraprec=400) if f.degree() > 0 else []
            mods = sorted((abs(r) for r in roots), reverse=True)
            top = mods[0] if mods else mp.mpf(0)
            print(f"  factor deg {f.degree():3d} mult {mult}: max |root| = {mp.nstr(top, 20)}   {f if f.degree() <= 24 else '(long)'}")
            if best is None or top > best[0]:
                best = (top, f)
        print(f"  dominant modulus {mp.nstr(best[0], 30)}  rate {mp.nstr(best[0] / 2, 20)}; minimal polynomial of the dominant eigenvalue:")
        print("   ", best[1])


if __name__ == "__main__":
    main()
