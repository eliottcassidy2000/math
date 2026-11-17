#!/usr/bin/env python3
"""lattices2.py (H-2, H-6): delayed-choice lattices.  At order n the diagonal polynomial P_q = q^n + sum_{m in M} a_m q^m
with M = {m = n (mod 2), 1 <= m <= n-2}.  Two tournaments with equal P_1..P_k have coefficient differences Delta a in the
integer kernel of the rows (c^m)_{m in M}, c = 1..k, whose Z-basis is the band vectors with entries the coefficients of
prod_{c=1}^k (X - c^2) placed on k+1 consecutive elements of M (triangular, leading entry 1).  Hence the lattice of
Delta P_{k+1} is (k+1)^{m_0} prod_c ((k+1)^2 - c^2) Z, m_0 = min M.  This script builds the band basis, computes the gcd
of its images, compares with the formula, and brute-forces small kernel vectors as an independent check."""
import itertools, math
from sympy import symbols, Poly, prod
X = symbols('X')
for k in (2, 3, 4):
    band = Poly(prod(X - c * c for c in range(1, k + 1)), X).all_coeffs()[::-1]   # increasing powers
    C = sum(v * (k + 1) ** (2 * j) for j, v in enumerate(band))
    print(f"k = {k}: band vector {band}, constant prod((k+1)^2 - c^2) = {C}")
    for n in range(5, 15):
        M = list(range(1 if n % 2 else 2, n - 1, 2))
        if len(M) < k + 1: print(f"   n = {n}: fewer than {k+1} free coefficients, no delayed choice at level {k+1}"); continue
        m0 = M[0]; g_basis = 0
        for j in range(len(M) - k):
            img = sum(v * (k + 1) ** M[j + i] for i, v in enumerate(band)); g_basis = math.gcd(g_basis, img)
        formula = (k + 1) ** m0 * C
        # brute force: small integer vectors in the kernel
        g_bf = 0; rng = range(-4, 5)
        if len(M) <= 6:
            for w in itertools.product(rng, repeat=len(M)):
                if all(sum(x * c ** m for x, m in zip(w, M)) == 0 for c in range(1, k + 1)):
                    g_bf = math.gcd(g_bf, sum(x * (k + 1) ** m for x, m in zip(w, M)))
        print(f"   n = {n}: M = {M}; band-basis lattice {g_basis} Z; formula {formula}; brute force (|entries| <= 4): {g_bf if len(M) <= 6 else 'skipped'}")
