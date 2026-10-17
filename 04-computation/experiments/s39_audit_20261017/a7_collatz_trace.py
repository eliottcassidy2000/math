#!/usr/bin/env python3
"""S39 audit, item 7: the Collatz affine digraph modulo p (local_systems_dictionary_20261017.md section 4).

Independent of the session script: the right side of the trace formula is computed by a dynamic programme over the
states (j, c_w mod p) (c_w evolves as c -> c for the letter x/2 and c -> 3c + 2^i for the letter (3x+1)/2 at position i),
not by enumerating the 2^k words; the left side by exact integer matrix powers.  Also: Perron-Frobenius data, the
equidistribution numbers at k = 20, the arborescence counts (matrix-tree theorem), and the factorisation of the
characteristic polynomial on the nonconstant functions.
Usage: python -B a7_collatz_trace.py
"""
import numpy as np
import sympy
from sympy import primerange


def digraph(p):
    i2 = pow(2, -1, p)
    A = [[0] * p for _ in range(p)]
    for x in range(p):
        A[x][x * i2 % p] += 1
        A[x][(3 * x + 1) * i2 % p] += 1
    return A


def rhs_dp(p, k):
    """sum over words w of length k with 3^j = 2^k (mod p) of (p [p | c_w] - 1), plus 2^k; and (B_k, N_k)."""
    # state counts: dict (j, c mod p) -> number of words; after i letters c_w uses den 2^i
    st = {(0, 0): 1}
    for i in range(k):
        nst = {}
        for (j, c), m in st.items():
            nst[(j, c)] = nst.get((j, c), 0) + m                                     # letter x/2: c unchanged
            key = (j + 1, (3 * c + pow(2, i, p)) % p)
            nst[key] = nst.get(key, 0) + m                                          # letter (3x+1)/2: c -> 3c + 2^i
        st = nst
    B = Nk = 0
    for (j, c), m in st.items():
        if (pow(3, j, p) - pow(2, k, p)) % p == 0:
            B += m
            if c == 0:
                Nk += m
    return 2 ** k + p * Nk - B, B, Nk


bad = []
for p in primerange(5, 62):
    A = sympy.Matrix(digraph(p))
    Ak = sympy.eye(p)
    for k in range(1, 21):
        Ak = Ak * A
        if Ak.trace() != rhs_dp(p, k)[0]:
            bad.append((p, k))
print('trace formula tr(A_p^k) = 2^k + sum_(3^j = 2^k) (p[p | c_w] - 1): every prime 5 <= p < 62, k <= 20; failures', bad)
assert not bad

print('\nPerron-Frobenius and the equidistribution numbers:')
for p in (5, 7, 11, 13):
    A = np.array(digraph(p), dtype=float)
    ev = sorted(np.abs(np.linalg.eigvals(A)), reverse=True)
    rho = ev[1]
    t, B, Nk = rhs_dp(p, 20)
    bound = (p - 1) * rho ** 20 / p
    print('  p = %2d: rho = %.4f; k = 20: N/B = %d/%d = %.4f vs 1/p = %.4f; |N - B/p| = %.1f <= (p-1) rho^k / p = %.1f'
          % (p, rho, Nk, B, Nk / B if B else float('nan'), 1 / p, abs(Nk - B / p), bound))
    assert abs(Nk - B / p) <= bound + 1e-9
for p in (11, 13):
    zero = [k for k in range(1, 25) if rhs_dp(p, k)[1] == 0]
    print('  p = %d: B_k = 0 (no word with 3^j = 2^k mod p) for k in %s' % (p, zero[:12]))

print('\narborescences and the characteristic polynomial on the nonconstant functions:')
z = sympy.symbols('z')
for p in (5, 7, 11, 13, 17, 19, 23):
    A = sympy.Matrix(digraph(p))
    Lap = 2 * sympy.eye(p) - A
    tau = Lap[1:, 1:].det()
    cp = A.charpoly(z).as_expr()
    q, r = sympy.div(cp, z - 2, z)
    fl = sympy.factor_list(q, z)[1]
    print('  p = %2d: tau = %s; factors of charpoly/(z - 2): %s' % (p, tau, [(sympy.factor(f), m) for f, m in fl if sympy.degree(f, z) <= 4] +
                                                                  [('deg %d' % sympy.degree(f, z), m) for f, m in fl if sympy.degree(f, z) > 4]))
print('ALL CHECKS PASSED')
