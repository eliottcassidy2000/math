#!/usr/bin/env python3
"""minrec.py -- minimal linear recurrences of the constant-letter run sequences H(a^n) of span 2 (Berlekamp-
Massey over Q on 60 terms, verified on 120 terms) and the dominant root of each; the recurrence's
characteristic polynomial must divide the letter's characteristic polynomial on the observable space."""
import numpy as np, sympy as sp
from fractions import Fraction
from ribbon import transfer_system, letters, count_with_matrices
w = 2
sysw = transfer_system(w)
def bm(seq):
    """Berlekamp-Massey over Q: returns connection polynomial coefficients c (seq[n] = sum c[i] seq[n-i-1])."""
    C = [Fraction(1)]; B = [Fraction(1)]; L = 0; m = 1; b = Fraction(1)
    for n in range(len(seq)):
        d = Fraction(seq[n]) + sum(C[i] * seq[n - i] for i in range(1, L + 1))
        if d == 0:
            m += 1; continue
        T = C[:]
        coef = d / b
        C = C + [Fraction(0)] * (len(B) + m - len(C))
        for i in range(len(B)):
            C[i + m] -= coef * B[i]
        if 2 * L <= n:
            L = n + 1 - L; B = T; b = d; m = 1
        else:
            m += 1
    return L, [-c for c in C[1:L + 1]]
names = {frozenset([1]): '{1}', frozenset([2]): '{2}', frozenset([1, 2]): '{1,2}'}
for a in letters(w):
    if not a: continue
    seq = [count_with_matrices(w, [a] * n, sysw) for n in range(0, 120)]
    L, c = bm(seq[:60])
    ok = all(seq[n] == sum(c[i] * seq[n - i - 1] for i in range(L)) for n in range(L, 120))
    x = sp.symbols('x')
    poly = sp.expand(x**L - sum(sp.Rational(c[i].numerator, c[i].denominator) * x**(L - i - 1) for i in range(L)))
    roots = np.roots([float(cc) for cc in sp.Poly(poly, x).all_coeffs()])
    rho = max(abs(roots))
    print(f"letter {names[a]}: minimal recurrence order {L} (valid for all n >= {L}, checked to n = 119: {ok}); characteristic polynomial {sp.factor(poly)}; dominant root {rho:.6f}")
