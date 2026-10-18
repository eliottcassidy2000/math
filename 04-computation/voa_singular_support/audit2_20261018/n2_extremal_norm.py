#!/usr/bin/env python3
"""<u_n,u_n> for u_n = G^-_{-n+1/2} ... G^-_{-3/2}|0> (n odd modes) as a polynomial in c: interpolate mod P from
values at many c (degree <= n), then reconstruct rational coefficients and factor with sympy."""
import sympy as sp
from fractions import Fraction
from aud_algebras import N2, GM
import math

P = 2147483647
def norm(c, n):
    va = N2(c, P)
    mono = tuple(sorted([(GM, -(2 * i + 3)) for i in range(n)]))
    vec = {mono: 1}
    for X in [va.omega(X) for X in mono]:
        vec = va.apply_vec(X, vec)
    return vec.get((), 0)

def ratrec(a):
    r0, r1, s0, s1 = P, a % P, 0, 1
    B = math.isqrt(P // 2)
    while r1 > B:
        qq = r0 // r1; r0, r1 = r1, r0 - qq * r1; s0, s1 = s1, s0 - qq * s1
    return Fraction(r1, s1)

C = sp.symbols('c')
for n in range(1, 6):
    pts = list(range(2, n + 4))
    vals = [norm(Fraction(t), n) for t in pts]
    # Lagrange interpolation mod P -> coefficients
    coeffs = [0] * (len(pts))
    poly = [0]
    import itertools
    def polymul(a, b):
        out = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i + j] = (out[i + j] + x * y) % P
        return out
    tot = [0] * len(pts)
    for i, xi in enumerate(pts):
        num = [1]; den = 1
        for j, xj in enumerate(pts):
            if j != i:
                num = polymul(num, [(-xj) % P, 1]); den = den * (xi - xj) % P
        f = vals[i] * pow(den, P - 2, P) % P
        for k_, a in enumerate(num):
            tot[k_] = (tot[k_] + a * f) % P
    expr = sum(sp.Rational(ratrec(a).numerator, ratrec(a).denominator) * C ** k_ for k_, a in enumerate(tot))
    print(f'n={n}: <u,u> =', sp.factor(expr))
