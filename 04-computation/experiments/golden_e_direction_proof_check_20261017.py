#!/usr/bin/env python3
"""Exact check of the last case in the proof that clock_duality's direction e is Cartesian on K + uK + vK + eK (S41).

Case v-coefficient zero: a unit a + u b + e g (a, b, g in K, g != 0) forces c0 = a + u b = t Z g with t real and
|g|^2 (1 + t^2 N) = 1, hence (a + u b)^2 (1 + conj(d) u) = kappa' u (d + u) with kappa' = kappa g^2 in K
(conj(u) = 1/u). Reduce modulo the minimal polynomial u^3 = m2 u^2 + m1 u + m0 of u over Q(j) and compare coefficients:
  [1]   a^2 + b^2 dbar m0 = 0
  [u]   2ab + a^2 dbar + b^2 dbar m1 - kappa' d = 0
  [u^2] b^2 + 2ab dbar + b^2 dbar m2 - kappa' = 0.
b = 0 forces a = 0, kappa' = 0, |g| = 1 (pure axis). For b != 0 put a = lam b; [u^2] gives kappa'. [u] is then
quadratic in lam, and becomes linear after substituting lam^2 = -dbar m0 from [1] (2 - 2|d|^2 != 0). Consistency
requires lam^2 = -dbar m0. This script computes m0, m1, m2 exactly over Q(j) and
tests that consistency condition exactly in Q(zeta15, j) (sympy, minimal polynomials).
Usage: python golden_e_direction_proof_check_20261017.py
Session opus-2026-10-17-S41.
"""
import sympy as sp

J, U, W, Zt = sp.symbols('j u w zeta')
# u = (w + 3 - j)/(w + 3 + j)  =>  w = (3 - j - (3 + j) u)/(u - 1)
w_of_u = (3 - J - (3 + J) * U) / (U - 1)
cubic = sp.together(w_of_u ** 3 + 2 * w_of_u ** 2 - 4 * w_of_u - 4)
num = sp.expand(sp.numer(cubic))
num = sp.expand(num.subs(J ** 2, -11))
num = sp.Poly(sp.expand(num), U)
# reduce powers of j: replace j^2 -> -11 repeatedly
def redj(e):
    e = sp.expand(e)
    for _ in range(5):
        e = sp.expand(e.subs(J ** 4, 121).subs(J ** 3, -11 * J).subs(J ** 2, -11))
    return e
coeffs = [redj(c) for c in num.all_coeffs()]
lead = coeffs[0]
print('cubic in u over Q(j): degree', len(coeffs) - 1)
# monic: u^3 + c2/lead u^2 + c1/lead u + c0/lead = 0  ->  u^3 = m2 u^2 + m1 u + m0
def div(x, y):
    # x / y in Q(j): multiply by conjugate
    yc = y.subs(J, -J)
    return redj(sp.expand(redj(x * yc)) / redj(sp.expand(y * yc)))
m2 = redj(-div(coeffs[1], lead)); m1 = redj(-div(coeffs[2], lead)); m0 = redj(-div(coeffs[3], lead))
print('m2 =', sp.nsimplify(m2), '\nm1 =', sp.nsimplify(m1), '\nm0 =', sp.nsimplify(m0))
# numeric sanity check at the intended embedding: j = i sqrt(11), w the real root in (1, 2)
wr = [r for r in sp.Poly(W ** 3 + 2 * W ** 2 - 4 * W - 4).nroots() if abs(sp.im(r)) < 1e-12 and 1 < sp.re(r) < 2][0]
jn = sp.I * sp.sqrt(11)
un = complex(sp.N(((wr + 3 - jn) / (wr + 3 + jn)), 30))
lhs = un ** 3; rhs = complex(sp.N((m2 * U ** 2 + m1 * U + m0).subs({J: jn, U: un}), 30))
print('|u^3 - (m2 u^2 + m1 u + m0)| at the embedding:', abs(lhs - rhs))
# d = 1 + zeta5, dbar = 1 + zeta5^-1; test the consistency condition numerically at all embeddings of zeta5 and both j
import cmath
worst = None
for k in (1, 2, 3, 4):
    z5 = cmath.exp(2j * cmath.pi * k / 5)
    for js in (jn, -jn):
        mm0 = complex(sp.N(m0.subs(J, js), 40)); mm1 = complex(sp.N(m1.subs(J, js), 40)); mm2 = complex(sp.N(m2.subs(J, js), 40))
        d = 1 + z5; db = 1 + 1 / z5; dd = d * db
        lam = (db ** 2 * mm0 - db * mm1 + d + dd * mm2) / (2 - 2 * dd)
        val = lam ** 2 + db * mm0
        print('zeta5^%d, j = %s: lam^2 + dbar m0 = %.6g %+.6gi' % (k, '+' if js == jn else '-', val.real, val.imag))
        worst = abs(val) if worst is None else min(worst, abs(val))
print('minimum |lam^2 + dbar m0| over embeddings: %.6g  (nonzero => no unit with b != 0; the case is closed)' % worst)

# Exact version: lam^2 + dbar m0 as an element of Q(zeta5, j), reduced modulo zeta5^4 + zeta5^3 + zeta5^2 + zeta5 + 1 and
# j^2 = -11; it is nonzero iff some reduced coefficient is nonzero.
z = sp.symbols('z')
dsym = 1 + z; dbsym = 1 + z ** 4                       # zeta5^-1 = zeta5^4
lam = (dbsym ** 2 * m0 - dbsym * m1 + dsym + dsym * dbsym * m2) / (2 - 2 * dsym * dbsym)
num_l, den_l = sp.fraction(sp.together(lam))
expr = sp.expand(num_l ** 2 + dbsym * m0 * den_l ** 2)     # = den^2 (lam^2 + dbar m0)
def reduce_zj(e):
    e = sp.expand(e)
    e = sp.rem(sp.Poly(e, z), sp.Poly(z ** 4 + z ** 3 + z ** 2 + z + 1, z)).as_expr()
    return redj(sp.expand(e))
red_expr = reduce_zj(expr)
den_red = reduce_zj(sp.expand(den_l ** 2))
P = sp.Poly(red_expr, z, J)
print('exact numerator of den^2 (lam^2 + dbar m0) in the basis zeta5^i j^k:', P.as_expr())
print('denominator den^2 reduced:', den_red)
print('EXACT: lam^2 + dbar m0 is nonzero:', P.as_expr() != 0 and den_red != 0)
