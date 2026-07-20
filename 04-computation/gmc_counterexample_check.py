#!/usr/bin/env python3
"""
Exact verification of a claimed counterexample to the Gaussian Moments
Conjecture GMC(n)  [Derksen - van den Essen - Zhao, arXiv:1506.05192,
Israel J. Math. 219 (2017) 917-928, Conjecture 1.1].

GMC(n):  P in C[x_1..x_n].  If  E[P(X)^m] = 0  for all m >= 1,
         then for every Q in C[x_1..x_n],  E[Q(X) P(X)^m] = 0  for m >> 0.
         (X = standard Gaussian vector in R^n.)

CLAIMED COUNTEREXAMPLE (n = 4 real Gaussians):
    P = (1 + Z2) * ( W1*(1 - Z1) + W2 )
    Q = Z2
with  E[P^m] = 0 for all m >= 1  but  E[Q P^m] = m!  for all m >= 1.

PACKAGING.  x1,x2,x3,x4 iid N(0,1).  Set
    Z1 = (x1 + i x2)/sqrt(2),  W1 = (x1 - i x2)/sqrt(2)
    Z2 = (x3 + i x4)/sqrt(2),  W2 = (x3 - i x4)/sqrt(2)
so (Z1,W1) and (Z2,W2) are two independent conjugate complex-Gaussian pairs.
Wick contractions:  E[Zi Wj] = delta_ij,  E[Zi Zj] = E[Wi Wj] = 0.
Hence the exact monomial moment rule

    E[ Z1^a1 Z2^a2 W1^b1 W2^b2 ] = a1! * a2!   if a1 == b1 and a2 == b2
                                 = 0            otherwise.

Everything below is exact integer arithmetic -- no floating point, no sampling.
"""

from math import factorial

# ---------------------------------------------------------------- polynomials
# monomial = (a1, a2, b1, b2) = exponents of (Z1, Z2, W1, W2)
# polynomial = dict monomial -> integer coefficient


def padd(p, q):
    r = dict(p)
    for m, c in q.items():
        r[m] = r.get(m, 0) + c
    return {m: c for m, c in r.items() if c != 0}


def pmul(p, q):
    r = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            m = (m1[0] + m2[0], m1[1] + m2[1], m1[2] + m2[2], m1[3] + m2[3])
            r[m] = r.get(m, 0) + c1 * c2
    return {m: c for m, c in r.items() if c != 0}


def ppow(p, k):
    r = {(0, 0, 0, 0): 1}
    for _ in range(k):
        r = pmul(r, p)
    return r


def expect(p):
    """Exact Gaussian expectation via the complex Wick monomial rule."""
    tot = 0
    for (a1, a2, b1, b2), c in p.items():
        if a1 == b1 and a2 == b2:
            tot += c * factorial(a1) * factorial(a2)
    return tot


# ---------------------------------------------------------------- the example
ONE = {(0, 0, 0, 0): 1}
Z1 = {(1, 0, 0, 0): 1}
Z2 = {(0, 1, 0, 0): 1}
W1 = {(0, 0, 1, 0): 1}
W2 = {(0, 0, 0, 1): 1}

# P = (1 + Z2) * ( W1*(1 - Z1) + W2 )
one_minus_Z1 = padd(ONE, {(1, 0, 0, 0): -1})
P = pmul(padd(ONE, Z2), padd(pmul(W1, one_minus_Z1), W2))
Q = Z2

print("P expanded (exponents are (Z1,Z2,W1,W2)):")
for m in sorted(P):
    print("   %+d * Z1^%d Z2^%d W1^%d W2^%d" % (P[m], m[0], m[1], m[2], m[3]))
print()

MMAX = 9
print(" m |        E[P^m] |      E[Q P^m] |          m! | match")
print("---+---------------+---------------+-------------+------")
ok_p = ok_q = True
for m in range(1, MMAX + 1):
    Pm = ppow(P, m)
    ep = expect(Pm)
    eq = expect(pmul(Q, Pm))
    fm = factorial(m)
    ok_p &= (ep == 0)
    ok_q &= (eq == fm)
    print(" %d | %13d | %13d | %11d | %s"
          % (m, ep, eq, fm, "yes" if eq == fm else "NO"))

print()
print("E[P^m] == 0        for 1 <= m <= %d : %s" % (MMAX, ok_p))
print("E[Q P^m] == m!     for 1 <= m <= %d : %s" % (MMAX, ok_q))

# ------------------------------------------------- sanity check of the machinery
# Known: E[Z^a W^a] = a!, E[Z^a W^b] = 0 for a != b, E[|Z|^2] = 1.
assert expect(pmul(ppow(Z1, 3), ppow(W1, 3))) == 6
assert expect(pmul(ppow(Z1, 3), ppow(W1, 2))) == 0
assert expect(pmul(Z1, W1)) == 1
assert expect(pmul(Z1, W2)) == 0
assert expect(pmul(Z1, Z1)) == 0
print("Wick machinery sanity checks: passed")
print()

# ---------------------------------------------------------------- CLOSED FORM
# PROOF (valid for every m >= 1).
#
# Write P^m = (1+Z2)^m * (W1(1-Z1) + W2)^m and binomially expand the second factor:
#     P^m = sum_k C(m,k) * [(1-Z1)^k W1^k] * [(1+Z2)^m W2^(m-k)].
# (Z1,W1) is independent of (Z2,W2), so the expectation factorises.
#
#   E[(1-Z1)^k W1^k]      = (-1)^k k!            (only the Z1^k term survives)
#   E[(1+Z2)^m W2^(m-k)]  = C(m,m-k)(m-k)! = m!/k!
# hence
#   E[P^m] = sum_k C(m,k) (-1)^k k! * m!/k! = m! * sum_k C(m,k)(-1)^k
#          = m! * (1-1)^m = 0     for every m >= 1.                        [A]
#
#   E[Z2 (1+Z2)^m W2^(m-k)] = C(m, m-k-1)(m-k)!   (needs k <= m-1)
# hence
#   E[Q P^m] = sum_{k=0}^{m-1} C(m,k)(-1)^k k! * C(m,m-k-1)(m-k)!
# Using C(m,k)k! = m!/(m-k)! and C(m,m-k-1)(m-k)! = m!(m-k)/(k+1)!, the k-th
# term is (-1)^k (m!)^2 / [(m-k-1)!(k+1)!].  Substituting j = k+1:
#   E[Q P^m] = -(m!)^2 sum_{j=1}^{m} (-1)^j /[(m-j)! j!]
#            = -m! * [ (1-1)^m - 1 ] = m!    for every m >= 1.             [B]


def closed_form_EPm(m):
    return sum(binom(m, k) * (-1) ** k * factorial(k) * (factorial(m) // factorial(k))
               for k in range(m + 1))


def closed_form_EQPm(m):
    return sum(binom(m, k) * (-1) ** k * factorial(k)
               * binom(m, m - k - 1) * factorial(m - k)
               for k in range(m))


def binom(a, b):
    if b < 0 or b > a:
        return 0
    return factorial(a) // (factorial(b) * factorial(a - b))


print("Closed-form identities [A] and [B], checked to m = 60:")
bad = []
for m in range(1, 61):
    if closed_form_EPm(m) != 0 or closed_form_EQPm(m) != factorial(m):
        bad.append(m)
print("   sum_k C(m,k)(-1)^k m!      == 0   for all 1<=m<=60 : %s" % (not bad))
print("   sum_k ... (the E[QP^m] sum) == m!  for all 1<=m<=60 : %s" % (not bad))
if bad:
    print("   FAILURES at m =", bad)

# cross-check closed form against the brute-force Wick expansion
for m in range(1, MMAX + 1):
    Pm = ppow(P, m)
    assert expect(Pm) == closed_form_EPm(m) == 0
    assert expect(pmul(Q, Pm)) == closed_form_EQPm(m) == factorial(m)
print("   closed form agrees with brute-force expansion for m <= %d" % MMAX)

# -------------------------------------------- independent Monte-Carlo cross-check
# This validates the PACKAGING (that Z1,Z2,W1,W2 really are built from four
# iid real N(0,1) variables) without using the Wick rule at all.
try:
    import numpy as np

    rng = np.random.default_rng(20260720)
    N = 4_000_000
    x1, x2, x3, x4 = (rng.standard_normal(N) for _ in range(4))
    z1 = (x1 + 1j * x2) / np.sqrt(2.0)
    w1 = (x1 - 1j * x2) / np.sqrt(2.0)
    z2 = (x3 + 1j * x4) / np.sqrt(2.0)
    w2 = (x3 - 1j * x4) / np.sqrt(2.0)
    Pv = (1 + z2) * (w1 * (1 - z1) + w2)
    Qv = z2
    print()
    print("Monte-Carlo cross-check with %d samples of FOUR REAL N(0,1):" % N)
    print("  m |     MC E[P^m] (exact 0) |   MC E[Q P^m] |    m!")
    for m in (1, 2, 3):
        Pm = Pv ** m
        a = Pm.mean()
        b = (Qv * Pm).mean()
        print("  %d | %10.4f%+10.4fi | %10.4f%+10.4fi | %5d"
              % (m, a.real, a.imag, b.real, b.imag, factorial(m)))
except ImportError:
    print("(numpy unavailable - skipping Monte-Carlo cross-check)")


# ============================================================================
# APPENDIX: independent check of the Alpoge (19 July 2026) counterexample to
# the Jacobian Conjecture in dimension 3, as published at
# https://jacobianfun.org/jacobian-explained  (no arXiv preprint as of
# 2026-07-20).  Relevant here because DEZ Theorem 1.6 says
#     GMC(n) for all n  ==>  JC(n) for all n,
# so JC(3) false already forces GMC(N) false for SOME finite N -- though the
# published chain (Prop 3.2: GMC(2n) => Ker(F_n) is MZ) does NOT localise
# that N to 4.  The P,Q above do localise it to N = 4.
# ============================================================================
def check_jc_counterexample():
    try:
        from sympy import symbols, Matrix, simplify, Rational, expand
    except ImportError:
        print("(sympy unavailable - skipping JC check)")
        return
    x, y, z = symbols('x y z')
    u = 1 + x * y
    F = Matrix([u**3 * z + y**2 * u * (4 + 3 * x * y),
                y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
                2 * x - 3 * x**2 * y - x**3 * z])
    print()
    print("Alpoge JC(3) counterexample:")
    print("   det J =", simplify(expand(F.jacobian([x, y, z]).det())),
          "  component degrees =",
          [expand(f).as_poly(x, y, z).total_degree() for f in F])
    for p in [(0, 0, Rational(-1, 4)),
              (1, Rational(-3, 2), Rational(13, 2)),
              (-1, Rational(3, 2), Rational(13, 2))]:
        val = [simplify(f.subs({x: p[0], y: p[1], z: p[2]})) for f in F]
        print("   F%-22s = %s" % (str(p), val))
    print("   => constant nonzero Jacobian determinant but not injective.")


check_jc_counterexample()
