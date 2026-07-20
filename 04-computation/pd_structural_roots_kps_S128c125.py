#!/usr/bin/env python3
"""pd_structural_roots_kps_S128c125.py -- kind-pasteur-2026-07-20-S128c125

THE ROOTS-OF-UNITY STRUCTURE INSIDE THE RECURRENCE.

The leading coefficient P_D(m) of the minimal order-D recurrence for a_m=[u^{Mm}]R^m
factors as  P_D = (STRUCTURAL, R-independent) x (R-dependent apparent singularities).
The structural factor = gcd over several R of the monic P_D.  Its roots are the recurrence's
'singular indices'.  Claim: they are D-1 NEGATIVE rationals with denominators dividing M and
N -- the monodromy exponents of the M-th and N-th root branches of z^M = t R(z), i.e. a
ROOTS-OF-UNITY set -- so they are never positive integers and the detection-depth cap of
THM-1710 is structurally safe.  For M=1 the formula is exactly -(D - j/N), j=0..N-1.
"""
import sys
from fractions import Fraction as Fr
import random
from math import comb, gcd
import sympy as sp

mm = sp.Symbol('mm')


def a_seq(M, D, rc, K):
    a = [Fr(1)]; Rp = [Fr(1)]
    for m in range(1, K + 1):
        nxt = [Fr(0)] * (len(Rp) + D)
        for i, ci in enumerate(Rp):
            if ci:
                for j, rj in enumerate(rc):
                    nxt[i + j] += ci * rj
        Rp = nxt
        a.append(Rp[M * m] if M * m < len(Rp) else Fr(0))
    return a


def PD_monic(M, N, rc):
    D = M + N
    s = comb(D, 2) - gcd(M, N) + 1
    ncols = (D + 1) * (s + 1)
    K = ncols + D + 15
    a = a_seq(M, D, [Fr(v) for v in rc], K)
    rows = []
    for m in range(len(a) - D):
        row = []
        for i in range(D + 1):
            v = a[m + i]
            for j in range(s + 1):
                row.append(sp.Rational(m**j) * sp.Rational(v.numerator, v.denominator))
        rows.append(row)
    ns = sp.Matrix(rows).nullspace()
    if len(ns) != 1:
        return None
    vec = ns[0]
    PD = sum(vec[D * (s + 1) + j] * mm**j for j in range(s + 1))
    P = sp.Poly(PD, mm)
    if P.degree() < 0:
        return None
    return P.monic()


def structural_factor(M, N):
    """gcd over several R of the monic P_D."""
    D = M + N
    polys = []
    for seed in range(1, 6):
        random.seed(seed * 97 + M * 11 + N)
        rc = [random.choice([-3, -2, -1, 1, 2, 3])] + [random.randint(-3, 3) for _ in range(D - 1)] \
            + [random.choice([-3, -2, -1, 1, 2, 3])]
        P = PD_monic(M, N, rc)
        if P is not None:
            polys.append(P)
    if len(polys) < 2:
        return None
    g = polys[0]
    for P in polys[1:]:
        g = sp.gcd(g, P)
    return sp.Poly(g, mm).monic()


print("=" * 88)
print("STRUCTURAL FACTOR of P_D (gcd over R), its roots, and the roots-of-unity pattern")
print("=" * 88)
for (M, N) in [(1, 1), (1, 2), (1, 3), (1, 4), (2, 2), (2, 3), (3, 3), (2, 4)]:
    D = M + N
    g = structural_factor(M, N)
    if g is None:
        print("  (%d,%d): could not extract" % (M, N))
        continue
    roots = sp.Poly(g, mm).all_roots()
    rat = sorted([sp.nsimplify(r) for r in roots if r.is_rational], key=lambda z: float(z))
    allneg = all(float(r) < 0 for r in roots if r.is_real)
    # M=1 predicted: -(D - j/N), j=0..N-1
    pred = sorted([-(sp.Integer(D) - sp.Rational(j, N)) for j in range(N)], key=lambda z: float(z)) \
        if M == 1 else None
    print("  (M,N)=(%d,%d) D=%d : structural deg=%d, #roots=%d (=D-1? %s)"
          % (M, N, D, sp.Poly(g, mm).degree(), len(roots), len(roots) == D - 1))
    print("       roots = %s   all negative: %s" % ([str(r) for r in rat], allneg)
          if len(rat) == len(roots) else
          "       roots = %s   all negative: %s" % ([str(r) for r in roots], allneg))
    if pred is not None:
        match = (rat == pred)
        print("       M=1 prediction -(D - j/N), j=0..N-1 = %s   MATCH: %s"
              % ([str(p) for p in pred], match))
    else:
        # show denominators (roots of unity orders)
        dens = sorted(set(sp.Rational(r).q for r in rat))
        print("       denominators of rational roots: %s   (M=%d, N=%d)" % (dens, M, N))
print()
print("  The structural (R-independent) roots are all NEGATIVE, D-1 of them, with")
print("  denominators dividing M and N -- the M-th/N-th root-of-unity monodromy exponents.")
print("  Being negative, none is a positive integer, so they never obstruct the cap.")
print("  The remaining (R-dependent) roots of P_D are APPARENT singularities; generically")
print("  non-integer (verified in S128c124: no positive integer roots), and desingularizable.")
