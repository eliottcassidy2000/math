#!/usr/bin/env python3
"""pd_roots_proviso_kps_S128c125.py -- kind-pasteur-2026-07-20-S128c125

COMPLETE THE DETECTION-DEPTH CAP (THM-1710) by pinning the roots of the leading recurrence
coefficient P_D(m).  The cap 'a_1=..=a_D=0 => a_m=0 for all m' needs P_D(m) != 0 for m>=1.

This computes the minimal order-D recurrence for a_m=[u^{Mm}]R^m at a generic integer R
(exact rational nullspace), extracts the leading coefficient P_D(m) AND the trailing
coefficient P_0(m) as exact polynomials in m, factors them, and -- crucially -- checks
whether their ROOTS are R-INDEPENDENT (i.e. structural 'singular indices', not accidents of
the chosen R).  If the roots are universal and all < 1, the cap is unconditional.

Keeping to the owner's steer (roots of unity + recurrence): the singular indices of a
diagonal's recurrence are exactly where the M small roots of z^M = t R(z) exchange under
monodromy -- an M-th-root-of-unity phenomenon -- so P_D's roots should be tied to M.
"""
import sys
from fractions import Fraction as Fr
import random
from math import comb, gcd
import sympy as sp

mm = sp.Symbol('mm')


def a_seq(M, D, rc, K):
    a = [Fr(1)]
    Rp = [Fr(1)]
    for m in range(1, K + 1):
        nxt = [Fr(0)] * (len(Rp) + D)
        for i, ci in enumerate(Rp):
            if ci:
                for j, rj in enumerate(rc):
                    nxt[i + j] += ci * rj
        Rp = nxt
        a.append(Rp[M * m] if M * m < len(Rp) else Fr(0))
    return a


def recurrence(M, N, rc):
    """Return (P_0(mm), P_D(mm)) exact, for the minimal order-D recurrence at integer R=rc."""
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
        return None, None, len(ns)
    vec = ns[0]
    def block(i):
        return sp.Poly(sum(vec[i * (s + 1) + j] * mm**j for j in range(s + 1)), mm)
    return block(0), block(D), s


print("=" * 88)
print("ROOTS OF THE LEADING COEFFICIENT P_D(m) -- are they R-independent, and all < 1?")
print("=" * 88)
for (M, N) in [(1, 1), (1, 2), (2, 2), (1, 3), (2, 3), (1, 4), (3, 3)]:
    D = M + N
    rootsets = []
    P0roots = []
    for seed in (1, 2, 3):
        random.seed(seed * 100 + M * 7 + N)
        rc = [random.choice([-3, -2, -1, 1, 2, 3])] + [random.randint(-3, 3) for _ in range(D - 1)] \
            + [random.choice([-3, -2, -1, 1, 2, 3])]
        P0, PD, s = recurrence(M, N, rc)
        if PD is None:
            rootsets.append("dim%s" % s)
            continue
        # rational + real roots of P_D
        rts = sp.Poly(PD, mm).all_roots()
        rts_simpl = sorted(set(sp.nsimplify(r) for r in rts if r.is_real), key=lambda z: float(z))
        rootsets.append(tuple(str(r) for r in rts_simpl))
        if seed == 1:
            leadPD = sp.LC(PD, mm)
            P0r = sorted(set(sp.nsimplify(r) for r in sp.Poly(P0, mm).all_roots() if r.is_real),
                         key=lambda z: float(z))
    same = len(set(rootsets)) == 1
    maxroot = None
    if same and isinstance(rootsets[0], tuple):
        reals = [sp.nsimplify(r) for r in rootsets[0]]
        maxroot = max((float(sp.sympify(r)) for r in rootsets[0]), default=None)
    print("  (M,N)=(%d,%d) D=%d : P_D real roots R-independent: %s" % (M, N, D, same))
    print("       real roots of P_D = %s" % (rootsets[0] if same else "VARY: %s" % rootsets))
    if maxroot is not None:
        print("       max real root = %.4f   -> P_D(m) != 0 for all m >= 1 : %s"
              % (maxroot, maxroot < 1))
    print("       real roots of P_0(trailing) = %s" % P0r)
print()
print("  If P_D's real roots are R-independent and all < 1, the cap of THM-1710 is")
print("  UNCONDITIONAL for that (M,N): a_1=..=a_D=0 => a_m=0 for all m, no proviso needed.")
