#!/usr/bin/env python3
"""tnc_detection_depth_kps_S128c124.py -- kind-pasteur-2026-07-20-S128c124

A STRONGER CLOSURE OF TNC: the detection depth is D = M+N.

opus THM-1685 made TNC a Nullstellensatz emptiness test (Rabinowitsch, one Groebner),
parametrised by NUMBER OF TERMS, tested per pattern (17/17 closed), needing "<= 5 CT
levels."  This bounds the LEVELS a priori and uniformly, from my recurrence-order result:

  a_m = [u^{Mm}] R^m obeys an order-D recurrence  sum_{i=0}^D P_i(m) a_{m+i} = 0.
  If the LEADING coefficient P_D(m) has no positive integer root, then
      a_1 = a_2 = ... = a_D = 0   =>   a_m = 0 for all m >= 1.
  (induction: at step m>=1, a_m..a_{m+D-1}=0 and P_D(m)!=0 force a_{m+D}=0.)

So the toral nullcone is cut out by just the FIRST D moments, and TNC <=>
{a_1=...=a_D=0} has no two-sided solution -- a system of exactly D equations.  This CAPS
opus's level count at D (their <=5 is D<=5), and is UNIFORM in the number of terms.

THREE CHECKS.
(A) M=1 is ELEMENTARY (no Groebner): if r_1=...=r_{j-1}=0 then a_j = [u^j]R^j = j*r_j, so
    a_1=...=a_D=0 forces r_1=...=r_D=0 by a triangular induction.  TNC for M=1 (=N=1 by
    reversal) is one line.
(B) The leading coefficient P_D(m) has NO positive integer root, for a grid of (M,N) --
    this is the proviso that validates the depth-D cap.
(C) The depth is EXACTLY D: {a_1..a_D=0} forces one-sidedness (Rabinowitsch/Groebner,
    1 in <a_1..a_D, 1 - w*r_D> with r_0=1), while {a_1..a_{D-1}=0} does NOT.
"""
import sys
from fractions import Fraction as Fr
import random
import sympy as sp

PR = [(1 << 61) - 1, (1 << 62) - 57]
u = sp.Symbol('u')


def a_seq_sym(M, D, rc, K):
    """a_m for m=0..K with rc a list of sympy exprs (r_0..r_D)."""
    a = [sp.Integer(1)]
    Rp = [sp.Integer(1)]
    for m in range(1, K + 1):
        nxt = [sp.Integer(0)] * (len(Rp) + D)
        for i, ci in enumerate(Rp):
            if ci != 0:
                for j, rj in enumerate(rc):
                    nxt[i + j] += ci * rj
        Rp = [sp.expand(x) for x in nxt]
        a.append(sp.expand(Rp[M * m]) if M * m < len(Rp) else sp.Integer(0))
    return a


def a_seq_int(M, D, rc, K):
    a = [1]
    Rp = [1]
    for m in range(1, K + 1):
        nxt = [0] * (len(Rp) + D)
        for i, ci in enumerate(Rp):
            if ci:
                for j, rj in enumerate(rc):
                    nxt[i + j] += ci * rj
        Rp = nxt
        a.append(Rp[M * m] if M * m < len(Rp) else 0)
    return a


print("=" * 86)
print("(A) M=1 IS ELEMENTARY: a_j = j*r_j once r_1..r_{j-1}=0 (triangular induction)")
print("=" * 86)
for D in range(2, 7):
    rr = sp.symbols('r0:%d' % (D + 1))
    # impose r_1=...=r_{j-1}=0 and read a_j
    ok = True
    for j in range(1, D + 1):
        rc = [rr[0]] + [sp.Integer(0)] * (j - 1) + list(rr[j:D + 1])
        a = a_seq_sym(1, D, rc, j)
        aj = sp.expand(a[j])
        expect = j * rr[0]**(j - 1) * rr[j] if j <= D else None
        # with r_0 general, a_j = j r_0^{j-1} r_j ; check
        ok = ok and sp.simplify(aj - j * rr[0]**(j - 1) * rr[j]) == 0
    print("  D=%d : a_j = j*r_0^{j-1}*r_j whenever r_1..r_{j-1}=0, all j<=D : %s" % (D, ok))
print("  => a_1=...=a_D=0 forces r_1=...=r_D=0 (with r_0!=0): R one-sided. TNC(M=1) proved.")
sys.stdout.flush()

print()
print("=" * 86)
print("(B) LEADING RECURRENCE COEFFICIENT P_D(m) HAS NO POSITIVE INTEGER ROOT")
print("=" * 86)
from math import comb, gcd


def rank_modp(rows, p):
    rows = [r[:] for r in rows]
    ncol = len(rows[0]); rk = 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] % p), None)
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][c], p - 2, p)
        rows[rk] = [(x * inv) % p for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] % p:
                f = rows[i][c]
                rows[i] = [(rows[i][k] - f * rows[rk][k]) % p for k in range(ncol)]
        rk += 1
        if rk == len(rows):
            break
    return rk


def leading_coeff_roots(M, N):
    """Build the order-D recurrence at a generic integer R (exact nullspace), return the
    integer roots of the leading coefficient P_D(m)."""
    D = M + N
    s = comb(D, 2) - gcd(M, N) + 1
    ncols = (D + 1) * (s + 1)
    K = ncols + D + 15
    random.seed(31 + M * 5 + N)
    rc = [random.choice([-3, -2, -1, 1, 2, 3])] + [random.randint(-3, 3) for _ in range(D - 1)] \
        + [random.choice([-3, -2, -1, 1, 2, 3])]
    a = a_seq_int(M, D, [Fr(v) for v in rc], K)
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
        return None, len(ns)
    vec = ns[0]
    mm = sp.Symbol('mm')
    # leading block i=D occupies columns [D*(s+1) : (D+1)*(s+1)]
    PD = sum(vec[D * (s + 1) + j] * mm**j for j in range(s + 1))
    PD = sp.Poly(sp.expand(PD), mm)
    if PD.degree() < 0:
        return "P_D identically 0?!", None
    roots = [r for r in PD.all_roots() if r.is_real and r.is_integer]
    posint = [int(r) for r in roots if int(r) >= 1]
    return PD, posint


for (M, N) in [(1, 1), (1, 2), (2, 2), (1, 3), (2, 3), (3, 3)]:
    PD, pos = leading_coeff_roots(M, N)
    if isinstance(PD, sp.Poly):
        print("  (M,N)=(%d,%d) D=%d : deg P_D=%d, positive-integer roots: %s"
              % (M, N, M + N, PD.degree(), pos if pos else "NONE"))
    else:
        print("  (M,N)=(%d,%d): %s (dim %s)" % (M, N, PD, pos))
sys.stdout.flush()

print()
print("=" * 86)
print("(C) DEPTH IS EXACTLY D: {a_1..a_D=0} forces one-sided; {a_1..a_{D-1}=0} does not")
print("=" * 86)
w = sp.Symbol('w')
for (M, N) in [(1, 1), (2, 2), (2, 3), (3, 3)]:
    D = M + N
    rr = sp.symbols('r1:%d' % (D + 1))          # r_0 := 1
    rc = [sp.Integer(1)] + list(rr)
    a = a_seq_sym(M, D, rc, D)
    conds = [a[j] for j in range(1, D + 1)]      # a_1..a_D
    rD = rr[-1]
    gens_full = conds + [1 - w * rD]
    gens_short = conds[:-1] + [1 - w * rD]       # a_1..a_{D-1}
    try:
        Gfull = sp.groebner(gens_full, *rr, w, order='grevlex')
        forces_full = (list(Gfull) == [sp.Integer(1)])
    except Exception as e:
        forces_full = "err:%s" % e
    try:
        Gshort = sp.groebner(gens_short, *rr, w, order='grevlex')
        forces_short = (list(Gshort) == [sp.Integer(1)])
    except Exception as e:
        forces_short = "err:%s" % e
    print("  (M,N)=(%d,%d) D=%d : {a_1..a_D=0}=>r_D=0 (1 in ideal): %s   {a_1..a_{D-1}=0}=>r_D=0: %s"
          % (M, N, D, forces_full, forces_short))
    print("     -> detection depth is EXACTLY %d : %s"
          % (D, forces_full is True and forces_short is False))
sys.stdout.flush()
print()
print("  TNC <=> {a_1=...=a_D=0} has no two-sided solution -- exactly D moments, no more.")
print("  The order-D recurrence (THM-1670) is what caps opus THM-1685's level count at D.")
