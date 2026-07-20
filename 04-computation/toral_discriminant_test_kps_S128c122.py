#!/usr/bin/env python3
"""toral_discriminant_test_kps_S128c122.py -- kind-pasteur-2026-07-20-S128c122  (part B, fixed)

Two fixes over the previous part B:
  * confirm r(1,2) = 3 over several GENERIC instances (the earlier [2,3] disagreement was a
    degenerate instance giving a spurious 2);
  * the coeff-degree s and all nullspaces use GENERIC random two-sided R, never the
    palindromic all-ones probe that made the nullspace dimension jump.

QUESTION (b), stated as a normalisation-independent ZERO-TEST (the nullspace is
1-dimensional, so "is the (i=0,j=0) entry zero" is well defined up to the free scalar):
  * generic R (disc != 0): is P_0(0) != 0?         [off-locus: the descent-at-0 works]
  * double-root R (disc = 0): is P_0(0) = 0?        [is {P_0(0)=0} = {disc=0}?]
At (1,1), P_0(0) = disc(R) by hand, so both hold.  Test (1,2) and (2,2) honestly.
"""
import sys
from fractions import Fraction as Fr
import random
import sympy as sp

PR = [(1 << 61) - 1, (1 << 62) - 57]
u = sp.Symbol('u')


def seq_toral(M, D, rc, K):
    a = [rc[0] * 0 + 1]
    Rpow = [rc[0] * 0 + 1]
    for m in range(1, K + 1):
        nxt = [rc[0] * 0] * (len(Rpow) + D)
        for i, ci in enumerate(Rpow):
            if ci:
                for j, rj in enumerate(rc):
                    nxt[i + j] = nxt[i + j] + ci * rj
        Rpow = nxt
        a.append(Rpow[M * m] if M * m < len(Rpow) else rc[0] * 0)
    return a


def rank_modp(rows, p):
    rows = [r[:] for r in rows]
    ncol = len(rows[0])
    rk = 0
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


def has_rec_p(a, r, s, p, margin=20):
    ncols = (r + 1) * (s + 1)
    maxm = len(a) - 1 - r
    if maxm + 1 < ncols + margin:
        return None
    rows = []
    for m in range(maxm + 1):
        row = []
        for i in range(r + 1):
            base = a[m + i] % p
            mj = 1
            for j in range(s + 1):
                row.append(mj * base % p)
                mj = mj * m % p
        rows.append(row)
    return rank_modp(rows, p) < ncols


def order_of(a, S_BIG, max_r):
    for r in range(1, max_r + 1):
        hits = [has_rec_p(a, r, S_BIG, p) for p in PR]
        if None in hits:
            return None
        if all(hits):
            return r
    return None


def min_s(a, r, S_BIG):
    for s in range(S_BIG + 1):
        if has_rec_p(a, r, s, PR[0]) and has_rec_p(a, r, s, PR[1]):
            return s
    return None


def P0_is_zero(M, D, rc, r, s):
    ncols = (r + 1) * (s + 1)
    K = ncols + r + 15
    a = seq_toral(M, D, [Fr(v) for v in rc], K)
    rows = []
    for m in range(len(a) - r):
        row = []
        for i in range(r + 1):
            v = a[m + i]
            for j in range(s + 1):
                row.append(sp.Rational(m**j) * sp.Rational(v.numerator, v.denominator))
        rows.append(row)
    ns = sp.Matrix(rows).nullspace()
    if len(ns) != 1:
        return None, len(ns)
    return (ns[0][0] == 0), 1


print("=" * 84)
print("CONFIRM r(1,2) over several generic instances")
print("=" * 84)
random.seed(3)
orders = []
for _ in range(6):
    rc = [random.choice([-3, -2, -1, 1, 2, 3]), random.randint(-3, 3),
          random.randint(-3, 3), random.choice([-3, -2, -1, 1, 2, 3])]
    a = seq_toral(1, 3, rc, 200)
    orders.append(order_of(a, 9, 5))
print("  orders over 6 generic R: %s   -> generic r(1,2) = %d (mode)"
      % (orders, max(set(orders), key=orders.count)))
sys.stdout.flush()

print()
print("=" * 84)
print("(b) DISCRIMINANT ZERO-TEST for (1,1), (1,2), (2,2)")
print("=" * 84)
for (M, N, r) in [(1, 1, 2), (1, 2, 3), (2, 2, 4)]:
    D = M + N
    S_BIG = D * (D - 1) // 2 + 6
    # find s on a GENERIC R
    random.seed(500 + M * 7 + N)
    rcg = [random.choice([-3, -2, -1, 1, 2, 3])] + [random.randint(-3, 3) for _ in range(D - 1)] \
        + [random.choice([-3, -2, -1, 1, 2, 3])]
    Kg = (D + 2) * (S_BIG + 1) + 30
    s = min_s(seq_toral(M, D, rcg, Kg), r, S_BIG)
    print()
    print("  (M,N)=(%d,%d) D=%d order r=%d coeff-deg s=%s" % (M, N, D, r, s))
    if s is None:
        print("     could not pin s; skip")
        continue
    # generic disc != 0
    gz = []
    random.seed(600 + M * 7 + N)
    while len(gz) < 4:
        rc = [random.randint(-4, 4) for _ in range(D + 1)]
        if rc[0] == 0 or rc[D] == 0:
            continue
        Rp = sum(sp.Integer(rc[j]) * u**j for j in range(D + 1))
        if sp.discriminant(Rp, u) == 0:
            continue
        z, nd = P0_is_zero(M, D, rc, r, s)
        if z is None:
            print("     [generic] nullspace dim %d at rc=%s (s off?)" % (nd, rc))
            continue
        gz.append(z)
    print("     generic (disc!=0): P_0(0)==0? %s  => P_0(0)!=0 off disc locus: %s"
          % (gz, all(x is False for x in gz)))
    # double-root disc = 0
    dz = []
    random.seed(700 + M * 7 + N)
    tries = 0
    while len(dz) < 4 and tries < 200:
        tries += 1
        rho = random.choice([-2, -1, 1, 2])
        deg_S = D - 2
        Sc = [random.randint(-3, 3) for _ in range(deg_S + 1)]
        if Sc[0] == 0:
            Sc[0] = 1
        if Sc[-1] == 0:
            Sc[-1] = 1
        Rp = sp.expand((u - rho)**2 * sum(sp.Integer(Sc[j]) * u**j for j in range(deg_S + 1)))
        pc = sp.Poly(Rp, u).all_coeffs()[::-1]
        rc = [int(x) for x in pc]
        if len(rc) != D + 1 or rc[0] == 0 or rc[D] == 0:
            continue
        z, nd = P0_is_zero(M, D, rc, r, s)
        if z is None:
            continue
        dz.append((z, rc, sp.discriminant(Rp, u)))
    print("     double-root (disc=0): P_0(0)==0? %s" % [z for z, _, _ in dz])
    allz = bool(dz) and all(z for z, _, _ in dz)
    print("     -> {disc=0} => {P_0(0)=0} on the sample: %s" % allz)
    if dz and not allz:
        bad = next(rc for z, rc, _ in dz if not z)
        print("        HONEST NEGATIVE beyond (1,1): disc=0 but P_0(0)!=0, e.g. R=%s" % bad)
sys.stdout.flush()
