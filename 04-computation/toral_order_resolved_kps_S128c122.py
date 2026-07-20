#!/usr/bin/env python3
"""toral_order_resolved_kps_S128c122.py -- kind-pasteur-2026-07-20-S128c122  (v2, fast+flush)

DECISIVE order r(M,N) of a_m = [u^{Mm}] R^m, resolving the (2,2) conflict, plus the honest
test of "is the descent's trailing-coefficient locus {P_0(0)=0} equal to {disc(R)=0}?"

(A) ORDER: test has_rec(r, S_BIG) for r = 1,2,... ASCENDING; the first True is the minimal
    order.  Justification: a genuinely order-r0 sequence admits NO recurrence of order < r0
    at ANY coefficient degree, so a single generous S_BIG suffices for every r.  Two primes
    must agree; >= 20 holdout rows.  Fast (r0 rank calls, not a full (r,s) scan).

(B) DISCRIMINANT, done as a NORMALISATION-INDEPENDENT ZERO-TEST.  The first pass compared
    primitive-normalised P_0(0) to disc(R) and got varying ratios -- but primitive
    normalisation divides by a per-R integer content, so that ratio is meaningless.  What IS
    well defined (the nullspace is 1-dimensional, so any entry is zero-or-not independent of
    scaling) is WHETHER the trailing coefficient at m=0 vanishes.  So:
      * generic R (disc != 0): is P_0(0) != 0?   [if yes, off-locus descent works]
      * double-root R (disc = 0): is P_0(0) = 0?  [if yes, {P_0(0)=0} contains {disc=0}]
    At (1,1) I proved by hand P_0(0) = disc, so both hold there.  This tests the rest
    HONESTLY -- if disc=0 does NOT force P_0(0)=0, the discriminant story is false beyond
    (1,1) and that is reported as a negative.
"""
import sys
from fractions import Fraction as Fr
import random
from math import comb
import sympy as sp

PR = [(1 << 61) - 1, (1 << 62) - 57]


def flush():
    sys.stdout.flush()


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


print("=" * 88)
print("(A) r(M,N), DECISIVE (ascending-r, one big S, two primes, >=20 holdout)")
print("=" * 88)
flush()
random.seed(11)
table = {}
for D in range(2, 7):
    for M in range(1, D // 2 + 1):
        N = D - M
        S_BIG = D * (D - 1) // 2 + 6
        max_r = D + 1
        K = (max_r + 1) * (S_BIG + 1) + max_r + 30
        got = []
        for _ in range(2):
            rc = [random.choice([-3, -2, -1, 1, 2, 3])] + \
                 [random.randint(-3, 3) for _ in range(D - 1)] + \
                 [random.choice([-3, -2, -1, 1, 2, 3])]
            got.append(order_of(seq_toral(M, D, [int(x) for x in rc], K), S_BIG, max_r))
        r = got[0] if got[0] == got[1] else "DISAGREE %s" % got
        table[(M, N)] = r
        print("  (M,N)=(%d,%d)  D=%d  ->  r = %-12s  [C(D,M)=%d, D=%d, D-[M==N]=%d]"
              % (M, N, D, str(r), comb(D, M), D, D - (1 if M == N else 0)))
        flush()

print()
print("  r(M,N) summary:  (order 2 <=> orthogonal family)")
for (M, N), r in sorted(table.items()):
    print("     (%d,%d) D=%d : r=%s" % (M, N, M + N, r))
flush()


# ---- exact rational nullspace, minimal, returns the (i=0,j=0) entry zero-test ----
def trailing0_zero(M, D, rc, r, s):
    """Return (is_P0_0_zero, nullspace_dim) using EXACT rational nullspace at integer R."""
    ncols = (r + 1) * (s + 1)
    K = ncols + r + 12
    a = seq_toral(M, D, [Fr(v) for v in rc], K)
    rows = []
    for m in range(len(a) - r):
        row = []
        for i in range(r + 1):
            val = a[m + i]
            for j in range(s + 1):
                row.append(sp.Rational(m ** j) * sp.Rational(val.numerator, val.denominator))
        rows.append(row)
    ns = sp.Matrix(rows).nullspace()
    if len(ns) != 1:
        return None, len(ns)
    return (ns[0][0] == 0), 1


def s_at_order(a, r, S_BIG):
    for s in range(S_BIG + 1):
        if has_rec_p(a, r, s, PR[0]) and has_rec_p(a, r, s, PR[1]):
            return s
    return None


print()
print("=" * 88)
print("(B) DOES {disc(R)=0} FORCE P_0(0)=0?  normalisation-independent zero-test")
print("=" * 88)
flush()
u = sp.Symbol('u')
for (M, N) in [(1, 1), (1, 2), (2, 2), (1, 3)]:
    D = M + N
    r = table[(M, N)]
    if not isinstance(r, int):
        print("  (%d,%d): order unresolved, skipping" % (M, N))
        continue
    S_BIG = D * (D - 1) // 2 + 6
    K = (D + 2) * (S_BIG + 1) + 30
    s = s_at_order(seq_toral(M, D, [1] + [1] * (D - 1) + [1], K), r, S_BIG)
    print()
    print("  (M,N)=(%d,%d) D=%d order r=%d coeff-deg s=%d" % (M, N, D, r, s))
    # generic R (disc != 0): expect P_0(0) != 0
    random.seed(50 + M * 7 + N)
    gz = []
    for _ in range(4):
        rc = [random.randint(-4, 4) for _ in range(D + 1)]
        if rc[0] == 0 or rc[D] == 0:
            continue
        Rp = sum(sp.Integer(rc[j]) * u**j for j in range(D + 1))
        if sp.discriminant(Rp, u) == 0:
            continue
        z, nd = trailing0_zero(M, D, rc, r, s)
        gz.append(z)
    print("     generic R (disc!=0): P_0(0)==0 ?  %s   -> so P_0(0)!=0 off the disc locus: %s"
          % (gz, all(zz is False for zz in gz)))
    # double-root R (disc = 0): R = (u - rho)^2 * S(u), integer rho, S two-sided
    dz = []
    random.seed(200 + M * 7 + N)
    tries = 0
    while len(dz) < 4 and tries < 60:
        tries += 1
        rho = random.choice([-2, -1, 1, 2])
        Scoef = [random.randint(-3, 3) for _ in range(D - 2 + 1)]
        if D - 2 >= 0:
            if Scoef[0] == 0:
                Scoef[0] = 1
            if Scoef[-1] == 0:
                Scoef[-1] = 1
        Spoly = sum(sp.Integer(Scoef[j]) * u**j for j in range(len(Scoef)))
        Rp = sp.expand((u - rho)**2 * Spoly)
        pc = sp.Poly(Rp, u).all_coeffs()[::-1]        # low -> high
        rc = [int(x) for x in pc]
        if len(rc) != D + 1 or rc[0] == 0 or rc[D] == 0:
            continue
        assert sp.discriminant(Rp, u) == 0
        z, nd = trailing0_zero(M, D, rc, r, s)
        if z is None:
            continue
        dz.append((z, rc))
    allzero = all(zz for zz, _ in dz)
    print("     double-root R (disc=0): P_0(0)==0 ?  %s" % [zz for zz, _ in dz])
    print("     -> {disc=0} forces P_0(0)=0 (discriminant story holds here): %s" % allzero)
    if not allzero:
        print("        HONEST NEGATIVE: the trailing-coefficient locus is NOT the discriminant")
        print("        locus beyond (1,1); example R with disc=0 but P_0(0)!=0: %s"
              % next(rc for zz, rc in dz if not zz))
flush()
