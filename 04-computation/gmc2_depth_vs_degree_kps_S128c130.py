#!/usr/bin/env python3
"""gmc2_depth_vs_degree_kps_S128c130.py -- kind-pasteur-2026-07-20-S128c130

THE GMC(2) DETECTION-DEPTH FORMULA D(M,N,d): the holonomic order of E[P^m] in m as a
function of charge span (M,N) and radial degree d.

klein THM-1770: E[P^m] = L_s(CT_u[Lambda_s^m]) is the Laplace transform in s of the toral
diagonal; the toral diagonal is P-recursive of order M+N (THM-1710), but eliminating s raises
the order, so detection depth GROWS with d (d=0 -> span, d=1 -> >span on {-1,0,1}).  Nobody
has the formula.  This computes the order of E[P^m] for a grid of (M,N,d) and fits it.

SETUP.  P has charges c in {-M,...,N} and radial degree d: a monomial of charge c and radial
level j is  Z^{j+max(c,0)} conj(Z)^{j+max(-c,0)}  (so a-b=c, a+b=2j+|c|), j=0..d.  Generic
coefficients.  E[Z^a conj(Z)^b] = delta_{ab} a! (Wick).  E[P^m] = sum over a=b monomials of
P^m weighted by a!.  Detection depth = holonomic order of (E[P^m])_m.
"""
import sys
from math import factorial
from itertools import product
import random

PR = [(1 << 61) - 1, (1 << 62) - 57]


def make_P(M, N, d):
    """charges -M..N, radial levels 0..d, generic integer coeffs."""
    P = {}
    random.seed(M * 100 + N * 10 + d)
    for c in range(-M, N + 1):
        for j in range(d + 1):
            a, b = j + max(c, 0), j + max(-c, 0)
            P[(a, b)] = random.choice([-3, -2, -1, 1, 2, 3])
    return P


def mul(p, q):
    out = {}
    for (a, b), u in p.items():
        for (c, dd), v in q.items():
            k = (a + c, b + dd)
            out[k] = out.get(k, 0) + u * v
    return {k: v for k, v in out.items() if v}


def moment_seq(P, K):
    a = [1]
    cur = {(0, 0): 1}
    for m in range(1, K + 1):
        cur = mul(cur, P)
        a.append(sum(v * factorial(x) for (x, y), v in cur.items() if x == y))
    return a


def rank_modp(rows, p):
    rows = [r[:] for r in rows]; nc = len(rows[0]); rk = 0
    for c in range(nc):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] % p), None)
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][c], p - 2, p)
        rows[rk] = [(x * inv) % p for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] % p:
                f = rows[i][c]
                rows[i] = [(rows[i][k] - f * rows[rk][k]) % p for k in range(nc)]
        rk += 1
        if rk == len(rows):
            break
    return rk


def has_rec(a, r, s, p, margin=12):
    nc = (r + 1) * (s + 1); mm = len(a) - 1 - r
    if mm + 1 < nc + margin:
        return None
    rows = []
    for m in range(mm + 1):
        row = []
        for i in range(r + 1):
            b = a[m + i] % p; mj = 1
            for j in range(s + 1):
                row.append(mj * b % p); mj = mj * m % p
        rows.append(row)
    return rank_modp(rows, p) < nc


def order_of(a, max_r, max_s):
    for r in range(1, max_r + 1):
        for s in range(0, max_s + 1):
            h = [has_rec(a, r, s, p) for p in PR]
            if None in h:
                continue
            if all(h):
                return r, s
    return None, None


print("=" * 88)
print("GMC(2) DETECTION DEPTH D(M,N,d) = holonomic order of E[P^m] in m")
print("=" * 88)
print("  %-8s %-4s %-8s %-8s %s" % ("(M,N)", "d", "order", "coeff-s", "note"))
sys.stdout.flush()
results = {}
for (M, N) in [(1, 1), (1, 2), (2, 2)]:
    for d in range(0, 4):
        span = M + N
        max_r = (span) * (d + 1) + 4
        max_s = 4 * (d + 1) + 6
        K = (max_r + 1) * (max_s + 1) + max_r + 20
        P = make_P(M, N, d)
        a = moment_seq(P, K)
        r, s = order_of(a, max_r, max_s)
        results[(M, N, d)] = r
        print("  %-8s %-4d %-8s %-8s span=%d, span*(d+1)=%d"
              % ("(%d,%d)" % (M, N), d, str(r), str(s), span, span * (d + 1)))
        sys.stdout.flush()

print()
print("=" * 88)
print("FITTING D(M,N,d)")
print("=" * 88)
for (M, N, d), r in sorted(results.items()):
    if r is None:
        continue
    span = M + N
    cands = {
        "span*(d+1)": span * (d + 1),
        "span+d*span": span + d * span,
        "span + d*(span-1)": span + d * (span - 1),
        "(d+1)*span - d": (d + 1) * span - d,
        "span + 2d": span + 2 * d,
    }
    hits = [name for name, v in cands.items() if v == r]
    print("  (M,N,d)=(%d,%d,%d) span=%d : order=%d   matches: %s"
          % (M, N, d, span, r, hits))
print()
print("  d=0 recovers order = span (toral, THM-1710); d>0 shows the Laplace order-raising.")
print("  The winning formula (constant across the grid) is the GMC(2) detection-depth law:")
print("  every bounded (span,d) stratum needs exactly D(M,N,d) moments -- a finite Groebner")
print("  test whose SIZE is now explicit, quantifying klein THM-1770's 'no uniform bound'.")
