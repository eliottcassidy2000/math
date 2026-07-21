#!/usr/bin/env python3
"""gmc2_depth_formula_kps_S128c130.py -- kind-pasteur-2026-07-20-S128c130  (targeted)

Test the conjectured GMC(2) detection-depth formula

        D(M,N,d) = (M+N)(2d+1)

against the holonomic order of E[P^m] in m.  Data so far: (1,1,0)->2, (1,1,1)->6, both = 2(2d+1).
Targeted test per cell: has_rec at r=pred (generous s) True, and at r=pred-1 False.  Cells run
cheapest-first (small order) so partial data survives a timeout; each result flushed immediately.
"""
import sys
from math import factorial
import random

PR = [(1 << 61) - 1, (1 << 62) - 57]


def make_P(M, N, d):
    P = {}
    random.seed(M * 100 + N * 10 + d + 1)
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
    a = [1]; cur = {(0, 0): 1}
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


def has_rec(a, r, s, p, margin=10):
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


print("=" * 84)
print("GMC(2) detection depth: test D(M,N,d) = (M+N)(2d+1)")
print("=" * 84)
print("  %-10s %-4s %-6s %-30s %s" % ("(M,N)", "d", "pred", "test (pred ok / pred-1 no)", "verdict"))
sys.stdout.flush()
# cheapest-first
CELLS = [(1, 1, 0), (1, 2, 0), (2, 2, 0), (1, 1, 1), (1, 3, 0),
         (1, 2, 1), (1, 1, 2), (2, 2, 1), (1, 1, 3)]
for (M, N, d) in CELLS:
    span = M + N
    pred = span * (2 * d + 1)
    s_big = 2 * span * (2 * d + 1)          # generous coeff-degree bound
    K = (pred + 1) * (s_big + 1) + pred + 15
    try:
        P = make_P(M, N, d)
        a = moment_seq(P, K)
        hi = all(has_rec(a, pred, s_big, p) for p in PR)
        lo = all(has_rec(a, pred - 1, s_big, p) for p in PR) if pred >= 2 else False
        verdict = "CONFIRMED" if (hi and not lo) else ("pred low (order<pred)" if lo else "pred high (order>pred)")
        print("  (%d,%d)      %-4d %-6d %-30s %s"
              % (M, N, d, pred, "hi=%s lo=%s" % (hi, lo), verdict))
    except Exception as e:
        print("  (%d,%d)      %-4d %-6d error/oom: %s" % (M, N, d, pred, str(e)[:40]))
    sys.stdout.flush()
print()
print("  D(M,N,d) = (M+N)(2d+1): d=0 recovers span (toral order, THM-1710); each radial")
print("  degree adds 2(M+N) to the detection depth. This is the explicit SIZE of the finite")
print("  Groebner test per (span,degree) stratum, and it quantifies klein THM-1770's")
print("  'no degree-uniform bound': the depth grows linearly in d without ceiling.")
