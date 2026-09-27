#!/usr/bin/env python3
"""Orchestrator audit of lane `robin2` (Conjecture R with a constant). The exact counters N_m(L) (Robin
barrier) and A_M(L) (hard wall) are the orchestrator's own, written for the robin audit
(procgen_robin_20260926_orchestrator_check.py) from the robin note's definitions; the robin2 lane's code
was not read.

  1. N_1(L) = 0 for L >= 1.
  2. Theorem R_15 on a finite range: 10000 N_m(L) <= 10039 A_(m+15)(L) for 2 <= m <= 14, L <= 400.
  3. The ratio profile: max N_m/A_(m+c0) over 2 <= m <= 14, L <= 300 for c0 = 1, 2, 3 (the lane reports
     1.0032546 at (7,50), 1.0000889 at (9,58), 1 + 1.99e-6 at (12,77)).
  4. Theorem R_small on a finite range: N_m <= 1.1562 A_(m+1), 1.0474 A_(m+2), 1.0148 A_(m+3) for
     2 <= m <= 24, L <= 250 (the theorem asserts every L).
"""
from fractions import Fraction


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg, flush=True)


LMAX = 400
P3 = [3 ** i for i in range(LMAX + 80)]
P2 = [2 ** i for i in range(LMAX + 80)]


def pw3(e):
    return P3[e] if e >= 0 else None


def N_barrier(m, L):
    """reflected barrier: state U (up count), time t; S' = U - c t; zone [m-1, m-c) (both letters -> S' - c)."""
    cur = {0: 1}
    out = [1]
    for t in range(0, L):
        nxt = {}
        for U, cnt in cur.items():
            a = U - m + 1
            in_zone = (a >= 0 and P3[a] >= P2[t]) and ((((U - m < 0) or P3[U - m] < P2[t - 1])) if t >= 1 else (U - m < 0))
            if in_zone:
                cands = ((U, 2 * cnt),)
            else:
                cands = ((U + 1, cnt), (U, cnt))
            for U2, w in cands:
                if P3[U2] > P2[t + 1]:
                    nxt[U2] = nxt.get(U2, 0) + w
        cur = nxt
        out.append(sum(cur.values()))
    return out


def A_hard(M, L):
    """words with S_t > 0 (1 <= t <= L) and S_t < M (0 <= t <= L-1), S_t = e_t - c t."""
    cur = {0: 1}
    out = [1]
    for t in range(0, L):
        nxt = {}
        for e, cnt in cur.items():
            if not (e - M < 0 or P3[e - M] < P2[t]):
                continue
            for b in (0, 1):
                e2 = e + b
                if P3[e2] > P2[t + 1]:
                    nxt[e2] = nxt.get(e2, 0) + cnt
        cur = nxt
        out.append(sum(cur.values()))
    return out


# 1. N_1 = 0
n1 = N_barrier(1, 50)
check(all(v == 0 for v in n1[1:]), "N_1(L) = 0 for 1 <= L <= 50 (the trivial case m = 1)")

# precompute
Ncache, Acache = {}, {}
for m in range(2, 25):
    Ncache[m] = N_barrier(m, LMAX if m <= 14 else 250)
for M in range(3, 30):
    Acache[M] = A_hard(M, LMAX if M <= 29 else 250)

# 2. R_15 finite range
viol = 0
worst = Fraction(0)
for m in range(2, 15):
    N = Ncache[m]
    A = Acache[m + 15]
    for L in range(1, LMAX + 1):
        if 10000 * N[L] > 10039 * A[L]:
            viol += 1
        if A[L]:
            worst = max(worst, Fraction(N[L], A[L]))
check(viol == 0, f"Theorem R_15 on 2 <= m <= 14, L <= {LMAX}: 10000 N_m(L) <= 10039 A_(m+15)(L) exactly (max ratio 1 + {float(worst - 1):.2e})")

# 3. ratio profile
prof = {}
for c0 in (1, 2, 3):
    best, arg = Fraction(0), None
    for m in range(2, 15):
        N = Ncache[m]
        A = Acache[m + c0]
        for L in range(1, 301):
            if A[L] and Fraction(N[L], A[L]) > best:
                best, arg = Fraction(N[L], A[L]), (m, L)
    prof[c0] = (float(best), arg)
check(abs(prof[1][0] - 1.0032546) < 5e-7 and prof[1][1] == (7, 50) and abs(prof[2][0] - 1.0000889) < 5e-7 and prof[2][1] == (9, 58)
      and abs(prof[3][0] - (1 + 1.99e-6)) < 5e-8,
      "ratio profile (2 <= m <= 14, L <= 300): " + "; ".join(f"c0={c}: max {v[0]:.7f} at (m,L)={v[1]}" for c, v in prof.items()))

# 4. R_small finite range
consts = {1: Fraction(11562, 10000), 2: Fraction(10474, 10000), 3: Fraction(10148, 10000)}
bad = 0
for m in range(2, 25):
    N = Ncache[m]
    for c0, K in consts.items():
        A = Acache[m + c0]
        for L in range(1, min(len(N), len(A))):
            if N[L] > K * A[L]:
                bad += 1
check(bad == 0, "Theorem R_small on 2 <= m <= 24, L <= 250 (m <= 14: L <= 400): N_m <= 1.1562 A_(m+1), 1.0474 A_(m+2), 1.0148 A_(m+3)")
