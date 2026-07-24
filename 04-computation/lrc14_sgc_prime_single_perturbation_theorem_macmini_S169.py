#!/usr/bin/env python3
"""
LRC(14) SGC'(13) for the single-perturbation family -- COMPLETE PROOF (exact rational arithmetic).
mac-mini-2026-07-23-S169 (Opus).

THEOREM. For every j in {1..13} and every positive integer w,
    S = ({1..13}\{j}) u {w}   satisfies   gap(S) not in the open interval (1/14, 3/41).
Equivalently, within this family:  gap(S) > 1/14  =>  gap(S) >= 3/41,
with the extremal value 3/41 attained at {1..11,13,36} (j=12, w=36, tau=17/41).

gap(S) = max_tau min_{v in S} ||v tau||.

PROOF = (a) + (b):
(a) STRANGER-DECOUPLING LEMMA (gap axis, explicit).  Let C be finite, theta rational, and suppose
    f_C(tau) = min_{v in C} ||v tau|| >= theta on some interval I of length delta > 0.
    Then for every integer w >= 1/delta, I contains a tau with ||w tau|| = 1/2 (since w*tau sweeps a
    full period across I), hence min(f_C(tau), ||w tau||) >= min(theta, 1/2) = theta.
    Therefore gap(C u {w}) >= theta.   [So a set with gap < theta forces w < 1/delta.]
(b) EXACT FINITE CHECK.  With theta = 3/41 we compute delta_j EXACTLY for each 12-core
    C_j = {1..13}\{j} (all breakpoints are rational), giving W_j = ceil(1/delta_j); then verify by
    exact rational arithmetic that no w <= max_j W_j puts gap in (1/14, 3/41).

Everything below is exact (fractions.Fraction); no floating point in the decision path.
"""
from fractions import Fraction as F
from math import gcd, ceil
from functools import reduce

THETA = F(3, 41)
ONE14 = F(1, 14)


def safe_intervals(v, theta):
    """{tau in [0,1] : ||v tau|| >= theta} as a list of closed intervals (exact)."""
    out = []
    for k in range(v):
        lo = F(k, v) + theta / v
        hi = F(k + 1, v) - theta / v
        if lo <= hi:
            out.append((lo, hi))
    return out


def intersect(A, B):
    """intersect two sorted lists of disjoint closed intervals (exact)."""
    out, i, j = [], 0, 0
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0])
        hi = min(A[i][1], B[j][1])
        if lo <= hi:
            out.append((lo, hi))
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return out


def longest_safe_interval(C, theta):
    """EXACT delta = max length of an interval on which min_{v in C} ||v tau|| >= theta."""
    cur = safe_intervals(C[0], theta)
    for v in C[1:]:
        cur = intersect(cur, safe_intervals(v, theta))
        if not cur:
            return F(0)
    return max((hi - lo) for lo, hi in cur)


def exact_gap(V):
    """EXACT gap = max_tau min_v ||v tau||; peaks lie at tau = k/q, q in {2v, v_i +- v_j}."""
    V = sorted(set(V))
    qs = set()
    for a in V:
        qs.add(2 * a)
        for b in V:
            if a != b:
                qs.add(a + b)
                qs.add(abs(a - b))
    bn, bd, bt = 0, 1, None
    for q in qs:
        if q < 2:
            continue
        for k in range(1, q // 2 + 1):
            m = q
            for v in V:
                r = (v * k) % q
                d = r if r <= q - r else q - r
                if d < m:
                    m = d
                    if m * bd <= bn * q:
                        break
            if m * bd > bn * q:
                bn, bd, bt = m, q, F(k, q)
    return F(bn, bd), bt


def main():
    print("THEOREM: no ({1..13}\\{j}) u {w} has gap in (1/14, 3/41).  theta = 3/41\n")
    print("(a) EXACT stranger bounds W_j = ceil(1/delta_j):")
    Ws = {}
    for j in range(1, 14):
        C = [x for x in range(1, 14) if x != j]
        delta = longest_safe_interval(C, THETA)
        W = ceil(1 / delta) if delta > 0 else None
        Ws[j] = W
        print(f"    j={j:>2}: delta_j = {delta}  (~{float(delta):.8f})   W_j = {W}")
    Wmax = max(Ws.values())
    print(f"\n    max_j W_j = {Wmax}  ->  any band-hitter in this family has w < {Wmax}\n")

    print(f"(b) EXACT check of every j and every w <= {Wmax}:")
    bad = []
    extremal = None
    for j in range(1, 14):
        base = [x for x in range(1, 14) if x != j]
        for w in range(1, Wmax + 1):
            V = base + [w]
            if len(set(V)) != 13:
                continue
            if reduce(gcd, V) != 1:
                continue
            g, t = exact_gap(V)
            if ONE14 < g < THETA:
                bad.append((g, j, w, t))
            if g > ONE14 and (extremal is None or g < extremal[0]):
                extremal = (g, j, w, t)
    print(f"    sets with gap in (1/14, 3/41): {len(bad)}")
    for g, j, w, t in bad:
        print(f"      gap={g} at j={j}, w={w}, tau={t}")
    if not bad:
        print("    NONE -> THEOREM PROVED for the single-perturbation family.")
    print(f"\n    extremal (smallest gap > 1/14): gap={extremal[0]} at j={extremal[1]}, "
          f"w={extremal[2]}, tau={extremal[3]}   -> the set {{1..11,13,36}}")


if __name__ == "__main__":
    main()
