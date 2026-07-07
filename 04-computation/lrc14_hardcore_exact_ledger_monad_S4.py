#!/usr/bin/env python3
"""
monad-explorer-2026-07-07-S4 -- EXACT G_P-intersected ledgers at the two hard cores
(HYP-4907; the floor's only load-bearing mixed shapes, per HYP-4917 window factoring).

Engine: general-E exact G2(P,E) = meas(G_P /\ {maxgap({e x}) > theta}) by Farey_{2S}-cell
piecewise-affine maxgap (per-cell domination asserted exactly; S = max(E)) intersected with
exact G_P intervals.  Extends my S2 AP engine (validated vs canon) to general E + G_P;
sanity row: G2(empty, {0..12}) must equal 477/1078.

Outputs per hard core P in {(9,10,11,12,13) [k=8], (10,11,12,13) [k=9]}:
  * EXHAUSTIVE exact min over co-offset sets E = {0} u (k-1 from {1..S}), S <= 14
  * exact G2, plain mu, quasi-independence R = G2/(meas(G_P)*mu)  [R >= 0.75 watch]
  * mu at theta = 3/14 (the corrected W-inequality threshold; "mu3 tables")
  * the balanced short relation vectors of the argmin shapes (d-perp probe: where the
    deficit sits for the floor-provers)
Cross-check row: kps-S60's worst small-p part P={1,2,3,4,5} (their ILedger weak spot;
my S3 mid-grid q=6 window covers it -- complementarity datum).

Tournament Analysis declaration:
  vertices: ledger rows (P, E-shape classes); pairwise observable: exact G2 vs m_P;
  switch/gauge: min row first; tie path: exact engine -> exhaustive sweep -> short vectors.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import sys

M_P = F(14249, 252252)

def farey(n):
    seq = [F(0, 1), F(1, n)]
    a, b, c, d = 0, 1, 1, n
    while c < n or (c == n and d == 1):
        k = (n + b) // d
        a, b, c, d = c, d, k * c - a, k * d - b
        seq.append(F(c, d))
        if c == 1 and d == 1:
            break
    return seq

_farey_cache = {}
def good_intervals(E, theta):
    """Exact maximal intervals of {x in (0,1): maxgap({e x mod 1 : e in E}) > theta}."""
    E = sorted(set(E)); m = len(E)
    S = max(max(E), 1)
    K = 2 * S if S > 1 else 3
    if K not in _farey_cache:
        _farey_cache[K] = farey(K)
    cells = _farey_cache[K]
    out = []
    for a, b in zip(cells, cells[1:]):
        mid = (a + b) / 2
        fl = [int(e * mid) if e else 0 for e in E]
        pts = sorted(range(m), key=lambda i: E[i] * mid - fl[i])
        gaps = []
        for s in range(m):
            i1, i2 = pts[s], pts[(s + 1) % m]
            sl = E[i2] - E[i1]
            ic = F(-(fl[i2] - fl[i1]) + (1 if s == m - 1 else 0))
            gaps.append((sl, ic))
        vm = [sl * mid + ic for sl, ic in gaps]
        ia = max(range(m), key=lambda i: vm[i])
        sa, ca = gaps[ia]
        va = [sl * a + ic for sl, ic in gaps]
        vb = [sl * b + ic for sl, ic in gaps]
        assert all(sa * a + ca >= va[i] for i in range(m)), ("dom-left", E, a, b)
        assert all(sa * b + ca >= vb[i] for i in range(m)), ("dom-right", E, a, b)
        ga, gb = sa * a + ca - theta, sa * b + ca - theta
        if ga <= 0 and gb <= 0:
            continue
        if ga > 0 and gb > 0:
            lo, hi = a, b
        else:
            xs = (theta - ca) / F(sa)
            lo, hi = (a, xs) if ga > 0 else (xs, b)
        if out and lo <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], hi))
        else:
            out.append((lo, hi))
    return out

def gp_intervals(P):
    if not P:
        return [(F(0), F(1))]
    bad = []
    for p in P:
        w = F(1, 14 * p)
        for j in range(p + 1):
            bad.append((max(F(j, p) - w, F(0)), min(F(j, p) + w, F(1))))
    bad.sort()
    merged = []
    for lo, hi in bad:
        if merged and lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    good, cur = [], F(0)
    for lo, hi in merged:
        if lo > cur:
            good.append((cur, lo))
        cur = max(cur, hi)
    if cur < 1:
        good.append((cur, F(1)))
    return good

def isect_measure(A, B):
    total = F(0)
    i = j = 0
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0]); hi = min(A[i][1], B[j][1])
        if lo < hi:
            total += hi - lo
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return total

def measure(A):
    return sum(hi - lo for lo, hi in A)

def balanced_short_vectors(E, l1max=6):
    """Relation-lattice vectors m: sum m_i e_i = 0, sum m_i = 0, 0 < |m|_1 <= l1max,
    support <= 4 (enumerate small supports)."""
    E = sorted(set(E)); out = []
    idx = range(len(E))
    for sup in range(3, 5):
        for T in combinations(idx, sup):
            # solve small integer combos
            from itertools import product
            rng = range(-3, 4)
            for m in product(rng, repeat=sup):
                if all(x == 0 for x in m) or sum(abs(x) for x in m) > l1max:
                    continue
                if sum(m) != 0:
                    continue
                if sum(mi * E[t] for mi, t in zip(m, T)) != 0:
                    continue
                if next(x for x in m if x != 0) < 0:
                    continue
                out.append(tuple((E[t], mi) for t, mi in zip(T, m) if mi))
    return sorted(set(out))[:12]

if __name__ == "__main__":
    THETA = F(1, 7)
    # sanity: G2(empty, {0..12}) == mu(AP_13) = 477/1078
    g = good_intervals(list(range(13)), THETA)
    assert measure(g) == F(477, 1078), measure(g)
    print("sanity: G2(empty,{0..12}) = 477/1078 MATCH")
    print()
    for P, k, S in [((9, 10, 11, 12, 13), 8, 14), ((10, 11, 12, 13), 9, 13),
                    ((1, 2, 3, 4, 5), 8, 12)]:
        gp = gp_intervals(P)
        mgp = measure(gp)
        kk = 13 - len(P)
        rows = []
        n_done = 0
        for extra in combinations(range(1, S + 1), kk - 1):
            E = (0,) + extra
            gi = good_intervals(list(E), THETA)
            g2 = isect_measure(gp, gi)
            mu = measure(gi)
            rows.append((g2, mu, E))
            n_done += 1
        rows.sort()
        g2min, mu_at, Emin = rows[0]
        R = g2min / (mgp * mu_at) if mu_at else F(0)
        mu3 = measure(good_intervals(list(Emin), F(3, 14)))
        print(f"P={P} (k={k}) meas(G_P)={mgp}~{float(mgp):.5f}  shapes={n_done} (S<={S})")
        print(f"  EXACT min G2 = {g2min} ~ {float(g2min):.5f}  ({float(g2min/M_P):.2f}x m_P)"
              f"  at E={Emin}")
        print(f"  plain mu there = {mu_at} ~ {float(mu_at):.5f};  R = G2/(mGP*mu) ~ {float(R):.4f}")
        print(f"  mu_(3/14) at argmin shape = {mu3} ~ {float(mu3):.5f}")
        worst5 = [(float(r[0]), r[2]) for r in rows[:5]]
        print(f"  five lowest G2 shapes: {worst5}")
        sv = balanced_short_vectors(Emin)
        print(f"  balanced short relation vectors of argmin (d-perp probe): {sv}")
        print()
    print("DONE.")
