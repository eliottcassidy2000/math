#!/usr/bin/env python3
"""
klein-2026-07-07-S155 (part 3) -- SIGNED Hunter-tree floor at the k=8 criticality
(fixes part-2's sign bug: for negative differences d the hit event is
frac(|d|x) in [1-theta, 1), not (0, theta]; interior vertices mix signs).

W_i = {x : no orbit point in (P_i, P_i + theta]} = INTERSECT_{d in D_i} {frac(d x) not in (0,theta]}
with D_i = {e_a - e_i : a != i} SIGNED. hit_d = {frac(|d|x) in (0,theta]} if d>0,
{frac(|d|x) in [1-theta,1)} if d<0. Each hit has measure exactly theta (one-sided sieve).
At k=8, theta=1/7: Bonferroni base 1 - 7theta = 0, so Hunter's spanning-tree bound gives
   meas W_i >= max_{spanning tree T} Sum_{(a,b) in T} m(d_a, d_b)   [RIGOROUS]
with m(d_a,d_b) = meas(hit_a ∩ hit_b) EXACT via interval intersections.

ALSO: extended universal pair-mass scan (signed, q2 <= 200) for the floor lemma candidate
      m(d_a, d_b) >= theta^2 = 1/49 (equality cases?), and same-sign vs mixed-sign minima.
"""
import numpy as np
from fractions import Fraction as F
from math import gcd
from itertools import combinations

TH = F(1, 7)

def m_signed_exact(d1, d2, theta: F):
    """EXACT meas(hit_{d1} ∩ hit_{d2}) for SIGNED nonzero integers d1, d2.
    hit_d = {frac(|d| x) in W_sign} with W_+ = (0,theta], W_- = [1-theta, 1).
    Reduce by g=gcd(|d1|,|d2|) (x -> x/g measure preserving): q_i = |d_i|/g.
    Intervals: A_j = ((j+o1)/q1, (j+o1+theta)/q1] where o1 = 0 for +, 1-theta for -."""
    g = gcd(abs(d1), abs(d2)); q1, q2 = abs(d1)//g, abs(d2)//g
    o1 = F(0) if d1 > 0 else 1 - theta
    o2 = F(0) if d2 > 0 else 1 - theta
    tot = F(0)
    for j in range(q1):
        a_lo, a_hi = (j + o1)/q1, (j + o1 + theta)/q1
        l_min = int((a_lo*q2 - o2 - theta)) - 1
        l_max = int((a_hi*q2 - o2)) + 1
        for l in range(l_min, l_max + 1):
            b_lo, b_hi = (l + o2)/q2, (l + o2 + theta)/q2
            lo, hi = max(a_lo, b_lo), min(a_hi, b_hi)
            if hi > lo:
                # clip to [0,1)
                lo2, hi2 = max(lo, F(0)), min(hi, F(1))
                if hi2 > lo2: tot += hi2 - lo2
    return tot

def hunter_floor_vertex(E, i_elem, theta: F):
    D = [e - i_elem for e in E if e != i_elem]
    n = len(D)
    edges = []
    for a in range(n):
        for b in range(a+1, n):
            edges.append((m_signed_exact(D[a], D[b], theta), a, b))
    edges.sort(key=lambda t: t[0], reverse=True)
    parent = list(range(n))
    def find(u):
        while parent[u] != u: parent[u] = parent[parent[u]]; u = parent[u]
        return u
    tot = F(0); used = 0
    for w, a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb; tot += w; used += 1
            if used == n-1: break
    return tot

def grid(NG): return (np.arange(NG) + 0.5)/NG

def true_W(E, i_elem, x, theta=1/7):
    ok = np.ones_like(x, bool)
    for e in E:
        if e == i_elem: continue
        d = e - i_elem
        v = (abs(d)*x) % 1.0
        if d > 0: ok &= ~((v > 0) & (v <= theta))
        else:     ok &= ~(v >= 1-theta)
    return ok.mean()

if __name__ == "__main__":
    rng = np.random.default_rng(31555)
    x = grid(60013)
    print("=== signed engine validation ===")
    for (d1, d2) in [(1, -1), (1, 2), (-3, 4), (2, -5), (6, 7), (-6, -7)]:
        ex = m_signed_exact(d1, d2, TH)
        v1 = (abs(d1)*x) % 1.0; v2 = (abs(d2)*x) % 1.0
        h1 = ((v1 > 0) & (v1 <= 1/7)) if d1 > 0 else (v1 >= 6/7)
        h2 = ((v2 > 0) & (v2 <= 1/7)) if d2 > 0 else (v2 >= 6/7)
        print(f"  (d1,d2)=({d1:>3},{d2:>3}): exact = {str(ex):>8} = {float(ex):.6f}   numeric = {(h1&h2).mean():.6f}")

    print("\n=== SIGNED Hunter-tree floors at k=8 (rigorous: meas W_i >= tree mass; base 1-7theta=0) ===")
    bank = {
        "AP {1..8}": [1,2,3,4,5,6,7,8],
        "perforated {0,2..8}": [0,2,3,4,5,6,7,8],
        "spread": [0,5,11,17,26,33,41,50],
        "two-cluster 4+4": [0,1,2,3,100,101,102,103],
        "adv-from-part2": [0,4,7,9,21,29,30,36],
    }
    for nm_, E in bank.items():
        best = (F(-1), None, 0.0)
        for i in E:
            tm = hunter_floor_vertex(E, i, TH)
            if tm > best[0]:
                best = (tm, i, true_W(E, i, x))
        viol = "  [VIOLATION!]" if float(best[0]) > best[2] + 2e-3 else ""
        print(f"  {nm_:>22}: best i={best[1]}: Hunter = {float(best[0]):.4f} (exact {best[0]})   true W_i = {best[2]:.4f}{viol}")

    print("\n=== adversarial min of best-vertex SIGNED Hunter floor (k=8, jump moves) ===")
    def hbest(EE):
        return max(hunter_floor_vertex(EE, i, TH) for i in EE)
    gmin = (F(10), None)
    for trial in range(18):
        H = int(rng.choice([9, 14, 22, 34]))
        E = sorted(rng.choice(np.arange(0, H+1), size=8, replace=False).tolist())
        cur = hbest(E)
        for step in range(22):
            ii = int(rng.integers(8)); new = int(rng.integers(0, int(rng.choice([12, 24, 36]))+1))
            if new in E: continue
            cand = sorted(set(E) - {E[ii]} | {new})
            if len(cand) != 8: continue
            c = hbest(cand)
            if c < cur: E, cur = cand, c
        if cur < gmin[0]: gmin = (cur, tuple(E))
    print(f"  min best-vertex Hunter floor = {float(gmin[0]):.4f} (exact {gmin[0]}) at {gmin[1]}")
    print(f"  => rigorous uniform-candidate: mu_1/7(E_8) >= meas W_i >= this (adversarial; sup proof open)")

    print("\n=== extended universal pair-mass scan (signed; is min m = theta^2 = 1/49 ~ 0.0204?) ===")
    worst_ss = (F(10), None); worst_ms = (F(10), None)
    for q1 in range(1, 11):
        for q2 in range(q1, 201):
            if gcd(q1, q2) != 1: continue
            if q1 == q2 and q1 != 1: continue
            vss = m_signed_exact(q1, q2, TH)   # same sign
            if q1 != q2 and vss < worst_ss[0]: worst_ss = (vss, (q1, q2))
            vms = m_signed_exact(q1, -q2, TH)  # mixed sign
            if vms < worst_ms[0]: worst_ms = (vms, (q1, -q2))
    print(f"  same-sign  min m = {worst_ss[0]} = {float(worst_ss[0]):.6f} at {worst_ss[1]}")
    print(f"  mixed-sign min m = {worst_ms[0]} = {float(worst_ms[0]):.6f} at {worst_ms[1]}")
    print(f"  (note (1,-1): m = 0 exactly -- antipodal windows are disjoint; excluded by q1!=q2 in same-sign scan;")
    print(f"   the Hunter tree must avoid +-same-|d| pairs, or use them knowingly as 0-weight edges)")
