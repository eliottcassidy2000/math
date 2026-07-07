#!/usr/bin/env python3
"""
klein-2026-07-07-S155 (part 2) -- the TOURNAMENT/PHASE-CLOCK bridge and the HUNTER-TREE
floor at the k=8 criticality.

FRAME (Tournament Analysis, THM-373/374/381 lineage): at time x the orbit {frac(e_i x)}
carries the theta-successor digraph D_theta(x): i -> j iff frac((e_j - e_i)x) in (0, theta].
 - Each edge indicator is a PAIR-DISTANCE event: frac(d x) in (0, theta], d = e_j - e_i,
   and each pair distance is EXACTLY uniform in x (the kps-S59 forgotten factoid).
   => E[#edges] = k(k-1)*theta EXACTLY for EVERY integer set E.  [proved, one line]
 - B(x) := #{i : outdeg_i = 0} = #gaps > theta  (a vertex heads a big gap iff no
   successor within theta).  mu_theta(E) = P(B >= 1) exactly.
 - At (k, theta) = (8, 1/7): E[outdeg_i] = 7*(1/7) = 1 -- the CRITICAL mean-degree-1
   regime. (k<=7: subcritical, and pigeonhole makes B>=1 deterministic; k=13: 12/7
   supercritical.)  The binding k=8 leg of THM-530 is the criticality.
 - The digraph movie x -> D_theta(x) is a THM-373-style comparator movie whose walls are
   at x = m/d (pair denominators): the wall count Sum |e_j - e_i| IS the piece count of
   the S154/S155 crossing lemma. Tournament lens = the crossing lemma's bookkeeping.

HUNTER-TREE FLOOR (rigorous per-shape, positive EXACTLY at k=8):
   W_i = {outdeg_i = 0} = intersection over d in D_i = {e_a - e_i} of {frac(dx) not in (0,theta]}.
   Hunter (1976) upper bound on the union of the 7 'hit' events with any spanning tree T:
     meas W_i >= 1 - 7*theta + Sum_{(a,b) in T} m_ab = Sum_{(a,b) in T} m_ab   at theta=1/7,
   where m_ab = meas{frac(d_a x) in (0,theta], frac(d_b x) in (0,theta]} -- an EXACT
   two-frequency lattice-line measure. Since 1 - 7*theta = 0, the tree mass IS the floor.

THIS SCRIPT:
 1. verify E[#edges] = 56/7 = 8 exactly (numeric, multiple shapes)
 2. exact m_ab engine (lattice line in a box) + verify against numeric
 3. Hunter-tree floors: max-spanning-tree mass per vertex; compare to true meas W_i;
    minimize the best-vertex Hunter floor adversarially over 8-sets (is it uniformly
    bounded below? conjecture: m_ab >= theta^2-ish universally)
 4. worst-case pair mass: min m_ab over coprime (q1,q2) grids (the universal-floor lemma
    candidate: m_ab as q2 -> infinity should -> theta^2; find the true min)
 5. PZ-on-B at k=8: E[B], E[B^2], ratio vs need_8 = 0.675 (bank + adversarial)
"""
import numpy as np
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations

TH = F(1, 7)
NEED8 = 1 - F(2243, 5880) + F(14249, 252252)

def grid(NG): return (np.arange(NG) + 0.5)/NG

def m_ab_exact(d1, d2, theta: F):
    """EXACT meas{frac(d1 x) in (0,theta], frac(d2 x) in (0,theta]}, d1,d2 >= 1.
    Reduce by g = gcd: x -> x/g measure-preserving => depends on (q1,q2) = (d1,d2)/g.
    Event: frac(q1 t) in (0,th], frac(q2 t) in (0,th], t uniform on [0,1).
    Solve exactly by interval intersection: for each j1 in 0..q1-1 the first condition is
    t in (j1/q1, j1/q1 + th/q1]; intersect with the q2-comb intervals."""
    g = gcd(d1, d2); q1, q2 = d1//g, d2//g
    tot = F(0)
    th = theta
    # intervals A: (j/q1, (j+th)/q1], j=0..q1-1 ; B: (l/q2, (l+th)/q2]
    for j in range(q1):
        a_lo, a_hi = F(j, q1), F(j, q1) + th/q1
        # overlapping l range: l/q2 < a_hi and (l+th)/q2 > a_lo
        l_min = int((a_lo*q2 - th)) - 1
        l_max = int((a_hi*q2)) + 1
        for l in range(max(0, l_min), min(q2-1, l_max) + 1):
            b_lo, b_hi = F(l, q2), F(l, q2) + th/q2
            lo, hi = max(a_lo, b_lo), min(a_hi, b_hi)
            if hi > lo: tot += hi - lo
    return tot

def w_and_pairs_numeric(E, i_elem, x, theta=1/7):
    D = [e - i_elem for e in E if e != i_elem]
    hits = [(((d*x) % 1.0) > 0) & (((d*x) % 1.0) <= theta) for d in D]
    H = np.stack(hits)                     # (7, NG)
    W = ~H.any(axis=0)
    return W.mean(), D, H

def max_tree_mass(D, theta: F):
    """max spanning tree on the 7 'hit' events with edge weight m_ab (exact); Kruskal."""
    n = len(D)
    edges = []
    for a in range(n):
        for b in range(a+1, n):
            edges.append((m_ab_exact(abs(D[a]), abs(D[b]), theta), a, b))
    edges.sort(reverse=True)
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

def gapcount_B(E, x, theta=1/7):
    A = np.asarray(E, float)
    P = np.sort((A[None, :] * x[:, None]) % 1.0, axis=1)
    G = np.diff(P, axis=1); wrap = P[:, :1] + 1.0 - P[:, -1:]
    Gall = np.concatenate([G, wrap], axis=1)
    return (Gall > theta).sum(axis=1)

if __name__ == "__main__":
    rng = np.random.default_rng(8155)
    x = grid(60013)

    print("=== 1. exact edge-count identity: E[#theta-edges] = k(k-1)theta (= 8 at k=8, theta=1/7) ===")
    for E in ([1,2,3,4,5,6,7,8], [0,3,7,12,19,27,44,60], [0,1,2,4,8,16,32,64]):
        cnt = 0.0
        for i in E:
            for j in E:
                if i == j: continue
                d = j - i
                v = (np.asarray(d, float)*x) % 1.0
                cnt += (((v > 0) & (v <= 1/7)).mean())
        print(f"  E={E}: measured E[#edges] = {cnt:.5f}   (exact 8)")

    print("\n=== 2. exact pair-mass engine m_ab vs numeric ===")
    for (d1, d2) in [(1, 2), (1, 6), (6, 7), (3, 11), (5, 12)]:
        ex = m_ab_exact(d1, d2, TH)
        v1 = ((d1*x) % 1.0); v2 = ((d2*x) % 1.0)
        nm = (((v1 > 0) & (v1 <= 1/7)) & ((v2 > 0) & (v2 <= 1/7))).mean()
        print(f"  (d1,d2)=({d1},{d2}): exact = {str(ex):>10} = {float(ex):.6f}   numeric = {nm:.6f}   theta^2 = {float(TH*TH):.6f}")

    print("\n=== 3. HUNTER-TREE floor at k=8 criticality: meas W_i >= max-tree pair mass (exact) ===")
    bank = {
        "AP {1..8}": [1,2,3,4,5,6,7,8],
        "perforated {0,2..8}": [0,2,3,4,5,6,7,8],
        "spread": [0,5,11,17,26,33,41,50],
        "two-cluster 4+4": [0,1,2,3,100,101,102,103],
    }
    for nm_, E in bank.items():
        best = (F(-1), None, 0.0)
        for i in E:
            Wm, D, _ = w_and_pairs_numeric(E, i, x)
            tm = max_tree_mass(D, TH)
            if tm > best[0]: best = (tm, i, Wm)
        print(f"  {nm_:>22}: best vertex i={best[1]}: Hunter floor = {float(best[0]):.4f} (exact {best[0]})   true meas W_i = {best[2]:.4f}")
    # adversarial: minimize the best-vertex Hunter floor over 8-sets
    xs = grid(4099)
    gmin = (F(10), None)
    for trial in range(20):
        H = int(rng.choice([9, 14, 22, 34]))
        E = sorted(rng.choice(np.arange(0, H+1), size=8, replace=False).tolist())
        def hfloor(EE):
            b = F(-1)
            for i in EE:
                D = [e - i for e in EE if e != i]
                t = max_tree_mass(D, TH)
                if t > b: b = t
            return b
        cur = hfloor(E)
        for step in range(25):
            ii = int(rng.integers(8)); new = int(rng.integers(0, int(rng.choice([12, 24, 36]))+1))
            if new in E: continue
            cand = sorted(set(E) - {E[ii]} | {new})
            if len(cand) != 8: continue
            c = hfloor(cand)
            if c < cur: E, cur = cand, c
        if cur < gmin[0]: gmin = (cur, tuple(E))
    print(f"  ADVERSARIAL MIN of best-vertex Hunter floor: {float(gmin[0]):.4f} (exact {gmin[0]}) at {gmin[1]}")
    print(f"  (positive uniformly? floor gives mu >= this rigorously via mu >= meas W_i)")

    print("\n=== 4. universal pair-mass floor: min m_ab over coprime pairs (q1<q2<=40) ===")
    worst = (F(10), None)
    for q1 in range(1, 13):
        for q2 in range(q1+1, 41):
            if gcd(q1, q2) != 1: continue
            v = m_ab_exact(q1, q2, TH)
            if v < worst[0]: worst = (v, (q1, q2))
    print(f"  min m_ab = {worst[0]} = {float(worst[0]):.6f} at (q1,q2)={worst[1]}   [theta^2 = {float(TH*TH):.6f}]")
    print(f"  => universal-floor lemma candidate: m_ab >= {worst[0]} for all coprime pairs (to be proved by monotonicity in q)")

    print("\n=== 5. PZ-on-B at k=8 vs need_8 = %.5f ===" % float(NEED8))
    for nm_, E in bank.items():
        B = gapcount_B(E, x)
        EB = B.mean(); EB2 = (B.astype(float)**2).mean()
        mu = (B >= 1).mean()
        print(f"  {nm_:>22}: E[B]={EB:.4f} E[B^2]={EB2:.4f} PZ={EB*EB/EB2:.4f}  mu={mu:.4f}")
    gmin2 = (2.0, None)
    for trial in range(24):
        H = int(rng.choice([9, 14, 22, 34]))
        E = sorted(rng.choice(np.arange(0, H+1), size=8, replace=False).tolist())
        def pz(EE, xg):
            B = gapcount_B(EE, xg); e1 = B.mean(); e2 = (B.astype(float)**2).mean()
            return e1*e1/e2 if e2 > 0 else 0.0
        cur = pz(E, xs)
        for step in range(40):
            ii = int(rng.integers(8)); new = int(rng.integers(0, int(rng.choice([12, 24, 36]))+1))
            if new in E: continue
            cand = sorted(set(E) - {E[ii]} | {new})
            if len(cand) != 8: continue
            c = pz(cand, xs)
            if c < cur - 1e-3: E, cur = cand, c
        v = pz(E, grid(15013))
        if v < gmin2[0]: gmin2 = (v, tuple(E))
    print(f"  ADVERSARIAL MIN PZ-on-B (k=8) = {gmin2[0]:.4f} at {gmin2[1]}   ({'>= need_8 OK' if gmin2[0] >= float(NEED8) else 'below need_8'})")
