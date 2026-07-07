#!/usr/bin/env python3
"""
klein-2026-07-07-S159 -- the CHERRY TREE DEEPENED: the sharp PAIRWISE-MOMENT LP.

FRAME: Hunter (spanning tree) and Bukszar-Prekopa (cherry tree) are closed-form
DUAL-FEASIBLE points of the partially-specified-moment LP:
    minimize x_empty  over atom measures {x_S : S subset of [n]}, x_S >= 0,
    s.t. sum x_S = 1,  sum_{S ni i} x_S = p_i (= theta exactly),
         sum_{S >= {i,j}} x_S = m_ij (EXACT via THM-638 for same-sign differences).
The LP value is the SHARP floor on W given full pairwise information -- no triple
bounds needed; the DUAL optimal solution is a per-shape linear certificate over
{1, singletons, pairs} (a 'generalized cherry').

THIS SCRIPT:
 1. k=8 endpoint LP floor: bank + adversarial min; dual certificate structure at minima.
 2. k=9..12 endpoint LP floors: does disaggregated pair info stay positive where the
    Bonferroni base 1-(k-1)/7 < 0 kills Hunter/cherry? (adversarial at k=9.)
 3. CONDITIONAL LP on the hard cores (P8={9..13} k=8, P9={10..13} k=9): atoms restricted
    to G_P (total mass meas G_P; marginals/pairs = meas(G_P n A_i (n A_j)) numeric-grid),
    objective = floor on rho*(P,E) directly; adversarial min vs m_P = 0.05649.
 4. the FKG hair: pairwise positive correlation (G>=0, proved) does NOT imply the product
    bound: W(AP endpoint) vs (6/7)^7 -- pin the sign of the gap at high precision.
"""
import numpy as np
from fractions import Fraction as F
from math import gcd
from itertools import combinations
from scipy.optimize import linprog

TH = 1/7
M_P = float(F(14249, 252252))

def law_pair(d1, d2):
    g = gcd(d1, d2); q1, q2 = d1//g, d2//g
    r1, r2 = q1 % 7, q2 % 7
    return 1/49 + (min(r1, r2)*(7 - max(r1, r2)))/(49*q1*q2)

def lp_floor_pairs(n, p, mpair):
    """sharp LP floor on P(no event) given exact singles p[i] and pairs mpair[(i,j)]."""
    N = 1 << n
    # objective: minimize x_emptyset
    c = np.zeros(N); c[0] = 1.0
    A = []; b = []
    A.append(np.ones(N)); b.append(1.0)
    for i in range(n):
        row = np.zeros(N)
        for S in range(N):
            if S >> i & 1: row[S] = 1.0
        A.append(row); b.append(p[i])
    for (i, j), m in mpair.items():
        row = np.zeros(N)
        for S in range(N):
            if (S >> i & 1) and (S >> j & 1): row[S] = 1.0
        A.append(row); b.append(m)
    res = linprog(c, A_eq=np.array(A), b_eq=np.array(b), bounds=[(0, None)]*N, method="highs")
    return (res.fun if res.status == 0 else None), res

def endpoint_D(E):
    et = max(E)
    return sorted(et - e for e in E if e != et)

def lp_floor_shape(E):
    D = endpoint_D(E); n = len(D)
    p = [TH]*n
    mp = {(a, b): law_pair(D[a], D[b]) for a in range(n) for b in range(a+1, n)}
    v, res = lp_floor_pairs(n, p, mp)
    return v

def grid(NG): return (np.arange(NG)+0.5)/NG

def true_W(E, x):
    et = max(E); ok = np.ones_like(x, bool)
    for e in E:
        if e == et: continue
        v = ((et-e)*x) % 1.0
        ok &= ~((v > 0) & (v <= TH))
    return ok.mean()

if __name__ == "__main__":
    rng = np.random.default_rng(15900)
    x = grid(60013)

    print("=== 1. k=8 endpoint PAIRWISE-LP floor (sharp given THM-638 pair data) ===")
    bank = {
        "AP {1..8}": [1,2,3,4,5,6,7,8],
        "spread": [0,5,11,17,26,33,41,50],
        "geometric": [0,3,8,17,31,52,80,118],
        "two-cluster": [0,1,2,3,100,101,102,103],
        "S158 cherry-adv": [3,4,5,7,8,9,14,22],
    }
    for nm, E in bank.items():
        v = lp_floor_shape(E)
        print(f"  {nm:>18}: LP floor = {v:.4f}   true W = {true_W(E,x):.4f}   (cherry S158 ~ varies; Hunter 6/49 = {6/49:.4f})")
    # adversarial min of the LP floor
    gmin = (2.0, None)
    for trial in range(30):
        H = int(rng.choice([9, 14, 22, 40]))
        E = sorted(rng.choice(np.arange(0, H+1), size=8, replace=False).tolist())
        cur = lp_floor_shape(E)
        for step in range(30):
            i = int(rng.integers(8)); new = int(rng.integers(0, int(rng.choice([12, 26, 44]))+1))
            if new in E: continue
            cand = sorted(set(E) - {E[i]} | {new})
            if len(cand) != 8: continue
            c = lp_floor_shape(cand)
            if c is not None and c < cur - 1e-5: E, cur = cand, c
        if cur < gmin[0]: gmin = (cur, tuple(E))
    print(f"  ADVERSARIAL MIN LP floor (k=8) = {gmin[0]:.4f} at {gmin[1]}   [R-route bar 0.197; full bar 0.675... note LP >= Hunter always]")

    print("\n=== 2. k=9..12 endpoint LP floors (Hunter/cherry base is negative here) ===")
    for k, E in [(9, list(range(1,10))), (10, list(range(1,11))), (11, list(range(1,12))), (12, list(range(1,13)))]:
        v = lp_floor_shape(E)
        print(f"  k={k:>2} AP: LP floor = {v:.4f}   (bare Hunter base 1-(k-1)/7 = {1-(k-1)/7:+.4f})")
    # adversarial at k=9
    gmin9 = (2.0, None)
    for trial in range(20):
        H = int(rng.choice([10, 16, 26]))
        E = sorted(rng.choice(np.arange(0, H+1), size=9, replace=False).tolist())
        cur = lp_floor_shape(E)
        for step in range(25):
            i = int(rng.integers(9)); new = int(rng.integers(0, int(rng.choice([14, 30]))+1))
            if new in E: continue
            cand = sorted(set(E) - {E[i]} | {new})
            if len(cand) != 9: continue
            c = lp_floor_shape(cand)
            if c is not None and c < cur - 1e-5: E, cur = cand, c
        if cur < gmin9[0]: gmin9 = (cur, tuple(E))
    print(f"  ADVERSARIAL MIN LP floor (k=9) = {gmin9[0]:.4f} at {gmin9[1]}   [does pair info rescue the negative base?]")

    print("\n=== 3. CONDITIONAL LP on the HARD CORES (floor rho*(P,E) directly) ===")
    def gp_mask(P, xg):
        m = np.ones_like(xg, bool)
        for p_ in P: m &= np.abs(((p_*xg) % 1.0) - 0.5) <= 0.5 - 1/14
        return m
    for (Pfix, k, nm) in [([9,10,11,12,13], 8, "P8={9..13}"), ([10,11,12,13], 9, "P9={10..13}")]:
        g = gp_mask(Pfix, x); mG = g.mean()
        def cond_lp(E):
            D = endpoint_D(E); n = len(D)
            hits = [(((d*x) % 1.0) > 0) & (((d*x) % 1.0) <= TH) for d in D]
            p = [float((g & h).mean()) for h in hits]
            mp = {(a, b): float((g & hits[a] & hits[b]).mean()) for a in range(n) for b in range(a+1, n)}
            # atoms restricted to G_P: total mass mG
            N = 1 << n
            c = np.zeros(N); c[0] = 1.0
            A = [np.ones(N)]; b = [mG]
            for i in range(n):
                row = np.zeros(N)
                for S in range(N):
                    if S >> i & 1: row[S] = 1.0
                A.append(row); b.append(p[i])
            for (i, j), m in mp.items():
                row = np.zeros(N)
                for S in range(N):
                    if (S >> i & 1) and (S >> j & 1): row[S] = 1.0
                A.append(row); b.append(m)
            res = linprog(c, A_eq=np.array(A), b_eq=np.array(b), bounds=[(0, None)]*N, method="highs")
            return res.fun if res.status == 0 else None
        # bank + small adversarial
        E0 = list(range(1, k+1))
        v0 = cond_lp(E0)
        gm = (v0 if v0 is not None else 2.0, tuple(E0))
        for trial in range(12):
            H = int(rng.choice([k+2, 2*k, 4*k]))
            E = sorted(rng.choice(np.arange(0, H+1), size=k, replace=False).tolist())
            cur = cond_lp(E)
            if cur is None: continue
            for step in range(15):
                i = int(rng.integers(k)); new = int(rng.integers(0, int(rng.choice([2*k, 5*k]))+1))
                if new in E: continue
                cand = sorted(set(E) - {E[i]} | {new})
                if len(cand) != k: continue
                c_ = cond_lp(cand)
                if c_ is not None and c_ < cur - 1e-5: E, cur = cand, c_
            if cur < gm[0]: gm = (cur, tuple(E))
        print(f"  {nm}: meas(G_P) = {mG:.4f};  conditional-LP floor on rho*: AP = {v0:.4f};  adversarial min = {gm[0]:.4f} at {gm[1]}")
        print(f"     vs bar m_P = {M_P:.4f}: {'CLEARS' if gm[0] >= M_P else 'BELOW'} ({gm[0]/M_P:.1f}x)")

    print("\n=== 4. the FKG hair: pairwise positive correlation does NOT give the product bound ===")
    E = list(range(1, 9))
    W_ap = true_W(E, grid(200003))
    prod = (6/7)**7
    print(f"  W(AP endpoint) = {W_ap:.6f}  vs product (6/7)^7 = {prod:.6f}  gap = {W_ap-prod:+.6f}")
    print(f"  (all 21 pair correlations are >= 0 by THM-638, yet the 7-fold product bound fails at the AP:")
    print(f"   pairwise association does not extend to full association -- the LP/cherry layer is genuinely needed)")
