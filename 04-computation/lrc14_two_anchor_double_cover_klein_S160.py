#!/usr/bin/env python3
"""
klein-2026-07-07-S160 -- the 2-ANCHOR TAIL through the Z/2 DOUBLE COVER + the
WEIGHTED-CHERRY dual extraction.

(1) DOUBLE-COVER IDENTITIES for boxeph's PA_2(E) = P(max(gap@0, gap@1/2) > 1/7):
    write r_a/l_a = right/left one-sided min distances at anchor a in {0, 1/2}.
    (i)  gap@0(2E, x) = 2*(min(r_0,r_1/2) + min(l_0,l_1/2))   [the cover collapses anchors]
         => gap@0(2E) <= 2*min(gap@0, gap@1/2) <= 2*max(...) pointwise.
    (ii) ALL-ODD family: orbit at x+1/2 = orbit at x shifted by 1/2 exactly
         => gap@{1/2}(E, x) = gap@0(E, x+1/2) => the two anchor-gap LAWS coincide.
    (iii) inclusion-exclusion: PA_2 = P(g0>t) + P(g1/2>t) - P(both>t).
    Verify (i),(ii) pointwise; measure the parity split of P(g1/2>t) for mixed families.

(2) INDEPENDENT CHECK of boxeph's 2-anchor inf table (adversarial, jump moves) vs
    T_k = {0.6185, 0.5057, 0.3956, 0.2747, 0.1429, 0.0565}.

(3) WEIGHTED-CHERRY DUALS: solve the k=8 pairwise LP at the AP and at the barrier shape
    {1,3,5,6,7,9,25,38}; extract dual weights (y0; y_i; y_ij); test the AP's dual as a
    FEASIBLE certificate on other shapes (feasibility is shape-independent — it is a
    property of the quadratic form; only the VALUE changes via the law masses).
"""
import numpy as np
from itertools import combinations
from math import gcd
from scipy.optimize import linprog

TH = 1/7
TK = {8: 0.6185, 9: 0.5057, 10: 0.3956, 11: 0.2747, 12: 0.1429, 13: 0.0565}

def grid(NG): return (np.arange(NG)+0.5)/NG

def anchor_gaps(E, x):
    """returns (gap@0, gap@1/2) arrays."""
    A = np.asarray(E, float)
    P = (A[None, :] * x[:, None]) % 1.0                 # (NG, k) raw fracs
    r0 = P.min(axis=1)                                   # right at 0 = min frac
    l0 = (1.0 - P).min(axis=1)                           # left at 0
    Q = (P - 0.5) % 1.0
    r1 = Q.min(axis=1)                                   # right at 1/2
    l1 = (1.0 - Q).min(axis=1)
    return r0 + l0, r1 + l1, (r0, l0, r1, l1)

def PA2(E, x, t=TH):
    g0, g1, _ = anchor_gaps(E, x)
    return float((np.maximum(g0, g1) > t).mean())

if __name__ == "__main__":
    rng = np.random.default_rng(16060)
    x = grid(60013)

    print("=== (1i) double-cover identity: gap@0(2E,x) = 2(min r + min l) ===")
    worst = 0.0
    for E in ([1,2,3,4,5,6,7,8], [1,3,8,15,22,40,41,55], [2,4,6,8,10,12,13,14,16,18,20,22,24]):
        g0_2E, _, _ = anchor_gaps([2*e for e in E], x)
        _, _, (r0, l0, r1, l1) = anchor_gaps(E, x)
        pred = 2*(np.minimum(r0, r1) + np.minimum(l0, l1))
        worst = max(worst, float(np.abs(g0_2E - pred).max()))
    print(f"  max pointwise |gap@0(2E) - 2(min r + min l)| over 3 shapes x 60013 pts: {worst:.2e}")

    print("\n=== (1ii) all-odd parity identity: law(gap@1/2) == law(gap@0) ===")
    for nm, E in [("all-odd {1,3,..,15}", [1,3,5,7,9,11,13,15]),
                  ("all-odd spread", [1,3,7,15,21,33,41,55]),
                  ("mixed {1..8}", [1,2,3,4,5,6,7,8]),
                  ("mostly-even", [2,4,6,8,10,12,14,1])]:
        g0, g1, _ = anchor_gaps(E, x)
        p0 = (g0 > TH).mean(); p1 = (g1 > TH).mean()
        print(f"  {nm:>22}: P(g0>1/7) = {p0:.4f}   P(g1/2>1/7) = {p1:.4f}   diff = {p1-p0:+.4f}")

    print("\n=== (2) independent adversarial check of boxeph's PA_2 inf table ===")
    for k in (8, 9, 13):
        gmin = (2.0, None)
        for trial in range(24):
            H = int(rng.choice([k+2, 2*k, 3*k, 5*k]))
            E = sorted(rng.choice(np.arange(1, H+1), size=k, replace=False).tolist())
            xs = grid(8009)
            cur = PA2(E, xs)
            for step in range(30):
                i = int(rng.integers(k)); new = int(rng.integers(1, int(rng.choice([2*k, 4*k, 7*k]))+1))
                if new in E: continue
                cand = sorted(set(E) - {E[i]} | {new})
                if len(cand) != k: continue
                c = PA2(cand, xs)
                if c < cur - 1e-4: E, cur = cand, c
            v = PA2(E, x)
            if v < gmin[0]: gmin = (v, tuple(E))
        note = "boxeph inf: " + {8:"0.766",9:"0.685",13:"0.360"}[k]
        print(f"  k={k:>2}: adversarial inf PA_2 = {gmin[0]:.4f} at {gmin[1]}   [{note}; T_k = {TK[k]}]"
              f"   {'DISCHARGES' if gmin[0] >= TK[k] else 'BELOW'}")

    print("\n=== (3) WEIGHTED-CHERRY DUALS at k=8 ===")
    def law_pair(d1, d2):
        g = gcd(d1, d2); q1, q2 = d1//g, d2//g
        r1, r2 = q1 % 7, q2 % 7
        return 1/49 + (min(r1, r2)*(7 - max(r1, r2)))/(49*q1*q2)
    def lp_with_dual(E):
        et = max(E); D = sorted(et - e for e in E if e != et); n = len(D)
        p = [TH]*n
        mp = {(a,b): law_pair(D[a], D[b]) for a in range(n) for b in range(a+1,n)}
        N = 1 << n
        c = np.zeros(N); c[0] = 1.0
        A = [np.ones(N)]; b = [1.0]
        for i in range(n):
            row = np.zeros(N)
            for S in range(N):
                if S >> i & 1: row[S] = 1.0
            A.append(row); b.append(p[i])
        keys = list(mp.keys())
        for (i, j) in keys:
            row = np.zeros(N)
            for S in range(N):
                if (S >> i & 1) and (S >> j & 1): row[S] = 1.0
            A.append(row); b.append(mp[(i,j)])
        res = linprog(c, A_eq=np.array(A), b_eq=np.array(b), bounds=[(0,None)]*N, method="highs")
        y = res.eqlin.marginals  # dual variables (sign convention: value = -y.b ?)
        return res.fun, y, D, keys
    for nm, E in [("AP {1..8}", [1,2,3,4,5,6,7,8]), ("barrier {1,3,5,6,7,9,25,38}", [1,3,5,6,7,9,25,38])]:
        v, y, D, keys = lp_with_dual(E)
        print(f"  {nm}: LP floor = {v:.4f}")
        print(f"    dual y0 = {-y[0]:+.4f}; singles: {[f'{-y[1+i]:+.3f}' for i in range(7)]}")
        pairw = {keys[t]: -y[8+t] for t in range(len(keys))}
        big = sorted(pairw.items(), key=lambda kv: -abs(kv[1]))[:8]
        print(f"    top pair weights (D = {D}):")
        for (a, b), w in big:
            print(f"      (d{a+1},d{b+1}) = ({D[a]},{D[b]}): weight {w:+.4f}")
