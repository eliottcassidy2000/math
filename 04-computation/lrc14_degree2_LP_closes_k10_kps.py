"""
lrc14_degree2_LP_closes_k10_kps.py  (kind-pasteur-2026-06-27-S31ag)

Does the EXACT degree-2 moment-LP close the cover bound at k=10 (the loosest
binding row, cap_10=55/91 = pair-Pascal, THM-576/HYP-3092)?

For a config with miss-count distribution q (N in {0..6}), the cover bound is
q_0 = meas(S7) <= cap_k. The degree-2 LP bound:
   max q_0  s.t.  sum q_t = 1, sum t q_t = S1, sum C(t,2) q_t = S2, q_t >= 0.
This is the BEST bound using only S0,S1,S2 (tighter than Paley-Zygmund).
If max_q0_LP(S1,S2) <= cap_k for ALL k-clusters, the cover bound CLOSES at
degree 2 (pure pairwise) for that k. Expect: closes k=10 (and k>=10), fails
k=8,9 (the higher-Pascal dip = mac-mini's Hankel obligation).
"""
import sys, itertools, random
from fractions import Fraction as F
from math import comb

def sector_of(p): return int((p % 1) * 7)

def q_dist(E):
    E = sorted(set(E)); b = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * e + 1): b.add(F(m, 7 * e))
    b = sorted(b); q = [F(0)] * 7
    for i in range(len(b) - 1):
        x0, x1 = b[i], b[i + 1]
        if x1 <= x0: continue
        cov = set(sector_of(e * ((x0 + x1) / 2)) for e in E)
        t = 7 - len(cov)
        if 0 <= t <= 6: q[t] += x1 - x0
    return q

def deg2_LP_max_q0(S1, S2):
    """max q0 s.t. sum q=1, sum t q=S1, sum C(t,2) q=2... =S2, q>=0, t in 0..6.
    Vertex method: optimum supported on {0} + at most 2 other points (3 constraints).
    Enumerate supports {0,a,b} (a<b in 1..6), solve 3x3, keep feasible max q0."""
    best = F(0)
    pts = list(range(7))
    for a, b in itertools.combinations(range(1, 7), 2):
        # q0 + qa + qb = 1 ; a qa + b qb = S1 ; C(a,2) qa + C(b,2) qb = S2
        # solve for qa,qb from last two, q0 = 1-qa-qb
        import numpy as np
        Amat = [[a, b], [comb(a,2), comb(b,2)]]
        det = Amat[0][0]*Amat[1][1] - Amat[0][1]*Amat[1][0]
        if det == 0: continue
        qa = (S1*Amat[1][1] - Amat[0][1]*S2) / det
        qb = (Amat[0][0]*S2 - S1*Amat[1][0]) / det
        q0 = 1 - qa - qb
        if qa >= 0 and qb >= 0 and q0 >= 0:
            if q0 > best: best = q0
    # also supports with a single extra point {0,a}: q0+qa=1, a qa=S1 => qa=S1/a, needs C(a,2)qa=S2
    for a in range(1, 7):
        qa = F(S1, a) if a else F(0)
        if comb(a,2)*qa == S2 and 0 <= qa <= 1:
            q0 = 1 - qa
            if q0 > best: best = q0
    return best

CAPS = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91), 11: F(66,91)}

if __name__ == "__main__":
    sys.stdout.reconfigure(line_buffering=True)
    random.seed(1010)
    for k in (8, 9, 10, 11):
        cap = CAPS[k]
        configs = [tuple(range(k)), tuple(2*i for i in range(k))]
        for _ in range(400):
            cfg = tuple(sorted([0] + random.sample(range(1, 26), k-1)))
            if len(set(cfg)) == k: configs.append(cfg)
        viol = 0; max_lp = F(0); max_q0 = F(0); worst = None
        for E in configs:
            q = q_dist(E)
            S1 = sum(t*q[t] for t in range(7))
            S2 = sum(comb(t,2)*q[t] for t in range(7))
            lp = deg2_LP_max_q0(S1, S2)
            if q[0] > max_q0: max_q0 = q[0]
            if lp > max_lp: max_lp = lp; worst = E
            if lp > cap + F(1,10**9): viol += 1
        status = "CLOSES at degree-2 (pairwise)" if viol == 0 else f"FAILS ({viol} viol) -> needs dip"
        print(f"k={k}: cap={float(cap):.5f}  max q0(actual)={float(max_q0):.5f}  "
              f"max deg2-LP bound={float(max_lp):.5f}  -> {status}")
