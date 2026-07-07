#!/usr/bin/env python3
"""
klein-2026-07-07-S154 (part 2) -- the BISECTION CHAIN quantities.

(a) EXACT AP bisection chain at theta=1/7 (pure arithmetic from the proved AP constants):
    Bis_end(AP_k) = mu(AP_{k-1}) - mu(AP_k); chain telescopes to 1 - mu(AP_13) = 601/1078.
(b) Three-distance insertion check: for the AP orbit, the new point k*x lands in a
    MAXIMAL gap of {1..k-1}x (classical three-distance insertion dynamics) -- verify the
    a.e. rate; also the >theta version (lands in the unique big gap when one exists).
(c) worst-case Delta_j over (E,j) with R3_j = 0 (is independence tight when no exact
    midpoint relation passes through j? -- the provability hook for the k=8 leg).
(d) per-k adversarial beta_k = max_E min_j Bis_j(1/7, E) for k=8..13 (jump descent),
    telescoped floors 1 - sum beta vs the per-k needs thr_k + m_P (THM-530 ledger).
"""
import numpy as np
from fractions import Fraction as F
from math import gcd
from functools import reduce

M_P = 14249/252252
THR = {8: 3637/5880, 9: None, 10: None, 11: None, 12: None, 13: 0.0}
# thr_k = 1 - min_P meas(G_P) at |P| = 13-k; minima from THM-530-A (psz=1..10 list)
GP_MIN = {1: F(6,7), 2: F(66,91), 3: F(55,91), 4: F(1979,4004), 5: F(2243,5880),
          6: F(3029,10780), 7: F(45107,229320), 8: F(2479,17640), 9: F(10601,114660),
          10: F(14249,252252)}
def thr_k(k):
    psz = 13 - k
    return 0.0 if psz == 0 else float(1 - GP_MIN[psz])

MU_AP = {7: F(1), 8: F(691,735), 9: F(247,294), 10: F(38,49), 11: F(1381,2205),
         12: F(13823,24255), 13: F(477,1078)}

def grid(NG): return (np.arange(NG) + 0.5)/NG

def gapmat(E, x):
    P = np.sort(np.outer(x, np.asarray(E, float)) % 1.0, axis=1)
    G = np.diff(P, axis=1, append=(P[:, :1] + 1.0))
    return P, G

def bis_j(E, j_elem, theta, x):
    Ej = [e for e in E if e != j_elem]
    Pj, Gj = gapmat(Ej, x)
    mgj = Gj.max(axis=1)
    big = Gj > theta
    nbig = big.sum(axis=1)
    sel = (nbig == 1) & (mgj <= 2*theta)
    idx = np.argmax(big, axis=1); rows = np.arange(len(x))
    gstart = Pj[rows, idx]; glen = Gj[rows, idx]
    win_lo = (gstart + glen - theta) % 1.0; win_len = 2*theta - glen
    pj = (float(j_elem)*x) % 1.0
    inwin = ((pj - win_lo) % 1.0) < win_len
    return float(np.mean(sel & inwin)), float(np.mean(np.where(sel, win_len, 0.0)))

def R3_through(E, j_elem):
    S = set(E); c = 0
    for a in E:
        if a == j_elem: continue
        b = 2*j_elem - a
        if b in S and b != j_elem and b > a: c += 1
    return c

def primitive(E):
    g = reduce(gcd, E); return tuple(sorted(e//g for e in E))

if __name__ == "__main__":
    rng = np.random.default_rng(4154)
    TH = 1/7
    print("=== (a) EXACT AP bisection chain, theta=1/7 (arithmetic from proved constants) ===")
    tot = F(0)
    for k in range(8, 14):
        b = MU_AP[k-1] - MU_AP[k]; tot += b
        print(f"  Bis_end(AP_{k}) = mu(AP_{k-1}) - mu(AP_{k}) = {b} = {float(b):.5f}")
    print(f"  chain total = {tot} = {float(tot):.5f}  (= 1 - mu(AP_13) = {F(1)-MU_AP[13]} ✓)"
          f"  {'OK' if tot == 1 - MU_AP[13] else 'MISMATCH'}")

    print("\n=== (b) three-distance insertion: does k*x land in a MAXIMAL gap of {1..k-1}x? ===")
    x = grid(20001)
    for k in (8, 10, 13):
        E1 = list(range(1, k))
        P, G = gapmat(E1, x)
        mg = G.max(axis=1)
        pj = (k*x) % 1.0
        # gap containing pj: vectorized left-neighbor count per row
        cnt = (P <= pj[:, None]).sum(axis=1)
        left = np.where(cnt > 0, cnt-1, k-2)  # gap index (wrap: if cnt=0, in wrap gap k-2)
        rows = np.arange(len(x))
        gap_of_pj = G[rows, left]
        in_max = np.isclose(gap_of_pj, mg, rtol=0, atol=1e-12)
        big = mg > TH
        onebig = (G > TH).sum(axis=1) == 1
        print(f"  k={k}: P(k*x in maximal gap) = {in_max.mean():.4f}   "
          f"P(in max gap | unique >1/7 gap) = {(in_max & onebig).sum()/max(1,onebig.sum()):.4f}")

    print("\n=== (c) worst-case Delta_j at R3_j = 0 (random + structured probes) ===")
    worst = (0.0, None, None)
    cases = 0
    for _ in range(150):
        E = sorted(rng.choice(np.arange(1, 46), size=13, replace=False).tolist())
        js = [j for j in E if R3_through(E, j) == 0]
        if not js: continue
        j = int(rng.choice(js)); cases += 1
        b, ind = bis_j(E, j, TH, x[:10001])
        d = b - ind
        if abs(d) > abs(worst[0]): worst = (d, tuple(E), j)
    # structured probes: near-AP with one outlier (approximate midpoints!)
    for out in (14, 15, 27, 40, 53):
        E = list(range(1, 13)) + [out]
        for j in E:
            if R3_through(E, j) != 0: continue
            b, ind = bis_j(E, j, TH, x[:10001])
            d = b - ind; cases += 1
            if abs(d) > abs(worst[0]): worst = (d, tuple(E), j)
    print(f"  cases with R3_j=0 probed: {cases}")
    print(f"  WORST |Delta_j| at R3=0: {worst[0]:+.4f}  at E={worst[1]} j={worst[2]}")
    print("  (if this stays <<0.238, the k=8 window-avoidance has provable-looking slack;"
          " approximate midpoints can still contribute)")

    print("\n=== (d) per-k adversarial beta_k = max_E min_j Bis_j (theta=1/7, jump descent) ===")
    xs = grid(4003); xa = grid(12007)
    def minbis(E, xg):
        return min(bis_j(E, j, TH, xg)[0] for j in E)
    needs, betas = {}, {}
    for k in range(8, 14):
        gworst = (0.0, None)
        for trial in range(16):
            H = int(rng.choice([k+2, k+6, 2*k, 3*k]))
            E = sorted(rng.choice(np.arange(1, H+1), size=k, replace=False).tolist())
            cur = minbis(E, xs)
            for step in range(40):
                i = int(rng.integers(k)); new = int(rng.integers(1, int(rng.choice([k+4, 2*k+4, 3*k+4]))+1))
                if new in E: continue
                cand = sorted(set(E) - {E[i]} | {new})
                if len(cand) != k: continue
                c = minbis(cand, xs)
                if c > cur + 1e-3: E, cur = cand, c
            v = minbis(E, xa)
            if v > gworst[0]: gworst = (v, primitive(E))
        # also seed with AP (the conjectured extremal)
        vAP = minbis(list(range(1, k+1)), xa)
        if vAP > gworst[0]: gworst = (vAP, tuple(range(1, k+1)))
        betas[k] = gworst
        needs[k] = thr_k(k) + M_P
        print(f"  k={k:>2}: beta_k(adv) = {gworst[0]:.4f} at {gworst[1]}   [AP value {vAP:.4f}]")
    print("\n  telescoped recursion floors vs per-k needs (need_k = thr_k + m_P):")
    run = 0.0
    for k in range(8, 14):
        run += betas[k][0]
        floor = 1.0 - run
        need = needs[k]
        print(f"  k={k:>2}: floor 1-sum(beta) = {floor:.4f}   need >= {need:.4f}   "
              f"margin {floor-need:+.4f}  {'HOLDS' if floor >= need else 'FAILS'}")
    print("\n  NOTE: beta_k here is the adversarial-empirical sup of min_j Bis_j; the telescoped")
    print("  floor is a CONJECTURAL uniform bound m_k >= 1 - sum_{i<=k} beta_i pending sup proofs.")
