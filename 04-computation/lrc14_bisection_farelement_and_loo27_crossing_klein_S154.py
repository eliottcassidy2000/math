#!/usr/bin/env python3
"""
klein-2026-07-07-S154 (part 3) -- two decisive quantities for the bisection program.

(A) EXACT 2/7 diameter crossing n2* via the Farey roof (opus-S134, verified 1e-15):
    maxgap(AP_n, x) on the Farey(n) cell [p/q, p'/q'] is LINEAR from 1/q to 1/q'.
    mu_theta(AP_n) = sum of per-cell superlevel lengths (exact rationals).
    n2* := largest n with mu_{2/7}(AP_n) >= m_P = 14249/252252.
    COMPOSITE COROLLARY (proved modulo the roof, this session's C2 + kps-S59's subset lemma):
      every 13-set E with SOME j such that E \ {e_j} has primitive diameter <= n2* - 1
      satisfies mu_{1/7}(E) >= mu_{2/7}(E\{e_j}) >= mu_{2/7}(AP_{D+1}) >= m_P.
    This covers "12 tame elements + ONE arbitrarily far element" -- a class OUTSIDE
    kps-S59's own D<=75 floor (the far element blows up diam(E)).
    Sanity anchors: mu_{2/7}(AP_12)=1313/6468, mu_{2/7}(AP_13)=829/4620 (opus-S134).

(B) FAR-ELEMENT LAW (target-2 evidence): for fixed E_j and e_j -> infinity,
    Bis_j = Ind_j + Delta_j with |Delta_j| <= C/e_j expected (crossing-count/Koksma:
    window endpoints have slopes <= ~2*max(E_j), P_j has slope e_j; on each linearity
    piece the sweep is uniform up to edge effects O(1/e_j) x #pieces).
    Measure Delta_j vs e_j for two bases; fit the decay.
"""
from fractions import Fraction as F
import numpy as np

M_P = F(14249, 252252)

def farey_pairs(n):
    """consecutive pairs (p/q, p'/q') in Farey(n) via the standard next-term walk."""
    a, b, c, d = 0, 1, 1, n
    while c <= n:
        yield (a, b, c, d)
        k = (n + b) // d
        a, b, c, d = c, d, k*c - a, k*d - b

def mu_theta_AP_exact(n, theta: F):
    """mu_theta(AP_n) exactly via the roof: per-cell linear 1/q -> 1/q'."""
    tot = F(0)
    for a, b, c, d in farey_pairs(n):
        # cell [a/b, c/d], width 1/(bd); roof endpoints 1/b -> 1/d, linear
        w = F(1, b*d)
        lo, hi = F(1, b), F(1, d)
        y0, y1 = lo, hi
        if y0 > y1: y0, y1 = y1, y0  # y0=min end, y1=max end
        if y0 >= theta:
            tot += w
        elif y1 <= theta:
            pass
        else:
            # linear from one end to other: superlevel fraction = (y1-theta)/(y1-y0)
            tot += w * (y1 - theta) / (y1 - y0)
    return tot

def gapmat(E, x):
    P = np.sort(np.outer(x, np.asarray(E, float)) % 1.0, axis=1)
    G = np.diff(P, axis=1, append=(P[:, :1] + 1.0))
    return P, G

def bis_ind(Ej, ej, theta, x):
    Pj, Gj = gapmat(Ej, x)
    mgj = Gj.max(axis=1)
    big = Gj > theta
    sel = (big.sum(axis=1) == 1) & (mgj <= 2*theta)
    idx = np.argmax(big, axis=1); rows = np.arange(len(x))
    gstart = Pj[rows, idx]; glen = Gj[rows, idx]
    win_lo = (gstart + glen - theta) % 1.0; win_len = 2*theta - glen
    p = (float(ej)*x) % 1.0
    inwin = ((p - win_lo) % 1.0) < win_len
    return float(np.mean(sel & inwin)), float(np.mean(np.where(sel, win_len, 0.0)))

if __name__ == "__main__":
    print("=== (A) exact mu_{2/7}(AP_n) via the Farey roof; crossing vs m_P = 14249/252252 ~ 0.05649 ===")
    th = F(2, 7)
    n2star = None
    anchors = {12: F(1313, 6468), 13: F(829, 4620)}
    for n in range(10, 42):
        mu = mu_theta_AP_exact(n, th)
        tag = ""
        if n in anchors:
            tag = "  ANCHOR " + ("MATCH" if mu == anchors[n] else f"MISMATCH (expect {anchors[n]})")
        flag = ">= m_P" if mu >= M_P else "< m_P"
        if mu >= M_P: n2star = n
        print(f"  n={n:>2}: mu_2/7(AP_n) = {str(mu):>22} = {float(mu):.6f}  {flag}{tag}")
    print(f"\n  n2* (largest n with mu_2/7(AP_n) >= m_P) = {n2star}")
    print(f"  => COMPOSITE (C2 + subset lemma, proved modulo roof): every 13-set with SOME")
    print(f"     leave-one-out of primitive diameter <= {n2star-1} has mu_1/7 >= m_P --")
    print(f"     covers '12 tame + one arbitrarily far element' (outside kps-S59's D<=75 zone).")

    print("\n=== (B) far-element decay: Delta_j = Bis_j - Ind_j vs e_j (theta=1/7) ===")
    x = (np.arange(120011) + 0.5)/120011
    for base_name, Ej in [("AP7 {1..7}", [1,2,3,4,5,6,7]),
                          ("AP12 {1..12}", list(range(1,13))),
                          ("mixed {1,2,4,7,9,12,15}", [1,2,4,7,9,12,15])]:
        M = max(Ej)
        print(f"  base E_j = {base_name} (M={M}):")
        for ej in [2*M+1, 5*M+1, 20*M+1, 100*M+1, 500*M+1]:
            b, ind = bis_ind(Ej, ej, 1/7, x)
            d = b - ind
            print(f"    e_j = {ej:>6}: Bis={b:.5f}  Ind={ind:.5f}  Delta={d:+.5f}   Delta*e_j/M = {d*ej/M:+.2f}")
    print("  (law check: Delta*e_j/M roughly constant => |Delta| <= C*M/e_j; C the crossing constant)")
