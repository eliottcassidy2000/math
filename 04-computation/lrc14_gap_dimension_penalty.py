#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) GAP DIMENSION PENALTY — prove that every proper d>=2 generalized AP (GAP)
E with |E|=k (k=8,9,10) has L_y(E) < cap_k.

A 2-dim GAP: E = {a*i + b*j : 0<=i<L1, 0<=j<L2}, L1*L2 = k (proper means L1,L2>=2),
   embedded with 0 in E (translate so min=0). Steps a,b chosen so the GAP is proper
   (all k elements distinct) -> Freiman dimension exactly 2.
A 3-dim GAP: E = {a*i+b*j+c*l : 0<=i<L1,...}, L1*L2*L3 = k.

We sweep many step vectors and dimensions, compute L_y exactly (Fraction), and report
the MAXIMUM L_y over all d>=2 GAPs (the WORST GAP). Compare to cap_k.

caps: cap_8=0.38153, cap_9=0.49426, cap_10=0.6044.

Builds on lrc14_empty_sector_distribution_kps.py (dist_p, L_y, g_poly).
kind-pasteur-2026-06-19  GAP dimension penalty.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd

sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

# ---- engine (copied from lrc14_empty_sector_distribution_kps.py) ----
def N_at(E, x):
    hit = set()
    for e in E:
        v = e * x
        v = v - (v.numerator // v.denominator)
        s = (v.numerator * 7) // v.denominator
        hit.add(s)
    return sum(1 for j in range(1, 7) if j not in hit)

def dist_p(E):
    E = sorted(set(E))
    bps = set([Fraction(0), Fraction(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7 * e + 1):
            bps.add(Fraction(a, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [Fraction(0)] * 7
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo + hi) / 2
        t = N_at(E, mid)
        p[t] += (hi - lo)
    return p

def g_poly(k):
    g = []
    for t in range(7):
        if k == 8:
            val = Fraction((t-1)*(t-2)*(t-4)*(t-5), 40)
        elif k in (9, 10):
            val = Fraction(-(t-2)*(t-3)*(t-6), 36)
        else:
            val = Fraction((t-3)*(t-4), 12)
        g.append(val)
    return g

def L_y(E, k):
    p = dist_p(E); g = g_poly(k)
    return sum(p[t] * g[t] for t in range(7)), p

CAP = {8: Fraction(38153,100000), 9: Fraction(49426,100000), 10: Fraction(6044,10000)}
# tighter cap reference values for reporting
CAPF = {8: 0.38153, 9: 0.49426, 10: 0.6044}

def gap2(L1, L2, a, b):
    """2-dim GAP with side lengths L1,L2 and steps a,b; translated so min=0."""
    elts = sorted(set(a*i + b*j for i in range(L1) for j in range(L2)))
    if len(elts) != L1*L2:
        return None  # not proper (elements collide -> lower dim)
    m = elts[0]
    return [e - m for e in elts]

def gap3(L1, L2, L3, a, b, c):
    elts = sorted(set(a*i + b*j + c*l for i in range(L1) for j in range(L2) for l in range(L3)))
    if len(elts) != L1*L2*L3:
        return None
    m = elts[0]
    return [e - m for e in elts]

def factor_pairs(k):
    return [(d, k//d) for d in range(2, k) if k % d == 0 and d <= k//d]

def factor_triples(k):
    out = []
    for d1 in range(2, k+1):
        if k % d1: continue
        for d2 in range(d1, k+1):
            if (k//d1) % d2: continue
            d3 = k // (d1*d2)
            if d3 >= d2:
                out.append((d1, d2, d3))
    return out

def primitive(E):
    g = reduce(gcd, [e for e in E if e != 0]) if any(e!=0 for e in E) else 1
    return g == 1

def is_full_AP(E):
    """E (sorted, min 0) is a full arithmetic progression?"""
    E = sorted(E)
    if len(E) <= 2:
        return True
    d = E[1] - E[0]
    return all(E[i+1]-E[i] == d for i in range(len(E)-1))

def excess(E):
    """|E+E| - (2k-1).  excess==0  <=>  E is a full AP (Freiman-Vosper)."""
    S = set(a+b for a in E for b in E)
    return len(S) - (2*len(E) - 1)

if __name__ == "__main__":
    print("=== LRC(14) GAP DIMENSION PENALTY: max L_y over proper d>=2 GAPs ===\n")
    STEP_MAX = 20  # range of step values to sweep
    for k in [8, 9, 10]:
        cap = CAPF[k]
        Lc, _ = L_y(list(range(k)), k)
        print(f"--- k={k}  cap={cap}  L_y(consec)={float(Lc):.5f} ---")
        best2 = (Fraction(-1), None)
        seen = set()
        # 2-dim GAPs
        for (L1, L2) in factor_pairs(k):
            for a in range(1, STEP_MAX+1):
                for b in range(1, STEP_MAX+1):
                    E = gap2(L1, L2, a, b)
                    if E is None: continue
                    if not primitive(E): continue
                    if is_full_AP(E): continue   # excess 0 -> dim 1, handled by AP case
                    key = tuple(E)
                    if key in seen: continue
                    seen.add(key)
                    Lv, _ = L_y(E, k)
                    if Lv > best2[0]:
                        best2 = (Lv, (L1, L2, a, b, E))
        Lv, info = best2
        L1,L2,a,b,E = info
        print(f"  WORST 2-dim GAP: sides({L1}x{L2}) steps({a},{b}) E={E}")
        print(f"    excess={excess(E)} (dim>=2)")
        print(f"    L_y={Lv} = {float(Lv):.5f}   margin to cap = {cap-float(Lv):.5f}   ratio L_y/consec = {float(Lv/Lc):.4f}")

        # 3-dim GAPs
        best3 = (Fraction(-1), None)
        seen3 = set()
        for (L1,L2,L3) in factor_triples(k):
            for a in range(1, STEP_MAX+1):
                for b in range(1, STEP_MAX+1):
                    for c in range(1, STEP_MAX+1):
                        E = gap3(L1,L2,L3,a,b,c)
                        if E is None: continue
                        if not primitive(E): continue
                        if is_full_AP(E): continue
                        key = tuple(E)
                        if key in seen3: continue
                        seen3.add(key)
                        Lv3, _ = L_y(E, k)
                        if Lv3 > best3[0]:
                            best3 = (Lv3, (L1,L2,L3,a,b,c,E))
        if best3[1] is not None:
            Lv3, info3 = best3
            L1,L2,L3,a,b,c,E = info3
            print(f"  WORST 3-dim GAP: sides({L1}x{L2}x{L3}) steps({a},{b},{c}) E={E}")
            print(f"    L_y={float(Lv3):.5f}   margin to cap = {cap-float(Lv3):.5f}   ratio L_y/consec = {float(Lv3/Lc):.4f}")
        print()
