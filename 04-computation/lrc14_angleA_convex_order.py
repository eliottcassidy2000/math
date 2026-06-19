#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE A: coupling / convex-order on N_E for LRC(14) endgame.
Goal: PROVE consec maximizes L_y = E[g(N)] for k=8 and k=9.

We test candidate stochastic orders on the distribution p=(p_0,...,p_6) of N:
 - Is p_0=meas(S7) itself maximized by consec? (exact, over many competitors)
 - Does N_consec dominate N_E in some order (usual SD, increasing-convex, etc.)?
 - Which SINGLE linear inequality on (p_0..p_6) holds for all E and gives L_y<=L_y(consec)?

Engine copied from lrc14_empty_sector_distribution_kps.py (exact, fractions.Fraction).
kind-pasteur-2026-06-19 ANGLE-A.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd

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

def L_y(p, k):
    g = g_poly(k)
    return sum(p[t] * g[t] for t in range(7))

def consec(k):
    return list(range(k))

def gen_competitors(k, maxspread, limit=None):
    """all sets {0=e0<e1<...<e_{k-1}<=maxspread}, dissociation not required."""
    out = []
    elts = list(range(1, maxspread+1))
    for combo in itertools.combinations(elts, k-1):
        E = [0]+list(combo)
        if reduce(gcd, E) != 1:
            continue
        out.append(E)
        if limit and len(out) >= limit:
            break
    return out

if __name__ == "__main__":
    for k in [8, 9]:
        print(f"\n{'='*70}\nk={k}\n{'='*70}")
        g = g_poly(k)
        supp = [t for t in range(7) if g[t]!=0]
        print(f"g support {supp}, g={[str(x) for x in g]}")
        C = consec(k)
        pc = dist_p(C)
        Lc = L_y(pc, k)
        print(f"consec={C}")
        print(f"  p_consec = {[str(x) for x in pc]}")
        print(f"  p_consec ~ {[f'{float(x):.5f}' for x in pc]}")
        print(f"  L_y(consec) = {Lc} ~ {float(Lc):.6f}")
        print(f"  p0(consec)=meas S7 = {float(pc[0]):.6f}")

        # cumulative CDFs of consec
        cdf_c = [sum(pc[:t+1]) for t in range(7)]   # F(t)=P(N<=t)
        tail_c = [sum(pc[t:]) for t in range(7)]     # P(N>=t)

        spread = {8:14, 9:13}[k]   # keep runtime modest; consec spread is k-1
        comps = gen_competitors(k, maxspread=spread)
        print(f"  testing {len(comps)} competitors with max <= {spread}")

        # collected violations
        p0_beats = 0          # competitor p0 > consec p0
        Ly_beats = 0          # competitor L_y > consec L_y (should be 0)
        usd_holds = 0         # usual stochastic dominance N_E <= N_consec everywhere
        icx_holds = 0         # increasing-convex order
        # the "tail at support" inequality candidate
        worst_p0_gap = Fraction(0)
        worst_Ly_gap = Fraction(0)

        for E in comps:
            if E == C: continue
            p = dist_p(E)
            L = L_y(p, k)
            if p[0] > pc[0]: p0_beats += 1
            if L > Lc:
                Ly_beats += 1
                print(f"    !! L_y BEATS consec: {E} L={float(L):.6f}")
            gap0 = pc[0]-p[0]
            if gap0 < worst_p0_gap: worst_p0_gap = gap0
            gapL = Lc-L
            if gapL < worst_Ly_gap: worst_Ly_gap = gapL
            # usual SD: N_E <=_st N_consec  <=> P(N_E>=t) <= P(N_consec>=t) all t
            tail_e = [sum(p[t:]) for t in range(7)]
            if all(tail_e[t] <= tail_c[t] for t in range(7)):
                usd_holds += 1
            # increasing-convex: sum_{s>=t} P(N_E>=s) <= same for consec, all t
            icx_e = [sum(tail_e[s] for s in range(t,7)) for t in range(7)]
            icx_c = [sum(tail_c[s] for s in range(t,7)) for t in range(7)]
            if all(icx_e[t] <= icx_c[t] for t in range(7)):
                icx_holds += 1

        ncomp = len([E for E in comps if E!=C])
        print(f"  --- results over {ncomp} competitors ---")
        print(f"  consec p0 is the LARGEST: {'YES' if p0_beats==0 else f'NO ({p0_beats} beat it)'}")
        print(f"  consec L_y is the LARGEST: {'YES' if Ly_beats==0 else f'NO ({Ly_beats} beat it)'}")
        print(f"  worst p0 gap (consec-comp, want >=0): {float(worst_p0_gap):.6f}")
        print(f"  worst L_y gap (consec-comp, want >=0): {float(worst_Ly_gap):.6f}")
        print(f"  usual stoch dom N_E<=N_consec holds for {usd_holds}/{ncomp}")
        print(f"  increasing-convex order holds for {icx_holds}/{ncomp}")
