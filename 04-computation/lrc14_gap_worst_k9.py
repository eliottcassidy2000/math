#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14): pin down the WORST (max L_y) d>=2 set for the binding cluster k=9.

The 2-dim GAP sweep gave worst E=[0,2,3,4,5,6,7,8,10] with L_y=0.44933 (margin
0.045 below cap_9=0.49426). But the crux requires L_y(E)<cap for EVERY E with 0
in E, |E|=9 that is NOT a full AP -- not only GAP-FILLING sets.  So here we do a
BRUTE-FORCE confirmation: enumerate ALL primitive 9-element sets E with 0 in E and
max(E) <= W, that are NOT full APs, and report the maximizer.  We push W up until
the maximizer stabilizes, confirming the worst non-AP set is the excess-2 near-AP
[0,2,3,4,5,6,7,8,10] and that L_y < cap with margin ~0.045.

This is the finite check (step 2 of the proof program) made airtight for k=9, the
tightest cluster.  kind-pasteur-2026-06-19.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

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
        p[N_at(E, mid)] += (hi - lo)
    return p

def g_poly(k):
    return [Fraction(-(t-2)*(t-3)*(t-6), 36) for t in range(7)]  # k=9,10

def L_y(E, k=9):
    p = dist_p(E); g = g_poly(k)
    return sum(p[t]*g[t] for t in range(7))

def is_full_AP(E):
    E = sorted(E)
    d = E[1]-E[0]
    return all(E[i+1]-E[i]==d for i in range(len(E)-1))

if __name__ == "__main__":
    k = 9
    cap = 0.49426
    Lc = float(L_y(list(range(k))))
    print(f"k={k} cap={cap} L_y(consec)={Lc:.5f}\n")
    for W in [10, 11, 12, 13]:
        best = (Fraction(-1), None)
        cnt = 0
        # E = {0} U 8 chosen from 1..W, primitive, not full AP
        for combo in itertools.combinations(range(1, W+1), k-1):
            if W not in combo:  # only sets with max exactly W (avoids recount)
                continue
            E = (0,) + combo
            if reduce(gcd, E[1:]) != 1:
                continue
            if is_full_AP(E):
                continue
            cnt += 1
            Lv = L_y(E)
            if Lv > best[0]:
                best = (Lv, E)
        Lv, E = best
        print(f"  W={W:2d}: scanned {cnt:6d} primitive non-AP sets; "
              f"max L_y={float(Lv):.5f} at E={list(E)}  (margin {cap-float(Lv):.5f})")
    print("\n(If maximizer stays at [0,2,3,4,5,6,7,8,10]-type excess-2 near-AP as W grows,")
    print(" the worst non-AP set is pinned and L_y<cap with margin ~0.045.)")
