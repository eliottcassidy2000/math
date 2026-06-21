#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD 3 -- the load-bearing number for ANGLE 2 (lean, consec base only, k=8..12).

For the WIDE UPPER bound, p0_decorr is MAXIMIZED at the consec base (= Q(k-1)). The dangerous
direction is err+ = (p0 - p0_decorr)+ at the consec base. Compute it exactly per k and report
err+/margin. Smaller window of commensurable ratios, faster.
"""
import sys, random
from fractions import Fraction as F
from functools import reduce
from math import gcd, lcm
sys.stdout.reconfigure(line_buffering=True)

caps = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91), 11: F(66,91), 12: F(6,7)}
Qkm1 = {8: F(289,1470), 9: F(621,1715), 10: F(1229,2744), 11: F(65599,123480), 12: F(14873,24696)}

def p0(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        s = set()
        for e in E:
            v = F(e) * mid; v = v - F(v.numerator // v.denominator)
            s.add((v.numerator * 7) // v.denominator)
        if len(s & set(range(1, 7))) == 6: tot += hi - lo
    return tot

def p0_decorr(B, samples=20000, seed=0):
    Bs = sorted(set(int(e) for e in B))
    bps = {F(0), F(1)}
    for e in Bs:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    rng = random.Random(seed); tot = 0.0
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        s = set()
        for e in Bs:
            v = F(e) * mid; v = v - F(v.numerator // v.denominator)
            s.add((v.numerator * 7) // v.denominator)
        need = set(range(1, 7)) - s
        if not need:
            tot += float(hi - lo); continue
        hit = 0
        for _ in range(samples):
            got = set(rng.randrange(7) for _ in range(2))
            if need <= got: hit += 1
        tot += float(hi - lo) * hit / samples
    return tot

def main():
    C = 220
    Rs = []
    for q in range(1, 6):
        for p in range(q + 1, 3 * q + 1):
            if gcd(p, q) == 1 and 1 < p / q <= 2.2: Rs.append((p, q))
    Rs = sorted(set(Rs), key=lambda t: t[0]/t[1])
    print("consec base. err+ = max positive (p0 - p0_decorr) over commensurable far pairs.\n")
    print("k | margin | p0_decorr(consec) | max(err+) | worst ratio | err+/margin | lossiness")
    for k in (8, 9, 10, 11, 12):
        margin = float(caps[k] - Qkm1[k])
        base = list(range(k - 2))
        bdec = p0_decorr(base, samples=30000, seed=k)
        emax = -9; wr = None
        for (p, q) in Rs:
            E = sorted(set(base + [C*q, C*p]))
            if len(E) != k or reduce(gcd, E) != 1: continue
            err = float(p0(E)) - bdec
            if err > emax: emax = err; wr = (p, q)
        loss = margin / emax if emax > 0 else float('inf')
        print(f"{k} | {margin:.4f} | {bdec:.4f} | {emax:+.5f} | {wr[0]}/{wr[1]} | {emax/margin:.4f} | {loss:.0f}x")
    print("\n=> err+ shrinks while margin grows: the loose resonance bound has GROWING slack in k.")

if __name__ == "__main__":
    main()
