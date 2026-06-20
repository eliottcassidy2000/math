#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THE ONE OPEN INPUT of the LRC(14) sector-route proof (HYP-2653): the uniform
decorrelation constant  C(k) = sup_{E', w}  w|Delta_w|  (sigma-INDEPENDENT).
kind-pasteur-2026-06-19-S14. EXACT rationals throughout (number theory, not floats).

For E = E' u {w}, w=max E, the far-element plateau error is EXACTLY
   Delta_w = sum_{N=1 cells c} integral_{a_c}^{b_c} [1_{sigma_{s_c}}(w x) - 1/7] dx
           = (1/w) sum_c [ G0(w b_c - s_c/7) - G0(w a_c - s_c/7) ],
where cell_c=[a_c,b_c] is a maximal interval on which the E'-orbit misses EXACTLY
inner sector s_c, and G0(y)=integral_0^{frac y}(1_{[0,1/7)}-1/7) is the tent:
   G0(f) = (6/7) f            for 0<=f<=1/7
   G0(f) = 6/49 - (f-1/7)/7   for 1/7<=f<=1.        |G0| <= 6/49.
GOAL: the EXACT value of sup w|Delta_w| and its worst (E',w) — the structured
breakpoint discrepancy that must be O(1) surviving resonances (w sharing factors
with the e's / mod 7).  Need C <= ~1.95 (dovetail HYP-2653b) for the proof to close.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def G0(y):
    """exact tent G0 at real y (Fraction); 1-periodic."""
    f = y - (y.numerator // y.denominator)        # frac in [0,1)
    if f <= Fraction(1,7):
        return Fraction(6,7) * f
    return Fraction(6,49) - (f - Fraction(1,7)) / 7

def missed_sector_cells(Ep):
    """Return list of (a,b,s): maximal intervals [a,b] where E'-orbit {frac(e x)}
       misses EXACTLY inner sector s in {1..6}.  Exact Fraction endpoints."""
    Ep = sorted(set(Ep))
    bps = {Fraction(0), Fraction(1)}
    for e in Ep:
        if e == 0: continue
        for a in range(0, 7*e + 1):
            bps.add(Fraction(a, 7*e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    cells = []
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo+hi)/2
        hit = set()
        for e in Ep:
            v = e*mid; v = v - (v.numerator//v.denominator)
            hit.add((v.numerator*7)//v.denominator)
        missed = [j for j in range(1,7) if j not in hit]
        if len(missed) == 1:
            cells.append((lo, hi, missed[0]))
    # merge adjacent cells with same missed sector (so runs telescope correctly)
    merged = []
    for (a,b,s) in cells:
        if merged and merged[-1][2]==s and merged[-1][1]==a:
            merged[-1] = (merged[-1][0], b, s)
        else:
            merged.append((a,b,s))
    return merged

def wDelta(Ep, w, cells=None):
    """exact w*Delta_w = sum_c [G0(w b - s/7) - G0(w a - s/7)]."""
    if cells is None: cells = missed_sector_cells(Ep)
    tot = Fraction(0)
    for (a,b,s) in cells:
        tot += G0(w*b - Fraction(s,7)) - G0(w*a - Fraction(s,7))
    return tot

def p0(E):
    E=sorted(set(E)); bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; hit=set()
        for e in E:
            v=e*mid; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        if len(hit)==7: tot+=hi-lo
    return tot
def dist_p1(Ep):  # p_1(E') = meas{E' misses exactly 1}
    return sum(b-a for (a,b,s) in missed_sector_cells(Ep))
def primitive(E): return reduce(gcd,E)==1

if __name__ == "__main__":
    print("=== (0) VERIFY the exact reduction: w*Delta_w(engine) == w*(p0(E)-plateau) ===")
    for Ep, w in [([0,1,2,3,4,5,6,7], 30), ([0,1,2,3,4,5,6,7], 23), ([0,1,2,3,4,5,7], 19)]:
        cells = missed_sector_cells(Ep)
        plat = p0(Ep) + Fraction(1,7)*dist_p1(Ep)
        E = Ep + [w]
        lhs = w*(p0(E) - plat)
        rhs = wDelta(Ep, w, cells)
        print(f"   E'={Ep} w={w}: w*(p0-plat)={lhs} == wDelta={rhs}  [{'MATCH' if lhs==rhs else 'MISMATCH!!'}]  (#cells={len(cells)})")

    print("\n=== (1) EXACT sup of w*|Delta_w| over consec cores + range of w ===")
    for kcore in [7, 8, 9]:   # E' = consec_{kcore} = {0,..,kcore-1}; full set k=kcore+1
        Ep = list(range(kcore))
        cells = missed_sector_cells(Ep)
        best = (Fraction(0), None)
        for w in range(kcore, 400):
            if not primitive(Ep + [w]): continue
            val = abs(wDelta(Ep, w, cells))
            if val > best[0]: best = (val, w)
        print(f"   E'=consec_{kcore} (#N=1 cells={len(cells)}): max w*|Delta_w| = {best[0]} = {float(best[0]):.5f} at w={best[1]}")

    print("\n=== (2) sup over MANY bounded cores (span<=12) for k=9 (E' has 8 elts) ===")
    kcore = 8; overall = (Fraction(0), None, None)
    cnt = 0
    for tail in itertools.combinations(range(1, 13), kcore-1):
        Ep = (0,) + tail
        if not primitive(Ep): continue
        cnt += 1
        cells = missed_sector_cells(Ep)
        for w in range(13, 200):
            if not primitive(Ep + (w,)): continue
            val = abs(wDelta(Ep, w, cells))
            if val > overall[0]: overall = (val, Ep, w)
    print(f"   scanned {cnt} cores x w-range: SUP w*|Delta_w| = {overall[0]} = {float(overall[0]):.5f}")
    print(f"   worst (E',w) = {overall[1]}, w={overall[2]}")
    print(f"   => empirical C(9) <= {float(overall[0]):.4f}  (dovetail needs C < (cap_9-Q(8))*16; cap_9=1979/4004)")
