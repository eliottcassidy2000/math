#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) — the WEYL LIMIT of meas(S7) under a far element, and the rate.
kind-pasteur-2026-06-19-S12.

CLAIM (far-element contraction limit):
  Let B0 be a fixed bounded set (with 0). Consider E = B0 u {w}.
  As w -> infinity through values coprime / dissociated from B0, frac(w x)
  equidistributes on [0,1) INDEPENDENTLY of the bounded orbit O(x)={frac(b x)}.
  Then meas(S7)(E) -> LIM(B0) where the far element contributes an independent
  uniform extra sample. Specifically, conditional on the bounded orbit missing
  a set of sectors, the far point fills at most ONE of them, each with prob 1/7.

  Concretely, with M_j(B0) = measure{x: bounded orbit O(x) misses exactly the
  sector set... }, the limit is computed by an INTEGRAL over x of the
  conditional probability that the single uniform far-sample lands in the unique
  remaining missed sector (if exactly one sector is missed by B0) and 0 if >=2
  missed (one extra point can fill at most one). If B0 already hits all 1..6,
  the far point is irrelevant -> contributes meas(S7)(B0).

  LIM(B0) = meas(S7)(B0) + (1/7) * measure{x: O(x) misses EXACTLY ONE of 1..6}.

We verify this Weyl limit by:
  (a) computing meas(S7)(B0) and the "exactly one missed" measure directly,
  (b) computing meas(S7)(B0 u {w}) for large prime-ish w and comparing.
"""
import sys
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def orbit_missed_sectors(E, x):
    hit = set()
    for e in E:
        v = e*x; v = v - (v.numerator//v.denominator)
        hit.add((v.numerator*7)//v.denominator)
    return [j for j in range(1,7) if j not in hit]

def dist_missed_count(E):
    """measure of {x: exactly t of sectors 1..6 missed} for t=0..6 (Fraction)."""
    E = sorted(set(E))
    bps = set([Fraction(0), Fraction(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7*e + 1):
            bps.add(Fraction(a, 7*e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [Fraction(0)]*7
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo+hi)/2
        t = len(orbit_missed_sectors(E, mid))
        p[t] += (hi-lo)
    return p

def meas_S7(E):
    return dist_missed_count(E)[0]

def primitive(E):
    return reduce(gcd, [e for e in E if e != 0]) == 1

if __name__ == "__main__":
    print("=== WEYL LIMIT of meas(S7) under one far element ===\n")
    print("Predicted: LIM(B0) = meas(S7)(B0) + (1/7)*P(B0 misses EXACTLY 1 sector)\n")
    bases = [
        [0,1,2,3,4,5,6],       # consec k-1=7  (so full k=8)
        [0,1,2,3,4,5],         # consec 6
        [0,1,2,3,4],           # consec 5
        [0,2,3,5,6],           # non-AP base
        [0,1,3,4,6,7],         # non-AP base
    ]
    for B0 in bases:
        p = dist_missed_count(B0)
        mS7 = float(p[0]); p1 = float(p[1])
        lim = mS7 + p1/7.0
        print(f"B0={B0}")
        print(f"  meas(S7)(B0)={mS7:.5f}  P(exactly1 missed)={p1:.5f}  =>  LIM={lim:.5f}")
        # test with several large dissociated w
        for w in [101, 211, 307, 401, 503, 601, 701, 809, 907, 1009]:
            E = B0 + [w]
            if not primitive(E):
                continue
            m = float(meas_S7(E))
            print(f"    w={w:4d}: meas(S7)={m:.5f}   (LIM-pred {lim:.5f}, diff {m-lim:+.5f})")
        print()
