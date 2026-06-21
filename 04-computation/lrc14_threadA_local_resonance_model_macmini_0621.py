#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A (mac-mini 2026-06-21) -- EXACT LOCAL MODEL of one resonance & the
per-resonance trade-off (why consec wins the SUM but loses a=1,5 individually).

From lrc14_threadA_resonant_perA: over the full-residue stratum (k=8), consec
is the UNIQUE aggregate max of measS7 = sum_{a=1..6} W_a(E), but it is NOT the
max of each individual W_a: it LOSES a=1 (to [0,2,3,4,5,6,7,8]) and a=5 (to
[0,1,2,4,5,6,7,10]), and WINS a=2,3,4,6. So consec-max is genuinely AGGREGATE
even inside the resonant decomposition. This script pins the exact local model
of W_a and the trade-off, to extract the precise obstruction.

LOCAL MODEL of W_a(E)  (the cover length in the cell [a/7-1/14, a/7+1/14]).
Put x = a/7 + y, y in (-1/14, 1/14). For e in E:
   7 frac(e x) = 7 frac(e a /7 + e y) = (e a mod 7) + frac-correction from 7 e y.
The sector of e at x is:
   s_e(y) = (e a + floor(7 e y - {boundary}) ) mod 7   -- piecewise constant in y,
   jumping when 7 e y crosses an integer offset from the residue boundary.
At y=0 the covered set is { e a mod 7 : e in E } = a*R where R = {e mod 7} (=Z/7).
As |y| grows each e's sector drifts at LOCAL SPEED 7e. The cell stays fully
covered until some residue class loses its LAST representative.

So W_a = length of y-interval (within the cell) on which every residue r in Z/7
is HIT by at least one e. Each e hits residue (e a + j) mod 7 on a y-subinterval;
the residue r is covered iff at least one e currently sits in r. This is a
1-D interval-cover problem; W_a is exact-rational.

KEY: the local SLOPE of e is 7e (how fast its sector drifts). Residue r=ea mod 7
is realized at y=0 by all e with that residue; the FASTEST-drifting (largest |7e|
i.e. largest |e|) leaves its residue first, the SLOWEST (smallest |e|) stays
longest. So the residue's coverage interval is governed by the SLOWEST e in that
class on each side -- i.e. by min |e| per residue class. consec has the
SMALLEST possible representatives 0..7 (one per residue, residue 0 doubled by
0 and 7), which keeps every residue covered longest -> per-resonance this should
be near-optimal but the DOUBLING of residue 0 (by 0 and 7) and the asymmetry
of which side is governed by which e is what breaks the per-resonance optimality
for a=1,5.

This script:
 (1) compute, for each shape & each a, the per-residue "left/right survival" y
     and verify W_a equals the interval where all residues survive.
 (2) show the trade-off: which residue's survival is the binding constraint per a.
 (3) test the AGGREGATE-Markoff hypothesis: sum_a W_a as a function of the
     residue->magnitude assignment; consec uses magnitudes {0,1,2,3,4,5,6,7}.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def covered_sectors(E, xm):
    secs = set()
    for e in E:
        v = e * xm; v = v - (v.numerator // v.denominator)
        secs.add((v.numerator * 7) // v.denominator)
    return secs

def measS7_arcs(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    arcs = []
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        if len(covered_sectors(E, (lo+hi)/2)) == 7: arcs.append((lo, hi))
    return arcs

def W_a(E, a):
    arcs = measS7_arcs(E)
    lo_b, hi_b = F(2*a-1, 14), F(2*a+1, 14)
    w = F(0)
    for lo, hi in arcs:
        l = max(lo, lo_b); h = min(hi, hi_b)
        if h > l: w += h - l
    return w

def binding_residue(E, a):
    """At the LEFT and RIGHT edge of the a-cover, which residue is about to die?
       Determined by examining covered sectors just past the cover boundary."""
    arcs = measS7_arcs(E)
    lo_b, hi_b = F(2*a-1, 14), F(2*a+1, 14)
    sub = [(max(lo,lo_b), min(hi,hi_b)) for lo,hi in arcs if min(hi,hi_b)>max(lo,lo_b)]
    if not sub: return None
    L = min(s[0] for s in sub); R = max(s[1] for s in sub)
    # just outside: which residue dropped?
    eps = F(1, 7*max(abs(e) for e in E)*100)
    inside = covered_sectors(E, (L+R)/2)
    out_left = covered_sectors(E, L - eps) if L-eps > F(0) else set()
    out_right = covered_sectors(E, R + eps) if R+eps < F(1) else set()
    missL = inside - out_left; missR = inside - out_right
    return (L, R, sorted(missL), sorted(missR))

def residues(E): return frozenset(e % 7 for e in E)
def is_full_residue(E): return residues(E) == frozenset(range(7))
def primitive(E): return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

if __name__ == "__main__":
    print("#"*78); print("# LOCAL RESONANCE MODEL & PER-RESONANCE TRADE-OFF (THREAD A)"); print("#"*78)
    k = 8
    C = consec(k)
    competitors = {
        "consec":   C,
        "a1-beater": [0,2,3,4,5,6,7,8],
        "a5-beater": [0,1,2,4,5,6,7,10],
    }
    print("\nPer-resonance widths W_a (a=1..6) and which residue binds each side:")
    for name, E in competitors.items():
        print(f"\n  {name} = {E}  (residues {sorted(residues(E))}):")
        tot = F(0)
        for a in range(1, 7):
            w = W_a(E, a); tot += w
            br = binding_residue(E, a)
            print(f"    a={a}: W={float(w):.6f}={w}  bind={br[2:] if br else None}")
        print(f"    SUM = {float(tot):.6f} = {tot}")

    # Aggregate-Markoff: each residue class r in Z/7 is realized by e's; the
    # per-resonance survival for residue r at resonance a depends on the e with
    # residue (r a^{-1} mod 7). The "leg" lengths are min|e| per residue. Print
    # consec's per-residue MIN magnitude and the two beaters'.
    print("\nPer-residue MIN |e| (the slowest representative = the survival leg):")
    for name, E in competitors.items():
        byres = defaultdict(list)
        for e in E: byres[e%7].append(abs(e))
        legs = {r: min(byres[r]) for r in range(7)}
        print(f"  {name}: { {r:legs[r] for r in range(7)} }")
