#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A (mac-mini 2026-06-21) -- ANGLE 2 FOSTER/EFFECTIVE-RESISTANCE SUMMARY.

VERDICT on the electrical-network hypothesis for the LRC(14) wall (LAYER 3):

CONFIRMED (electrical dictionary, exact):
  D1. g(t)=2t(7-t) = 14 * R_eff(0,t) on the cycle graph C_7  (R_eff = t(7-t)/7).
  D2. FOSTER CONSERVATION: the residue legs are a TRANSVERSAL of Z/7, so
        sum_{r in Z/7} g(r) = sum_{t=0}^{6} 2t(7-t) = 112  is CONSTANT for every
        full-residue shape. This is the exact analogue of Foster's theorem
        (sum of edge effective-resistances on C_7 = N-1 = 6).
  D3. sumW(E) is EXACTLY a FUNCTION of the per-residue CONDUCTANCE vector
        c_r = sum_{e in E, e!=0, e mod7=r} 1/|e|  (0/2260 collisions). Parallel
        clocks in a residue class ADD conductance (1/|e|) exactly as parallel
        resistors. THIS is the sharp invariant the leg-profile (min|e|) missed.
  D4. consec UNIQUELY MAXIMIZES the TOTAL conductance T = sum_r c_r = sum 1/|e|
        over the stratum (trivially: {1..7} is the max harmonic sum of 7 distinct
        nonzero integers).
  D5. The e=0 clock is a PERFECT SHORT on sector 0: frozen at sector 0 for ALL x
        (speed 0). It permanently rails sector 0, freeing the SECOND residue-0
        clock (e=7) to handle the drift. This is exactly why doubling residue 0
        (topology (2,1,1,1,1,1,1)) wins: 0.327 vs 0.208 for any other doubling.

FAILED (the clean monotone-reduction does NOT exist):
  F1. sumW is NOT globally monotone in T (concordance ~0.51): T-max and sumW-max
        coincide AT consec but the functions do not track. So "T-max => sumW-max"
        is NOT a valid reduction.
  F2. Rayleigh contraction (decreasing a magnitude raises sumW) FAILS even within
        a FIXED topology (39562 inversions at span<=20). The electrical-monotone
        principle is FALSE here -- adding conductance does not monotonically widen
        survival, because the survival window is a MIN over sectors (a bottleneck),
        not a sum/flow. The bottleneck (min) structure breaks Rayleigh.
  F3. consec is NOT conductance-dominant coordinatewise (1184 shapes beat it in
        some residue by adding a parallel clock).

NET ELECTRICAL PICTURE (the honest synthesis):
  sum_a W_a is a MIN-COVERAGE (bottleneck) functional of a conductance vector that
  lives on C_7 with a Foster-conserved residue budget (sum g(r)=112). consec is
  the unique simultaneous (i) total-conductance maximizer AND (ii) perfect-short
  on the identity sector. The wall is NOT a Rayleigh-monotone extremality (min
  breaks monotonicity) -- it is a BOTTLENECK extremality: consec maximizes the
  WORST-COVERED sector's survival across all resonances by (a) shortest clocks
  (max conductance) and (b) the e=0 perfect short removing sector 0 from the min.

This REFRAMES the wall: LAYER 3 is a min-of-conductances (Chebyshev/bottleneck)
extremality on C_7 under a Foster budget, NOT a Foster-additive one. The right
tool is a MINIMAX (not a convex/Schur additive) argument: the AP equalizes the
per-sector survival bottlenecks (the standard signature of an extremal minimax).
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
    E = sorted(set(int(e) for e in E)); bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1); arcs = []
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        if len(covered_sectors(E, (lo+hi)/2)) == 7: arcs.append((lo, hi))
    return arcs
def W_a(E, a):
    lo_b, hi_b = F(2*a-1, 14), F(2*a+1, 14); w = F(0)
    for lo, hi in measS7_arcs(E):
        l = max(lo, lo_b); h = min(hi, hi_b)
        if h > l: w += h - l
    return w
def sumW(E): return sum(W_a(E, a) for a in range(1, 7))
def g(t, N=7): return 2*t*(N-t)

if __name__ == "__main__":
    print("ANGLE 2 SUMMARY: electrical model of the LRC(14) wall")
    print(f"  D2 Foster conservation sum_t g(t) = {sum(g(t) for t in range(7))} (=112)")
    print(f"  D1 g(t)/14 = R_eff(0,t) on C_7: {[ (t, F(g(t),14)) for t in range(1,7)]}")
    C = [0,1,2,3,4,5,6,7]
    print(f"  consec sumW = {sumW(C)} = {float(sumW(C)):.6f}")
    print(f"  D5 perfect short: sumW(consec)={float(sumW(C)):.4f} vs drop e=0 "
          f"-> {float(sumW([1,2,3,4,5,6,7])):.4f}")
    print("  See docstring for full CONFIRMED / FAILED verdict.")
