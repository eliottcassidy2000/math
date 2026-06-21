#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_angle1_survival_scheduling_opus_0621.py  (opus-2026-06-21, ANGLE 1)

THE SURVIVAL-WIDTH SCHEDULING VIEW of W_a (the LRC(14) WALL, LAYER 3).

GOAL: model W_a EXACTLY as a coverage/scheduling problem on Z/7 with drifting
clocks, verify the model against the brute (cell-decomposition) W_a, then TEST a
falsifiable SCHEDULING-EXTREMAL hypothesis: that the AP (consec) is the schedule
maximizing sum_a W_a because its speeds {0,7,14,...,7(k-1)} are maximally
"staggered".

THE EXACT LOCAL MODEL (derived from threadA + verified here):
  At resonance a (a=1..6), put x = a/7 + y, y near 0.  For clock e:
     7 frac(e x) = 7 frac(e a/7 + e y) = (e a mod 7) + (7 e y) mod 7  (continuous lift)
  The SECTOR of e is floor of that, i.e. e sits at sector
     s_e(y) = ( e a + floor(7 e y + theta_e) ) mod 7
  where at y=0 it sits at residue (e a mod 7) -- the START of sector (e a mod 7).
  Actually the boundary: at y=0, 7 frac(e x) = (e a mod 7) exactly (an INTEGER),
  so e sits AT the lower boundary of sector (e a mod 7).  As y increases, e moves
  UP through sectors at rate 7e (per unit y); as y decreases, DOWN.

  So in the local (y) coordinate, each clock e is a point moving on the circle
  R/7Z at signed velocity 7e, starting (at y=0) exactly at the integer point
  (e a mod 7).  Sector r in {0..6} is COVERED iff some clock lies in [r, r+1).
  W_a = the length of the maximal y-interval around 0 on which all 7 sectors
  are covered.  (The cell is [-1/14, 1/14] in y; cover may end inside or at cell
  boundary; W_a is the intersection.)

This is EXACTLY a 1-D moving-points-cover-the-circle problem:
  - 7 unit sectors on a circle of circumference 7.
  - k points (the clocks), point e at position p_e(y) = e*a + 7*e*y  (mod 7).
  - the cover survives while every sector [r,r+1) contains >=1 point.
  - W_a = first-failure time, symmetric-ish around y=0.

REDUCTION to PER-SECTOR SURVIVAL ARCS (the scheduling core):
  Sector r empties (going right, y>0) when its LAST point leaves, i.e. the point
  in r with the SMALLEST velocity-to-exit.  A point at integer start q=e*a mod 7
  occupies sector q at y=0 and sector q for y in (something).  Actually since all
  starts are integers, at y=0 a point sits at the LEFT edge of its start sector.
  For y>0 (moving up if e>0) it's inside sector q for y in (0, 1/(7e)).  For y<0
  it's in sector q-1.

  This makes each sector's right-survival a MIN over a discrete schedule.  W_a^+
  (right survival) = the y>0 time until the first sector goes empty.  Similarly
  W_a^- for y<0.  W_a = W_a^+ + W_a^-  (clipped to the half-cell 1/14 each side).

TESTS:
 (V) VERIFY the moving-points model reproduces the brute W_a exactly (consec +
     several adversaries, k=8).
 (S1) SCHEDULING DECOMPOSITION: compute W_a^+ and W_a^- and the binding sector
      for each, exact.  Show consec's structure.
 (S2) THE STAGGERING HYPOTHESIS (falsifiable): the AP speeds are an arithmetic
      progression 0,7,...,7(k-1).  Conjecture: among all speed-multisets {7e},
      the AP maximizes sum_a (W_a^+ + W_a^-).  We will test a SHARPER local
      claim: at each resonance a, the residue-to-startposition map e -> e*a mod 7
      is a BIJECTION-with-doubling (residue 0 doubled).  For the AP the start
      positions a*{0..k-1} mod 7 are themselves an AP mod 7, i.e. EQUISPACED /
      round-robin.  Hypothesis H1: "round-robin start + smallest speeds = max
      survival".  Test by perturbing.
 (S3) ROUND-ROBIN EXTREMAL TEST: replace consec by a shape with the SAME residue
      multiset (full residues, residue 0 doubled) but DIFFERENT magnitudes, so
      same start-positions per a but different speeds.  Does increasing any speed
      (making a clock faster) strictly DECREASE sum_a W_a?  (monotonicity in
      speed = the scheduling-extremal mechanism).
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

# ---------------------------------------------------------------------------
# BRUTE W_a via cell decomposition (ground truth, from threadA)
# ---------------------------------------------------------------------------
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

def W_a_brute(E, a):
    arcs = measS7_arcs(E)
    lo_b, hi_b = F(2*a-1, 14), F(2*a+1, 14)
    w = F(0)
    for lo, hi in arcs:
        l = max(lo, lo_b); h = min(hi, hi_b)
        if h > l: w += h - l
    return w

# ---------------------------------------------------------------------------
# SCHEDULING MODEL of W_a
# ---------------------------------------------------------------------------
# Local coordinate y in [-1/14, 1/14].  Point e at position p_e(y) = (e*a + 7*e*y) mod 7
# on circle of circumference 7.  Sector r=0..6 covered iff some point in [r, r+1) mod 7.
# We compute the maximal interval [y_lo, y_hi] containing 0 on which all 7 covered,
# then clip to [-1/14, 1/14].
#
# Each point e contributes coverage of sector r over a y-interval. Since starts are
# integers q = e*a mod 7 and speed 7e, the point sits in sector ((q + floor(7 e y)) mod 7).
# For e>0, y>0: occupies sector (q + j) mod 7 for y in (j/(7e), (j+1)/(7e)), j=0,1,2,...
# For e>0, y<0: occupies sector (q - 1 + j ... ) -- handle via the same with negative.
# Easiest EXACT approach: collect all breakpoints (y where any point crosses a sector
# boundary) within the cell, evaluate coverage on each sub-cell midpoint. This is the
# SAME as brute but in local coords -- so instead we directly model survival arcs.

def sector_of_point(e, a, y):
    """sector in {0..6} of clock e at local time y (circle circumference 7)."""
    pos = (F(e*a) + F(7*e)*y)  # real position on R, mod 7
    pos = pos - (pos.numerator // pos.denominator) * 1  # not mod7 yet
    # reduce mod 7
    q = pos / 7
    pos = pos - (q.numerator // q.denominator) * 7
    if pos < 0: pos += 7
    # sector = floor(pos)
    return (pos.numerator // pos.denominator) % 7

def covered_all_at(E, a, y):
    secs = set()
    for e in E:
        if e == 0:
            secs.add((0 * a) % 7)  # stationary at sector 0
            continue
        secs.add(sector_of_point(e, a, y))
    return len(secs) == 7

def W_a_schedule(E, a):
    """Exact W_a from the local moving-points model. Breakpoints = sector crossings."""
    half = F(1, 14)
    bps = {F(0)}
    for e in E:
        if e == 0: continue
        ae = abs(e)
        # crossings happen when 7*e*y + e*a is integer => y = (m - e*a)/(7*e)
        # within (-1/14, 1/14): 7*e*y in (-e/2, e/2) -> integer offsets near e*a
        # enumerate integer values of (7*e*y + e*a) in the cell
        lo_val = F(7*e)*(-half) + F(e*a)
        hi_val = F(7*e)*(half) + F(e*a)
        lo_i = min(lo_val, hi_val); hi_i = max(lo_val, hi_val)
        m = lo_i.numerator // lo_i.denominator
        while m <= hi_i.numerator // hi_i.denominator + 1:
            # y where 7*e*y + e*a = m
            y = F(m - e*a, 7*e)
            if -half <= y <= half:
                bps.add(y)
            m += 1
    bps = sorted(bps | {-half, half})
    # find maximal run of covered sub-cells containing y=0
    # build covered flags per sub-cell
    subs = []
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        mid = (lo+hi)/2
        subs.append((lo, hi, covered_all_at(E, a, mid)))
    # the sub-cell containing 0:
    total = F(0)
    # find index of cell containing 0 (0 is a breakpoint; take cells adjacent that are covered)
    # accumulate contiguous covered cells spanning across 0
    # left from 0:
    n = len(subs)
    # locate cell with lo<=0<hi  OR cells adjacent to 0 boundary
    # Since 0 is a breakpoint, two cells touch 0: the one ending at 0 and starting at 0.
    # accumulate rightward
    i = 0
    while i < n and subs[i][1] <= 0:
        i += 1
    # subs[i] starts at 0 (or contains 0); accumulate covered cells from i upward
    j = i
    while j < n and subs[j][2]:
        total += subs[j][1] - subs[j][0]
        j += 1
    # accumulate covered cells from i-1 downward
    j = i - 1
    while j >= 0 and subs[j][2]:
        total += subs[j][1] - subs[j][0]
        j -= 1
    return total

# ---------------------------------------------------------------------------
def residues(E): return frozenset(e % 7 for e in E)
def is_full_residue(E): return residues(E) == frozenset(range(7))
def primitive(E): return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

def measS7_sum(E):
    return sum(W_a_brute(E, a) for a in range(1, 7))

if __name__ == "__main__":
    print("="*78)
    print("ANGLE 1: SURVIVAL-WIDTH SCHEDULING MODEL of W_a")
    print("="*78)

    k = 8
    tests = {
        "consec":   consec(k),
        "a1-beater":[0,2,3,4,5,6,7,8],
        "a5-beater":[0,1,2,4,5,6,7,10],
        "rand1":    [0,1,2,3,4,5,7,9],
        "rand2":    [0,1,3,4,5,6,8,11],
    }

    print("\n[V] VERIFY scheduling model == brute W_a (per a, and SUM):")
    all_ok = True
    for name, E in tests.items():
        print(f"\n  {name} = {E}")
        tot_b = F(0); tot_s = F(0); ok = True
        for a in range(1, 7):
            wb = W_a_brute(E, a)
            ws = W_a_schedule(E, a)
            tot_b += wb; tot_s += ws
            match = (wb == ws)
            ok = ok and match
            print(f"    a={a}: brute={float(wb):.6f}  sched={float(ws):.6f}  {'OK' if match else 'MISMATCH!'}")
        print(f"    SUM brute={float(tot_b):.6f}  sched={float(tot_s):.6f}  {'OK' if tot_b==tot_s else 'MISMATCH!'}")
        all_ok = all_ok and ok and (tot_b == tot_s)
    print(f"\n  => scheduling model {'EXACTLY reproduces' if all_ok else 'DEVIATES from'} brute W_a.")
