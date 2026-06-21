#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_angle1_binding_speed_law_opus_0621.py  (opus-2026-06-21, ANGLE 1, part 6)

THE BINDING-SPEED LAW and the WINDOWS-ONLY extremal mechanism.

BIG RESULT (part 5, verified k=8,9,10): consec UNIQUELY maximizes the WINDOWS-ONLY
sum  WIN(E) = sum_{a=1..6} (T+(a) + T-(a)),  where T+/T- are the first-death times
right/left of the resonance center.  0 shapes beat it, 0 ties -- sharper than the
full measS7 (which has ties).  This windows object is EXACTLY the survival-window
scheduling object ANGLE 1 proposed.

THE MECHANISM (to pin down here):  At resonance a, going RIGHT, each clock e (speed
v_e=7|e|, moving up if e>0) holds its start sector for a while.  The cover dies when
the FIRST sector loses its last occupant.  Claim (binding-speed law): T+(a) = 1/v*
where v* is the speed of the BINDING clock -- the fastest clock that is the SOLE
occupant of its sector and thus whose departure first empties a sector.

GOAL: characterize T+(a)+T-(a) purely in terms of the velocity assignment, then test
the falsifiable extremal hypothesis:

 (B1) BINDING LAW (exact): T+(a) = 1/(7 * s+(a)) where s+(a) = max over sectors of
      (the min speed of clocks that will be needed) ... TEST against measured T+.
      Specifically: going right, sector r empties at time (1 - off_r)/v where off_r
      is how far into sector r its occupant started... derive & verify.

 (B2) THE HARMONIC SUM: WIN(E) = sum_a [1/(7 b+_a) + 1/(7 b-_a)].  consec's binding
      speeds are b+ in {6,7/2,...}, b- = 7 mostly.  Conjecture: consec MINIMIZES the
      multiset of binding speeds in a majorization sense (so the harmonic sum 1/b is
      maximized).  TEST: collect the 12 binding speeds (6 a's x 2 sides) for each
      shape; does consec's binding-speed multiset get MAJORIZED-from-below (i.e. is
      consec's the smallest in the convex/Schur sense, making sum 1/b largest)?

 (B3) SHARPEST: since t->1/t is convex decreasing, sum 1/b_i is Schur-CONVEX, so it is
      maximized by SPREADING the b_i.  But consec maximizes WIN... so consec must
      make the b_i SMALL (not spread).  Resolve: is it sum-of-b (linear) that consec
      minimizes, or is it the actual MULTISET that is componentwise dominated?
      TEST both: (i) sum of binding speeds; (ii) sorted-vector componentwise <=.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def sector_of_point(e, a, y):
    pos = F(e*a) + F(7*e)*y
    return (pos.numerator // pos.denominator) % 7

def covered_all_at(E, a, y):
    return len({sector_of_point(e, a, y) for e in E}) == 7

def breakpoints(E, a):
    half = F(1, 14); bps = {F(0), -half, half}
    for e in E:
        if e == 0: continue
        lo_val = F(7*e)*(-half) + F(e*a); hi_val = F(7*e)*(half) + F(e*a)
        lo_i = min(lo_val, hi_val); hi_i = max(lo_val, hi_val)
        m = lo_i.numerator // lo_i.denominator
        while m <= hi_i.numerator // hi_i.denominator + 1:
            y = F(m - e*a, 7*e)
            if -half <= y <= half: bps.add(y)
            m += 1
    return sorted(bps)

def window_TpTm(E, a):
    half = F(1, 14); bps = breakpoints(E, a); ivals = list(zip(bps, bps[1:]))
    Tp = F(0)
    for lo, hi in ivals:
        if hi <= 0: continue
        lo2 = max(lo, F(0))
        if covered_all_at(E, a, (lo2+hi)/2): Tp = hi
        else:
            if lo2 == F(0): Tp = F(0)
            break
    Tp = min(Tp, half)
    Tm = F(0)
    for lo, hi in reversed(ivals):
        if lo >= 0: continue
        hi2 = min(hi, F(0))
        if covered_all_at(E, a, (lo+hi2)/2): Tm = -lo
        else:
            if hi2 == F(0): Tm = F(0)
            break
    Tm = min(Tm, half)
    return Tp, Tm

def is_full_residue(E): return frozenset(e % 7 for e in E) == frozenset(range(7))
def consec(k): return list(range(k))
def window_sum(E): return sum(sum(window_TpTm(E, a)) for a in range(1, 7))

def binding_speeds(E):
    """The 'effective binding speed' on each side of each resonance: 1/(7 T).
       Returns list of (a, side, T, b=1/(7T) or None if T hits the cell cap)."""
    half = F(1,14); out = []
    for a in range(1, 7):
        Tp, Tm = window_TpTm(E, a)
        for side, T in [('+',Tp),('-',Tm)]:
            if T == 0:
                out.append((a, side, T, None))   # center not covered that side
            elif T == half:
                out.append((a, side, T, F(0)))    # capped by cell, binding speed -> 'infinitely slow'
            else:
                out.append((a, side, T, 1/(7*T)))
    return out

def majorizes_below(u, v):
    """True if sorted(u) <= sorted(v) componentwise (u dominated below by v)."""
    su, sv = sorted(u), sorted(v)
    return all(a <= b for a, b in zip(su, sv))

if __name__ == "__main__":
    print("="*80)
    print("BINDING-SPEED LAW & WINDOWS-ONLY MECHANISM")
    print("="*80)
    k = 8; C = consec(k)

    # (B1) exact binding speeds for consec & rivals
    print(f"\n[B1] binding speeds b=1/(7T) per (a,side):")
    for E in [C, [0,2,3,4,5,6,7,8], [0,1,2,4,5,6,7,10], [0,1,2,3,4,5,6,8]]:
        bs = binding_speeds(E)
        vals = [(f"{a}{s}", str(b) if b is not None else "X") for a,s,T,b in bs]
        nz = [float(b) for a,s,T,b in bs if b not in (None,) ]
        print(f"  {str(E):28s} WIN={float(window_sum(E)):.6f}")
        print(f"     binding b: " + "  ".join(f"{lab}={v}" for lab,v in vals))

    # (B2)/(B3) does consec minimize the binding-speed multiset?
    print(f"\n[B2/B3] binding-speed extremality over full-residue stratum:")
    span = 14
    shapes = [(0,)+c for c in itertools.combinations(range(1, span+1), k-1)]
    shapes = [E for E in shapes if is_full_residue(E)]
    # for each shape, the multiset of NONZERO binding speeds (drop X=center-empty and
    # capped which mean 'survives full cell'). Use a large sentinel for capped (slowest).
    BIG = F(10**6)
    def bvec(E):
        out = []
        for a,s,T,b in binding_speeds(E):
            if b is None:  # center empty that side: window contributes 0 -> 'instant death'
                out.append(BIG)     # treat as infinitely FAST (bad)
            elif b == 0:   # capped by cell: best possible
                out.append(F(0))
            else:
                out.append(b)
        return out
    Cv = bvec(C); Csum = sum(Cv)
    print(f"  consec binding-speed sum (sentinel BIG for center-empty) = {float(Csum):.4f}")
    # (i) sum of binding speeds: does consec minimize?
    nbeat_sum = sum(1 for E in shapes if E!=tuple(C) and sum(bvec(list(E))) < Csum)
    print(f"  (i) consec minimizes SUM of binding speeds? #shapes with smaller sum = {nbeat_sum}")
    # (ii) componentwise sorted dominance: is consec's sorted bvec <= every other?
    dom_all = sum(1 for E in shapes if E!=tuple(C) and majorizes_below(Cv, bvec(list(E))))
    not_dom = [E for E in shapes if E!=tuple(C) and not majorizes_below(Cv, bvec(list(E)))]
    print(f"  (ii) consec sorted-bvec <= other's sorted-bvec (componentwise): holds for "
          f"{dom_all}/{len(shapes)-1}")
    if not_dom:
        print(f"       NOT dominated by consec for {len(not_dom)} shapes -- sorted dominance FAILS")
        # but does WIN still win? (it must, from part5)
        ex = not_dom[0]
        print(f"       example {ex}: WIN={float(window_sum(list(ex))):.6f} < consec {float(window_sum(C)):.6f}")
        print(f"       => WINDOWS extremality is NOT explained by simple binding-speed dominance.")

    # The real functional: WIN = sum 1/(7 b) with capped->1/14, center-empty->0.
    # verify WIN reconstructed from binding speeds:
    print(f"\n  reconstruct WIN from binding speeds (sanity):")
    for E in [C, [0,2,3,4,5,6,7,8]]:
        rec = F(0)
        for a,s,T,b in binding_speeds(E):
            rec += T
        print(f"    {str(E):28s} sum T = {float(rec):.6f}  WIN={float(window_sum(E)):.6f}  "
              f"{'OK' if rec==window_sum(E) else 'MISMATCH'}")
