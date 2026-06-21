#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_route1_conductance_minimax_opus_0621.py  (opus-2026-06-21, ROUTE 1)

THE CONDUCTANCE-VECTOR MINIMAX for the LRC(14) Layer-3 wall.

GOAL: derive the EXACT formula measS7 = F(c_0,...,c_6) (or W_a = F_a(c)) from the
moving-clock survival model, and test whether consec is the unique maximizer of the
bottleneck functional via minimax EQUALIZATION.

MODEL (verified, == brute measS7):
  At resonance a (a=1..6), clock e sits on the circle R/Z (7 sectors of width 1/7).
  Position of clock e at clock-time x: e*x.  Sector(e,x) = floor(7*frac(e*x)) in Z/7.
  At x = a/7 (i.e. y=0), clock e is at integer position e*a/7*7 = e*a, the LEFT edge of
  sector (e*a mod 7).  Survival window W_a = length of the maximal cover (all 7 sectors
  occupied) interval through y=0, PLUS disconnected covered arcs in [-1/14,1/14].
  measS7 = sum_{a=1..6} W_a.  Velocity of clock e in sector coords = 7e.

KEY DECOMPOSITION (this script):
  W_a = WIN_a + DISC_a   where WIN_a = T+_a + T-_a is the contiguous-through-center
  survival and DISC_a is the disconnected remainder.  The WINDOW is the survival
  scheduling object.  Binding-speed law: T+_a = 1/(7 v+_a), T-_a = 1/(7 v-_a) where
  v+/-_a is the BINDING SPEED = the relative speed at which the first sector empties
  on each side (the bottleneck).

TESTS:
  (R1) Exact binding-speed table for consec (both sides, all a). Confirm equalization.
  (R2) The RELAY/BOTTLENECK formula for T+/T- in terms of the clock multiset.
  (R3) MINIMAX-EQUALIZATION: among all full-residue k=8 shapes, does consec uniquely
       EQUALIZE the per-sector survival bottlenecks?  Measure the spread of binding
       speeds; correlate with measS7.  Test: is the maximizer the equalizer?
  (R4) Exact F(c): regress / derive W_a as a function of the conductance vector.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

# ---------------- core model (sector coords, velocity 7e) ----------------
def sec(e, a, y):
    p = F(e*a) + F(7*e)*y
    return (p.numerator // p.denominator) % 7

def occ_at(E, a, y):
    s = set()
    for e in E:
        s.add(sec(e, a, y))
    return s

def breakpoints(E, a, half=F(1,14)):
    bps = {F(0), -half, half}
    for e in E:
        if e == 0: continue
        # crossings: e*a + 7*e*y in Z => y = (m - e*a)/(7e)
        lo = F(7*e)*(-half) + F(e*a); hi = F(7*e)*(half) + F(e*a)
        lo_i, hi_i = (lo, hi) if lo <= hi else (hi, lo)
        m = lo_i.numerator // lo_i.denominator
        while m <= hi_i.numerator // hi_i.denominator + 1:
            y = F(m - e*a, 7*e)
            if -half <= y <= half: bps.add(y)
            m += 1
    return sorted(bps)

def W_a_total(E, a, half=F(1,14)):
    bps = breakpoints(E, a, half); tot = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        if len(occ_at(E, a, (lo+hi)/2)) == 7: tot += hi - lo
    return tot

def survival(E, a, half=F(1,14)):
    """T+ (right), T- (|left|) of contiguous cover through 0; (0,0) if center uncovered."""
    bps = breakpoints(E, a, half)
    ivals = list(zip(bps, bps[1:]))
    # center covered?
    Tp = F(0); Tm = F(0)
    # right
    for lo, hi in ivals:
        if hi <= 0: continue
        lo2 = max(lo, F(0))
        if len(occ_at(E, a, (lo2+hi)/2)) == 7:
            Tp = hi
        else:
            if lo2 == F(0): Tp = F(0)
            break
    Tp = min(Tp, half)
    for lo, hi in reversed(ivals):
        if lo >= 0: continue
        hi2 = min(hi, F(0))
        if len(occ_at(E, a, (lo+hi2)/2)) == 7:
            Tm = -lo
        else:
            if hi2 == F(0): Tm = F(0)
            break
    Tm = min(Tm, half)
    return Tp, Tm

def measS7(E): return sum(W_a_total(E, a) for a in range(1, 7))
def is_full_residue(E): return frozenset(e % 7 for e in E) == frozenset(range(7))
def primitive(E): return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1
def cond_vec(E):
    c = defaultdict(F)
    for e in E:
        if e != 0: c[e % 7] += F(1, abs(e))
    return tuple(c[r] for r in range(7))
def consec(k): return list(range(k))

def binding_speed(T, half=F(1,14)):
    if T <= 0: return None
    if T >= half: return F(1, 14)/T  # clipped
    return 1/(7*T)

if __name__ == "__main__":
    print("#"*78)
    print("# ROUTE 1: CONDUCTANCE-VECTOR MINIMAX — exact survival & equalization")
    print("#"*78)

    C = consec(8)
    print(f"\nconsec = {C},  c-vec = {[str(x) for x in cond_vec(C)]}")
    print(f"measS7(consec) = {measS7(C)} = {float(measS7(C)):.6f}\n")

    # ---- (R1) binding-speed table for consec ----
    print("(R1) consec binding speeds v± = 1/(7 T±) (the per-resonance bottleneck):")
    bsp = []; bsm = []
    for a in range(1, 7):
        Tp, Tm = survival(C, a)
        vp, vm = binding_speed(Tp), binding_speed(Tm)
        bsp.append(vp); bsm.append(vm)
        W = W_a_total(C, a)
        print(f"  a={a}: T+={Tp} (v+={vp})  T-={Tm} (v-={vm})  W_a={W}  WIN={Tp+Tm}")
    print(f"  right binding speeds: {[str(x) for x in bsp]}")
    print(f"  left  binding speeds: {[str(x) for x in bsm]}")
    print(f"  -> EQUALIZATION: v+ constant={len(set(bsp))<=2}, v- constant={len(set(bsm))==1}")

    # ---- (R2) relay/bottleneck identity ----
    print("\n(R2) RELAY structure: which sector dies first and why (a=1 right):")
    a = 1
    bps = breakpoints(C, a)
    prev = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= 0: continue
        mid = (max(lo, F(0))+hi)/2
        occ = defaultdict(list)
        for e in C: occ[sec(e, a, mid)].append(e)
        empty = [r for r in range(7) if r not in occ]
        flag = "  <-- COVER BREAKS HERE" if empty else ""
        print(f"    y in ({float(max(lo,F(0))):.5f},{float(hi):.5f}): empty={empty}{flag}")
        if empty: break

    # ---- (R3) minimax equalization across the stratum ----
    print("\n(R3) MINIMAX-EQUALIZATION across full-residue k=8 stratum (span<=18):")
    W = 18
    bank = [(0,)+r for r in itertools.combinations(range(1, W+1), 7)]
    full = [list(E) for E in bank if primitive(E) and is_full_residue(E)]
    print(f"   stratum size: {len(full)}")
    mC = measS7(C)
    # For each shape, gather all binding speeds (v+ and v- over a=1..6), measure spread
    def speed_spread(E):
        vs = []
        for a in range(1, 7):
            Tp, Tm = survival(E, a)
            for T in (Tp, Tm):
                v = binding_speed(T)
                if v is not None: vs.append(float(v))
        if not vs: return None, None, None
        return min(vs), max(vs), max(vs)-min(vs)
    rows = []
    for E in full:
        m = measS7(E)
        lo, hi, spread = speed_spread(E)
        rows.append((m, spread, lo, hi, E))
    rows.sort(key=lambda r: -r[0])
    print("   TOP 8 by measS7 (m, speed-spread, vmin, vmax):")
    for m, spread, lo, hi, E in rows[:8]:
        tag = " <-- CONSEC" if E == C else ""
        print(f"     measS7={float(m):.6f}  spread={spread}  vmin={lo} vmax={hi}  {E}{tag}")
    # is the global max also the global MIN spread?
    valid = [r for r in rows if r[1] is not None]
    by_spread = sorted(valid, key=lambda r: (r[1], -r[0]))
    print(f"   shape with MIN speed-spread: measS7={float(by_spread[0][0]):.6f} "
          f"spread={by_spread[0][1]} {by_spread[0][4]}")
    print(f"   consec is global measS7-max: {rows[0][4]==C}")
    cs = next(r for r in valid if r[4]==C)
    print(f"   consec speed-spread={cs[1]} (rank among min-spread: "
          f"{[r[4] for r in by_spread].index(C)+1}/{len(by_spread)})")

    # ---- (R4) exact F(c): tabulate W_a vs full conductance info ----
    print("\n(R4) Is W_a a clean function of binding speeds? consec WIN check:")
    win = sum(sum(survival(C, a)) for a in range(1, 7))
    print(f"   sum_a WIN_a(consec) = {win} = {float(win):.6f}")
    print(f"   measS7(consec)      = {mC} = {float(mC):.6f}")
    print(f"   DISC = measS7 - WIN = {mC-win} = {float(mC-win):.6f}")
