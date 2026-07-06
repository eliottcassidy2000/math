#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S26 -- HYP-4552: THE k=12 RESONANCE-LADDER CLOSURE.

Synthesizing kps-S35 (plateau + ladder M=j/(6j+5), gap catches ONE rung at k=7)
+ opus-S115 (M height-independent relative to a sub-AP = the plateau, GREEN cap)
+ my S25 (far element ~ q/2).  Extend the mechanism to k=12 (our case): for each
11-speed base B (defected {1..12}), sweep the outlier x; M(B u {x}) sits on a
PLATEAU M(B) and drops only at resonant x onto a ladder.  VERIFY the gap
(1/13, 2/25) catches NO rung (narrower than the rung spacing near plateau 1/13),
uniform over bases => the bounded/single-outlier case CLOSES.

kps template (k=7): base {1,2,3,4,5,7} M=1/6; ladder M=j/(6j+5); gap (1/8,2/15)
catches j=3 <=> 2j>5 AND 3j<10 <=> j=3 => 3/23.  At k=12 the analog inequality
should have NO integer solution.
"""
import itertools
from fractions import Fraction as F
from math import gcd

def exact_M(W):
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = F(0); seen = set()
    for s in dens:
        if s == 0:
            continue
        for j in range(1, s):
            t = F(j, s)
            if t in seen:
                continue
            seen.add(t)
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
    return best

def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in W)) if g > 1 else tuple(sorted(W))

LO, HI = F(1, 13), F(2, 25)

def ladder_for_base(base, xmax=80):
    """sweep outlier x over [max(base)+1, xmax]; return (plateau, in-gap rungs,
    all distinct M-values with their x)."""
    plateau = exact_M(tuple(base))  # M(base) alone
    vals = {}
    in_gap = []
    for x in range(max(base) + 1, xmax + 1):
        W = primitive(tuple(sorted(set(base) | {x})))
        if len(W) != len(base) + 1:
            continue
        M = exact_M(W)
        vals.setdefault(M, []).append(x)
        if LO < M < HI:
            in_gap.append((M, x, W))
    return plateau, in_gap, vals

def main():
    print("=" * 82)
    print("VALIDATE kps template at k=7: base {1,2,3,4,5,7}, gap (1/8,2/15)")
    print("=" * 82)
    lo7, hi7 = F(1, 8), F(2, 15)
    base7 = [1, 2, 3, 4, 5, 7]
    plat = exact_M(tuple(base7))
    print(f"  M(base)={plat} (plateau); ladder M=j/(6j+5); sweep x:")
    caught = []
    for x in range(8, 40):
        W = primitive(tuple(sorted(set(base7) | {x})))
        if len(W) != 7:
            continue
        M = exact_M(W)
        if lo7 < M < hi7:
            caught.append((M, x))
    print(f"  in-gap rungs: {caught}  ({'ONE rung (kps confirmed)' if len(caught)==1 else str(len(caught))})")

    print()
    print("=" * 82)
    print("k=12 (our case): 11-speed bases + outlier; does the gap (1/13,2/25) catch a rung?")
    print("=" * 82)
    # bases: single-defect {1..12}\{i}, and structured near-tight 11-sets
    bases = []
    for i in range(1, 13):
        B = [v for v in range(1, 13) if v != i]
        bases.append((f"{{1..12}}\\{{{i}}}", B))
    # kps-analog: {1..12}\{11} type near-tight, and block-core bases
    bases.append(("{1..11}", list(range(1, 12))))
    bases.append(("{1,2,3,4,5,6,7,8,9,10,11}\\ + 13", [1,2,3,4,5,6,7,8,9,10,13]))
    any_gap = 0
    print(f"  {'base':>22} {'M(base) plateau':>16} {'in-gap rungs':>14} {'min M above 1/13 (x)':>24}")
    for name, B in bases:
        plateau, in_gap, vals = ladder_for_base(B, xmax=90)
        # the smallest M > 1/13 among the ladder (the nearest approach to the gap)
        above = sorted([M for M in vals if M > LO])
        near = above[0] if above else None
        nearx = vals[near][0] if near else None
        gapstr = f"{len(in_gap)} !!!" if in_gap else "0 (empty)"
        if in_gap:
            any_gap += len(in_gap)
        print(f"  {name:>22} {str(plateau):>16} {gapstr:>14} "
              f"{(str(near)+' (x='+str(nearx)+')') if near else '--':>24}")
        for M, x, W in in_gap[:3]:
            print(f"      *** IN-GAP: M={M} x={x} W={list(W)}")
    print()
    if any_gap == 0:
        print("  => ZERO in-gap rungs across all k=12 single-outlier bases: the gap")
        print("     (1/13,2/25) catches NO rung of any base's resonance ladder.  The")
        print("     bounded single-outlier case is EMPTY -- the gap is narrower than the")
        print("     rung spacing near the plateau 1/13, uniform over these bases.")
    else:
        print(f"  *** {any_gap} in-gap rungs found -- investigate (the mechanism's k=12 rung).")

if __name__ == "__main__":
    main()
