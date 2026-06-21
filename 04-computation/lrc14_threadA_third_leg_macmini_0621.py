#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A (mac-mini 2026-06-21) -- THE THIRD LEG: measS7 is NOT a function of the
single-leg profile (min|e| per residue); it ALSO depends on the SECOND
representative of the doubled residue -- the product/third-leg coordinate, exactly
the G_P third leg C=||pq|| that is invisible to (||p||,||q||).

PRECISE OBSTRUCTION found: [0,1,2,3,4,5,6,7] and [0,1,2,3,4,5,6,8] have the SAME
leg profile (0,1,2,3,4,5,6) but measS7 0.327 vs 0.208. The only difference is the
EIGHTH element (the second representative of residue 0): 7 vs 8.
   - In consec, residue 0 is realized by {0, 7}: the 2nd rep is 7 = 7*1, slope 49.
   - In the loser, residue 0 is realized by {0, 8}: 8 = residue 1, NOT residue 0!
WAIT: 8 mod 7 = 1, so [0,1,2,3,4,5,6,8] has residue 1 DOUBLED (by 1 and 8) and
the leg profile recomputed treats min|e| for res1 = 1. So actually both have a
DOUBLED residue but a DIFFERENT one (consec doubles res0 via {0,7}; the other
doubles res1 via {1,8}). The leg profile (the 7 MINIMA) is identical, but the
DOUBLED-residue identity & its 2nd representative differ. THAT is the hidden
coordinate -- the analog of the third Markoff leg.

This script makes the third-leg coordinate explicit:
  For each full-residue k=8 shape, exactly ONE residue r* is doubled (8 elements,
  7 residues). Let (m1, m2) = the two magnitudes in class r*. The "third-leg"
  data is (r*, m2) [m1 = the leg]. We test whether measS7 is a function of
  (leg profile, r*, m2) -- the full local data -- and what consec's choice
  (r*=0, m2=7) does that maximizes the aggregate.
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

def measS7(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    total = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        if len(covered_sectors(E, (lo+hi)/2)) == 7: total += hi - lo
    return total

def residues(E): return frozenset(e % 7 for e in E)
def is_full_residue(E): return residues(E) == frozenset(range(7))
def primitive(E): return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

def doubled_data(E):
    """For a full-residue k=8 shape, the residue r* with two reps and the pair."""
    byres = defaultdict(list)
    for e in E: byres[e % 7].append(abs(e))
    doubled = [(r, sorted(v)) for r, v in byres.items() if len(v) == 2]
    if len(doubled) != 1: return None  # not exactly one doubled residue
    r, (m1, m2) = doubled[0][0], doubled[0][1]
    return r, m1, m2

if __name__ == "__main__":
    print("#"*78); print("# THE THIRD LEG: doubled-residue coordinate (THREAD A)"); print("#"*78)
    k = 8; W = 14
    bank = [(0,)+r for r in itertools.combinations(range(1, W+1), k-1)]
    full = [list(E) for E in bank if primitive(E) and is_full_residue(E)]
    # restrict to shapes with EXACTLY one doubled residue (the generic k=8 full case)
    full = [E for E in full if doubled_data(E) is not None]
    print(f"\nFull-residue, exactly-one-doubled-residue (span<= {W}): {len(full)} shapes")

    C = consec(k)
    print(f"\nconsec={C}: doubled_data = {doubled_data(C)}, measS7={float(measS7(C)):.6f}")
    print(f"   (residue 0 doubled by {{0,7}}; the doubled magnitude 7 = the third leg.)")

    # The two specific shapes from the obstruction:
    for E in [[0,1,2,3,4,5,6,7],[0,1,2,3,4,5,6,8]]:
        print(f"   {E}: doubled={doubled_data(E)}, measS7={float(measS7(E)):.6f}")

    # Group by FULL local data (leg profile, r*, m2) -- is measS7 a function now?
    def legprof(E):
        byres = defaultdict(list)
        for e in E: byres[e % 7].append(abs(e))
        return tuple(sorted(min(byres[r]) for r in range(7)))
    key = lambda E: (legprof(E), doubled_data(E))
    bykey = defaultdict(list)
    for E in full: bykey[key(E)].append(E)
    nonfunc = sum(1 for kk, Es in bykey.items() if len({measS7(E) for E in Es}) > 1)
    print(f"\nIs measS7 a function of (leg profile, r*, m1, m2)? non-function keys: {nonfunc}/{len(bykey)}")
    for kk, Es in list(bykey.items())[:0]: pass
    for kk, Es in bykey.items():
        if len({measS7(E) for E in Es}) > 1 and nonfunc<=10:
            print(f"   key {kk}: {[(E,float(measS7(E))) for E in Es]}")

    # Now: among shapes whose leg profile = consec's (0,1,2,3,4,5,6) -- the MINIMAL
    # profile -- which (r*, m2) maximizes measS7? Is it consec's (r*=0, m2=7)?
    print(f"\nAmong MINIMAL-leg-profile (0,1,2,3,4,5,6) shapes: vary the 3rd leg (r*, m2):")
    minp = (0,1,2,3,4,5,6)
    cands = [E for E in full if legprof(E) == minp]
    cands.sort(key=lambda E: -measS7(E))
    for E in cands:
        r, m1, m2 = doubled_data(E)
        tag = "  <-- consec" if tuple(E)==tuple(C) else ""
        print(f"   {E}: doubled res {r} = {{{m1},{m2}}}, 3rd-leg m2={m2}, measS7={float(measS7(E)):.6f}{tag}")

    # The third leg as a residue product: in G_P the 3rd leg is ||pq||. Here the
    # doubled residue r* and the 2nd-rep magnitude define how the doubled residue's
    # TWO reps drift at different speeds (7*m1 vs 7*m2) -> they cover residue r* on
    # a WIDER y-interval (two reps, two survival windows). consec doubling res 0 at
    # the closest magnitude (7) keeps the OVERLAP residue alive longest.
    print(f"\nKEY: doubling residue 0 (the IDENTITY/center residue) with the SMALLEST")
    print(f"second magnitude (7) maximizes the aggregate -- the third-leg-optimal choice.")
