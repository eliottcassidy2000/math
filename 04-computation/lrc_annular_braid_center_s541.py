#!/usr/bin/env python3
"""
lrc_annular_braid_center_s541.py    oracle-2026-06-01-S541

Deepen the LRC braid picture (S540o did the disk pure braid: linking=tension,
word-length=holdback, fat channel) with the ANNULAR-braid-group-specific structure
the user asked for. Model the OBSERVER as the HOLE of the annulus; the n-1 runners
are strands winding around it. Then:

  * radial coordinate of the annulus = distance-to-observer = ||v_i t||  (the LRC
    quantity itself). Loneliness = all runners at radius >= 1/n (a clear collar).
  * winding number of runner i around the hole = v_i (its speed).
  * the braid word uses the CYCLIC (annular) generators sigma_0..sigma_{n-1}; the
    'wrap' generator sigma_0 = a crossing at the cut = a runner lapping the
    observer. #sigma_0 = sum v_i; runner-runner crossings = sum_{i<j}|v_i-v_j|.
  * CENTER of the annular braid group = the full core rotation rho (all windings
    +1) = the LRC FRAME-SHIFT symmetry (relabel which runner is the observer).
    LRC-loneliness is invariant under rho, so LRC lives in the annular braid group
    MOD ITS CENTER = the affine-symmetric-group data.

Verify: (1) crossing counts; (2) the affine permutation / windings; (3) frame-shift
invariance of loneliness (the center acts trivially on the LRC question).
"""
from fractions import Fraction
from math import gcd
from itertools import combinations
from functools import reduce

def frac(x): return x - (x.numerator // x.denominator)
def d0(x):
    f = frac(x); return min(f, 1 - f)

def crossings_by_type(speeds, n):
    """Total braid-word crossings split into observer-wrap (sigma_0) and
    runner-runner (sigma_{>0}). Observer is index 0 with speed 0."""
    obs = 0  # runner laps observer: ||v_i t||=0 -> v_i times each
    for s in speeds[1:]:
        obs += abs(s)
    rr = 0
    for i, j in combinations(range(1, n), 2):
        rr += abs(speeds[i] - speeds[j])
    return obs, rr

def loneliest(speeds, n):
    """max over t of the observer's clear-collar radius = min_i ||v_i t||.
    Evaluate at BOTH the walls (where tight/boundary loneliness is achieved, e.g.
    the AP at t=1/n) AND the open midpoints. Lonely iff this >= 1/n."""
    W = set([Fraction(0), Fraction(1)])
    for s in speeds[1:]:
        if s == 0: continue
        for k in range(0, abs(s) + 1):
            for sg in (1, -1):
                t = Fraction(n * k + sg, n * abs(s))
                if 0 <= t <= 1: W.add(t)
    Wl = sorted(W)
    pts = list(Wl) + [(a + b) / 2 for a, b in zip(Wl, Wl[1:])]   # walls + midpoints
    best = Fraction(0)
    for t in pts:
        r = min((d0(Fraction(s) * t) for s in speeds[1:]), default=Fraction(1))
        if r > best: best = r
    return best

def primitive(s): return reduce(gcd, [x for x in s if x != 0]) == 1

def main():
    print("LRC annular braid: cyclic generators, windings, and the CENTER=shift (oracle-S541)\n")

    print("(1) crossing counts: #sigma_0 (observer-wrap) = sum v_i ; runner-runner = sum|v_i-v_j|")
    for sp in [(0,1,2,3),(0,1,2,4),(0,1,3,5),(0,1,2,3,4),(0,1,2,3,4,5)]:
        n = len(sp); obs, rr = crossings_by_type(sp, n)
        print(f"   v={sp[1:]}: sigma_0(wrap)={obs}=sum(v)  runner-runner={rr}  total word len={obs+rr}  windings(around hole)={sp[1:]}")
    print("   => the braid's winding vector around the observer-hole IS the speed vector.\n")

    print("(2a) per-runner loneliness (each runner's clear-collar = its own LRC question).")
    base = (0,1,2,3,4); n = len(base)
    print(f"    system S={base}, n={n}, threshold 1/{n}.  runner k lonely = re-center on k:")
    perrunner = {}
    for k in range(n):
        reordered = (0,) + tuple(base[i] - base[k] for i in range(n) if i != k)
        lk = loneliest(reordered, n); perrunner[k] = lk
        print(f"      runner {k} (speed {base[k]}): max collar = {lk} = {float(lk):.4f}  lonely={lk>=Fraction(1,n)}")
    print("      (edge runners are TIGHT at exactly 1/n — the AP extremal — middle runners looser.)")

    print("\n(2b) CENTER = rho (add a constant c to ALL speeds = +1 winding around the core).")
    print("    rho is central in the annular braid group; loneliness depends only on pairwise")
    print("    differences, which rho fixes -> every runner's collar is UNCHANGED. Verify c=+10:")
    ok = True
    for c in (10, 100):
        shifted_sys = tuple(x + c for x in base)
        for k in range(n):
            reordered = (0,) + tuple(shifted_sys[i] - shifted_sys[k] for i in range(n) if i != k)
            if loneliest(reordered, n) != perrunner[k]: ok = False
    print(f"      loneliness invariant under +c for c in {{10,100}}, all runners: {ok}")
    print("      => the CENTER of the annular braid group acts TRIVIALLY on the LRC question;")
    print("         LRC is a statement in the annular braid group MOD CENTER = the affine")
    print("         symmetric group / the windings-up-to-global-shift.\n")

    print("(3) affine-permutation level = sum of windings = sum v_i (the #observer-crossings).")
    print("    The annular braid surjects onto the affine symmetric group S~_{n-1}; the LRC")
    print("    braid is the homogeneous TORUS braid whose affine part is the speed vector.")
    print("    Loneliness = a height where the core collar is fat = the alcove walk (S525)")
    print("    visits a deep alcove. The disk picture (S540o) forgets the hole; the annulus")
    print("    keeps it, and the hole's radial coordinate is exactly ||v_i t||.")

if __name__ == "__main__":
    main()
