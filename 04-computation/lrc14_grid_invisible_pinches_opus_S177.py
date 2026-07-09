"""
lrc14_grid_invisible_pinches_opus_S177.py   (opus-2026-07-09-S177)

GRID-INVISIBLE PINCHES (mac-mini MISTAKE-130) and the lemniscate node.  The good set
  G = {x in [0,1): maxgap{frac(e_i x)} > 1/7}
is PINCHED at the measure-zero rationals where maxgap(x) = 1/7 EXACTLY (the strict boundary): there G
touches its complement, splitting into arcs that meet at a point.  These pinches are INVISIBLE to any
uniform grid (measure zero + the exact-boundary value), so a grid OVER-MERGES adjacent arcs =>
grid-measured maxIntG OVER-estimates the true widest arc (this is why the 117k widest-arc search failed).

The pinch = the LEMNISCATE NODE: maxgap=1/7 is the self-crossing where G pinches (r=0 on r^2=cos2theta);
it is also the TIGHT LRC point (maxgap=1/7 <=> M(S)=1/14 exactly, the extremal), the measure-zero
exact-resonance (opus-S168 broken clock), and the Eisenstein/cyclotomic certificate denominator
(OPEN-Q-110).  DEMONSTRATION: at increasing grid resolution, #arcs of G INCREASES and maxIntG DECREASES
(pinches resolve, over-merged arcs split) -- the signature of grid-invisible pinches.  Also: the
non-strict good set {maxgap >= 1/7} INCLUDES the pinches (mac-mini's surviving criterion, M>=1/14).
"""
import numpy as np
from math import gcd
from functools import reduce

INV7 = 1.0 / 7


def maxgap_series(E, x):
    E = np.asarray(E, float)
    Ph = np.mod(np.outer(x, E), 1.0); Ph.sort(axis=1)
    g = np.diff(Ph, axis=1)
    g = np.concatenate([g, (1.0 - Ph[:, -1] + Ph[:, 0])[:, None]], axis=1)
    return g.max(axis=1)


def arcs_and_maxint(E, NG, strict=True):
    """#arcs of G={maxgap>1/7} (or >=) and maxIntG (widest arc fraction) at grid resolution NG."""
    x = (np.arange(NG) + 0.5) / NG
    mg = maxgap_series(E, x)
    good = (mg > INV7) if strict else (mg >= INV7 - 1e-12)
    gi = good.astype(int)
    up = ((gi - np.roll(gi, 1)) == 1)
    narcs = int(up.sum())
    if gi.all():
        return 1, 1.0
    if narcs == 0:
        return 0, 0.0
    # widest arc (circular): run lengths of 1s
    idx = np.where(up)[0]
    maxrun = 0
    for s in idx:
        t = s; ln = 0
        while gi[t % NG] == 1:
            ln += 1; t += 1
            if ln > NG: break
        maxrun = max(maxrun, ln)
    return narcs, maxrun / NG


def prim(E):
    E = sorted(E)
    return reduce(gcd, [E[i + 1] - E[i] for i in range(len(E) - 1)]) == 1


print("=" * 100)
print("GRID-INVISIBLE PINCHES: #arcs of G and maxIntG*spread vs grid resolution (should DRIFT as pinches")
print("resolve).  main claim: coarse grid OVER-MERGES => over-estimates maxIntG (mac-mini MISTAKE-130).")
print("=" * 100)
clusters = {
    "knife-edge-ish {0,7,14,...}": [0, 7, 14, 21, 26, 29, 37, 44, 51, 58, 67, 75, 82],   # 7-structured
    "dissociated Sidon           ": [0, 1, 3, 7, 12, 20, 30, 44, 60, 78, 95, 115, 140],
    "near-AP 3*{0..9}+{1,2,3}    ": sorted(set([3 * i for i in range(10)] + [1, 2, 40, 41])),
}
for name, E in clusters.items():
    E = sorted(set(E))
    if len(E) != 13:
        print(f"  {name}: |E|={len(E)} skip"); continue
    spread = max(E)
    print(f"\n  {name}  (spread={spread})")
    print(f"    {'NG':>9} {'#arcs(strict)':>13} {'maxIntG*spread':>15} {'#arcs(nonstrict)':>16}")
    prev_mi = None
    for NG in [2000, 20000, 200000, 2000000]:
        na, mi = arcs_and_maxint(E, NG, strict=True)
        na_ns, _ = arcs_and_maxint(E, NG, strict=False)
        drift = "" if prev_mi is None else f"  (d={mi*spread - prev_mi:+.3f})"
        print(f"    {NG:>9} {na:>13} {mi*spread:>15.4f} {na_ns:>16}{drift}")
        prev_mi = mi * spread
print()
print("  READING: if #arcs(strict) RISES and maxIntG*spread FALLS with resolution, the coarse grid was")
print("  OVER-MERGING across grid-invisible pinches (maxgap dipping to exactly 1/7) -- confirming that no")
print("  uniform grid measures the true widest arc.  The pinch = the LEMNISCATE NODE = the tight LRC point")
print("  (M=1/14).  The non-strict good set (>=1/7) INCLUDES the pinches (mac-mini's surviving M>=1/14).")
print("=" * 100)
