#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_p0_vs_rhostar_bridge_kpswf10.py   (kind-pasteur 2026-06-21, THREAD B)

THE highest-value THREAD B question: pin the EXACT relationship between the two
LRC(14) phi-phase objects, resolving the 1/7 vs 2/7 factor-of-2.

  SECTOR ROUTE (L7):  p0(E)  = meas{ x in [0,1) : {sector(e_i x): i} = Z/7 }
                              = meas{ all 7 of the 1/7-sectors are HIT by {frac(e_i x)} }.
                       sector(y) = floor(7 frac(y)).  p0 <= cap_k => survival.

  THM-527 (rho*):     g*(E)  = meas{ x in [0,1) : maxgap{frac(e_i x)} > 2/7 }   (ignore G_P)
                       rho*   = meas{ x in G_P : maxgap{frac(e_i x)} > 2/7 }.
                       g*(E) > 0 (with G_P) => M(S) >= 1/14   (the LRC bound).

THRESHOLDS:
  - p0 is about a 1/7-scale event (each 1/7-sector hit).
  - g* is about a 2/7-scale event (a >2/7 GAP in the phases).

CLEAN LOGICAL FACTS we test EXACTLY:
  (F1)  maxgap > 2/7  ==>  some 1/7-sector lies ENTIRELY inside the gap, hence is MISSED.
        So  {maxgap>2/7}  subset  {not all 7 sectors hit}  =>  g*  <=  1 - p0.
  (F2)  Is the reverse anywhere near?  Define the family of intermediate thresholds:
            gapT(E; t)   = meas{ maxgap{frac(e_i x)} > t },
            sect1(E)     = meas{ >= 1 sector entirely missed } = 1 - p0   (sectors all hit complement)
        and the "two-adjacent-sectors missed" object
            sect2adj(E)  = meas{ two CYCLICALLY-ADJACENT 1/7-sectors are BOTH missed }.
        CLAIM to test:  {two adjacent sectors missed}  ==>  maxgap > 1/7  (a 1/7 gap, NOT 2/7),
        and conversely {maxgap > 2/7} ==> two adjacent (in fact, an open 2/7) missed.
  (F3)  The SECTOR ROUTE's natural gap is 1/7 (one missed sector = a >= 1/7 hole boundary),
        which is HALF of THM-527's 2/7 need.  Does the closed p0<=cap give the 2/7 event,
        or only the (weaker) 1/7 event?  -> resolve whether the routes COINCIDE or are PARALLEL.

This script computes, EXACTLY (rational), for a battery of clusters E:
   p0(E), 1-p0, g*(E)=gap>2/7, gap>1/7, sect2adj, and the inclusion chain,
   verifying F1 (g* <= 1-p0) and quantifying the gap between g* and 1-p0.
"""
import itertools
from fractions import Fraction as Fr

P = 7                       # number of sectors / the modulus 7
THRESH = Fr(2, 7)           # THM-527 good-gap threshold (the '2/7')


# ----------------------------------------------------------------- helpers
def sector(yfrac):
    """floor(7*frac(y)) for yfrac a Fraction in [0,1)."""
    return int(P * yfrac)


def phase_breakpoints(E, extra_dens=None):
    """All x-breakpoints where any frac(e_i x) crosses a multiple of 1/7
       (sector boundaries) -- the right grid for BOTH p0 and the >2/7 gap test,
       because the gap endpoints are cluster phases (handled per-cell at the mid)."""
    E = [int(e) for e in E]
    bp = {Fr(0), Fr(1)}
    for e in E:
        if e == 0:
            # frac(0*x)=0 always; contributes the fixed phase 0 (no breakpoints)
            continue
        for t in range(0, P * e + 1):
            bp.add(Fr(t, P * e))
    if extra_dens:
        for d in extra_dens:
            for t in range(0, d + 1):
                bp.add(Fr(t, d))
    return sorted(b for b in bp if Fr(0) <= b <= Fr(1))


def phases_at(E, x):
    """The sorted multiset {frac(e_i x)} as Fractions in [0,1)."""
    ph = []
    for e in E:
        ph.append((int(e) * x) % 1)
    return sorted(ph)


def maxgap(phases):
    """Circular max gap of a sorted list of points in [0,1)."""
    if not phases:
        return Fr(1)
    if len(phases) == 1:
        return Fr(1)
    g = Fr(0)
    for a, b in zip(phases, phases[1:]):
        if b - a > g:
            g = b - a
    # wrap gap
    wrap = (phases[0] + 1) - phases[-1]
    if wrap > g:
        g = wrap
    return g


def sectors_hit(E, x):
    return set(sector((int(e) * x) % 1) for e in E)


# ------------------------------------------------ exact measures over a fine grid
def measure_objects(E):
    """Returns dict of exact-rational measures:
        p0          = meas{all 7 sectors hit}
        not_p0      = 1 - p0  (>=1 sector missed)
        gap_gt_2_7  = meas{maxgap > 2/7}   == g*(E)  (THM-527 good, ignoring G_P)
        gap_gt_1_7  = meas{maxgap > 1/7}
        sect2adj    = meas{two cyclically-adjacent sectors both missed}
       The grid: sector boundaries (j/(7e)) for the sector tests, which ALSO
       capture the >2/7 and >1/7 gap transitions because a gap edge moving past a
       1/7 multiple is exactly a sector-content change.  We additionally refine the
       gap tests on the cluster-phase coincidence grid to be safe."""
    E = [int(e) for e in E]
    # primary grid = sector boundaries
    grid = phase_breakpoints(E)
    # refine for the gap-threshold tests: also where a PAIR of phases is THRESH apart,
    # i.e. frac(e_i x) - frac(e_j x) = +-2/7 (mod 1).  Breakpoints at x = (m +- 2/7)/(e_i-e_j).
    refine = set(grid)
    for i in range(len(E)):
        for j in range(len(E)):
            d = E[i] - E[j]
            if d == 0:
                continue
            ad = abs(d)
            # frac(d x) = +-2/7  -> d x = m +- 2/7  -> x = (m +- 2/7)/d
            for m in range(0, ad + 1):
                for s in (THRESH, -THRESH, Fr(1, 7), -Fr(1, 7)):
                    val = Fr(m, ad) + s / ad
                    if Fr(0) <= val <= Fr(1):
                        refine.add(val)
    pts = sorted(refine)

    p0 = Fr(0)
    gap2 = Fr(0)
    gap1 = Fr(0)
    s2adj = Fr(0)
    for a, b in zip(pts, pts[1:]):
        if b <= a:
            continue
        w = b - a
        mid = (a + b) / 2
        ph = phases_at(E, mid)
        S = set(sector(p) for p in ph)
        if len(S) == P:
            p0 += w
        g = maxgap(ph)
        if g > THRESH:
            gap2 += w
        if g > Fr(1, 7):
            gap1 += w
        # two adjacent sectors both missed?
        missed = set(range(P)) - S
        adj = any(((s + 1) % P) in missed for s in missed)
        if adj:
            s2adj += w
    return {
        "p0": p0,
        "not_p0": Fr(1) - p0,
        "gap_gt_2_7": gap2,
        "gap_gt_1_7": gap1,
        "sect2adj": s2adj,
    }


def main():
    print("=" * 92)
    print("THREAD B: p0 (sector route) vs g*=rho*-without-G_P (THM-527).  Resolving 1/7 vs 2/7.")
    print("=" * 92)
    print("Objects (EXACT rational), all as meas over x in [0,1):")
    print("  p0          = all 7 sectors hit by {frac(e_i x)}            (sector route; p0<=cap => survive)")
    print("  1-p0        = >=1 sector missed")
    print("  g*=gap>2/7  = THM-527 'good' (maxgap{frac(e_i x)} > 2/7)    (g*>0 => M>=1/14)")
    print("  gap>1/7     = a >1/7 hole (the sector route's NATURAL gap)")
    print("  sect2adj    = two cyclically-adjacent sectors BOTH missed")
    print()
    print("KEY INCLUSIONS to verify:")
    print("  (A) {gap>2/7} subset {>=1 sector missed}  =>  g* <= 1-p0")
    print("  (B) {gap>2/7} subset {gap>1/7}            =>  g* <= gap>1/7")
    print("  (C) {gap>2/7} ?= {two adjacent sectors missed}  (the 2/7 <-> adjacency link)")
    print()

    battery = {
        "consec k=3 {0,1,2}":        [0, 1, 2],
        "consec k=4 {0..3}":         [0, 1, 2, 3],
        "consec k=5 {0..4}":         [0, 1, 2, 3, 4],
        "consec k=6 {0..5}":         [0, 1, 2, 3, 4, 5],
        "consec k=7 {0..6}":         [0, 1, 2, 3, 4, 5, 6],
        "consec k=8 {0..7}":         [0, 1, 2, 3, 4, 5, 6, 7],
        "consec k=9 {0..8}":         list(range(9)),
        "perforated k=7 {0..8}\\{1,7}": [0, 2, 3, 4, 5, 6, 8],   # part F min
        "Sidon k=5 {0,1,3,7,12}":    [0, 1, 3, 7, 12],
        "two-block k=6":             [0, 1, 2, 40, 41, 42],
        "wide k=5":                  [0, 5, 11, 18, 27],
    }

    hdr = f"{'cluster E':<30}{'p0':>9}{'1-p0':>9}{'g*=gap>2/7':>12}{'gap>1/7':>10}{'sect2adj':>10}"
    print(hdr)
    print("-" * len(hdr))
    all_A = True
    all_B = True
    all_C = True
    rows = []
    for name, E in battery.items():
        m = measure_objects(E)
        rows.append((name, E, m))
        A_ok = m["gap_gt_2_7"] <= m["not_p0"]
        B_ok = m["gap_gt_2_7"] <= m["gap_gt_1_7"]
        C_ok = m["gap_gt_2_7"] == m["sect2adj"]
        all_A &= A_ok
        all_B &= B_ok
        all_C &= C_ok
        print(f"{name:<30}{float(m['p0']):>9.4f}{float(m['not_p0']):>9.4f}"
              f"{float(m['gap_gt_2_7']):>12.4f}{float(m['gap_gt_1_7']):>10.4f}"
              f"{float(m['sect2adj']):>10.4f}")
    print("-" * len(hdr))
    print(f"\n(A) g* <= 1-p0 for ALL tested:        {all_A}")
    print(f"(B) g* <= gap>1/7 for ALL tested:     {all_B}")
    print(f"(C) g* == sect2adj for ALL tested:    {all_C}  (exact equality of the 2/7 event and adjacency)")

    # ---- quantify the slack 1-p0 - g* (how far the sector route over-shoots rho*) ----
    print("\n" + "=" * 92)
    print("SLACK: (1-p0) - g*  [how much the sector route's '1 missed sector' EXCEEDS the 2/7 gap]")
    print("=" * 92)
    print(f"{'cluster E':<30}{'1-p0':>10}{'g*':>10}{'slack=(1-p0)-g*':>18}{'g*/ (1-p0)':>12}")
    for name, E, m in rows:
        slack = m["not_p0"] - m["gap_gt_2_7"]
        ratio = float(m["gap_gt_2_7"] / m["not_p0"]) if m["not_p0"] > 0 else float("nan")
        print(f"{name:<30}{float(m['not_p0']):>10.4f}{float(m['gap_gt_2_7']):>10.4f}"
              f"{float(slack):>18.4f}{ratio:>12.3f}")

    print("""
READING:
  * If (A) holds universally: g* <= 1-p0, so a sector-route UPPER bound p0<=cap gives a
    LOWER bound g* ... NO: p0<=cap => 1-p0 >= 1-cap, which bounds 1-p0 BELOW; but g*<=1-p0
    is an UPPER bound on g*.  So p0<=cap does NOT directly lower-bound g*.  (The inequality
    points the WRONG way for transfer.)  This is the crux of THREAD B.
  * The sector route certifies via p0<=cap a DIFFERENT (1/7-scale) survival event; THM-527
    needs the 2/7-scale gap.  The slack column shows g* is STRICTLY less than 1-p0 (the 2/7
    event is rarer than the 1/7 'one missed sector' event).
  * (C) tests whether the 2/7 gap is EXACTLY the 'two adjacent sectors missed' event.
""")


if __name__ == "__main__":
    main()
