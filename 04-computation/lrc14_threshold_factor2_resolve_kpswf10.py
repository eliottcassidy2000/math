#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_threshold_factor2_resolve_kpswf10.py   (kind-pasteur 2026-06-21, THREAD B Q3/Q4)

RESOLVE THE FACTOR-2.  THM-527 needs maxgap{frac(e_i x)} > 2/7.  The sector route
(p0) is a 1/7-scale object.  Three candidate "sector-style" events at the 2/7 scale,
all EXACT-rational, to see which (if any) EQUALS the 2/7-gap event g*:

  EV0 := {>=1 of the 7 (1/7-)sectors entirely missed}           = 1-p0       (1/7 scale)
  EV1 := {some 1/7-sector lies ENTIRELY inside a phase-gap}      (a fully-empty 1/7 cell)
  EV2 := {some 2/7-LONG arc [a, a+2/7) is entirely phase-free}   (a fully-empty 2/7 window)
  EV3 := {two CYCLICALLY-ADJACENT empty 1/7-sectors}             (adjacency)
  G2  := {maxgap{frac(e_i x)} > 2/7}                              = g* (THM-527 good)
  G1  := {maxgap{frac(e_i x)} > 1/7}

PROVABLE CONTAINMENTS (test EXACTLY that they hold as measures, and find which are EQ):
  maxgap > 2/7  <=>  exists open arc of length > 2/7 that is phase-free
                 ==>  exists a closed 2/7 arc phase-free  (EV2),  and conversely up to bdry
  maxgap > 2/7  ==>  >=1 full 1/7-sector empty (EV1)  [an open >2/7 interval contains a full 1/7 cell]
  But EV1 (one empty cell) only needs maxgap > 1/7-ish, NOT 2/7.  So EV1 != G2 in general.

GOAL: determine whether g*=G2 equals a CLEAN sector-style event at the 2/7 scale (EV2),
so a "2/7 sector route" p0' := 1 - meas(EV2) would have  p0' <= cap'  <=>  g* >= 1-cap' > 0.
THAT would let a (rescaled) sector route CLOSE THM-527.  We test EV2 == G2 exactly.

Also Q4: the '3.5-sector' / 7-window-at-2/7-scale cover idea: meas{the phases hit a
2/7-net} -- we compute the natural dual and compare to g*.
"""
import itertools
from fractions import Fraction as Fr

P = 7
T2 = Fr(2, 7)
T1 = Fr(1, 7)


def phases_at(E, x):
    return sorted((int(e) * x) % 1 for e in E)


def maxgap(phases):
    if len(phases) <= 1:
        return Fr(1)
    g = max(b - a for a, b in zip(phases, phases[1:]))
    wrap = (phases[0] + 1) - phases[-1]
    return max(g, wrap)


def empty_cell_exists(phases, cell):
    """Does some [j*cell,(j+1)*cell) (j=0..1/cell-1) contain NO phase? cell divides 1."""
    ncell = int(1 / cell)
    occupied = set()
    for p in phases:
        occupied.add(int(p / cell))   # which cell (floor)
    return len(occupied) < ncell


def empty_arc_of_length(phases, L):
    """Does there exist an arc [a,a+L) (a anywhere) with NO phase strictly inside?
       Equivalent to: maxgap (circular) >= L when L<=maxgap, but as a *closed* empty
       window of EXACT length L it exists iff maxgap >= L.  (We test >= L and > L.)"""
    g = maxgap(phases)
    return g >= L, g > L


def gen_grid(E):
    E = [int(e) for e in E]
    bp = {Fr(0), Fr(1)}
    for e in E:
        if e == 0:
            continue
        # sector (1/7) boundaries
        for t in range(0, P * e + 1):
            bp.add(Fr(t, P * e))
    # pairwise threshold crossings at +-1/7 and +-2/7
    for i in range(len(E)):
        for j in range(len(E)):
            d = E[i] - E[j]
            if d == 0:
                continue
            ad = abs(d)
            for m in range(0, ad + 1):
                for s in (T1, -T1, T2, -T2):
                    val = Fr(m, ad) + s / ad
                    if Fr(0) <= val <= Fr(1):
                        bp.add(val)
    return sorted(b for b in bp if Fr(0) <= b <= Fr(1))


def compute(E):
    E = [int(e) for e in E]
    pts = gen_grid(E)
    M = {k: Fr(0) for k in
         ["p0", "EV1_emptycell", "EV2_emptyarc2_7_ge", "EV2_emptyarc2_7_gt",
          "EV3_adj", "G2", "G1"]}
    for a, b in zip(pts, pts[1:]):
        if b <= a:
            continue
        w = b - a
        mid = (a + b) / 2
        ph = phases_at(E, mid)
        S = set(int(P * p) for p in ph)
        if len(S) == P:
            M["p0"] += w
        # EV1: a fully empty 1/7 sector
        if empty_cell_exists(ph, T1):
            M["EV1_emptycell"] += w
        # EV2: an empty closed 2/7 arc (>= and >)
        ge, gt = empty_arc_of_length(ph, T2)
        if ge:
            M["EV2_emptyarc2_7_ge"] += w
        if gt:
            M["EV2_emptyarc2_7_gt"] += w
        # EV3: two adjacent empty sectors
        missed = set(range(P)) - S
        if any(((s + 1) % P) in missed for s in missed):
            M["EV3_adj"] += w
        g = maxgap(ph)
        if g > T2:
            M["G2"] += w
        if g > T1:
            M["G1"] += w
    return M


def main():
    print("=" * 100)
    print("FACTOR-2 RESOLUTION: which clean 'sector-style' event equals g* = {maxgap > 2/7} (THM-527)?")
    print("=" * 100)
    print("EV1_emptycell = a fully empty 1/7-sector        (>= 1/7 hole, the route's natural scale)")
    print("EV2_arc>2/7   = an empty OPEN 2/7 arc = {maxgap > 2/7} EXACTLY (def of circular maxgap)")
    print("EV3_adj       = two cyclically-adjacent empty 1/7-sectors")
    print("G2            = maxgap > 2/7  (THM-527 good).   G1 = maxgap > 1/7.")
    print()
    print("EXPECTED IDENTITIES:")
    print("  G2 == EV2_emptyarc2_7_gt   (tautology: maxgap>2/7 <=> an open empty 2/7 arc exists)")
    print("  G2  => EV1_emptycell       (open >2/7 interval contains a full closed 1/7 cell)")
    print("  G1 == EV1-ish but at 1/7   (one missed sector ~ a >1/7 hole)")
    print()

    battery = {
        "consec k=4":  [0, 1, 2, 3],
        "consec k=5":  [0, 1, 2, 3, 4],
        "consec k=6":  list(range(6)),
        "consec k=7":  list(range(7)),
        "consec k=8":  list(range(8)),
        "consec k=9":  list(range(9)),
        "consec k=11": list(range(11)),
        "consec k=13": list(range(13)),
        "perf k=7 {0..8}\\{1,7}": [0, 2, 3, 4, 5, 6, 8],
        "Sidon k=5":   [0, 1, 3, 7, 12],
    }
    hdr = (f"{'E':<22}{'p0':>8}{'EV1cell':>9}{'EV2>2/7':>9}{'EV3adj':>9}"
           f"{'G2':>8}{'G1':>8}{'G2==EV2?':>10}{'G2<=EV1?':>10}")
    print(hdr)
    print("-" * len(hdr))
    G2_eq_EV2 = True
    G2_sub_EV1 = True
    G2_eq_G1minus = []
    for name, E in battery.items():
        M = compute(E)
        eqEV2 = (M["G2"] == M["EV2_emptyarc2_7_gt"])
        subEV1 = (M["G2"] <= M["EV1_emptycell"])
        G2_eq_EV2 &= eqEV2
        G2_sub_EV1 &= subEV1
        print(f"{name:<22}{float(M['p0']):>8.4f}{float(M['EV1_emptycell']):>9.4f}"
              f"{float(M['EV2_emptyarc2_7_gt']):>9.4f}{float(M['EV3_adj']):>9.4f}"
              f"{float(M['G2']):>8.4f}{float(M['G1']):>8.4f}{str(eqEV2):>10}{str(subEV1):>10}")
    print("-" * len(hdr))
    print(f"\nG2 == EV2(open 2/7 empty arc) for ALL:  {G2_eq_EV2}   (tautological identity, sanity)")
    print(f"G2 <= EV1(one empty 1/7 cell) for ALL:  {G2_sub_EV1}   (the 2/7 event is RARER than 1 empty cell)")

    print("""
CONCLUSION on the factor-2:
  * g* = G2 = {an OPEN empty 2/7 arc exists}.  This is a 2/7-scale event, NOT the
    1/7-scale event EV1 (one empty sector).  EV1 (>= one empty 1/7 cell) is the dual of
    the sector-route p0 only at the 1/7 scale: 1-p0 = meas{>=1 sector missed}.  But
    'sector missed' (no phase in a 1/7 cell) is EV1, not G2.  So 1-p0 = EV1 >= G2, strictly.
  * THEREFORE: a sector route AT THE 1/7 SCALE controls 1-p0 = EV1, which OVER-counts g*.
    A route that EQUALS g* must be at the 2/7 scale (EV2).  The natural '2/7 sector route'
    is p0' = 1 - meas(EV2) = 1 - g*, and 'all 2/7-windows hit' = NOT g*.  So closing g*>0
    is the SAME as a 2/7-scale 'not-all-windows-covered' statement -- a DIFFERENT cap than
    the 1/7 cap_k.  The 1/7 sector route (p0<=cap_k) and THM-527 (g*>0) are PARALLEL: same
    phi-phase, DIFFERENT threshold (1/7 vs 2/7), and the 1/7 route's inequality direction
    does not transfer to the 2/7 object.
""")


if __name__ == "__main__":
    main()
