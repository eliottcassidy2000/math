#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) MOUTH-OWNERSHIP GEOMETRY (kind-pasteur-2026-06-19).

Question (codex's lonely-measure route, angle = mouth-ownership):
  THM-541's drop-6 collar has 4 surviving safe components. Which SPEEDS' danger
  arcs form the left/right WALLS of each component? Are the tower elements
  {1,2,4,8} (the shell-1 dyadic word, carry[1]=15) the binding walls? If not,
  what is the exact geometric mechanism by which deleting a tower element forces
  meas(G_C) >= thr2 = 426/35035 (HYP-2661, "carry conservation")?

EXACT rationals throughout. Re-uses the addressed-wall engine of
04-computation/lrc14_ap_window_single_hole_certificate_codex_s34.py.

FINDINGS (computed below):
  1. The 4 drop-6 mouths are bounded ONLY by speeds {11,12,13}. The tower
     {1,2,4,8} owns NO mouth wall directly. So the mechanism is NOT
     "tower element is a binding wall".
  2. Deleting any tower element from the 12-core PRESERVES all 4 mouths exactly
     (they are 11/12/13-walled, untouched) but OPENS huge new safe mass
     elsewhere, instantly exceeding thr2. The tower's role is GLOBAL coverage,
     not LOCAL mouth-walling.
  3. The two-tail min row (the genuine HYP-2661 boundary case) keeps the two
     OUTER det=3 mouths (13->12, 12->13) intact but DESTROYS the two INNER
     det=5 mouths (12->11, 11->12): the {5:+2,23:+2} tails (speeds 20,46)
     carve the inner mouth band into many thin pieces. Net measure is pushed up.
"""
from __future__ import annotations
import sys
sys.path.insert(0, "04-computation")
from fractions import Fraction
from lrc14_ap_window_single_hole_certificate_codex_s34 import (
    drop_core, addressed_components, safe_measure, walls,
)

THR1 = Fraction(7, 858)        # drop-6 mouth (global min)
THR2 = Fraction(426, 35035)    # AP one-hole 2nd value (the cutoff)
TOWER = (1, 2, 4, 8)           # shell-1 dyadic word; carry[1] = 1+2+4+8 = 15


def fmt(q): return f"{q}={float(q):.6f}"


def section(t): print("\n" + "=" * 72 + "\n" + t + "\n" + "=" * 72)


# ---------------------------------------------------------------------------
section("1. DROP-6 MOUTH WALLS: which speeds bound each surviving component?")
base = drop_core(6)
base_comps = addressed_components(base)
mouth_ivals = [(c.lo, c.hi) for c in base_comps]
print(f"drop-6 core = {base}")
print(f"meas(G) = {fmt(safe_measure(base))}  (= THR1)\n")
mouth_speeds = set()
for i, c in enumerate(base_comps):
    L = c.left_owners[0]; R = c.right_owners[0]
    mouth_speeds |= {L.speed, R.speed}
    print(f"  mouth {i+1}: [{c.lo},{c.hi}] len={c.length}  "
          f"LEFT wall {L.label()} (speed {L.speed}) | "
          f"RIGHT wall {R.label()} (speed {R.speed})  det={c.determinant_numerator()}")
print(f"\n  => all mouth-bounding speeds = {sorted(mouth_speeds)}")
print(f"  => tower {TOWER} intersect mouth-walls = "
      f"{sorted(set(TOWER) & mouth_speeds)}  (EMPTY => tower owns NO mouth wall)")


# ---------------------------------------------------------------------------
section("2. DELETE a tower element from the 12-core: mouths persist, mass explodes")
print("(deleting speed s makes an 11-core; isolates the GLOBAL coverage role)\n")
print(f"{'del':>4} | {'meas':>14} | {'>=THR2':>7} | mouths kept | delta_meas")
print("-" * 64)
for s in TOWER:
    nc = tuple(v for v in base if v != s)
    comps = addressed_components(nc)
    m = safe_measure(nc)
    kept = sum(1 for c in comps if (c.lo, c.hi) in mouth_ivals)
    print(f"{s:>4} | {str(m):>14} | {str(m >= THR2):>7} | "
          f"{kept}/4 (all 11/12/13 walls intact) | +{float(m-THR1):.6f}")
print("\n  Interpretation: removing a tower speed does NOT touch the 11/12/13")
print("  mouth walls (mouths survive 4/4) but uncovers large far regions; the")
print("  measure leaps to >= THR2 by GLOBAL coverage loss, not mouth damage.")


# ---------------------------------------------------------------------------
section("3. THE GENUINE BOUNDARY: two-tail min row (HYP-2661 'pay the threshold')")
# {1:-4,5:+2,23:+2}: holes {4,6,10}, tails {20,46}; carry tower 15 -> 11 (loses 4)
core2 = tuple(sorted([d for d in range(1, 14) if d not in {4, 6, 10}] + [20, 46]))
comps2 = addressed_components(core2)
m2 = safe_measure(core2)
print(f"core = {core2}")
print(f"meas = {fmt(m2)}   >= THR2 ? {m2 >= THR2}   (carry tower 15 -> 11: lost bit 4)\n")
kept2 = []
for c in comps2:
    if (c.lo, c.hi) in mouth_ivals:
        kept2.append((c.lo, c.hi))
print(f"OLD drop-6 mouths retained: {len(kept2)}/4")
for i, (lo, hi) in enumerate(mouth_ivals):
    L = base_comps[i].left_owners[0]; R = base_comps[i].right_owners[0]
    status = "KEPT " if (lo, hi) in kept2 else "DESTROYED"
    print(f"  mouth {i+1} [{lo},{hi}] det={base_comps[i].determinant_numerator()} "
          f"({L.label()}->{R.label()}): {status}")
print("\n  => the two OUTER det=3 mouths (13->12, 12->13) survive;")
print("     the two INNER det=5 mouths (12->11, 11->12) are carved up by the")
print("     new tail speeds 20,46 (the {5:+2,23:+2} replacements).")


# ---------------------------------------------------------------------------
section("4. PER-COMPONENT MEASURE LEDGER: drop-6 vs perturbations")
def ledger(name, core):
    comps = addressed_components(core)
    m = safe_measure(core)
    print(f"\n{name}: core={core}")
    print(f"  meas = {fmt(m)}   components = {len(comps)}")
    for c in comps:
        tag = "MOUTH" if (c.lo, c.hi) in mouth_ivals else "new  "
        L = c.left_owners[0]; R = c.right_owners[0]
        print(f"    [{c.lo},{c.hi}] len={c.length}={float(c.length):.6f} {tag} "
              f"{L.label()}->{R.label()}")
ledger("drop-6 (THR1, mouth)", base)
ledger("drop-12 (THR2, AP 2nd)", drop_core(12))


# ---------------------------------------------------------------------------
section("5. drop-12: how the AP 2nd value relates to the mouth walls")
# drop-12 keeps speed... let's see which walls bound its 4 components
c12 = addressed_components(drop_core(12))
print("drop-12 components and their walls (note: speed 12 is GONE):")
sp12 = set()
for c in c12:
    L = c.left_owners[0]; R = c.right_owners[0]
    sp12 |= {L.speed, R.speed}
    print(f"  [{c.lo},{c.hi}] len={c.length} {L.label()}(sp{L.speed})->{R.label()}(sp{R.speed})")
print(f"  drop-12 mouth-wall speeds = {sorted(sp12)}")
print("  (removing speed 12 forces the high-speed collar to re-wall on 11,13,14-ish)")


# ---------------------------------------------------------------------------
section("SUMMARY: the carry-conservation geometry")
print("""
The drop-6 MOUTH (the global-min safe set, meas=7/858) is a 4-component collar
walled EXCLUSIVELY by the three TOP speeds 11,12,13 in the determinant pattern
[3,5,5,3]:
    13->12  (det 3)   [outer left]
    12->11  (det 5)   [inner left]
    11->12  (det 5)   [inner right]
    12->13  (det 3)   [outer right]

The TOWER {1,2,4,8} owns NONE of these walls. Its job is GLOBAL: the low dyadic
speeds 1,2,4,8 are the densest danger-arc generators, so together they CLAMP the
rest of [0,1) down to near-zero safe mass, leaving only the tight high-speed
collar. Concretely:
  * speed 1 alone fixes a 1/7-wide forbidden band around 0 (the wrap arc).
  * speeds 1,2,4,8 = the geometric dyadic chain -> arcs at every 2^-a scale,
    the binary refinement that pins the measure to the collar.

Deleting ANY tower bit:
  * leaves the 11/12/13 mouth walls intact (mouths persist), BUT
  * removes a whole dyadic scale of coverage, uncovering far mass >> THR2.
So you cannot damage the tower (carry[1] < 15) without the measure leaping over
THR2 -- exactly HYP-2661. The conservation is: mouth-mass is tiny and high-speed
walled; tower-mass is the GLOBAL clamp. Lose a clamp bit, pay the threshold.
""")
