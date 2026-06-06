#!/usr/bin/env python3
"""
S637 / HYP-2315 — The chromatic number of the plane: what is genuinely provable.

Three honest computations:
  (A) The upper bound 7 IS the Eisenstein prime above 7: 7 = N(3+w), and the
      hexagonal 7-coloring = reduction mod (3+w). The 6 units map onto F_7*.
  (B) The 7-coloring's robustness: one fixed hexagonal 7-coloring simultaneously
      forbids ALL monochromatic distances in an interval; compute its ratio.
  (C) The density barrier: why the simple measurable-density bound chi_m >= 1/m1
      can prove chi >= 5 but PROVABLY can never reach 6.

No external libs. Pure arithmetic / exact where possible.
"""

import math

print("=" * 70)
print("(A) The upper-bound 7 = the Eisenstein prime above 7")
print("=" * 70)
# Eisenstein integers a + b*w, w = e^{2pi i/3}, w^2 + w + 1 = 0.
# Norm N(a+bw) = a^2 - a*b + b^2.
def norm(a, b):
    return a*a - a*b + b*b

print(f"  N(3 + w) = 3^2 - 3*1 + 1^2 = {norm(3,1)}   (the prime above 7)")
print(f"  7 splits in Z[w] iff 7 = 1 mod 3:  7 mod 3 = {7 % 3}  -> splits")
print(f"  smallest Eisenstein norms (Loeschian numbers) up to 13:")
loesch = sorted({norm(a,b) for a in range(-5,6) for b in range(-5,6)} - {0})
print("   ", [n for n in loesch if n <= 13])
print("    -> 3 (ramified, triangular-lattice 3-coloring), then 4, then 7.")
print("    7 is the first SPLIT prime; its quotient Z[w]/(3+w) = F_7.")

# In F_7, w = -3 = 4. The 6 units {+-1,+-w,+-w^2} mod 7:
w = 4
units = sorted({1 % 7, (-1) % 7, w % 7, (-w) % 7, (w*w) % 7, (-(w*w)) % 7})
print(f"\n  w == -3 == {w} (mod 7);  check w^2+w+1 = {(w*w + w + 1) % 7} (mod 7)")
print(f"  6 Eisenstein units mod 7 = {units}  == F_7* = {list(range(1,7))}")
print(f"  closed hexagon nbhd (center 0 + 6 units) = {[0]+units} == all of F_7")
print(f"  evaluation point 2: 2^2+2+1 = {(4+2+1)%7} (mod 7) -> 2 is ALSO a cube root")
print(f"  the two primitive cube roots of unity mod 7 are 2 and 4 = 2^2; 2*4={2*4%7}")

print("\n" + "=" * 70)
print("(B) Robustness of ONE hexagonal 7-coloring (forbidden-distance interval)")
print("=" * 70)
# Hexagon circumradius R (center->vertex). Centers form a triangular lattice with
# nearest-center spacing a = sqrt(3)*R (across-flats width). The 7-coloring's
# same-color centers are the index-7 Eisenstein sublattice, min distance a*sqrt(7).
# Distances REALIZED inside a color class: within one hexagon (0 .. 2R), and
# between nearest like-colored hexagons (>= closest approach). So distances in the
# open gap are FORBIDDEN monochromatically -> safe to be "the unit distance".
#   in-hexagon max diameter:        2R
#   nearest like-color center dist:  sqrt(7) * a = sqrt(7)*sqrt(3)*R = sqrt(21) R
#   closest approach of two same-color hexagons (translates): center - 2R
# Forbidden interval = ( 2R , (sqrt(21) - 2) R ).
import fractions
R = 1.0
lo = 2 * R
hi = (math.sqrt(21) - 2) * R
print(f"  hexagon circumradius R = {R}")
print(f"  forbidden monochromatic-distance interval: ( {lo:.4f} , {hi:.4f} )")
print(f"  robustness ratio (hi/lo) = (sqrt(21)-2)/2 = {(math.sqrt(21)-2)/2:.5f}")
print(f"  -> setting unit distance =1 works for R in (1/{hi:.4f}, 1/{lo:.4f})")
print(f"     = ( {1/hi:.4f} , {1/lo:.4f} )   (diameter 2R in ({2/hi:.4f}, {2/lo:.4f}))")
print("  So a single 7-coloring kills a whole BAND of distances, ratio ~1.29:")
print("  this slack is exactly why 7 is not known to be tight (room to merge colors).")

print("\n" + "=" * 70)
print("(C) The density barrier: why 1/m1 caps at 5, never 6")
print("=" * 70)
# m1 = sup density of a measurable 1-avoiding (independent) set in R^2.
# A measurable k-coloring: some class has density >= 1/k, and is 1-avoiding,
# so 1/k <= m1, i.e. k >= 1/m1.  This is the LRC p0 / 1-avoiding density bound.
m1_lo  = 0.2293   # Croft: a 1-avoiding set with this density EXISTS (lower bound on m1)
m1_cl  = 0.2598   # classical upper bound on m1 (fleet HYP-2278 uses this)
m1_rec = 0.2470   # recent LP upper bound (Ambrus et al.)
print(f"  known: {m1_lo} <= m1 <= {m1_cl}  (classical) / <= {m1_rec} (recent LP)")
print(f"  density bound chi >= ceil(1/m1):")
print(f"    1/m1 ranges in [1/{m1_cl}, 1/{m1_lo}] = [{1/m1_cl:.3f}, {1/m1_lo:.3f}]")
print(f"    reaching 5 is THRESHOLD-DEPENDENT: 1/{m1_cl}={1/m1_cl:.2f}<4 (only >=4);")
print(f"                                       1/{m1_rec}={1/m1_rec:.2f}>4 (gives >=5).")
print("    (the PROVEN chi_m>=5 is Falconer 1981, a refined, not single-class, argument)")
need = 1/6
print(f"  THE ROBUST BARRIER (m1-bound-INDEPENDENT): chi >= 6 is IMPOSSIBLE here.")
print(f"    need m1 < 1/6 = {need:.4f}, but the Croft set has density {m1_lo} > {need:.4f}.")
print(f"    so 1/m1 <= {1/m1_lo:.3f} < 6 ALWAYS, regardless of the m1 upper bound.")
print("  => the single-class density (= LRC lonely-density) bound is PROVABLY capped")
print("     below 6; the {5,6,7} distinction is irreducibly combinatorial (de Grey")
print("     graphs / multi-class LP), not single-class density. = the LRC Vitali wall.")
print("\n  This is the plane mirror of LRC(14): the first-moment/density bound is")
print("  vacuous past a point; the content lives in arc CORRELATIONS (HYP-2195).")
