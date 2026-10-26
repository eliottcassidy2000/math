---
id: THM-4773
title: "Fixed-twist rational lifts on the hourglass boundary: a base point y of P^2 lifts to the twisted seven-root surface X_d over a field k (char k != 2) iff the nonzero values d_i L_i(y) share one squareclass, with 2^(6-r) lifts; over an infinite field, open central-line E = 0 lifts exist iff the three parallel-pair port characters are trivial on the twist (-d_1 d_2, -d_3 d_4, -d_5 d_6 squares) and the conic w^2 = d_1 d_5 d_3^2 s^2 + d_3 d_5 t^2 has a k-point, which over Q is decided by Hilbert symbols; the Pythagorean twist (1,1,-1,1,-1,1,-1) makes this conic the unit circle and gives the rational point (0:3:3:4:4:5:5) from (3,4,5); on each outer line the lift locus is a twisted genus-five curve pair, hence finite over a number field; over Q the untwisted surface has no boundary lifts. Over number fields the dichotomy is: central line none or infinitely many, outer lines finitely many."
status: >
  PROVED (klein, 2026-10-25) for the pointwise criterion and the central-line theorem over infinite fields;
  finite fields require an actual nondegenerate conic point, as the F3 correction demonstrates; CITED for the genus-five
  components (six_seven_gale_cube_20261025 GENUS_FIVE_ELLIPTIC_FACTORS / BRANCH_CUBE_AND_DOUBLE_COVER; Session60)
  and Faltings' theorem. FINITE-EXACT: Hilbert-symbol decision agrees with brute-force existence (height <= 25) on
  300 random twists and eight displayed ones; untwisted surface has no lifts on E = 0 up to height 60; Pythagorean
  twist lifts exactly 4 of 6521 base points of height <= 12. Not independently audited.
session: klein-2026-10-25 (codex/session-tournament-twist-unify)
source: 05-knowledge/results/fixed_twists_conductors_reverse_potential_20261025/twists/FIXED_TWIST_LIFTS.md
scripts:
  - 05-knowledge/results/fixed_twists_conductors_reverse_potential_20261025/twists/check_central_conic.py (+ conic-normal.txt / conic-optimized.txt)
  - 05-knowledge/results/fixed_twists_conductors_reverse_potential_20261025/twists/check_boundary_lifts.py (+ lifts-normal.txt / lifts-optimized.txt)
related:
  - 05-knowledge/results/normseven_boundary_20261025/BOUNDARY.md (section 5: residual torsor; 5a: untwisted emptiness)
  - 01-canon/theorems/THM-4765-quadratic-twists-of-a-square-root-cover-average-to-the-branch-complement.md
  - 05-knowledge/results/session63_integral_observers_root_arrangements_20261025/twists/RESULT.md (fixed versus fresh twists)
  - 05-knowledge/results/session57_square_sum_magic_paths_20261025/SYNTHESIS.md (square-sum theme)
  - 01-canon/theorems/THM-4772-twisted-root-algebras-over-normal-bases-square-root-ideals-mixed-reduced-roots-port-classes-and-the-conductor-valuation-matrix.md
---

# THM-4773 — Fixed-twist rational lifts on the hourglass boundary

## Statement and proof

Forms `L = (E, E+x, E-x, E+z, E-z, E+x+z, E-x-z)`, twist `d`, surface `X_d : s_i^2 = d_i L_i`.

1. **Pointwise.** `y in P^2(k)` lifts iff the nonzero `d_i L_i(y)` lie in one class of `k^x/k^{x2}`;
   then there are `2^{6-r}` lifts, `r` the number of vanishing forms (the fibre count in the proof of
   THM-4765).
2. **Central line (infinite k for the nonemptiness equivalence).** For `(0:x:z)` with `xz(x+z) != 0` the six values are
   `d_1 x, -d_2 x, d_3 z, -d_4 z, d_5(x+z), -d_6(x+z)`; one squareclass forces `-d_1d_2, -d_3d_4, -d_5d_6`
   squares and `d_1d_3 xz`, `d_1d_5 x(x+z)` squares; `x = d_1d_3 s^2`, `z = t^2` reduces the last to the
   conic `w^2 = d_1d_5d_3^2 s^2 + d_3d_5 t^2`. A conic with one rational point has infinitely many, so
   the degenerate points are avoidable. Over `Q` the conic's solubility is Legendre's criterion via
   Hilbert symbols.
   Over finite k require a point avoiding st(d_1d_3 s^2+t^2)=0. The Pythagorean conic
   over F3 has four points, all excluded; the old whole-conic criterion was too broad.
3. **Examples over Q.** `d = 1`: `-1` is not a square, no lifts. Pythagorean `d`: unit circle, the point
   `(0:3:3:4:4:5:5)` satisfies all quadric relations. `d = (1,1,-1,1,-1,7,-7)`: `w^2 = 7s^2 + 7t^2`,
   `(7,7)_7 = -1`, no lifts.
4. **Outer lines and points over number fields.** The lift locus over an outer line is the `d`-twist of the two
   genus-five branch components (CITED), finite over number fields by Faltings (CITED); the fifteen
   special points are finite squareclass tests (`(0:0:1)` lifts iff `d_3 ~ -d_4 ~ d_5 ~ -d_6`).

## Mechanism and boundary

The three pair conditions are the parallel-pair Walsh characters of THM-4765 evaluated on the twist
(ports preserved); the conic carries the two remaining even characters on `E = 0`. Nothing is claimed
about interior rational points, finite-field counts beyond THM-4765, or integer magic squares.
