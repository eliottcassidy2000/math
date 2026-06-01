---
id: HYP-2033
status: PARTIALLY-TRUE
source: oracle-2026-06-01-S541
related:
  - HYP-2003
  - HYP-2012
  - THM-381
---

# HYP-2033: LRC is a property of the affine-symmetric-group shadow of a torus braid; the annular braid group's CENTER is the speed-shift symmetry

**Setup.** Put the observer at the HOLE of an annulus; the n-1 runners are strands
winding around it. The annulus' radial coordinate is exactly `||v_i t||`; runner i
winds `v_i` times around the hole. The LRC system's image in the annular braid group
`Z_{n-1}` (= affine braid group `Ã_{n-2}`) is the homogeneous TORUS braid with
winding vector `v`.

**VERIFIED (`lrc_annular_braid_center_s541.py`):**
1. Crossing counts: `#σ_0` (cut/observer-wrap generator) `= Σ v_i` (the affine
   level); runner-runner crossings `= Σ_{i<j}|v_i-v_j|` (holdback/tension, S540o).
   The winding vector around the hole IS the speed vector.
2. The CENTER of the annular braid group is the full core rotation `ρ` = add a
   constant `c` to every speed. Loneliness depends only on pairwise differences,
   which `ρ` fixes -> every runner's clear-collar radius is UNCHANGED under
   `v -> v + c·1` (checked c=10,100). So the center acts TRIVIALLY on the LRC
   question: LRC lives in the annular braid group MOD CENTER = the affine symmetric
   group `S̃_{n-1}`. This is the braid reason LRC depends only on differences.
3. Per-runner loneliness (AP `{0,1,2,3,4}`): edge runners tight at exactly `1/n`,
   middle runners looser (the AP extremal, recovered).

**CLAIM (open):** "the hole has a fat collar (`≥1/n`) at some height" is a property
of the `S̃_{n-1}` image of `v` ALONE — an alcove-geometry condition on the affine
permutation, independent of the central core rotation. If so, LRC(n) is a pure
alcove-geometry statement about the affine symmetric group element `v`, tying the
annular-braid picture to the S525 alcove walk and the S532 independent-pair
(= affine root) channels.

**TO DO:** (a) make the surjection `Z_{n-1} -> S̃_{n-1}` and the affine permutation
of the LRC torus braid explicit; (b) express "fat collar" as a condition on the
alcove the affine permutation lands in; (c) check it is `ρ`-invariant (must be, by
2) and whether it is the SAME condition as the S526 covering / S532 channel-state.

**Files:** `04-computation/lrc_annular_braid_center_s541.py` (+.out). Reflection:
`07-reflections/lrc-annular-braid-the-center-is-the-shift-s541.md`. Builds on the
disk pure-braid `oracle-S540o`.
