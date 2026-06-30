---
id: HYP-3717
title: THE THREE-GAP THEOREM gives the covering-min DIRECTLY (and the >1/14 margin), and the TORUS LIFT shows the speeds are a 2-row hexagonal patch -- (A) CF(n/Phi_6(n))=[0;n-1,n]; (B) the densest core {1..n-2} is a clean 2-GAP AP (gaps n, ~2*Phi_6/(n-1)) with deep hole 1/(n-1); the minimal killer lcm(n-1,n) inserts ONE point splitting the big gap into two (THREE-gap config), deep hole drops to 14/183=n/Phi_6(n)=M; M>1/n because a single insertion can only split 29->1+28 (can't reach the <=26 that 1/14 needs; margin 13/2562); (C) Eisenstein-embedding the speeds {1..12,182} into Z[w]/(14-w) gives a 2-ROW hexagonal patch (b in {0,1}) -- the covering-min is genuinely 2D-hexagonal, and the 1D three gaps are the patch's view of the hexagonal directions
status: VERIFIED for n=14 (binding config gaps {1,14,28}, M=28/2/183=14/183; core {1..12} gaps {14,29}; killer splits 29->1+28; CF=[0;13,14]; speeds embed in 2 hexagonal rows), structurally general. The three-gap gives M and M>1/n for the densest-core+lcm-killer FAMILY (rigorous mechanism); the GLOBAL covering-min + the precise 1D<->2D Kershner metric bridge remain OPEN.
source: klein-2026-06-29-S28
depends_on:
  - HYP-3715   # the covering-min is the zeta_6-line; M=covering radius; densest core + killer
  - HYP-3716   # the 2D covering optimality is Kershner (the side this reduces toward)
related:
  - HYP-3706   # the hexagonal/Eisenstein bridge
  - HYP-+2913  # prior three-gap / Steinhaus tight-locus work
  - HYP-3551   # M = j/D; 14/183; densest core + minimal killer
results:
  - 04-computation/torus_lift_threegap_klein.py
  - 05-knowledge/results/torus_lift_threegap_klein.out
---

# HYP-3717 — the three-gap theorem gives the covering-min; the torus lift gives a hexagonal patch

## (A) The continued fraction
`CF(n/Phi_6(n)) = [0; n-1, n]` (since `Phi_6(n) = n^2-n+1 = (n-1)n + 1`). Verified n=3,5,7,14. The two partial
quotients are the core size `n-1` and `n`; the convergent denominators are `1, n-1, Phi_6(n)`.

## (B) The three-gap theorem gives the covering-min DIRECTLY
At the binding rotation `a = 14 = zeta_6`:
- **the densest core `{1..12}`** maps to `{14k : k=1..12}` -- a clean ARITHMETIC PROGRESSION (step 14), with
  exactly **two** gaps `{14 (x11), 29 (x1)}` (the small step and the one big wrap-gap); deep hole `= 29/2`,
  i.e. `M_core = 1/13 = 1/(n-1)` (the core's lonely gap);
- **the minimal killer `182 = lcm(13,14)`** maps to `169`, which lands INSIDE the big `29`-gap and SPLITS it
  `29 -> 1 + 28` -- now a **three-gap** config `{1, 14, 28}` (the three-distance theorem); the deep hole
  drops to `28/2 = 14`, so `M = 14/183 = n/Phi_6(n)`.
- **the margin `M > 1/n`**: a single inserted point can only split the big gap as `29 = 1 + 28`, so the deep
  hole is `>= 28/2 = 14`; to reach `M = 1/14` would need the big gap `<= 2*183/14 = 26.1`, unattainable with
  one killer. Hence `M = 14/183 > 1/14`, margin `13/2562`. **The three-gap theorem is the mechanism behind
  the covering floor's positivity for the densest-core + minimal-killer family.**

## (C) The torus lift: the speeds are a 2-row hexagonal patch
Eisenstein-embedding the speeds `{1..12,182} mod (14 - w)` (`w = zeta_6`, norm `a^2+ab+b^2`):
```
1..6, 182 -> (1,0),(2,0),...,(6,0),(-1,0)     [row b=0]
7..12     -> (-7,1),(-6,1),...,(-2,1)          [row b=1]
```
the speeds occupy **two adjacent rows** (`b in {0,1}`) of the hexagonal lattice -- a compact 2D hexagonal
patch, NOT a single line. So the covering-min is genuinely 2D-hexagonal; the 1D three gaps `{1,14,28}` are
the patch's projection onto the time-line (the line's view of the hexagonal nearest-neighbour directions,
the `A2` 3-fold). The torus lift (2D) and the three-gap (1D) are two views of the same hexagonal object.

## The reduction status (honest)
- **RIGOROUS (this family):** the three-gap theorem computes `M = n/Phi_6(n)` for the densest-core +
  `lcm(n-1,n)`-killer construction and proves `M > 1/n` (the single-killer split cannot shave the big gap
  enough). So the conjecture holds, with an explicit three-gap mechanism, on this family.
- **OPEN:** that `14/183` is the GLOBAL covering-min (no exotic, non-densest-core covering beats it); and the
  precise 1D-three-gap `<->` 2D-Kershner metric bridge (the lift is the nonlinear Eisenstein reduction, not
  a linear rescaling, so "the 1D deep hole = the 2D hexagonal covering radius" is structural, not yet a clean
  metric equality).

## Net
The three-gap theorem gives the covering-min and its `>1/n` margin directly and rigorously for the
densest-core+killer family (the killer splits the core's big gap `29 -> 1 + 28`, deep hole `28/2`,
`M = n/Phi_6(n)`); the torus lift shows the construction is a 2-row hexagonal patch, confirming the `A2`
ambient. The two views unify the LRC covering-min as a hexagonal object; the global optimality and the exact
1D<->2D bridge remain the open core.
