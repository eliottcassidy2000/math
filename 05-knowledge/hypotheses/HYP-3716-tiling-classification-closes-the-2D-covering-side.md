---
id: HYP-3716
title: PURSUING THE GEOMETRIC CLAIM -- the 2D covering side of the LRC->hexagonal bridge is CLOSED: the rhombic-lattice covering density is minimized EXACTLY at theta=60deg (rhombus = 2 EQUILATERAL triangles = the hexagonal lattice, delta=2pi/sqrt27=1.2092 = Kershner), beating theta=90deg (square = 2 right triangles, delta=pi/2=1.571); the convex-tiling classification (Reinhardt 1918 + Rao 2017: triangles, all quads, 3 hexagons, 15 pentagons, no >=7-gons) makes the LATTICE-covering search FINITE -- centrally-symmetric convex Voronoi cells = hexagon or rectangle, hexagon wins; aperiodic (Socolar-Taylor = decorated hexagon) gives no covering gain (Kershner is general, all coverings). So 'the optimal covering is hexagonal' is a THEOREM; the bridge gap reduces to the LRC->2D-hexagonal reduction alone
status: VERIFIED (angle scan: min density at 60deg = 2pi/sqrt27 = Kershner; square = pi/2) + the finite-classification framing (Reinhardt/Rao) + Kershner/Fejes Toth (general). The 2D covering optimality is settled (theorem); the LRC->2D reduction (does the zeta_6-line inherit it) is the remaining OPEN step.
source: klein-2026-06-29-S27
depends_on:
  - HYP-3715   # the covering-min is the zeta_6-line in the hexagonal lattice
  - HYP-3706   # the hexagonal wallpaper bridge
related:
  - HYP-3710   # the A1->A2 Coxeter-Catalan ladder (A2 = the triangular lattice; the tridiagonal/Jacobi lead)
  - HYP-3705   # the covering-min reframing
results:
  - 04-computation/hexagonal_covering_angle_optimization_klein.py
  - 05-knowledge/results/hexagonal_covering_angle_optimization_klein.out
---

# HYP-3716 — the convex-tiling classification closes the 2D covering side

## The geometric core (verified): density minimized at the equilateral triangle
A 2D LATTICE covering's fundamental cell is a rhombus = TWO triangles. Scanning the rhombus angle `theta`
(unit basis vectors), covering density `delta(theta) = pi R(theta)^2 / sin(theta)` (`R` = covering radius =
deep-hole / Delaunay circumradius):
```
theta:  50    55    60      65    70    75    80    85    90
delta: 1.247 1.210 1.2092  1.210 1.240 1.287 1.356 1.449 1.571
```
**MINIMUM at `theta = 60deg`**, `delta = 2pi/sqrt(27) = 1.20920` (Kershner) -- the rhombus is then TWO
EQUILATERAL triangles (the hexagonal / triangular lattice). The right-angle case `theta=90deg` (square = two
right 45-45-90 triangles) gives `delta = pi/2 = 1.571`, far worse. So the equilateral triangle wins: "the
quadrilateral at the right angle is two right triangles; the optimum is the rhombus = two EQUILATERAL
triangles."

## The convex-tiling classification makes the search FINITE
Convex polygons that monohedrally tile the plane (Reinhardt 1918 + Rao 2017, the field now CLOSED): ALL
triangles, ALL quadrilaterals, exactly **3 hexagon** types, exactly **15 pentagon** types; **none with >= 7
sides**. For a LATTICE covering the Voronoi cell is **centrally symmetric**, convex, `<= 6` sides => only the
**HEXAGON or the RECTANGLE** (parallelogram). The angle scan settles the two: hexagon (60deg) beats
rectangle/square (90deg). So among lattice coverings the optimum is the EQUILATERAL-TRIANGULAR (hexagonal)
lattice -- a **finite** check, won by the hexagon. The triangular lattice is the `A2` "god": the 3-fold
(`p6m`) structure, the ambient lattice of the LRC covering-min's `zeta_6`-line (HYP-3715).

## Aperiodic gives no gain
Kershner 1939 is general (the hexagonal lattice is the thinnest covering among ALL coverings, not just
lattices). The Socolar-Taylor aperiodic monotile is a DECORATED HEXAGON (matching rules force aperiodicity)
-- still hexagonal; aperiodicity gives no covering-density improvement. So the optimum stays hexagonal/`A2`
even off the lattice. (Aperiodic order is a red herring for covering DENSITY; it may matter for the apex-7
column, a 7-fold story, not this 3-fold covering column -- a lead.)

## What this closes, and the one remaining gap
- **CLOSED (theorem):** the optimal planar covering is the hexagonal lattice (Kershner/Fejes Toth, all
  coverings), and the convex-tiling classification makes the lattice case a finite check (hexagon vs
  rectangle), won by the equilateral-triangular lattice at `theta=60`. The 2D covering side of the bridge is
  settled.
- **OPEN (the only remaining step):** the **LRC -> 2D-hexagonal reduction** -- that the LRC covering-min
  (the `zeta_6`-line in `Z/Phi_6(n) = Z[zeta_6]/(n-zeta_6)`, HYP-3715) corresponds to / inherits the 2D
  hexagonal covering optimality. The 1D cyclic `zeta_6`-line lives in the 2D hexagonal lattice; the gap is
  the cyclic-line-covering-radius `<->` 2D-Kershner step.

## A lead: the tridiagonal/A2 link to the Catalan thread
"The triangular tiling is the god, just tridiagonalized." The `A2`/triangular lattice's 3-fold structure
is the tridiagonal (3-term / Jacobi) form; the Catalan family (HYP-3710) are the moments of the free Jacobi
operator, and the `A2` Coxeter-Catalan `1,5,42,462` counts the hexagonal-chamber walks. So the covering
(this tiling work) and the sequence thread (HYP-3710) meet at the tridiagonal Jacobi / `A2` operator -- a
lead for computing the `zeta_6`-line covering radius spectrally.

## Net
The geometric claim's 2D half is now a theorem: the optimal covering is the equilateral-triangular
(hexagonal) lattice, by Kershner + the finite convex-tiling classification (the angle optimum is exactly
60deg, two equilateral triangles), with aperiodic giving no gain. The LRC->hexagonal bridge is reduced to
the single remaining LRC->2D reduction (the `zeta_6`-line inheriting Kershner).
