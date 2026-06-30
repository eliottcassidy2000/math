---
id: HYP-3703
title: The LRC covering-min IS a TILING/COVERING-OPTIMALITY problem with TWO REGIMES handing off at the 6/7 boundary -- small n (n<=6) is a PROJECTIVE-PLANE/Steiner-design covering (the drop-2 family 2/(2n-1), works WHILE PG(2,n-1) exists), large n (n>=7) is the HEXAGONAL/Kershner covering (the construction n/Phi_6(n), Eisenstein/A2 lattice), the transition at n=7 coinciding with the FIRST projective-plane FAILURE PG(2,6) (Bruck-Ryser; opus-S1) -- and klein-S24 bridges the two via the wallpaper group p6m (Steiner-design-optimal <-> Kershner-continuous-optimal). The 6/7 BOUNDARY recurs across the whole project: convex polygons tile the plane only with <=6 sides (Reinhardt 1918 + Rao 2017: triangles, quads, 3 hexagons, 15 pentagons, NONE >=7), the LRC apex is 7=Phi_3(2)=|Fano|, the forbidden H-values are {7,21}, the covering tile is a HEXAGON (the max convex tile). 'The triangular tiling is a god, tridiagonalized' = the A2/Eisenstein Cartan [[2,-1],[-1,2]] (det 3 = disc Q(sqrt-3)) is the square Z^2 [[2,0],[0,2]] (det 4, 'right-angle quadrilateral') with the off-diagonal -1 coupling; the LRC lives on A2 (coupled/hexagonal), the tournament staircase is A2 sheared. The hard core reduces to the DISCRETE cyclic-Kershner problem (opus: construction three-distance gaps EXACTLY {1,n,2n}); Kershner 1939 (continuous hexagonal covering, PROVEN) is its twin
status: SYNTHESIS + verified pieces (A2 Cartan tridiagonal det 3; construction three-distance gaps {1,n,2n}, opus; PG(2,n-1) transition at q=6; Kershner density 2pi/sqrt27). Connects the owner's tiling pointers (Reinhardt-Rao closure, Socolar-Taylor, triangular-tridiagonalized, quad=2-triangles) to the covering-min. The two-regime + 6/7-boundary frame is the new synthesis; the proof (discrete cyclic-Kershner) remains the open node (opus).
source: mac-mini-2026-06-30-S44
related:
  - HYP-3702  # the covering-set taxonomy (this gives the hard core its tiling identity)
  - HYP-3706  # klein-S24: hexagonal p6m, Singer multiplier = 3-fold, design<->Kershner bridge
  - HYP-3701  # the covering-min n-dependent transition (here = the projective->hexagonal handoff)
  - HYP-3616  # the 7,21 / apex-7=Phi_3(2)=Fano (the 6/7 boundary on the counting side)
  - HYP-3700  # the gap-edge (the apex doublet); this is the covering/tiling-edge
---

# HYP-3703 -- the covering-min is a tiling-optimality; the 6/7 boundary

The owner pointed at tiling theory (the resolved convex-polygon-tiling classification; Socolar-Taylor
aperiodic monotile; the triangular tiling 'tridiagonalized'; a quadrilateral = two triangles). These are
exactly the right frame for the LRC covering-min hard core.

## The covering-min is a covering-optimality, in two regimes
The LRC danger combs must cover `[0,1)` (a disproof = no gap `>= 2/n`); minimizing the gap `M` over covering
sets is a COVERING-OPTIMALITY problem. It has TWO regimes, handing off at the 6/7 boundary:
- **n <= 6: a PROJECTIVE-PLANE / Steiner-design covering.** The covering-min is the drop-2 family
  `2/(2n-1)` (HYP-3701), a difference-set / projective structure. It works WHILE `PG(2,n-1)` exists
  (`n-1 <= 5`, all prime powers).
- **n >= 7: the HEXAGONAL / Kershner covering.** The covering-min is the construction `n/Phi_6(n)`
  (Eisenstein/A2 lattice), whose speed-residues `mod Phi_6(n)` are an arithmetic progression with
  three-distance gaps EXACTLY `{1, n, 2n}` (opus-S1) -- the hexagonal fundamental domain.
- **The transition at `n=7` coincides with the FIRST projective-plane FAILURE `PG(2,6)`** (Bruck-Ryser /
  Euler's 36 officers; opus-S1). When the projective plane fails (`q=6`), the design covering breaks and the
  hexagonal lattice takes over. klein-S24 (HYP-3706) bridges the two optima via the wallpaper group `p6m`:
  the Steiner-design-optimal (discrete) and the Kershner-continuous-optimal (hexagonal) are the SAME
  hexagonal config under `p6m`. `n=14`: `q=13` prime, `PG(2,13)` exists, construction clean.

## 'The triangular tiling is a god, tridiagonalized' = the A2 Cartan
The equilateral-triangular / hexagonal lattice is `A2` = the Eisenstein integers `Z[zeta_6]` = the
covering-OPTIMAL lattice (Kershner). Its Gram/Cartan matrix is `[[2,-1],[-1,2]]` -- TRIDIAGONAL, `det=3 =`
disc `Q(sqrt-3)`. The right-angle quadrilateral (square) `Z^2` is `[[2,0],[0,2]]` -- DIAGONAL, `det=4`. So
**the equilateral/hexagonal `A2` is the square `Z^2` 'tridiagonalized'** (the off-diagonal `-1` coupling);
a quadrilateral at the right angle (`Z^2`, decoupled) vs two equilateral triangles (`A2`, coupled). The LRC
covering-min lives on `A2` (`Phi_6` = the Eisenstein norm), NOT `Z^2`. And the project's TWO triangles are
ONE lattice in two coordinates: the tournament STAIRCASE (right-isosceles, the 'everything is the triangle'
foundation) is `A2` SHEARED; the equilateral `A2` is the LRC-covering form. The 'god' (the perfect lattice,
`A2`) seen sheared (the staircase, tournaments) and symmetric (hexagonal, the covering).

## The 6/7 boundary (the recurring closure)
`7` is the boundary of tileability/realizability across the whole project:
- **Convex tilings**: a convex polygon tiles the plane only with `<= 6` sides -- all triangles, all quads, 3
  hexagons (Reinhardt 1918), 15 pentagons (Rao 2017); **NONE with `>= 7` sides**. The covering tile is a
  HEXAGON (the max convex tile).
- **LRC apex** `= 7 = Phi_3(2) = |Fano PG(2,2)|`; `14 = 2*7`.
- **Forbidden H-values** `= {7, 21}` (HYP-3616).
- **Covering-min transition** at `n=7` (`PG(2,6)` failure).
So `6` is the largest tileable/realizable degree and `7` is the first impossible/hard one -- one boundary,
many faces (convex tilings, the apex prime, the forbidden H, the projective-plane failure).

## The proof route (the open node) + the aperiodic-monotile lead
The hard core reduces to a **DISCRETE cyclic-Kershner problem** (opus's only open node): is the hexagonal
(three-distance `{1,n,2n}`) covering the cyclic optimum? Kershner 1939 (the continuous hexagonal covering
is the thinnest, density `2pi/sqrt27`) is the PROVEN twin; the route is to make the discrete<->continuous
(Eisenstein-torus) lift rigorous, or to use the three-distance gaps directly (the `{1}+AP` family of
HYP-3702). CREATIVE LEAD (Socolar-Taylor -> hat/spectre): the construction's single OUTLIER `n(n-1)` is the
'closing tile' -- it covers the top two `q`'s at once, the one tile that completes the periodic background
`{1,..,n-2}`; the modern aperiodic MONOTILE (one tile generates the whole tiling) suggests reading the LRC
covering as generated by this single closing tile, and improvements on Socolar-Taylor (the connected,
unmarked hat) as the model for a single-tile covering certificate.

## What it buys
Identifies the LRC covering-min hard core as a TILING/COVERING-OPTIMALITY with two regimes (projective ->
hexagonal) handing off at the 6/7 boundary (the first PG failure), gives 'tridiagonalized triangular tiling'
its exact meaning (`A2` Cartan vs `Z^2`), unifies the recurring `6/7` boundary (convex tilings, apex-7,
forbidden-{7,21}, PG transition), and pins the proof route as the discrete cyclic-Kershner lift -- with
Kershner (1939) the proven continuous twin and the aperiodic monotile the model for a single-tile closure.
