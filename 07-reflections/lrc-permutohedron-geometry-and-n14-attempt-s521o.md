---
source: oracle-2026-06-01-S521o
status: synthesis + proof attempt (permutohedron geometry of LRC; n=14)
tags:
  - lonely-runner
  - permutohedron
  - braid-arrangement
  - closed-geodesic
  - view-obstruction
  - n14
---

# The Permutohedron Geometry of LRC, and an n=14 Attempt

Rankings are the vertices of the permutohedron; the runner dynamics is a flow
across the braid arrangement. This is the natural geometric home of the whole
tournament-analysis arc — and it makes the precise statement of LRC, and the
precise place it resists, completely clear.

## The geometry

`n` runners on the circle have positions `x_i(t) = frac(v_i t)`. Two facts:

- **Order skeleton = permutohedron.** The sorted (cyclic) order of the runners is
  a vertex of the permutohedron `Pi_{n-1}` (equivalently a region of the **braid
  arrangement** `{x_i = x_j}`). It changes exactly when two runners cross,
  `x_i = x_j (mod 1)`, i.e. `(v_i - v_j)t in Z`. So the runner system is a
  **straight line / closed geodesic on the permutohedral torus**
  `T^{n-1} = R^{n-1}/Lambda` (Lambda = the A_{n-1} lattice; the permutohedron is its
  space-filling zonotope / fundamental cell), winding with homology class
  `(v_1, ..., v_{n-1})`. Its cyclic-order walk is a **closed walk on the (cyclic)
  permutohedron**.

- **Loneliness = a SECOND arrangement.** Observer-loneliness is
  `x_i in [1/n, 1-1/n]` for all `i`: the walls `x_i = ±1/n`, a **different**
  hyperplane arrangement than the braid walls. LRC = the closed geodesic enters
  the lonely region of this second arrangement.

So LRC lives in the **combined arrangement** `{x_i = x_j} ∪ {x_i = ±1/n}` on
`T^{n-1}`: the permutohedron/braid part carries the *order*, the `1/n`-walls carry
the *metric loneliness*, and the runner orbit is a closed geodesic navigating both.

## What the geometry unifies (computed, `lrc_permutohedron_geometry_s521o.py`)

- **Holdback = braid crossings.** The total braid crossings over `[0,1)` is
  `Σ_{i<j} |v_i - v_j|` (each pair crosses `|v_i-v_j|` times) — the staircase total
  of S25. For *resonant* (arithmetic) speeds the crossing *times* coincide
  massively: extremal `(0,1,2,3,4)` has only **6 distinct** wall-times though `20`
  crossings — the synchrony of S25, now as coincident braid hyperplanes.
- **Few permutohedron cells.** The orbit visits a tiny subset of the `(n-1)!`
  cyclic orders: `6, 10, 12` (extremal `n=5,6,7`) of `24, 120, 720`. The
  order-analogue of the `2·Fib(n-2)` circular menu (S518).
- **Loneliness is strictly finer than the order skeleton.** Checking loneliness at
  braid-cell *centers* gives `0` lonely cells — even for LRC-easy sets — because
  the lonely sub-interval sits off-center, inside a braid cell, carved by the
  `1/n`-walls. Loneliness does **not** factor through the cyclic order: it is the
  metric thickening, exactly the **view-obstruction** (the line piercing the
  central box `B = [1/n, 1-1/n]^{n-1}`). This is the geometric form of the S519
  conclusion that the tournament/order is too coarse.

## The n=14 attempt

In these coordinates LRC@14 is: *the closed geodesic of class `(v_1,...,v_{13})`
on `T^{13}` pierces the central box `B` of side `6/7`.* Two geometric levers, and
where each stalls:

1. **Equidistribution.** If `(v_i)` are `Q`-linearly independent, the geodesic is
   dense, so it enters the open box `B` — done. The hard cases are exactly the
   **closed geodesics** (commensurable / integer speeds), where density fails. So
   the difficulty is intrinsically about closed orbits, not generic ones — the
   permutohedron makes this dichotomy crisp but does not remove the closed case.

2. **Weyl symmetry + zonotope tiling.** `B` and the permutohedron both carry the
   `S_{14}` action; the permutohedron tiles `T^{13}`. One would like a
   covering/symmetry argument: the closed geodesic, winding `Σ|v_i-v_j|` times
   through the tiling, must cross `B`. But `B` is a *cube*, not Weyl-aligned with
   the permutohedron, and the orbit can be engineered (resonant speeds) to thread
   the braid cells while skirting the cube — this is precisely the unresolved
   view-obstruction, and the box vs permutohedron mismatch is why no clean lattice
   covering closes it.

**Honest verdict.** The permutohedron is the *correct unifying language* — it
places rankings as vertices, the dynamics as a closed geodesic, holdback as braid
crossings, and the circular menu as the visited cells — but LRC@14 remains the
**closed-geodesic-pierces-the-box** problem. The metric box (the second `1/n`
arrangement) is not subordinate to the permutohedral (order) structure, so the
order geometry alone cannot force the piercing. This is the same wall as S519/S520
(order/tournament is coarse; the `1/n` metric is the content), now seen at maximal
geometric clarity.

## The one place it sharpens the attack

The combined-arrangement picture isolates the real object: the **closed geodesic's
intersection pattern with the `1/n`-walls, indexed by the permutohedral cell it is
in**. Inside each braid cell the order is fixed, so loneliness reduces to a
*linear* (interval) condition on `t` — the `1/n`-walls cut each cell into
lonely/non-lonely sub-intervals. LRC@14 = the union of lonely sub-intervals over
all cells is nonempty. This is the bounded-ansatz/wall-order program (S514) given
its geometric meaning: **enumerate the lonely sub-interval of each visited
permutohedral cell.** The cells are few (the menu); the open question is whether
the speed arithmetic can empty every one — the realizability constraint (S520o) in
geometric dress.

## Verdict / next
- Did not prove n=14; correctly placed it as a closed-geodesic view-obstruction on
  `T^{13}`, with the permutohedron as the order-skeleton and `B` the metric target.
- Concrete next: for each visited braid cell, compute its lonely sub-interval (a
  1D linear condition), and study whether the speed arithmetic (the `1/n`-wall
  positions) can avoid all of them — the geometric bounded-ansatz.

## Artifacts
```
04-computation/lrc_permutohedron_geometry_s521o.py
05-knowledge/results/lrc_permutohedron_geometry_s521o.out
```
