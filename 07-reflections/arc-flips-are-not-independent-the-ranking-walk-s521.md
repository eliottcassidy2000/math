# Arc-flips are not independent: the LRC walk lives on the ranking variety, not the cube (S521)

*claudebox-2026-06-01-S521, thinking seed. On why the `2^{C(n,2)}` free-switch model
misses the essence, made quantitative for the LRC walk. Connects the tiling
model (cut ⊕ cycle), Rédei parity, the OCF conflict graph, and the round-tournament
menu (THM-387/HYP-1993).*

## The thesis

The labeled count `2^{C(n,2)}` treats each arc as an independent switch.  It is
the wrong coordinate system.  A tournament is a *ranking* (Hamiltonian path)
recursively composed of sub-rankings, some aligned with the main path
(transitive) and some against it (the cyclic strongly-connected blocks).  The
arcs are bound by hidden dependence at three levels:

1. **Rédei parity.** `H(T)` (the number of rankings) is always odd; no arc
   assignment yields an even ranking-count.  One global law over all `C(n,2)`
   switches.
2. **Cut ⊕ cycle (the tiling).** Over GF(2) the arcs split into an `(n-1)`-dim
   cut space (the score hierarchy = the base-path ranking) and a `C(n-1,2)`-dim
   cycle space (the tiles = deviations).  The base-path arcs *are* the ranking;
   tiles are the "against-the-current" sub-rankings.
3. **OCF.** `H = I(Ω(T),2)` is a function of the odd-cycle collection, not of
   individual arcs; one flip can create/kill many odd cycles at once.

The recursive "tree of sub-rankings" is, precisely, the **SCC condensation
recursion**: a transitive backbone of strongly-connected blocks, each block a
smaller cyclic tournament, recursing.  That is exactly the shape of the
round-tournament menu (transitive backbone + reversed up-set).

## The dependence is measurable in the LRC walk

For a speed set, let `t` run over `[0,1)`; the runner sub-tournament changes only
at walls (a pair crosses).  Computed (`m = n-1` movers):

| speeds | cells visited | hypercube `2^{C(m,2)}` | fraction |
|--------|---------------|------------------------|----------|
| (1,2,3,4,5)    | 12 | 1024  | 1.2% |
| (1,2,4,7,11)   | 56 | 1024  | 5.5% |
| (2,3,5,7,11,13)| 72 | 32768 | 0.2% |

The walk never enters the free cube; it is confined to the **round tournaments**
(`Circ(m) = A000016(m)` iso-classes — vanishing inside `A000568(m)`, itself
vanishing inside the cube).  Confinement deepens with `m`.  This is the hidden
dependence made into a number: the "independent axes" are a fiction; the reachable
set is a thin recursive-ranking variety.

## The two moves (single-arc, but not free)

Every generic wall is a **single arc flip** (Hamming distance 2 on the labeled
matrix), of exactly one of two kinds as `t` increases:

- **Order move** — two runners coincide (separation 0): an *adjacent
  transposition of the circular ranking*.  Changes the ranking.  Lives in the
  cut/score hierarchy.
- **Horizon move** — two runners go antipodal (separation 1/2): the pair crosses
  the half-turn boundary with the circular order fixed.  Flips the arc without
  reordering.  Lives in the cycle/half-turn cut.

Coincident walls (e.g. arithmetic-progression speeds) fire several of these at
once — the arithmetic of the speeds *couples* the flips (Hamming 4, 8, 20, 30 in
the data).  Even the simultaneity is dependence.

So the LRC walk is a path on the **affine permutohedron of circular rankings**
(order moves = adjacent transpositions, the Coxeter/type-A structure of my arc-menu
up-sets) decorated by **horizon bits** (which antipodal side each pair is on) —
not a walk on a Boolean cube.  Ranking and horizon are the cut and cycle halves of
the tiling, now seen dynamically.

## Why this is the leverage, not a nuisance

LRC is "almost always true" precisely because the walk *cannot* wander freely:
it is pinned to the round-tournament variety, and the lonely cell sits inside that
thin set (S521: the safe menu = the whole round set, up to the observer mark).
A would-be counterexample cannot exploit independent arcs; it must thread the
recursive-ranking variety while the observer mark stays blocked — a far smaller
search than `2^{C(n,2)}` suggests.  The right state space for an LRC proof is the
ranking walk (permutohedron + horizon bits), where moves are local and the
forbidden (observer-unsafe) region is a union of faces.

## Seed

Model LRC(n) explicitly as: a closed walk on the **affine permutohedron of
circular rankings** of the `n-1` movers (order moves) with horizon decorations
(horizon moves), the step sequence dictated by the speeds; the observer mark adds
`n-1` safety bits.  LRC(n) = every such walk hits the all-safe face.  Because the
permutohedron is `(n-2)`-dimensional and the walk is the projection of a straight
torus line, this is the genuine geometry: a line crossing a zonotopal fan, asked
to enter a central cell.  Counting/forcing on the permutohedron (not the cube) is
where the dependence becomes a tool.
