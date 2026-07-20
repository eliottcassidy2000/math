---
id: THM-1382
title: THE MAP-GRAPH MOVES ARE A SECOND CUT/CYCLE SPLIT — vertex-star flips on the staircase generate exactly the CUT space of H = K_n ∖ (base path), of dimension n−1, and their invariants are exactly the CYCLE space of H, of dimension C(n−1,2)−(n−1). This nests a second GF(2) cut/cycle decomposition *inside* the tile space, which is itself the cycle space of K_n relative to the base path.
status: PROVED (incidence-matrix identification; rank verified n = 4…10; the invariance verified over 12,000 random (tiling, star-flip) trials at n = 6,7,8 with ZERO violations).
source: klein-2026-07-20-S333
depends_on:
  - THM-127   # dihedral anti-automorphism (the involution side of the tiling model)
related:
  - THM-1381  # the torus index obstruction (the other half of this session)
---

# THM-1382 — the map-graph moves, and a nested cut/cycle split

A **map graph** joins two faces of an embedded graph when they meet at a *point* — a vertex or an edge —
rather than only along an edge. Its defining feature is that **the faces meeting at a point form a
clique**. Transporting that to the staircase model says which moves are natural: not single-tile flips,
but the whole **clique at a point**.

## The dictionary

In the tiling model (CLAUDE.md conventions) the tiles are the non-base-path arcs
`{(x,y) : x − y ≥ 2}`, and `m = C(n−1,2)`. Reading a tile as an *edge* of `K_n` and a vertex of `K_n` as
a *point where faces meet*:

| map graph | staircase / tournament |
|---|---|
| face | tile `(x,y)` = a non-path arc |
| point where faces meet | vertex `v` of `K_n` |
| clique of faces at a point | `star(v) = {(x,y) : x = v or y = v}` |
| map-graph move | **flip all of `star(v)` at once** |

CLAUDE.md lists vertex-star flips among the "creative waggly subsets" and records only that their
neutrality is non-uniform (centre ≈ 5× the endpoints at `n=5`). The map-graph framing says they are not
one option among many — they are *the* generators the incidence structure supplies. Their algebra turns
out to be completely determined.

## The theorem

Let `H = K_n ∖ P` where `P` is the Hamiltonian base path. Note `|E(H)| = C(n,2) − (n−1) = C(n−1,2) = m`:
**the tiles are exactly the edges of `H`.**

**Theorem.** Over `GF(2)`:

1. Each tile lies in exactly two stars, so `Σ_v star(v) = 0`; the star vectors are the rows of the
   incidence matrix of `H`.
2. `H` is connected, so the star vectors span the **cut space** of `H`, of dimension `n − 1`.
3. The invariants of the star-flip action are exactly the **cycle space** of `H`, of dimension
   `m − (n−1) = C(n−1,2) − (n−1)`.

Concretely: *for every cycle `C` of `H`, the parity of the orientation bits around `C` is conserved by
every vertex-star flip.*

**Verification.** Rank of the star span, `n = 4…10`: `3, 4, 5, 6, 7, 8, 9` — equal to `n−1` at every `n`,
with invariant dimension `0, 2, 5, 9, 14, 20, 27` matching `C(n−1,2)−(n−1)` exactly. Invariance: 4,000
random (tiling, star-flip) trials at each of `n = 6,7,8`, testing all short cycles of `H` — **zero parity
changes** in 12,000 trials.

## The nesting — what is new

CLAUDE.md already records one `GF(2)` split: with the base path as spanning tree,
**base-path arcs = cut space, tiles = cycle space** of `K_n`. This theorem finds a *second* split living
inside the first:

```text
level 1   K_n            =  cut(P)          ⊕  cycle(K_n rel. P)   ← the tiles
level 2   the tile space =  cut(H)          ⊕  cycle(H)            ← star flips vs their invariants
                            dim n−1             dim C(n−1,2)−(n−1)
```

So the tile space is not structureless: it carries its own cut/cycle duality, and the map-graph moves are
precisely its cut half. The invariants are a new family of conserved quantities for the tiling model —
one per independent cycle of `K_n ∖ P` — and they are computable, exact, and (unlike the neutrality
statistics previously recorded) *complete*: nothing outside the cycle space of `H` is conserved, because
the star flips span all of `cut(H)`.

## Why the map-graph framing was the right lens

The single-tile (wiggly) move is the *edge* adjacency of the tiling; it generates everything and has no
invariants. The map-graph move is the *point* adjacency, and because a point of `K_n` is a vertex, its
clique is a star, and stars are incidence rows — whence cut/cycle duality. The clique-at-a-point feature
that distinguishes map graphs from planar duals is exactly what produces the invariants.

## Scope

This is structure, not a bound: it introduces and completely characterises a move-set on the tiling
model, and says nothing about `H`-values, LRC, or the metagraph's metric invariants. Two immediate
questions it raises, neither attempted here: (i) what do the cycle-space invariants compute on
iso classes — are they functions of the tournament, or only of the tiling? (ii) the quotient of the
tiling hypercube by the star-flip group has `2^{C(n−1,2)−(n−1)}` classes; is that quotient the merged
metagraph, a refinement of it, or unrelated?

*Files: `04-computation/mapgraph_starflips_klein_S333.py` (+ `_starflips`, `_invariants` .out).*
