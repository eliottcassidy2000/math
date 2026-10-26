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

## Both open questions ANSWERED (klein-S334) — and both negatively

The two questions this theorem raised are settled by comparing the two partitions of the `2^m` tilings
directly: by tournament **iso class**, and by **star orbit**.

| `n` | tilings | iso classes | star orbits | iso refines orbits? | orbits refine iso? |
|---|---|---|---|---|---|
| 4 | 8 | 4 | 1 (`=2^0`) | *(vacuous)* | no |
| 5 | 64 | 12 | 4 (`=2^2`) | **NO** | no |
| 6 | 1024 | 56 | 32 (`=2^5`) | **NO** | no |

**(i) The cycle invariants do NOT descend to iso classes.** They are functions of the *tiling*, not of the
tournament. Explicit witness at `n = 5` (tiles `(3,1),(4,1),(5,1),(4,2),(5,2),(5,3)`): one iso class
contains 5 tilings spread over **4 distinct star orbits** — e.g. `(0,0,0,0,1,0)` and `(0,0,0,1,1,0)` are
isomorphic tournaments whose star-orbit representatives are `(0,0,0,0,1,0)` and `(0,0,0,0,0,0)`. A second
iso class spreads 9 tilings over 3 orbits.

**(ii) The star quotient is UNRELATED to the merged metagraph.** Neither partition refines the other in
either direction at `n = 5, 6`; the two are transverse. So it is not the merged metagraph, not a
refinement of it, and not a coarsening.

**What this means.** The algebra of §"The theorem" is exact and complete — but it is complete *for the
tiling model*, which is a tournament **together with a distinguished Hamiltonian path**. The star-flip
invariants are precisely quantities that see the path. This is a sharp instance of the CLAUDE.md warning
that the tiling model "has lower symmetry than the arc model" and "breaks the `S_n` isotropy": the cycle
space of `H = K_n ∖ P` is *defined from* `P`, so relabelling moves `P`, moves `H`, and moves the
invariants. Nothing in the construction was ever `S_n`-equivariant, and the census shows that is not a
repairable defect but the content.

The honest reading: THM-1382 characterises a move-set and its conserved quantities completely, at the
level below iso classes. Anyone wanting tournament invariants from this must first symmetrise over the
choice of Hamiltonian path — and the transversality in the table above suggests that symmetrisation
destroys the information rather than descending it.

*Method note.* Both answers are **witness-effective**: the negative is not deduced from an abstract
obstruction but exhibited as a concrete isomorphic pair in different orbits, checkable in seconds. That is
the standard this repo's tournament analysis should hold itself to — tournaments *are* binary relations,
and a claim about them should come with the relation that witnesses it.

## Scope

This is structure, not a bound: it introduces and completely characterises a move-set on the tiling
model, and says nothing about `H`-values, LRC, or the metagraph's metric invariants.

*Files: `04-computation/mapgraph_starflips_klein_S333.py` (+ `_starflips`, `_invariants` .out).*
