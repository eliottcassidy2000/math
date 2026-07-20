---
id: THM-1405
title: "THE STAR GROUP: point-adjacency does produce an invariant, and it does NOT descend to tournaments. (0) THE DUALITY IS EXACT. Tiles are the cycle-space coordinates of K_n with the base path as spanning tree; the vertex stars S_v span the CUT space, so 'stars are incidence rows' is literally that Γ = ⟨S_v⟩ is the row space of the vertex–edge incidence matrix restricted to non-tree edges. Σ_v S_v = 0 since every tile has two endpoints, and for n ≥ 4 no non-zero cut is supported on the base path (δ(v) has n−1 edges but at most 2 path arcs), so the restriction is injective on the cut space and rank Γ = n−1 EXACTLY. Hence the point-adjacency quotient has 2^{C(n−1,2)−n+1} cosets: 1, 4, 32, 512 at n = 4,5,6,7 — where the edge-adjacency (single-tile) move generates everything and has no invariants at all. (I) BUT THE INVARIANT IS A TILING FUNCTION, NOT A TOURNAMENT FUNCTION: 8 of 12 iso classes at n=5 and 50 of 56 at n=6 carry more than one coset. The n=4 case appears to descend and is VACUOUS — the quotient there is a single coset. (II) AND THE QUOTIENT IS UNRELATED TO THE MERGED METAGRAPH: no coset lies inside an iso class at any n (0 of 1, 0 of 4, 0 of 32), so the two partitions CROSS — neither refines the other — and the coset counts 1, 4, 32 do not match the merged-class counts 3, 10, 34. Star moves do not even preserve c₃ (1018 of 1024 tilings at n=6). (III) THE CGP QUESTION, in the reading the star sentence forces, CLOSES: at tile level the witness is FORCED to be the tile–vertex incidence graph of K_n minus the base path, which is planar iff that graph is, i.e. iff n ≤ 6 — and non-planar for ALL n ≥ 9 with no computation, since vertices 0,2,4,… are pairwise non-consecutive and give a clique of size ⌈n/2⌉ ≥ 5"
status: >
  (0) PROVED (the rank argument is complete and elementary) + VERIFIED at n = 4,5,6.
  (I),(II) PROVED-BY-EXHAUSTION over all 2^m tilings, n ≤ 6, with explicit mixed-class
  witnesses.  These are negative answers to the two questions asked, and the n=4
  vacuity is flagged because it is a trap: the quotient is trivial there, so "it
  descends" carries no information.
  (III) PROVED for the forced witness, with the alternating-clique argument giving
  non-planarity for all n ≥ 9 unconditionally.  SCOPE, stated rather than buried: CGP
  is an EXISTENTIAL over bipartite witnesses, so a non-planar natural witness is
  evidence and not a proof, and this says nothing about G^(≤k) on ISO CLASSES, which
  is a different graph.  That question remains open.
  Nothing here advances a named open problem.
source: kind-pasteur-2026-07-20-S128c109 (owner: the single-tile flip is edge adjacency and generates everything with no invariants; the map-graph move is point adjacency, a point of K_n is a vertex, its clique is a star, stars are incidence rows, whence duality — do the cycle invariants descend to iso classes, and is the hypercube quotient by the star group the merged metagraph, a refinement, or unrelated?)
depends_on: []
related: [THM-1390, THM-1400, HYP-8235]
script: 04-computation/star_group_quotient_kps_S128c109.py, cgp_planar_witness_kps_S128c109.py (+ .out)
---

# THM-1405 — the star group, and what it does and does not see

## 0. The duality is exact, and the rank is exactly n−1

In the `GF(2)` edge space of `K_n` with the base path as spanning tree, the **tiles are
the non-tree edges** (cycle-space coordinates) and the **vertex stars `δ(v)` span the cut
space**. So `Γ = ⟨S_v⟩`, the star group acting on tilings, is precisely the row space of
the vertex–edge incidence matrix restricted to non-tree edges — the owner's "stars are
incidence rows", literally.

> **`Σ_v S_v = 0`** — every tile has exactly two endpoints. So `rank Γ ≤ n−1`.
>
> **`rank Γ = n−1` exactly, for `n ≥ 4`.** A cut in the kernel of the restriction is
> supported on the base path. But `δ(S)` has `|S|·|S^c|` edges while the path has `n−1`,
> forcing `|S| = 1`; and `δ(v)` has `n−1` edges of which at most 2 are path arcs, so
> `n−1 ≤ 2`, impossible for `n ≥ 4`. ∎

Verified at `n = 4,5,6`. Hence the point-adjacency quotient has

> `|F₂^m / Γ| = 2^{C(n−1,2) − n + 1}` = **1, 4, 32, 512** at `n = 4,5,6,7`.

This is the owner's contrast made precise: the single-tile flip generates all of `F₂^m`
and so has **no** invariants; the star move generates only `Γ`, and the coset `b + Γ` is
the invariant it produces.

## I. The invariant does not descend

**Do the cycle invariants descend to iso classes — tournament functions, or only tiling
functions?** Only tiling functions.

| n | iso classes | classes carrying >1 coset |
|---|---|---|
| 4 | 4 | 0 — but **vacuous**, the quotient is a single coset |
| 5 | 12 | **8** |
| 6 | 56 | **50** |

The `n = 4` row is a trap and is flagged as such: with one coset everything trivially
descends. From `n = 5` on, most classes are mixed — witness class at `n=5` carries cosets
`{0,1,2}`. So `b + Γ` is **not** a tournament function.

Corroborating, at the level of the most basic tournament invariant: a star move fails to
preserve the 3-cycle count `c₃` for **1018 of 1024** tilings at `n = 6`. The star move is
not an isomorphism and does not act through the iso quotient at all.

## II. The quotient is unrelated to the merged metagraph

**Is the hypercube quotient by the star group the merged metagraph, a refinement, or
unrelated?** Unrelated — the partitions **cross**.

- Cosets lying entirely inside one iso class: **0 of 1, 0 of 4, 0 of 32** at `n = 4,5,6`.
  So the coset partition does not refine the iso partition.
- Iso classes are not unions of cosets either (that is §I).
- Coset counts `1, 4, 32` against merged-class counts `3, 10, 34`. Close at `n=6`
  (`32` vs `34`) and **not** equal — a near-miss worth not over-reading.

Neither partition refines the other at any `n ≥ 5`.

## III. The CGP question, in the reading the star sentence forces

Take the sentence at **tile** level: faces are tiles, two tiles adjacent when they share a
vertex of `K_n`, and the clique at the point `v` is exactly the star `S_v`. Then the
bipartite witness is **forced**:

> `B` = the tile–vertex incidence graph of `K_n` minus the base path.

The half-square of `B` on the tile side *is* the tile-adjacency graph, and the incidence
graph of a graph `H` is planar exactly when `H` is (it is a subdivision of `H`). So the
CGP criterion reduces to one planarity question, and it resolves:

| n | 3 | 4 | 5 | 6 | **7** | 8 | 9+ |
|---|---|---|---|---|---|---|---|
| `K_n` − path planar? | Y | Y | Y | Y | **N** | N | N |

with `n = 7` sitting exactly at `|E| = 3n−6 = 15`. And the failure is unconditional for
large `n` with no planarity routine at all: **vertices `0,2,4,…` are pairwise
non-consecutive**, so every edge among them survives removal of the base path, giving a
clique of size `⌈n/2⌉` — a `K₅` as soon as `n ≥ 9`.

**Scope, stated at its true strength.** CGP is an *existential* over bipartite witnesses,
so a non-planar natural witness is evidence, not a proof; and this is the tile-level
graph, not `G^(≤k)` on iso classes. **The question I flagged as open last session — is
`G^(≤k)` the half-square of a planar bipartite graph — is still open.** What closes here
is the tile reading, and it closes negatively past `n = 6`.

## Why this is a useful negative

The owner's diagnosis was right in its positive half: point adjacency *does* produce an
invariant where edge adjacency produces none, and the invariant has exactly the size the
incidence-matrix duality predicts. What fails is the transfer to tournaments. The star
group acts on the *cube*, not on the *folding* — the same distinction recorded in
`07-reflections/a-bijection-of-sets-is-not-a-bijection-of-structure.md` and in THM-1400's
half-cube correction. Two independent constructions have now produced invariants that
live on tilings and die under `S_n`, which suggests the pattern is structural rather than
accidental: **cut-space constructions are base-path-relative, and the base path is exactly
what the iso quotient forgets.**

## Named next

- The genuinely open form: is `G^(≤k)` on **iso classes** the half-square of a planar
  bipartite graph? The forced-witness trick does not apply there, since no witness is
  canonical.
- If a star-type invariant is wanted that *does* descend, it must be built from a
  base-path-**independent** subgroup. The natural candidate is the intersection of `Γ`
  over all choices of spanning path — worth one computation, and cheap at `n ≤ 6`.
