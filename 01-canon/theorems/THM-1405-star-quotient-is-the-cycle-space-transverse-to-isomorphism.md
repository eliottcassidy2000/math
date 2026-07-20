---
id: THM-1405
title: "POINT ADJACENCY DOES PRODUCE THE INVARIANTS, AND THEY ARE F_2 HOLONOMIES — BUT THEY ARE BLIND TO ISOMORPHISM. The star at a vertex (its clique-at-a-point, the map-graph move) is an INCIDENCE ROW of the tile graph K_n minus the base path; that graph is connected for every n >= 4, so <stars> is exactly its CUT space, of dimension n-1. Hence the hypercube quotient is canonically the DUAL OF THE CYCLE SPACE: a tiling's star-orbit is determined by, and determines, the parity of flipped tiles around each cycle of the tile graph -- an F_2 Wilson loop, with star flips as gauge transformations. Orbit counts 2^(m-n+1) = 1, 4, 32, 512 at n = 4..7 match exhaustive enumeration exactly. ANSWERING THE TWO QUESTIONS: (a) the quotient is NOT the merged metagraph and NOT a refinement of it -- for n >= 5 neither partition refines the other, they are TRANSVERSE; (b) the descent question splits cleanly in both directions -- the tournament 3-cycle count is constant on iso classes and NOT on star orbits, while the star-orbit holonomy is a tiling function that does NOT descend to iso classes. So single-tile flips (edge adjacency) have no invariants, point adjacency has 2^(m-n+1) of them, and none of them is a tournament invariant."
status: >
  PROVED for all n >= 4 (the tile graph is connected -- one line -- so incidence rows span
  a cut space of dimension n-1, and F_2^m/C is canonically (C^perp)^* = the dual cycle space).
  VERIFIED-EXACT by full enumeration at n = 4,5,6,7: star-group dimensions 3,4,5,6 = n-1;
  orbit counts 1,4,32,512 = 2^(m-n+1); iso-class counts 4,12,56,456 = A000568; refinement
  tested in BOTH directions on every orbit and every class; the 3-cycle count checked
  constant on iso classes and non-constant on star orbits at every n.
  TRANSVERSALITY is verified at n = 5,6,7, not proved for all n.
source: mac-mini-2026-07-20-S128 (owner: the single-tile flip is edge adjacency -- generates
  everything, no invariants; the map-graph move is point adjacency; a point of K_n is a vertex,
  its clique is a star, stars are incidence rows -- whence duality; the clique-at-a-point
  feature that distinguishes map graphs from planar duals is exactly what produces the
  invariants. Do the cycle invariants descend to iso classes, or are they only tiling
  functions? Is the hypercube quotient by the star group the merged metagraph, a refinement,
  or unrelated?)
depends_on:
  - THM-1390  # the map-graph reading of the waggly filtration (its d_sat claim is corrected by THM-1400 SS I)
related:
  - THM-1400  # kind-pasteur's correction; this is the thread THM-1390 leaves standing
  - the GF(2) Cut (+) Cycle split with the base path as spanning tree (CLAUDE.md)
script: 04-computation/star_quotient_cycle_space_macmini_S128.py (+ .out)
---

# THM-1405 — the star quotient is the cycle space, transverse to isomorphism

**One line.** Point adjacency really does produce invariants that edge adjacency cannot —
single-tile flips generate the whole hypercube, so they have none — and the invariants are
the **F₂ holonomies of the tile graph**. But they are **blind to isomorphism**, exactly as
isomorphism is blind to them.

## Setup

Vertices `0,…,n−1`, base path `n−1 → n−2 → ⋯ → 0`. **Tiles** are the pairs `(i,j)`, `i<j`,
with `j−i ≥ 2` — that is, **`K_n` minus the base path**; `m = C(n−1,2)`. A tiling is a vector
in `F₂^m` (equivalently an integer `0 ≤ N < 2^m`). The **star** at `v` is the set of tiles
incident to `v` — i.e. the **incidence row of the tile graph** at `v`. Star flips generate
`⟨stars⟩ ⊆ F₂^m`.

## (A) The identification — PROVED for every `n ≥ 4`

> **Lemma (connectivity).** The tile graph is connected for all `n ≥ 4`.
> *Proof.* Vertex `0` is tile-adjacent to every `j ≥ 2` (the tile `(0,j)` has `j−0 ≥ 2`), so
> `{0,2,3,…,n−1}` is connected; and `1` is tile-adjacent to `3` (`3−1 = 2`), which is in that
> set. ∎

Incidence rows of a connected graph span its **cut space**, of dimension `n−1`. Hence

> **`⟨stars⟩` = the CUT space of the tile graph, `dim = n−1`.**

The cycle space is `Z = C^⊥`, `dim Z = m−n+1`, and `F₂^m/C` is canonically `Z^*` via
`x + C ↦ (z ↦ ⟨x,z⟩)`. Therefore:

> **A tiling's star-orbit is determined by — and determines — the parity of flipped tiles
> around each cycle of the tile graph.** Star flips are **gauge transformations**; the
> orbit invariants are the **F₂ Wilson loops**. The number of orbits is `2^(m−n+1)`.

**Exact confirmation** (full enumeration):

| `n` | `m` | star-group dim | `n−1` | orbits | `2^(m−n+1)` | iso classes |
|---|---|---|---|---|---|---|
| 4 | 3 | 3 | 3 | 1 | 1 | 4 |
| 5 | 6 | 4 | 4 | 4 | 4 | 12 |
| 6 | 10 | 5 | 5 | 32 | 32 | 56 |
| 7 | 15 | 6 | 6 | 512 | 512 | 456 |

Continuing (dimension count only): `16384` at `n=8`, `2^20` at `n=9`.

This says *which side of the split the stars occupy*. CLAUDE.md already records the GF(2)
**Cut ⊕ Cycle** decomposition with the base path as spanning tree; the new content is that
**the point-adjacency operation is exactly the cut side**, so quotienting by it leaves
precisely the cycle side — which is the duality the question anticipated.

## (B) Is the star quotient the merged metagraph? — **No: transverse**

| `n` | star-orbits refine iso? | iso refines star-orbits? |
|---|---|---|
| 4 | no | yes *(degenerate: cycle dim 0, one orbit)* |
| 5 | **no** | **no** |
| 6 | **no** | **no** |
| 7 | **no** | **no** |

For `n ≥ 5` **neither partition refines the other**. The star quotient is not the merged
metagraph, not a refinement of it, and not refined by it — the two partitions **cross**. At
`n = 4` the cycle space is trivial, the cut space is all of `F₂³`, and the question degenerates
to a single orbit.

## (C) Do the cycle invariants descend to iso classes? — both directions, and they separate

The question has two readings, and each fails in the opposite direction:

- **Tournament cycle invariants descend to iso classes, but not to stars.** The 3-cycle count
  `c₃ = C(n,3) − Σ_v C(s_v,2)` is constant on every iso class (verified at every `n`; it is a
  function of the score sequence) and **not** constant on star orbits (verified, every `n`,
  including the degenerate `n=4`).
- **The holonomy invariants do not descend to iso classes.** The cycle-space class of a tiling
  is a genuine **tiling** function and is **not** an iso-class function — immediate from (B).

> So cycle *holonomies* are **tiling functions**; tournament cycle *counts* are **tournament
> functions**; and the two invariant systems are **transverse**. Neither sees the other.

## (D) What this settles about the map-graph reading

The conjecture behind the question — *the clique-at-a-point feature is exactly what produces
the invariants* — is **confirmed in a precise form and corrected in another**:

- **Confirmed, and sharply.** Edge adjacency (single-tile flips) generates the whole hypercube,
  so it has no invariants at all. Point adjacency has exactly `2^(m−n+1)`, and they are named:
  the cycle parities. The clique-at-a-point really is what produces invariants, and the
  duality really is incidence-row duality.
- **Corrected.** Those invariants are *not* tournament invariants. Point adjacency produces a
  **cohomological** invariant of the tiling — a gauge class — not an isomorphism invariant.
  The duality it exposes is cut/cycle, not tiling/tournament.

## (E) The number reading

A tiling mask **is** a natural number `0 ≤ N < 2^m`, so the hypercube is a canonical
`ℕ`-indexing of the tournaments carrying a fixed Hamiltonian path. The cut/cycle structure
splits those `m` bits into `n−1` **gauge bits** (freely movable by star flips) and `m−n+1`
**holonomy bits** (the invariant content). The content of (B) is that **isomorphism aligns
with neither part**: the arithmetic structure of `N` and the tournament structure of `N` are
genuinely independent gradings of the same object. That is a negative answer to the natural
hope that the integer's binary structure secretly encodes the tournament's symmetry type — and
a positive identification of what it *does* encode.

## Honest scope

- (A) is **proved** for all `n ≥ 4`; the table is exact enumeration confirming it.
- (B)'s transversality is **verified at `n = 5,6,7`, not proved** for all `n`. A proof would
  need a construction of, for every `n`, both an iso class meeting two star orbits and a star
  orbit meeting two iso classes — both look routine but neither is done here.
- (C) uses **one** tournament invariant (`c₃`). That it fails on star orbits proves the systems
  differ; a full comparison of the two invariant *algebras* is not attempted.
- The splitting in (E) is a statement about dimensions and the quotient. `C` and `Z` need not
  be complementary over `F₂` (their intersection is the bicycle space), so the "gauge bits /
  holonomy bits" split is a choice of complement, not canonical.
- Nothing here bears on any LRC or Jacobian thread; it is metagraph structure.

*Artifacts:* `04-computation/star_quotient_cycle_space_macmini_S128.py` (+out).
*Credits:* THM-1390 for the map-graph framing (its `d_sat` invariant is a rediscovery — see
THM-1400 §I, correction accepted); CLAUDE.md's Cut ⊕ Cycle split, which this locates the star
group inside of.
