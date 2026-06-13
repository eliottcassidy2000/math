---
theorem: THM-333
name: Apex Tile Strong Connectivity Theorem
status: PROVED
session: opus-2026-05-27-S1
verified: computationally n=3..6 (all apex-upward tilings are SC)
depends_on: THM-330
---

## Statement

In the tiling model on n vertices (vertices 0,...,n-1, base path n-1→n-2→...→0), let the **apex tile** be the tile (n-1, 0) — the arc between the highest-indexed vertex (n-1) and the lowest-indexed vertex (0).

If the apex tile is in the **upward orientation** (vertex 0 beats vertex n-1, bit=1), then the tournament is **strongly connected**.

## Proof

By THM-330, it suffices to show that the apex tile (n-1, 0) in upward orientation satisfies ALL cuts k ∈ {1,...,n-1}.

The apex tile (n-1, 0) crosses cut k iff y=0 < k and x=n-1 ≥ k, i.e., 1 ≤ k ≤ n-1. This holds for ALL n-1 cuts simultaneously.

When bit=1 (upward: 0→(n-1)), every cut k is satisfied by this tile. By THM-330, the tournament is strongly connected. □

## Remarks

**Converse fails:** Not every strongly connected tiling has the apex tile upward. At n=5: 32 of 64 tilings have apex downward (bit=0), and 18 of these 32 are still strongly connected (other tiles satisfy all cuts).

**Apex as Hamiltonian cycle closer:** The base path n-1→n-2→...→0 is a directed Hamiltonian path. The arc 0→(n-1) closes this path into the Hamiltonian cycle n-1→n-2→...→0→(n-1). By the theorem of Camion (1959), a tournament with a directed Hamiltonian cycle is strongly connected. This gives an independent proof.

**Geometric position:** In the staircase Young diagram, the apex tile (n-1, 0) occupies the **top-left corner** — the intersection of the vertical leg (source column) and horizontal leg (sink row). It is the unique tile that simultaneously connects the natural source (vertex n-1) and natural sink (vertex 0) of the base path.

## SC Fraction Change

When the apex tile is forced to be upward (bit=1), the tiling is ALWAYS SC. When apex is downward (bit=0): some tilings are SC, some are not. Specifically:

| n | SC with apex=1 | SC with apex=0 | Total SC |
|---|----------------|----------------|----------|
| 3 | 1/1 = 100%    | 0/1 = 0%       | 1/2      |
| 4 | 4/4 = 100%    | 1/4 = 25%      | 5/8      |
| 5 | 32/32 = 100%  | 18/32 = 56.25% | 50/64    |
| 6 | 512/512 = 100%| 391/512 = 76.4%| 903/1024 |

## Effect on H

The apex tile flip (apex=0 → apex=1) does NOT always increase H. In particular, it can **decrease** H for strongly connected tournaments that already have the "correct" local structure:

| n | Increase | Same | Decrease |
|---|----------|------|----------|
| 3 | 1 | 0 | 0 |
| 4 | 3 | 1 | 0 |
| 5 | 26 | 4 | **2** (both from SC, one is the regular tournament H=15→13) |
| 6 | 385 | 62 | **65** (ALL from SC tournaments) |

**Key observation:** ALL H-decreasing cases arise from SC tournaments. Non-SC tournaments uniformly benefit from the apex flip (it both creates SC and increases H). For SC tournaments near the H-maximum (regular, near-regular), flipping the apex can disrupt the optimal cycle structure.
