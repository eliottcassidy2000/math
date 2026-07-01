---
id: HYP-3811
title: THE TWIN SEPARATOR resolves the 5 mixed twin-pairs at n=6 -- 3 twins separate by tiling count / H, and the remaining 2 are INVARIANT-COSPECTRAL (identical H, |Aut|, c3, score sequence) SPECTRAL TWINS {12,22}(grid-sym 5) and {43,44}(grid-sym 7), separated ONLY by the adjacency eigenvalue spectrum. The twin pairing is NUMERICAL (by grid-sym count), NOT the flip.sigma involution (which does not pair cleanly). The 2 spectral twins are the tournament analog of DUAL polytopes (same combinatorial data, geometrically/spectrally distinct) -- there are exactly 2 dual pairs among the 5 Platonic solids (cube/octahedron, dodecahedron/icosahedron; tetrahedron self-dual), and n=6 is the first HYPERBOLIC size past the icosahedral n=5 (G_5 f-vector = (12,30,20) = icosahedron, kps S20bh) -- so the twins are the hyperbolic-onset refinement of the spherical Platonic structure.
status: CONFIRMED (n=6 exhaustive: the 5 mixed twins, their invariants, the tiling-count vs spectrum separator, the 2 invariant-cospectral spectral twins). The Platonic/2D-tiling connection is a STRUCTURAL synthesis (via kps S20bh's G_n<->Platonic-levels + the spherical/hyperbolic frame + the dual-polytope analogy), not a forced numerical identity.
source: mac-mini-2026-07-01-S86
related:
  - HYP-3809   # S84 the twin-pairing conjecture (grid-sym counts all-even multiplicities) -- this runs its separator
  - HYP-3810   # S85 T-join obstruction (the SC nodes being separated here)
  - HYP-3808   # S83 blue/black parity; n=6 transition
results:
  - 04-computation/twin_separator_n6_platonic_macmini_20260701.py
  - 05-knowledge/results/twin_separator_n6_platonic_macmini_20260701.out
---

# HYP-3811 -- the twin separator and the spectral twins (n=6, Platonic)

Owner's task: run the twin separator against the 5 remaining twins, and connect to 2D plane tilings /
Platonic solids. The "5 remaining twins" = the 5 twin-pairs among the 10 MIXED SC merged nodes at n=6
(the pure-blue pair is the trivial 6th); from HYP-3809's twin-pairing (grid-sym counts pair up: `{3:2, 5:2,
7:4, 9:2}`).

## The 5 twins and the separator hierarchy (n=6, exhaustive)
| twin | grid-sym | members (node: tiling count) | separated by |
|---|---|---|---|
| 1 | 3 | 5:5, 55:15 | tiling count (H) |
| 2 | 5 | 12:17, 22:17 | **adjacency SPECTRUM** (invariant-cospectral) |
| 3 | 7 | 30:37, 31:33 | tiling count (H) |
| 4 | 7 | 43:29, 44:29 | **adjacency SPECTRUM** (invariant-cospectral) |
| 5 | 9 | 46:41, 54:15 | tiling count (H) |

- **3 twins** (1,3,5) separate by **tiling count = H/|Aut|** (H = 5/45, 37/33, 41/45).
- **2 twins** (2,4) are **INVARIANT-COSPECTRAL**: the two members have IDENTICAL grid-sym count, tiling
  count, `H`, `|Aut|=1`, `c3`, and score sequence -- yet are distinct iso classes. They separate ONLY by
  the **adjacency eigenvalue spectrum**. These are the **SPECTRAL TWINS**: `{12,22}` (`scores (1,1,2,3,4,4)`,
  `H=17`, `c3=4`) and `{43,44}` (`scores (1,2,2,3,3,4)`, `H=29`, `c3=6`).

So the twin separator is a two-tier statistic: **tiling count first, adjacency spectrum for the residue.**
The twin pairing itself is by grid-sym count (numerical); `flip.sigma` does NOT induce it (it lands in many
nodes), so the twinning is a parity coincidence, not a structural involution -- an open sub-question.

## The Platonic / 2D-tiling connection (structural)
Via kps S20bh (`the-five-platonic-tournaments.md`): the 5 Platonic solids sit at 5 levels of tournament
theory, and `G_5` has the ICOSAHEDRON's f-vector `(12,30,20)` with `|S_5| = 2|A_5| = 2|Aut(icosahedron)|`
(the `x2` = complement). `n=5` is the pentagon -- the spherical/hyperbolic boundary; `n<=4` Platonic
(spherical), `n>=6` hyperbolic. So:
- **The twins are an `n=6` (hyperbolic-onset) phenomenon** -- the first size PAST the icosahedral `n=5`. The
  merged metagraph's spherical Platonic regularity (n<=5) fractures into twin-pairs at n=6, the same `n=6`
  where the black sea onsets, pure-black self-loops appear, the blue T-join goes non-bipartite (HYP-3808/10),
  and the minimal-flip gauge breaks (HYP-3798). The twinning is a 5th face of the `n=6` transition.
- **The 2 spectral twins mirror the 2 dual Platonic pairs.** Dual polytopes (cube/octahedron,
  dodecahedron/icosahedron) share symmetry group and combinatorial data but are geometrically distinct; the
  2 invariant-cospectral spectral twins share every combinatorial invariant but are spectrally distinct --
  the same "same-data / dual-realization" signature. There are exactly 2 dual pairs (tetrahedron self-dual),
  matching the 2 spectral twins.
- **2D-tiling frame:** the regular EUCLIDEAN plane tilings (curvature 0) are the flat boundary between
  spherical (Platonic, `n<=5`) and hyperbolic (`n>=6`); the `n=6` twins are the first tournament structures
  in the hyperbolic (negative-curvature) regime, the analog of hyperbolic tilings.

## New sub-conjectures
1. The spectral twins persist for all `n>=6`: some SC twins are invariant-cospectral but spectrally
   distinct, and their count relates to the `2` dual Platonic pairs / the hyperbolic genus of `G_n`.
2. The twin pairing (grid-sym count all-even) has a structural involution not yet identified (`flip.sigma`
   ruled out) -- perhaps a genuine cospectral-mate map on the fold.
3. The near-regular SC nodes (scores `(2,2,2,3,3,3)`: 46,54,55; and the high-`|Aut|=9` pure-blue node 23,
   scores `(1,1,1,4,4,4)`) are the "most Platonic" (spherical-symmetric) survivors -- the residual
   icosahedral regularity inside `n=6`.

## Honest scope
The 5 twins, the separator hierarchy (tiling-count then spectrum), and the 2 invariant-cospectral spectral
twins are EXACT (n=6 exhaustive). The Platonic/2D-tiling connection is a structural synthesis on top of kps
S20bh, honest and suggestive but not a numerical identity; the dual-pair <-> spectral-twin match (2 = 2) is
the sharpest concrete point. The twin involution is unidentified.
