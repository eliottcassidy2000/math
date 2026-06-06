---
id: THM-421
name: unit-distance-3n-common-neighbour-floor
status: PROVED (lower bound 17) + VERIFIED construction (upper bound 32)
date: 2026-06-06
session: monad-explorer-2026-06-06-S709
depends_on:
  - HYP-2262   # kissing cap is not the unit-distance cap (S702): two norm-1 layers
  - HYP-2267   # triangular sqrt(7) DISK patch beats 3N at N=43 (S703), now sharpened to 32
  - THM-412    # density quantization (the quantized 2D unit-distance density ladder)
  - HYP-2170   # UD graph = Cay(Z[zeta6],U6); the "3" = kissing/2 = Eisenstein rosette
---

# THM-421: The common-neighbour floor for beating 3N unit distances

## Setup

For a finite planar point set `P` (|P| = N), let `U(P)` be the number of
**unit distances**: pairs `{x,y} ⊆ P` with `‖x−y‖ = 1`. The unit-distance graph
has these pairs as edges. S702/S703 established that, unlike the **penny graph**
(min-distance ≥ 1, capped by Harborth `⌊3n−√(12n−3)⌋ < 3n` via the kissing bound
κ ≤ 6), a *general* unit-distance graph **can exceed 3N** by using a non-minimal
radius (the additive norm-R layer, HYP-2262): the triangular (Eisenstein) lattice
at unit² = 7 has every interior point of degree 12 (density → 6 > 3).

Open handoff **S702(b):** *"is triangular N = 43 provably the 2D-optimal finite
construction [beating 3N]?"* HYP-2267 gave only a constructive **upper** bound (43)
**over lattice disk patches**. This theorem supplies the missing rigorous **lower**
bound over **all** finite planar sets, and sharpens the upper bound to 39.

## Statement

**(A) Lower bound [PROVED].** Every finite planar point set with `U(P) > 3N`
satisfies `N ≥ 17`. More generally, in *any* graph in which every two vertices have
**at most 2 common neighbours**, average degree `> D` forces

```
        N  ≥  C(D,2) + 2 .
```

For the unit-distance threshold `U > 3N` (average degree `> 6 = κ`) this gives
`N ≥ C(6,2)+2 = 17`. **No configuration of 16 or fewer points beats 3N.**

**(B) Upper bound [VERIFIED].** Two improvements over HYP-2267's `N = 43` (a √7
Eisenstein disk centred at a lattice point), all in the same unit² = 7 Eisenstein
graph (layer `r_Q(7)=12`):
 - **best disk:** the √7 disk centred at the edge-midpoint `(½,0)` gives `N = 39`,
   `U = 118 > 117` (off-lattice centre balances the boundary deficit);
 - **best overall:** a 32-point **non-disk** subset gives `U = 97 > 96 = 3·32`
   (found by multistart anneal/shrink-repair, re-verified by independent exact
   integer recount; vertex list in `unit_distance_3n_minsearch_s709e.out`).

**(C) Consequence.** The true minimum `N*` over *all* planar point sets satisfies
`17 ≤ N* ≤ 32`.

## Proof of (A) [PROVED]

Two facts about the planar unit-distance graph `G`:

- **(CN) ≤ 2 common neighbours per pair.** A common neighbour of `x ≠ y` lies on
  both unit circles centred at `x` and at `y`; two distinct equal-radius circles
  meet in ≤ 2 points. (This is exactly the property behind the Kővári–Sós–Turán /
  Erdős `O(N^{3/2})` *upper* bound on unit distances — see the reflection.)

Count **cherries** (paths of length 2, i.e. pairs of edges sharing a centre):
the number of cherries is `Σ_v C(d_v, 2)` (choose the centre `v`, then 2 of its
neighbours). Each *unordered endpoint pair* `{x,y}` is the endpoint set of exactly
(`# common neighbours of x,y`) ≤ 2 cherries. Summing over the `C(N,2)` pairs:

```
        Σ_v C(d_v, 2)  ≤  2·C(N,2)  =  N(N−1).                         (★)
```

Now maximise `U = ½ Σ_v d_v` subject to (★). For fixed `Σ d_v`, the quantity
`Σ d_v(d_v−1) = Σ d_v² − Σ d_v` is minimised when the degrees are equal
(convexity / Cauchy–Schwarz), so the constraint (★) is loosest — hence `U` is
largest — at the regular degree `d`. Equal degrees: `N·d(d−1)/2 ≤ N(N−1)`, i.e.
`d(d−1) ≤ 2(N−1)`, i.e.

```
        d  ≤  (1 + √(8N−7)) / 2 ,        U ≤ N(1 + √(8N−7)) / 4 .       (★★)
```

`U > 3N` requires average degree `d > 6`, i.e. `√(8N−7) > 11`, i.e. `8N−7 > 121`,
i.e. `N > 16`. Hence `N ≥ 17`. The general-D statement is identical:
`d > D ⟺ √(8N−7) > 2D−1 ⟺ N > (D²−D+2)/2 = C(D,2)+1 ⟺ N ≥ C(D,2)+2`. ∎

*Exact integer check (`unit_distance_3n_floor_s709.py`):* maximising `U` over
**integer** degree sequences under (★) gives max-U(N) = …,43,**48**,**52**,… for
N = 15,16,17; first `> 3N` exactly at N = 17 (max-U(16)=48=3·16, max-U(17)=52>51).

## Proof of (B) [VERIFIED, exact integer arithmetic]

Eisenstein indices `(x,y)` with squared-distance form `Q(x,y)=x²+xy+y²`,
unit² = 7 (layer `r_Q(7)=12`). Order the integer box by squared distance to the
**off-lattice centre `(½,0)`** and take the first 39 points. Counting pairs at
`Q(Δ)=7` (exact, via the 12-element offset shell) gives `U = 118 > 117`. The
edge-midpoint centre balances the boundary deficit better than a lattice-point
centre, shaving 4 points off HYP-2267's disk. A fine center sweep (1/6 grid) and a
boundary add/remove + swap local search find nothing below 39
(`unit_distance_3n_{search_s709b, verify_s709c}.py`). ∎

## Local structure (why small balls fail; `unit_distance_3n_anneal_s709d.py`)

In the √7 Eisenstein graph **every adjacent pair has *exactly* 2 common neighbours**
(saturating the CN bound — each edge sits in 2 triangles), yet the 12 neighbours of
a vertex span only **12 internal edges**, so a "flower" (centre + 12 neighbours,
N=13) has just `U = 24`, average degree `3.69 < 6`. Density must be *built up over a
chunk*, not found in a single ball — which is why the minimiser is a 32-vertex blob,
not a tiny cluster, and why the floor (17) is not met.

## Why the gap 17 vs 32 (interpretation)

The floor `C(κ,2)+2 = 17` is a pure **design-theoretic** bound: it would be met by a
near-regular graph in which essentially *every* pair realises its full 2 common
neighbours (a 2-(N,…) design flavour). Euclidean geometry cannot pack common
neighbours that efficiently — most pairs of a planar set are too far apart to share
*any* unit-neighbour — so realizability inflates the floor up to the high-30s. The
gap **17 → 32 is the cost of Euclidean realizability**. The kissing cap `κ = 6`
appears in BOTH ends: the threshold to beat is `(κ/2)·N = 3N`, and the floor is
`C(κ,2)+2`.

## Open (handed to HYP-2285)

Exact value of `N*` in `[17, 32]`; whether a non-lattice config beats 39; whether
the floor can be pushed above 17 using the additional **K₄-free** property
(no 4 mutually unit-distant points → max clique 3) together with (CN).

## Relation to canon

Consistent with **THM-412** (density quantization): triangular is forced to density
6 and reaches it at the smallest popular norm D*=7 — that is *why* the best
construction is the √7 layer. Does **not** contradict HYP-2267's proved content;
it **sharpens** HYP-2267's unproven "N=43 minimizer" claim (43 → 39, and adds the
first rigorous floor). Credit S702 (HYP-2262) and S703 (HYP-2267, THM-412).
