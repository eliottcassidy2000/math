# Square-sum paths: growing prefixes and orientation guards

**Scope.** Let `G_n` be the simple graph on `1,...,n` with `{a,b}` an edge exactly when `a != b` and `a+b` is a square. A path is unoriented unless otherwise stated. The whole graph and its Hamilton paths are distinct from the user's A/B/C induced-component example; all counts below concern this precise prefix graph. The proofs in the insertion and distance sections hold for arbitrary finite simple graphs or the explicitly specified positive-integer graph. All fixed-order counts are finite exact computations, not eventual monotonicity or a universal spectrum claim.

## Exact insertion recursion

For a finite simple graph `G=(V,E)`, a nonempty `S⊆V`, and `x∈S`, write `F(S,x)` for the number of directed paths with vertex set exactly `S` and last vertex `x`. Thus `F({x},x)=1`, and otherwise

    F(S,x) = sum_{y∈S\{x}, {x,y}∈E} F(S\{x},y).

When adding a vertex `v` to an old graph on `W`, with old-neighbor set `N`, the number `H(W∪{v})` of directed Hamilton paths is exactly

    2 * sum_{x∈N} F(W,x)
    + sum_{A⊔B=W, A,B nonempty} (sum_{x∈A∩N} F(A,x))
                                    * (sum_{y∈B∩N} F(B,y)).

The second sum is over ordered partitions. To prove it, cut a directed path at `v`. If `v` is an endpoint, the old path's corresponding endpoint lies in `N`, and either end is possible. Otherwise its two nonempty sides form the ordered partition `(A,B)` with incident endpoints `x,y` in `N`. Reverse the second side to obtain its path ending at `y`. This is a bijection and there is no division by two. The formula proves why the previous total alone is insufficient: interior insertions require the entire two-path-cover endpoint interface. In square-sum prefixes, `N_n={k²−n: n<k²≤2n−1}`; distinctness and positivity are built in. This is not a tournament recurrence. The independent checker validates the identity on all 1,098 labelled graphs of orders 2 through 5.

## Three genuinely different orientation models

Suppose `n>1`, `G` has `m` actual edges and `U` unoriented Hamilton paths. In the bidirected graph each undirected edge is replaced by both arcs. Its directed count is deterministically `2U`. In the independent-orientation model every *actual edge* receives one independently chosen direction, producing `2^m` orientations, and

    E[H] = 2U / 2^(n−1).

This follows by summing the indicators of all `2U` directed paths, each fixing `n−1` edges. For a unique unoriented path, exactly `2 * 2^(m−n+1)` orientations contain a directed Hamilton path, and all have `H=1`; its reverse is incompatible with it. A graph edge absent from every Hamilton path is free in all these orientation choices. An edge forced in an undirected Hamilton path need not have a fixed orientation across different paths or orientations.

For any bank of directed paths, assign each path the arcs it requires. A set of path events is jointly feasible exactly when its union never requires both directions of the same actual edge. If the union fixes `s` edges, exactly `2^(m−s)` graph orientations realize that event intersection. Thus

    Z_G(z) = sum_{S compatible} 2^(m−|union of edges in S|) (z−1)^|S|.

The empty set contributes `2^m`. This inclusion–exclusion/binomial-moment expansion and the G25 polynomial were already proved in Session141: see pinned source `05-knowledge/results/session141_metric_aliquot_square_paths_20261102/square_paths/PROOF.md`, blob `3f98b8c8692ba70ee013bf226f5d20661357b18a` at repository commit `74e81ed37c9175d94e8fad57fc5809cf78cee572`. Our implementation independently recovers them; novelty is not claimed for that formula or census.

In a third model, independently orient *all* `N=binom(n,2)` pairs to obtain a full labelled tournament, but only count Hamilton paths using actual graph edges. Each graph orientation has `2^(N−m)` tournament extensions with the same **guarded** count. Therefore its generating polynomial is `2^(N−m) Z_G(z)`. The tournament also has many paths using nonedges; those are excluded and its ordinary Hamilton-path count is a different observable. No parity theorem for ordinary tournament Hamilton paths is transferred to the guarded count.

## Exact small-prefix cases

S147 proves that G15, G16, G17 are connected unicyclic and have exactly one undirected Hamilton path, and that G18 has none. In the first three, edge `{1,3}` is the unique edge absent from every Hamilton path, and every other graph edge is forced; each has `m=n`. Thus among their `2^n` actual-edge orientations precisely four contain a directed Hamilton path, each with exactly one. In G18 there is no such orientation. When there are no paths, saying every edge is forced is vacuous; we report no edge classification instead.

| n | m | U | Actual edges unused by all paths |
|---|---:|---:|---|
| 15,16,17 | n | 1 | `{1,3}` |
| 18 | 18 | 0 | no paths |
| 23 | 28 | 3 | `{1,15}`, `{3,6}`, `{3,13}`, `{4,5}` |
| 24 | 30 | 0 | no paths |
| 25 | 32 | 10 | `{1,15}` |
| 26 | 34 | 12 | `{1,3}` |
| 27 | 36 | 35 | none |

Exact guard orientation histograms are:

* G23: 268,435,072 with H=0; 384 with H=1.
* G25: 4,294,962,432 with H=0; 4,608 with H=1; 256 with H=2.
* G26: 17,179,857,952 with H=0; 10,240 with H=1; 928 with H=2; 64 with H=3.
* G27: 68,719,413,784 with H=0; 54,976 with H=1; 7,360 with H=2; 488 with H=3; 120 with H=4; 8 with H=5.

For G25, `N=300` and `m=32`; the full-tournament guarded polynomial is exactly

    2^268 * (4294962432 + 4608 z + 256 z²).

The factor counts free directions of the 268 *nonedges*, and does not convert `G25` into a tournament or count unrestricted tournament paths. Its first graph-orientation moment is `4608+2*256=5120=20*2^8`. Session141 also records the explicit directed-path compatibility graph `4 K2 ⊔ 12 K1` and the compatible pairs' 26-edge union.

The `orientation.py` path-enumeration and compatible-event calculation are checked against an independent endpoint subset DP through n=27; at n=23 and 25, every edge's incidence count is additionally checked through independent edge-deletion DP. The histogram routine is tested against direct orientation enumeration for all labelled graphs of orders 2–4 and K5. These checks independently reproduce prior Session141 values rather than claiming their discovery.

## The pair 8, 10 and the 3–4–5 labels

`8=9−1` and `10=9+1` are not adjacent since `8+10=18` is nonsquare, and cannot have a common neighbor: two squares cannot differ by 2 modulo 4. In every prefix G10 through G14, 8 has only neighbor 1 and 10 only neighbor 6. The path `8–1–3–6–10` exists, and no shorter path can connect their forced first/last neighbors, so their distance is exactly 4. From G15 onward, distance is exactly 3 by

    8 — 1 — 15 — 10     edge sums: 9, 16, 25 = 3², 4², 5².

The symmetric pair 0–9–16–0 has the same edge-sum labels `{9,16,25}`, but 0 lies outside the positive graph. This is an edge-label correspondence only: a three-edge path and triangle are not isomorphic as graphs. The edge 9–16 lies in the prior C group. The pair 8 and 10 are not the endpoints of the unique G15 Hamilton path: its endpoints are 8 and 9; 10 is internal.

A general positive-integer distance theorem: for every integer `t≥2`, put `x=t²−1` and `z=t²+1`. They are distinct positive nonneighbors, since `2t²` is never square. They have no common neighbor, again because their difference is 2 modulo 4. Their distance is exactly 3 via

    x — (8t²+1) — (8t²−1) — z,

with the three square sums `9t², 16t², 9t²`. These four positive vertices are distinct for `t≥2`. For t=3 the witness is 8–73–71–10; the smaller 8–1–15–10 has different square labels. This elementary construction gives an all-t existence bound, not shortest-label optimality.

## Replay

From this directory run `python orientation.py`, `python -O orientation.py`, `python test_independent.py`, `python -O test_independent.py`, `python test_geometry_and_edges.py`, and `python -O test_geometry_and_edges.py`. All checks use explicit exceptions, including under optimization. Geometry tests check prefix distances n=10,...,40 and the t-family through 100; the text gives the all-n and all-t proofs. No production repository files or publishing actions are needed to run these scripts.
