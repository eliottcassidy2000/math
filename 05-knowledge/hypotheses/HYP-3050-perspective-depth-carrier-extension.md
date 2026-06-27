---
id: HYP-3050
title: Exact non-node perspective-depth carriers after the first A000568 defect
status: EVIDENCE / exact first-failure scout and proof-interface synthesis; not a proof
source: codex-2026-06-26-S214
tangent: T1132
script: 04-computation/perspective_depth_ladder_codex_s214.py
result: 05-knowledge/results/perspective_depth_ladder_codex_s214.out
related:
  - HYP-3049
  - HYP-3047
  - HYP-3048
  - HYP-3040
  - HYP-3039
  - HYP-2210
  - HYP-2120
  - HYP-1978
  - HYP-1977
  - THM-381
  - THM-385
  - THM-260
  - THM-409
  - OPEN-Q-108
---

# HYP-3050: Exact Non-Node Perspective-Depth Carriers After The First A000568 Defect

## Claim

Extending HYP-3047, the old perspective coincidence

```text
P(m) = A000568(m+1)
```

where `P(m)` is the sum of vertex-orbits over all tournament isomorphism
classes on `m` vertices, is a controlled-forgetting warning.  It first fails at
base size `m=5`, equivalently next-level index `n=6`:

```text
P(5)=48,  A000568(6)=56,  defect=8.
```

The new point is that this defect is not repaired by simply looking deeper
inside each existing node perspective.  A directed k-depth node perspective
already recovers exact rooted node orbits by depth `3` at `m=5`; the missing
`8` are observer-extension / cut-vector data, not deeper node-neighborhood
data.

The HYP-3049 edge-perspective lift verifies the exact ordered-pair extension
identity at `m=5 -> 6`; this hypothesis records the wider exact edge/triple
carrier table and local depth ladder around that window.  The concurrent
HYP-3048 matrix atlas says how to make these carriers proof-facing: edge
sectors, triple types, and observer cuts should be rows, columns, or
sidecar-observability coordinates before any scalar matrix invariant is
trusted.

## Computation

Script:

```text
04-computation/perspective_depth_ladder_codex_s214.py
```

Result:

```text
05-knowledge/results/perspective_depth_ladder_codex_s214.out
```

Exact shift table:

```text
m classes(m) P(m) A000568(m+1) defect
1 1          1    1            0
2 1          2    2            0
3 2          4    4            0
4 4          12   12           0
5 12         48   56           8
```

Depth table for directed node perspectives:

```text
m exact d1 d2 d3 d4
1 1     1  1  1  1
2 2     2  2  2  2
3 4     4  4  4  4
4 12    10 12 12 12
5 48    36 47 48 48
```

Thus:

```text
m=4 needs depth 2 to recover exact rooted node perspectives.
m=5 needs depth 3.
depth 3 recovers P(5)=48, but still not A000568(6)=56.
```

## Edge, Clique, Cycle, Conflict Carriers

Exact non-node rooted carrier totals through the first-failure window:

```text
m arc_orbits triple_orbits transitive_triples cyclic_triples conflict_pairs
1 0          0             0                  0              0
2 1          0             0                  0              0
3 4          2             1                  1              0
4 16         12            8                  4              0
5 88         88            64                 24             0
```

Local depth totals:

```text
m edge_d1 edge_d2 edge_d3 all_triple_d1 all_triple_d2 cyclic_d1 cyclic_d2
1 0       0       0       0             0             0         0
2 1       1       1       0             0             0         0
3 4       4       4       2             2             1         1
4 15      16      16      12            12            4         4
5 75      88      88      86            88            24        24
```

The edge carrier is naturally dualistic: an edge has a tail and tip, and every
outside vertex lies in one of four sectors:

```text
tail -> w, tip -> w
tail -> w, w -> tip
w -> tail, tip -> w
w -> tail, w -> tip
```

This is the tournament analogue of a two-plate capacitor or pair-good switch:
the edge perspective records how the rest of the structure sees both sides of
a directed obligation.

The triple carrier separates transitive-clique perspectives from directed
cycle perspectives.  At `m=5`, exact triple perspectives split as
`64` transitive and `24` cyclic.  The equality

```text
arc_orbits(m=5) = triple_orbits(m=5) = 88
```

is a useful first-failure coincidence, but it should be treated as a carrier
hint, not as a bijection.  Conflict-pair perspectives require two disjoint
cycles and therefore need the `m=6` extension.

## Controlled-Forgetting Lesson

The controlled-forgetting ladder from HYP-3039 says a quotient is legal only
after the next hidden coordinate is exposed, killed, cut, or routed to named
debt.  Perspective depth gives a tournament version of the same rule:

```text
unrooted class
-> k-depth node perspective
-> exact rooted node perspective
-> observer-extension / cut perspective
```

At `m=5`, the node ladder is already exact by depth `3`, but the next-level
observer-coupled class count still has eight extra states.  Those eight states
are the first place where "node perspective" and "add one observer" stop being
interchangeable.

For LRC this matches the older correction: the source-perspective slice is
exact by THM-381 because deleting a source is canonical, while the full
observer-coupled problem needs marked cut data.  Edge perspectives, cycle
perspectives, and conflict perspectives are sidecar candidates for the
non-source residuals.

## Tournament Analysis

Vertices are perspective carriers, not runners:

```text
unrooted_A000568_class
node_depth_1_score_view
node_depth_2_neighbor_view
node_depth_3_exact_root
source_perspective
observer_extension_cut
edge_tail_tip_perspective
triple_clique_perspective
cycle_omega_perspective
conflict_pair_perspective
```

Pairwise observable:

```text
observer_coupling_retention,
root_payload_retention,
source_deletion_legality,
edge_duality,
cycle_conflict_visibility,
controlled_forgetting_stage,
proof_cost
```

Switch/gauge:

Orient toward the carrier that preserves the LRC predicate and names the
least hidden debt.  Tie Hamiltonian path:

```text
observer_extension_cut
> source_perspective
> node_depth_3_exact_root
> edge_tail_tip_perspective
> triple_clique_perspective
> cycle_omega_perspective
> conflict_pair_perspective
> node_depth_2_neighbor_view
> node_depth_1_score_view
> unrooted_A000568_class
```

The synthesis tournament is transitive under this gauge:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
```

## Assumption Challenge

Alternate vertices considered: tournament classes, nodes, node neighborhoods,
directed edges, unordered edges, transitive triples, directed 3-cycles, cliques,
cycle-conflict pairs, observer cuts, source roots, sink roots, and LRC proof
obligations.

The chosen carriers are perspective stages because the preserved predicate is
observer-coupled LRC schedulability.  Pure A000568 classes destroy the observer;
exact rooted node perspectives still destroy the observer-extension cut; edge
and cycle perspectives destroy less pair/cycle-local data but cannot replace
source-marked LRC payload.

The challenged assumption is that a deeper node neighborhood alone repairs the
first rooted defect.  The `m=5` depth table says no: depth `3` recovers exact
node roots, and the remaining gap is a different coordinate.

## Next Pull

1. Extend the exact non-node carrier table to `m=6` with a faster class
   generator, so conflict-pair perspectives become visible.
2. Define an observer-extension cut perspective: root the deleted observer
   together with its in/out cut into the `m`-tournament.
3. Compare edge perspectives with pair-good blocker teeth and residual
   capacitor cuts; both are two-ended obligations with outside-sector data.
4. Compare cycle and conflict perspectives with OCF `Omega(T)` vertices and
   conflict edges.
