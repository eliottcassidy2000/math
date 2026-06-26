---
id: HYP-3047
title: A000568 k-depth perspective ladder and incident carrier lift
status: SYNTHESIS / small exact tournament audit and proof-carrier index; not a proof
source: codex-2026-06-26-S211
tangent: T1129
script: 04-computation/a000568_k_depth_perspective_ladder_codex_s211.py
result: 05-knowledge/results/a000568_k_depth_perspective_ladder_codex_s211.out
related:
  - HYP-2120
  - HYP-2121
  - HYP-3046
  - HYP-3043
  - HYP-3042
  - HYP-3039
  - HYP-3018
  - HYP-3015
  - HYP-1824
  - HYP-1825
  - THM-381
  - THM-385
  - T1128
  - OPEN-Q-108
---

# HYP-3047: A000568 k-depth perspective ladder and incident carrier lift

## Claim

The old A000568/rooted-perspective count coincidence is the first small model
of the controlled-forgetting ladder.

Let `U(n)` be the number of unlabelled tournaments on `n` vertices and let
`P(m)` be the number of exact node perspectives, i.e. rooted tournament classes,
on all `m`-vertex tournament classes.  The prompt's shifted comparison is

```text
A000568(n) = U(n)
node perspectives on (n-1)-classes = P(n-1).
```

The first failure is therefore shifted `n=6`, equivalently root-size `m=5`:

```text
P(5) = 48, while U(6) = 56.
```

The k-depth node ladder explains why this is a controlled-forgetting alarm
rather than a request for deeper rooted-node memory.  At `m=5`, the iterative
node-depth refinement reaches the full exact rooted count by depth `2`:

```text
depth counts = [5, 41, 48, 48, 48].
```

But exact rooted nodes are still eight short of `U(6)`.  Thus the missing
payload is not another node-neighborhood layer.  It is the incident
word/cross-coupling coordinate needed to add the observer and distinguish the
six-vertex unrooted extension classes.

For LRC, this repeats the HYP-3039/HYP-3043 lesson:

```text
raw A000568 class
-> rooted/node perspective
-> k-depth rooted cache
-> directed-edge sector carrier
-> cycle/clique/conflict carrier
-> endpoint-owner / gap-pressure packet sheaf
-> proof-obligation automaton
```

The legal forgetful map is not "forget labels until only a scalar remains."
It is "forget a coordinate only after a sidecar retains it, reconstructs it,
annihilates it by a dual certificate, descends it, or routes it to named
residual debt."

## Evidence

The S211 script recomputes the exact small counts through the first failure:

```text
m  shifted_n  U(m)  P(m)  U(m+1)  U(m+1)-P(m)  source_roots  U(m-1)
1      2         1     1        1              0             1       1
2      3         1     2        2              0             1       1
3      4         2     4        4              0             1       1
4      5         4    12       12              0             2       2
5      6        12    48       56              8             4       4
```

This also reaffirms HYP-2120: the source slice is exact and gapless.  Source
roots on `m` vertices delete to `U(m-1)`.  That is the tournament toy model for
why LRC can use observer-source semantics safely, while arbitrary rooted
perspectives require sidecars.

The k-depth node ladder, where depth `k+1` sees in/out multisets of depth-`k`
node types, gives:

```text
m=3: [3, 4, 4, 4, 4], exact_rooted=4,  U(m+1)-exact=0
m=4: [4, 12, 12, 12, 12], exact_rooted=12, U(m+1)-exact=0
m=5: [5, 41, 48, 48, 48], exact_rooted=48, U(m+1)-exact=8
```

So depth helps reconstruct rooted perspectives but does not bridge the shifted
A000568 defect.  The defect is above rooted-node memory.

The first incident lifts have enough room to hold the missing data:

```text
m=5:
  exact directed-edge perspectives = 88
  exact directed-cycle perspectives = 24
  exact transitive-chain perspectives = 64
  shifted target U(6) = 56
```

This is not a bijection claim.  It is a capacity and carrier-selection signal:
directed edges and transitive cliques retain incident/cut data that rooted
nodes have forgotten.  Directed cycles retain chirality/conflict data that may
explain the old eight-class bridge to HYP-1824/HYP-1825.

## Burnside Rootless-Sidecar Check

A standalone companion audit,
`04-computation/tournament_perspective_ladder_codex_s211.py`, recomputes
`U(n)` and exact rooted counts by Burnside through `n=10`.  It gives the same
shift failure:

```text
R(5)=48
U(6)=56
defect=8
```

The nonzero Burnside terms for `U(6)` include:

```text
cycle_type=[5,1]         fixed_tournaments=8    fixed_vertices=1
cycle_type=[3,3]         fixed_tournaments=32   fixed_vertices=0
cycle_type=[3,1,1,1]     fixed_tournaments=128  fixed_vertices=3
cycle_type=[1,1,1,1,1,1] fixed_tournaments=32768 fixed_vertices=6
```

The `[3,3]` term is the rootless warning: two rotating triples can be fixed
with no fixed vertex.  A rooted node perspective cannot name that coordinate,
even after the k-depth ladder has recovered every rooted node type.  This is
the combinatorial analogue of the S210/T1128 matrix-atlas rule: a hidden mode
that is not observable from the chosen quotient must survive as a kernel,
observability, or Schur-complement sidecar rather than as another scalar.

## Proof-Carrier Index

1. **k-depth node ladder.**  Use as a controlled local cache.  It diagnoses
   when rooted type has been recovered.  Stop deepening once exact rooted type
   is reached unless a separate sidecar is added.
2. **Directed-edge perspective.**  For an edge `tail -> tip`, every other
   vertex lies in one of four sectors:
   `(tail beats x, tip beats x)`.  The cross-sector orientation matrix is the
   missing incident/coupling payload.  This is the most direct next lift.
3. **Directed-cycle conflict perspective.**  A rooted directed 3-cycle keeps
   chirality and frustration.  Use this when edge sectors fail to distinguish
   cyclic debt.
4. **Transitive-clique insertion perspective.**  A rooted chain classifies
   other vertices by cut position.  This models the add-a-new-observer problem
   more naturally than a lone root does.
5. **Edge-cycle incidence conflict.**  Use a bipartite carrier whose vertices
   are directed edges and cyclic witnesses.  This can test whether local
   sector data and chirality agree.
6. **Endpoint-owner packet sheaf.**  The LRC-facing lift: source predicate,
   incident threshold word, endpoint owners, gap-pressure fibers, route labels,
   and proof-obligation states.

## Tournament Analysis

S211 runs a qualitative Tournament Analysis over proof carriers, not runners.

Pairwise observable:

```text
retained source / incident / pair / cycle / insertion / owner / automaton payload
minus proof cost.
```

Switch:

```text
orient toward the carrier retaining more missing observer-coupling payload;
name order breaks exact ties.
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3_cycles=0
scc_sizes=[1,1,1,1,1,1,1,1]
hamiltonian_paths=1
```

One Hamiltonian path:

```text
endpoint_owner_packet_sheaf
-> edge_cycle_incidence_conflict
-> transitive_clique_insertion
-> directed_edge_perspective
-> directed_cycle_conflict
-> k_depth_node_ladder
-> exact_rooted_node
-> raw_A000568_class
```

## LRC Translation

Do not use raw A000568 classes or exact rooted nodes as final proof states.
They are cache layers.  The proof state must live at least at one of:

```text
observer-source threshold arcs,
incident threshold words,
directed-edge sector cross-couplings,
cycle/chirality conflict cells,
transitive-clique insertion cuts,
endpoint-owner / gap-pressure packets,
proof-obligation automata.
```

This explains how the current controlled-forgetting stack should absorb the
old perspective curiosity: HYP-2120 supplies the exact source slice; HYP-2121
supplies the observer-coupling defect; HYP-3043 supplies the lens map; and
HYP-3047 says the next lift after rooted node memory is incident/coupling
carriers, not deeper node neighborhoods.

## Incoming Matrix-Atlas Link

Incoming codex-S210 adds a tournament matrix atlas.  It connects directly to
this perspective ladder.  The directed-edge sector carrier should be emitted
both as a combinatorial sector deck and as a small block matrix indexed by the
four tip/tail sectors.  The directed-cycle conflict carrier is naturally read
through the skew sign matrix `S=A-A^T` and its signed walk/chirality traces.
The transitive-clique insertion carrier is a block-elimination or Schur
complement question: delete or insert the observer only after recording the
correction term.  Endpoint-owner packet sheaves are the LRC version of a
boundary/low-rank-update matrix with rows indexed by proof carriers rather
than runners.

So the next edge-perspective experiment should not emit only counts.  It
should emit sector matrices, skew-cycle traces, and low-rank update fields
that can be joined to endpoint-owner sidecars.

## Next Pulls

1. Build the exact extension map from `m=5` directed-edge perspectives to
   `U(6)` classes and isolate which sector-coupling patterns account for the
   eight-class defect.
2. Build an edge-cycle incidence graph over the `m=5` representatives and
   ask whether the defect splits by chirality, self-converse behavior, or the
   HYP-1824/HYP-1825 bridge.
3. Add a transitive-clique insertion-cut carrier to the LRC packet ledger:
   `observer_cut_position_word`, `incident_sector_deck`, and
   `cross_sector_orientation_word`.
4. Test the same perspective ladder on actual LRC threshold tournaments, not
   just abstract tournaments, with endpoint-owner labels retained.
5. Stress whether `depth=2` node color is enough to recover exact rooted type
   beyond `m=5`, and record any higher-depth collisions as guardrails rather
   than proof states.

## Assumption Challenge

Tournament vertices need not be runners, and in this pass runners were not the
right vertex set.  The considered vertices were tournament nodes, rooted
perspectives, directed edges, directed cycles, transitive clique insertion
cuts, edge-cycle conflicts, endpoint owners, gap-pressure fibers, automaton
states, and proof obligations.

The preserved LRC predicate is whether the carrier keeps observer-source and
incident-coupling data needed for a safe-box hit.  The destroyed data is exact
runner identity, full labelled extension rows, and scalar ordering.  That loss
is acceptable only with a retained/reconstructed/annihilated/descended sidecar
or named residual debt.
