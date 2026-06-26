---
id: HYP-3049
title: A000568 edge-perspective extension lift at the first rooted-count failure
status: EVIDENCE / exact first-failure lift and sidecar compression; not a proof
source: codex-2026-06-26-S213
tangent: T1131
script: 04-computation/a000568_edge_perspective_extension_codex_s213.py
result: 05-knowledge/results/a000568_edge_perspective_extension_codex_s213.out
related:
  - HYP-3048
  - HYP-3047
  - HYP-3043
  - HYP-3039
  - HYP-2120
  - HYP-2121
  - HYP-1824
  - HYP-1825
  - THM-381
  - THM-385
---

# HYP-3049: A000568 Edge-Perspective Extension Lift

## Claim

At the first shifted A000568/rooted-perspective failure,

```text
P(5)=48 < U(6)=56,
```

the missing coordinate is exactly exposed by the old-root/new-observer
incident word.

A rooted 5-tournament plus the new observer's incident word is equivalent to a
6-tournament with an ordered pair of marked vertices: old root first, new
observer second.  The S213 computation verifies the exact identity

```text
rooted 5-perspective + incident word states
  = exact ordered-pair perspectives on U(6)
  = 1408.
```

Forgetting the old/new endpoint role gives the directed-edge perspective:

```text
exact directed-edge perspectives on U(6) = 704
exact unordered-pair perspectives on U(6) = 704.
```

Thus directed-edge perspective is not just a metaphor.  It is the first
natural quotient of the fully retained extension state.  It is also already a
controlled-forgetting step, because it forgets which endpoint was the old root
and which endpoint was the new observer unless a sidecar records that role.

## Evidence

S213 generates all 56 six-vertex tournament classes by extending the 12
five-vertex classes with all `32` incident words:

```text
U(5)=12  P(5)=48  U(6)=56  gap=8
U(6) generated from U(5)+all incident words: 56
rooted 5-perspectives plus incident word states: 1408
exact ordered-pair perspectives on U(6): 1408
extension states equal ordered-pair perspectives: True
exact directed-edge perspectives on U(6): 704
```

The extension fibers are broad, confirming that the old failure was not an
orphan-class problem:

```text
target classes reached per rooted 5-parent hist:
{7: 1, 13: 12, 21: 15, 23: 15, 24: 5}

rooted 5-parents per target U(6) class hist:
{5: 3, 6: 3, 8: 2, 10: 1, 11: 5, 13: 2, 14: 4,
 15: 6, 16: 2, 18: 2, 19: 2, 20: 5, 21: 3,
 23: 10, 25: 2, 26: 2, 28: 2}
```

## Sector-Deck Compression

For an ordered pair `(r,o)`, classify every other vertex by the four-sector
key

```text
(r beats x, o beats x) in {0,1}^2.
```

The S213 deck of all ordered-pair sector signatures on each six-class gives:

```text
size    : individual_sigs= 35  unique_class_decks=55/56
internal: individual_sigs= 59  unique_class_decks=55/56
cross   : individual_sigs=362  unique_class_decks=56/56
full    : individual_sigs=422  unique_class_decks=56/56
```

So sector sizes and internal sector tournaments almost recover all six-classes,
but they miss exactly one pair.  That pair is a converse/chirality pair:

```text
rep=344 score=(2,2,2,3,3,3) c3=8 H=43 aut=1 self_converse=False
rep=345 score=(2,2,2,3,3,3) c3=8 H=43 aut=1 self_converse=False
canonical opposite(344)=345.
```

Cross-sector orientation separates it.  This is a clean toy model for the
controlled-forgetting repair: the hidden coordinate is the cross-sector
orientation, not another node-depth layer.

## LRC Translation

The safe proof-state ladder should be read as:

```text
raw A000568 class
-> rooted node perspective
-> incident word / ordered-pair extension
-> directed-edge sector deck
-> cross-sector orientation / chirality repair
-> endpoint-owner packet sheaf
-> proof-obligation automaton.
```

The directed edge is dualistic exactly as the prompt suggested, but the
old/new observer role is not optional in LRC.  A quotient may forget it only if
an endpoint-owner, threshold-arc, gap-pressure, or proof-obligation sidecar
retains or reconstructs it.

## Matrix-Atlas Link

This is the first concrete instantiation of HYP-3048's matrix-sidecar pull.
HYP-3048 asks for four-sector edge block matrices and sidecar observability
matrices after HYP-3047 showed node depth was insufficient.  S213 supplies the
first exact row for that matrix program:

```text
sector size/internal decks: 55/56
sector cross-orientation deck: 56/56
```

So the observability column that repairs the first A000568 rooted-count defect
is `cross_sector_orientation_word`.  In matrix language, the sidecar is a
four-block edge-sector matrix around the old-root/new-observer ordered pair,
and the one failed size/internal fiber is a converse/chirality pair.

## Tournament Analysis

S213 runs a carrier tournament whose vertices are sidecars and proof carriers,
not runners.  The observable rewards retained observer role, incident word,
cross-sector coupling, chirality, class separation, and LRC owner
compatibility, with proof cost breaking ties.

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3_cycles=0
scc_sizes=[1,1,1,1,1,1,1,1]
hamiltonian_paths=1
```

One Hamiltonian path:

```text
endpoint_owner_packet
-> sector_cross_deck
-> sector_full_deck
-> ordered_pair_extension
-> directed_edge_perspective
-> sector_size_deck
-> sector_internal_deck
-> rooted_5_node_cache
```

## Next Pulls

1. Fold `cross_sector_orientation_word` into HYP-3048's observability-matrix
   template as the first worked A000568 column.
2. Compare the sector-cross deck against HYP-1824/HYP-1825's eight-class
   chirality bridge.
3. Repeat the ordered-pair sector-deck audit at `n=7` using generated classes
   or canonical augmentation rather than full labelled enumeration.
4. Add `observer_endpoint_role`, `incident_sector_deck`, and
   `cross_sector_orientation_word` as candidate sidecars in LRC threshold
   packet experiments.
5. Test whether the converse-pair collision pattern is systematic for
   size/internal sector decks, or a first-failure accident.

## Assumption Challenge

The considered vertices were rooted nodes, ordered pairs, directed edges,
sector cells, converse/chirality pairs, endpoint owners, packet sidecars, and
proof obligations.  The selected Tournament Analysis vertices are carriers.

The preserved predicate is observer-source/incident-coupling data needed for a
safe-box hit.  The destroyed information is old/new endpoint role, labelled
runner identity, and full extension rows.  Those losses are acceptable only
with a sidecar or named residual debt.
