---
id: HYP-3478
title: LRC14 small-touch zero-edge geometry atlas
status: EVIDENCE / finite singleton-pocket atlas; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3477, HYP-3476 pair-current, HYP-3476 exception-router, HYP-3475, and S319 gate unit-delta split
tangent: T1438
technique: LTI-438
tournament_technique: LTT-338
script: 04-computation/lrc14_small_touch_geometry_atlas_codex_20260629.py
result: 05-knowledge/results/lrc14_small_touch_geometry_atlas_codex_20260629.out
reflection: 07-reflections/lrc14-small-touch-geometry-atlas-codex-20260629.md
related:
  - HYP-3477
  - HYP-3476
  - HYP-3475
  - HYP-3474
  - HYP-3473
  - HYP-3472
  - HYP-3471
  - HYP-3455
  - HYP-3453
  - HYP-3451
  - HYP-3438
  - THM-523
  - OPEN-Q-108
---

# HYP-3478: LRC14 Small-Touch Zero-Edge Geometry Atlas

## Claim

After HYP-3476 pair-current and HYP-3477 hard-mirror discharge, the six
small-touch/no-hard rows are not failed graph-current rows.  Their dead-cover
projections are edgeless because the dead components are isolated singleton
pockets: each dead component has a unique rank-`2` pair of blocker labels, one
`B0` owner and one `B1` owner, and every such pocket has an exact mirror
partner.

The six rows are:

```text
random_covering_001
random_covering_039
random_covering_062
random_covering_074
random_covering_086
random_covering_101
```

So the proof target should be a singleton-pocket theorem, not a larger
E/branch gate tuple search.

## Exact Readout

Executable scout:

```text
04-computation/lrc14_small_touch_geometry_atlas_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_small_touch_geometry_atlas_codex_20260629.out
```

Aggregate geometry:

```text
rows=6
dead_component_count_hist={2:5,4:1}
mirror_pair_count_hist={1:5,2:1}
projection_edge_hist={0:6}
edge_support_label_hist={0:6}
all_dead_components_singleton=True
all_hard_orbit_count_hist={0:6}
s319_min_gate_kind_hist={'branch_unit_delta':4,'delta_sidecar_packet':2}
```

The dead-pocket owner pairs are mirror-doubled:

```text
(165,179): 2
(81,153): 2
(63,129): 2
(9,81): 2
(15,99): 2
(133,169): 2
(7,175): 2
```

Every component has at least two complete E/branch gate touches:

```text
component_complete_gate_count_hist={2:2,3:6,4:2,6:2,9:2}
```

This means the obstruction is not lack of local gate contact.  It is that the
dead-cover projection has no shared owner label, hence no edge to cut.

## Row Grammar

Five rows have a single mirror pair of singleton dead pockets.  One row has two
mirror pairs:

```text
random_covering_001:
  (165,179) and (81,153)

random_covering_039:
  (63,129)

random_covering_062:
  (9,81)

random_covering_074:
  (15,99)

random_covering_086:
  (133,169)

random_covering_101:
  (7,175)
```

The incoming S319 gate unit-delta split separates the six:

```text
branch_unit_delta rows:
  random_covering_001
  random_covering_062
  random_covering_086
  random_covering_101

delta_sidecar rows:
  random_covering_039
  random_covering_074
```

Thus the cleanest next proof target is the four branch-unit singleton rows;
the two S319 delta-sidecar rows should carry their extra cover-delta payload.

## Best-Touch Refinement

A companion audit refines the branch-unit packet without replacing this atlas:

```text
04-computation/lrc14_small_touch_singleton_current_codex_20260629.py
05-knowledge/results/lrc14_small_touch_singleton_current_codex_20260629.out
```

It keeps the same zero-edge singleton-pocket carrier and then asks which
touching gate gives the smallest local current payload.  The split is:

```text
clean best-touch singleton current:
  random_covering_062
  random_covering_086
  random_covering_101

delta-sidecar singleton current:
  random_covering_039
  random_covering_074

asymmetric best-touch singleton current:
  random_covering_001
```

Aggregate refinement:

```text
best_delta_hist={(1,1):5,(2,1):1}
best_current_hist={(1,1):5,(2,1):1}
```

This does not contradict the S319 split.  `random_covering_001` is still a
branch-unit singleton row by the minimum gate-kind atlas, but it is not one of
the three clean best-touch rows: it has two mirror pairs, four isolated dead
components, and best touching current `(2,1)`.

## Proof Pull

Replace graph edge-current by a component-local lemma:

```text
mirror pair of rank-2 isolated dead pockets
+ unique owner labels
+ complete/touching E-branch gate sidecar
+ owner-pair residue/span word
=> terminal singleton-current discharge.
```

This is the shape HYP-3476 pair-current pointed toward but did not expose: a
zero-edge projection is not missing a cut; it is asking for a different
current carrier.

## Tournament Analysis

Vertices are geometric singleton-pocket proof carriers, not runners or raw row
names.

```text
pairwise_observable =
  predicate retention + component geometry + mirror payload
  + owner-pair payload + gate-touch payload + S319 split + scalar firewall
score_hist={7:1,38:1,55:1,59:2,62:1,65:1}
directed_3cycles=0
hamiltonian_path =
  G00_zero_edge_singleton_pocket_geometry
  -> G01_mirror_pair_dead_component_atlas
  -> G02_owner_pair_residue_span_word
  -> G03_gate_touch_complete_pocket_sidecar
  -> G04_s319_unit_delta_vs_delta_sidecar_split
  -> G05_single_best_gate_shadow
  -> G06_raw_row_name_exception_list
```

Assumption challenge: considered runners, gaps, fixed circle sections, section
boundaries, wall-crossing events, residues, cover arcs, dead-cover projection
nodes, singleton dead components, mirror component pairs, owner pairs,
E/branch gates, and terminal proof obligations.  The chosen carrier preserves
the singleton-current terminal-discharge predicate and destroys exact interval
location, mirror partner, owner-pair residue/span, S319 gate-kind class, and
complete gate touches if scalarized.
