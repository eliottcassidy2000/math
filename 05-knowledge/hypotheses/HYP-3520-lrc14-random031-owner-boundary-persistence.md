---
id: HYP-3520
title: LRC14 random031 owner-boundary persistence
status: EVIDENCE / finite owner-cobordism and quotient-price certificate; not an LRC14 proof
source: codex-2026-06-29 owner-boundary persistence pass after HYP-3512/HYP-3494/HYP-3493/HYP-3510/HYP-3511
tangent: T1520
technique: LTI-520
tournament_technique: LTT-420
script: 04-computation/lrc14_random031_owner_boundary_persistence_codex_20260629.py
result: 05-knowledge/results/lrc14_random031_owner_boundary_persistence_codex_20260629.out
reflection: 07-reflections/lrc14-random031-owner-boundary-persistence-codex-20260629.md
related:
  - HYP-3512
  - HYP-3494
  - HYP-3493
  - HYP-3511
  - HYP-3510
  - HYP-3490
  - HYP-3486
  - HYP-3485
  - HYP-3483
  - HYP-3482
  - HYP-3481
  - HYP-3402
  - HYP-3034
  - OPEN-Q-108
---

# HYP-3520: Random031 Owner-Boundary Persistence

## Claim

`random_covering_031` has a finite owner-boundary persistence certificate at
the current exploratory level.  The hard pair should be read as a forbidden
seam, not as a wall for phase flow:

```text
hard_components = (43,54)
forbidden_seam_owners = (23,45,93,113,147,169,173)
pure_bypass_flow_owners = (23,93,113)
owner_boundary_debt = (45,147,169,173)
pure_bypass_cells = 12
```

The lower-delta bypass is a moving flow carrier on the same two hard
components.  The four missing labels are stationary owner-boundary charge.  The
proof object is therefore an owner-current cobordism:

```text
forbidden_seam_owner_word - pure_bypass_flow_owner_word
  = (45,147,169,173).
```

This does not prove the LRC14 random031 terminal clause.  It makes the
remaining terminal clause sharper: any proof may delete the max-delta seam only
if it keeps an equivalent sidecar for this boundary word.

## Exact Readout

Executable:

```text
04-computation/lrc14_random031_owner_boundary_persistence_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_random031_owner_boundary_persistence_codex_20260629.out
```

The script imports HYP-3486/HYP-3481 and recomputes:

```text
witness_cells=282
cell_class_hist={'bypass':12,'free_hole':40,'ordinary':230}
legal_component_count=79
hard_seam_gate_count=2, mirror_pair=True
bypass_gate_count=2, mirror_pair=True
pure_bypass_rel_id=R01
pure_bypass_nodes=((527,1)..(532,1),(679,0)..(684,0))
pure_bypass_branch_hist={0:6,1:6}
pure_bypass_hit_components=(43,54)
pure_bypass_endpoint_ranks=(2,)
pure_bypass_external_horizontal_ports=0
```

The forbidden seam gates are:

```text
component=43 mask=branch0 delta=7 owners=(23,45,93,113,147,169,173)
component=54 mask=branch1 delta=7 owners=(23,45,93,113,147,169,173)
```

The lower-delta bypass gates are:

```text
component=43 mask=branch1 delta=2 owners=(23,93,113)
component=54 mask=branch0 delta=2 owners=(23,93,113)
```

The owner boundary matrix is:

```text
owner  seam  bypass  debt  dead_island  global_flow  seam-minus-bypass
23       1      1      0       1            1              0
45       1      0      1       1            1              1
55       0      0      0       0            1              0
93       1      1      0       1            1              0
113      1      1      0       1            1              0
147      1      0      1       0            1              1
169      1      0      1       0            1              1
173      1      0      1       0            1              1
```

So the finite boundary signature is:

```text
45:+1,147:+1,169:+1,173:+1.
```

The ordinary/global flow owner union is
`(23,45,55,93,113,147,169,173)`.  This is deliberately not the owner-boundary
certificate: subtracting the bypass owners from the global flow union leaves an
extra contaminant `55`.  The certificate must be local to the hard seam or to
an equivalent owner-current/relative-H1 sidecar.

## Quotient-Price Matrix

The exact sidecars are:

```text
seam_and_bypass_owner_words
owner_current_cobordism_matrix
hard_component_owner_map
observer_cut_payload_orbit
symbolic_relative_H1_boundary_class
```

Each reconstructs `(45,147,169,173)`.

The failing quotients explain the proof risk:

```text
bypass_owner_word_only              misses (45,147,169,173)
delta_route_phase_blocks            misses (45,147,169,173)
dead_island_owner_union_minus_bypass reconstructs only (45)
global_flow_owner_union_minus_bypass adds contaminant (55)
component_pair_only                 misses all owner labels
raw_bypass_count_12                 misses all owner labels
raw_282_witness_count               misses all terminal owner fields
```

This is the coordinate-resurrection discipline from HYP-3512 made local: a
quotient is legal only after it names the first destroyed coordinate and its
repair sidecar.

## Reframe

The pure bypass is not a hidden wall and not a failed rank-2 route.  It is an
isolated local stalk on the same mirror-paired components as the seam:
`external_horizontal_ports=0`.  That makes the seam-only owners a boundary
sidecar rather than an ordinary-neighbor charge.

The `n*2` side is the address/phase carrier `u=2t`: `282` witnesses, `12`
bypass hits, branch split `6/6`.  The `n+2` side is the owner-boundary seam:
seven owners with four persistent seam-only labels.  The useful theorem shape
is a span lemma retaining both sides, not a scalar choice between them.

## Proof Pull

The next proof-facing lemma can now be stated more narrowly:

```text
If a legal random031 seam-complement packet has the HYP-3486/HYP-3510
phase-flow receiver, the HYP-3511 free-hole bracket discharge, and the
HYP-3520 owner-current boundary word, then the pure bypass terminal is
discharged by boundary persistence rather than by projection-current deletion.
```

The remaining formal gap is not "find more bypass witnesses"; it is to turn the
owner-current matrix into a Lean-facing terminal lemma compatible with the
HYP-3490 private-label firewall and the HYP-3493 relative seam-sheaf table.

## Tournament Analysis

Vertices are quotient sidecars/proof interfaces, not runners, arcs, or raw
witnesses.  Pairwise observable: exact boundary-debt reconstruction plus local
`12`-cell bypass retention.  Switch: higher owner-boundary score, with a fixed
sidecar-cost tie path.

Script readout:

```text
score_hist={11:1,12:2,16:1,26:1,31:2,91:1,92:2,93:1,94:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1,1,1,1]
hamiltonian_path_count_under_fixed_tie_path=1
hamiltonian_path=
  seam_and_bypass_owner_words
  -> owner_current_cobordism_matrix
  -> hard_component_owner_map
  -> observer_cut_payload_orbit
  -> symbolic_relative_H1_boundary_class
  -> bypass_owner_word_only
  -> delta_route_phase_blocks
  -> dead_island_owner_union_minus_bypass
  -> global_flow_owner_union_minus_bypass
  -> raw_bypass_count_12
  -> raw_282_witness_count
  -> component_pair_only
```

## Assumption Challenge

Considered vertices: runners, raw gates, arcs, u-fibers, fixed circle sections,
section boundaries, wall-crossing events, residues, cover arcs, free-hole
packets, owner rows, hard-component owner maps, relative-H1 boundary classes,
and proof obligations.  Chosen vertices are quotient sidecars.  This preserves
the LRC predicate needed here--pure-bypass discharge with exact
owner-boundary debt--and intentionally destroys raw runner order, raw phase
count, and ordinary flow owner noise only after recording the repair sidecar.
