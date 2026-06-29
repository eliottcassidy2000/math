---
id: HYP-3520
title: LRC14 random031 owner-boundary persistence
status: EVIDENCE / finite owner-cobordism plus all-component persistence certificate; not an LRC14 proof
source: codex-2026-06-29 synthesis of the owner-boundary certificate, HYP-3513 firewall Nerode context, and the HYP-3493/HYP-3511 persistence audit
tangent: T1520
technique: LTI-520
tournament_technique: LTT-420
script: 04-computation/lrc14_random031_owner_boundary_persistence_codex_20260629.py
result: 05-knowledge/results/lrc14_random031_owner_boundary_persistence_codex_20260629.out
reflection: 07-reflections/lrc14-random031-owner-boundary-persistence-codex-20260629.md
related:
  - HYP-3523
  - HYP-3522
  - HYP-3513
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
the current exploratory level.  The hard pair is not a wall for phase flow.  It
is a forbidden seam, and the complement carries the phase-flow receiver:

```text
hard_components = (43,54)
forbidden_seam_owners = (23,45,93,113,147,169,173)
pure_bypass_flow_owners = (23,93,113)
owner_boundary_debt = (45,147,169,173)
pure_bypass_cells = 12
```

The lower-delta bypass is a moving flow carrier on the same two hard
components.  The four missing labels are stationary owner-boundary charge.  The
local proof object is therefore an owner-current cobordism:

```text
forbidden_seam_owner_word - pure_bypass_flow_owner_word
  = (45,147,169,173).
```

The global proof object is a three-exit right boundary over HYP-3493/HYP-3511:

```text
rank-2 route
or free-hole bracket
or pure-bypass owner-boundary.
```

This does not prove the LRC14 random031 terminal clause.  It makes the
remaining terminal clause sharper: any proof may delete or quotient the
max-delta seam only if it keeps an equivalent sidecar for the seam owner word
or has already emitted one of the three named exits.

## Exact Readout

Computed by:

```text
04-computation/lrc14_random031_owner_boundary_persistence_codex_20260629.py
05-knowledge/results/lrc14_random031_owner_boundary_persistence_codex_20260629.out
```

The script imports HYP-3486/HYP-3481 for the seam-local certificate and
HYP-3493/HYP-3511 for the all-component persistence audit:

```text
witness_cells=282
cell_class_hist={'bypass':12,'free_hole':40,'ordinary':230}
legal_component_count=79
hard_seam_gate_count=2, mirror_pair=True
bypass_gate_count=2, mirror_pair=True
pure_bypass_rel_id=R01
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
45:+1,147:+1,169:+1,173:+1
```

The all-component persistence partition is exact:

```text
component_count=79
flow_class_hist={'free_hole':14,'pure_bypass':1,'rank2_routed':64}
persistence_class_hist={'free_hole_bracket_persistent':14,
                        'pure_bypass_owner_boundary':1,
                        'rank2_owner_persistent':64}
unresolved_owner_boundary_components=[]
```

Free-hole full seam debt is no longer unexplained:

```text
full_seam_debt_components=14 bracketed free-hole packets
ordinary_bracketed_free_packets=10
half_open_free_packets=4
half_open_clusters=[(2,3),(8,13)]
half_open_cluster_rel_ids=[('R77','R76'),('R53','R05')]
half_open_cluster_boundary_components=[(4,7,90,93),(35,38,59,62)]
```

The only seam-only boundary debtor is the pure bypass:

```text
seam_only_boundary_debt_components=['R01']
R01 flow charge=(23,93,113)
R01 boundary charge=(45,147,169,173)
```

After the free-hole caps are bracketed separately, all nonfree components with
seam-owner support lie in one overlap carrier:

```text
owner_overlap_vertex_count=65
owner_overlap_component_count=1
owner_overlap_component_sizes=[65]
```

## Quotient Guardrails

The seam-local exact sidecars are:

```text
seam_and_bypass_owner_words
owner_current_cobordism_matrix
hard_component_owner_map
observer_cut_payload_orbit
symbolic_relative_H1_boundary_class
```

Each reconstructs `(45,147,169,173)`.

The failing seam-local quotients explain the proof risk:

```text
bypass_owner_word_only              misses (45,147,169,173)
delta_route_phase_blocks            misses (45,147,169,173)
dead_island_owner_union_minus_bypass reconstructs only (45)
global_flow_owner_union_minus_bypass adds contaminant (55)
component_pair_only                 misses all owner labels
raw_bypass_count_12                 misses all owner labels
raw_282_witness_count               misses all terminal owner fields
```

The all-component persistence guardrail agrees:

```text
clean_persistence_quotients =
  flow_class
  owner_boundary_persistence_class
  owner_presence_word
  seam_debt_word

lossy_persistence_quotients =
  component_size
  endpoint_rank_word
  mirror_closed_shadow
  raw_owner_count
  seam_owner_count
```

Raw owner count `3`, seam-owner count `3`, and endpoint rank `(2,)` all mix
the pure bypass with ordinary rank-2 discharge.  Mirror closure is useless as a
terminal quotient because all `79` legal components are mirror closed.  The
legal owner object is the seam word, not the scalar owner count.

## Seam-Sheaf Compression Canary

The script also imports the HYP-3493 relative seam-sheaf table and tests the
owner-boundary target against local discharge/free-hole context.  This is a
different compression question from the owner-current matrix: it asks which
component-level summaries accidentally identify the pure bypass stalk `R01`
with ordinary discharged components.

Exact readout:

```text
sheaf_component_count=79
sheaf_flow_class_hist={'free_hole':14,'pure_bypass':1,'rank2_routed':64}
sheaf_pure_bypass_rel_id=R01
sheaf_pure_bypass_owner_word=(23,93,113)
owner_boundary_persistence_word=PDPPOOO
free_hole_owner_union_hist={():14}
```

Safe component quotients:

```text
flow_class
allowed_exit
owner_union
sheet_pgf_bucket
```

Unsafe component quotients:

```text
owner_union_size: mixed at key 3
endpoint_ranks: mixed at key (2,)
branch_hist: mixed at key ((0,6),(1,6))
size: mixed at key 12
mirror_closed: mixed at key True
```

So the owner-current certificate and the seam-sheaf canary agree on the same
rule: the proof may compress the random031 terminal packet only after keeping
an owner word, pure-bypass exit label, legal sheet-PGF bucket, or equivalent
named repair sidecar.

HYP-3522 refines this owner-boundary word rather than replacing it.  Its
filtration keeps the transport word `(23,93,113)`, uses adjacent ordinary
branch-boundary cells to lift `(147,169)`, and leaves a residual pair
`(45,173)`.  HYP-3520 is therefore the persistence/quotient guardrail for the
whole boundary word; HYP-3522 is the next lemma shape for splitting that word.
HYP-3523 then packages the split as a bounded owner-carry transducer: stable
emissions `(23,93,113)` and `(147,169)`, residual carry `(45,173)`, and no
permission for scalar shadows to emit terminal proof digits.

## Reframe

The mirror-punctured cylinder should be read as two coupled but distinct
coordinates:

```text
n*2 phase flow in the seam complement
+ n+2 owner boundary on the deleted seam
=> rank-2 route or free-hole bracket or pure-bypass owner-boundary.
```

The pure bypass is an isolated local stalk on the same mirror-paired components
as the seam: `external_horizontal_ports=0`.  That makes the seam-only owners a
boundary sidecar rather than an ordinary-neighbor charge.

The free-hole packets are a cap system.  They carry full seam debt because they
have no owner support, but HYP-3511 brackets every cap by ordinary rank-2
boundary components.  The two half-open doublets are the only place where the
puncture picture remains visible.

HYP-3513 now supplies the firewall-side warning: existing colored axes decide
the private-firewall bit, but the full route still needs a sidecar.  HYP-3520
identifies the random031 route sidecar that cannot be compressed away: the
owner-current/seam word.

## Proof Pull

The next proof-facing lemma can be stated narrowly:

```text
If a legal random031 seam-complement packet has the HYP-3486/HYP-3510
phase-flow receiver, the HYP-3511 free-hole bracket discharge, and the
HYP-3520 owner-current boundary word, then the pure bypass terminal is
discharged by boundary persistence rather than by projection-current deletion.
```

The formal target is a finite right-boundary theorem:

1. Rank-2 components discharge by endpoint-rank-2 route.
2. Free-hole components discharge by HYP-3511 bracketed caps.
3. The pure bypass component discharges as flow owners `(23,93,113)` plus
   boundary debt `(45,147,169,173)`.

Any Lean, relative-H1, PGF, or owner-current step that uses owner count,
endpoint rank, component size, or mirror closure must carry the seam owner word
as a sidecar, or it mixes terminal classes.

## Tournament Analysis

Two compatible tournaments are now recorded.

The seam-local sidecar tournament uses quotient sidecars/proof interfaces as
vertices.  Pairwise observable: exact boundary-debt reconstruction plus local
`12`-cell bypass retention.  Switch: higher owner-boundary score, with a fixed
sidecar-cost tie path.

```text
score_hist={11:1,12:2,16:1,26:1,31:2,91:1,92:2,93:1,94:1}
directed_3cycles=0
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

The persistence tournament uses owner-boundary proof carriers as vertices.
Pairwise observable: terminal partition retention, seam owner word,
bypass/free-hole localization, and quotient legality.

```text
score_hist={17:1,28:1,29:1,58:1,61:1,78:1,98:1}
directed_3cycles=0
hamiltonian_path=
  owner_boundary_persistence_word
  -> seam_debt_coboundary_word
  -> free_hole_bracketed_debt
  -> pure_bypass_flow_charge
  -> owner_overlap_support_component
  -> endpoint_rank_shadow
  -> raw_owner_count_shadow
```

## Assumption Challenge

Considered vertices: runners, gaps, raw gates, arcs, u-fibers, fixed circle
sections, section boundaries, wall-crossing events, residues, cover arcs,
Fourier modes, dead islands, free-hole packets, owner rows, hard-component
owner maps, relative-H1 boundary classes, and proof obligations.

Chosen vertices: quotient sidecars and owner-boundary persistence carriers.
This preserves the LRC predicate needed here, pure-bypass discharge with exact
owner-boundary debt and a three-exit right boundary.  It intentionally destroys
raw runner order, raw phase counts, scalar owner counts, and ordinary flow
owner noise only after recording the repair sidecar.
