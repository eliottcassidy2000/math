---
id: HYP-3472
title: LRC14 dead-cover boundary-current audit
status: EVIDENCE / exact boundary-current audit; not a proof
source: codex-2026-06-29 continuation of HYP-3451/HYP-3453/HYP-3471 after the HYP-3462 AP84 corridor splice
tangent: T1432
technique: LTI-432
tournament_technique: LTT-332
script: 04-computation/lrc14_dead_cover_boundary_current_codex_20260629.py
result: 05-knowledge/results/lrc14_dead_cover_boundary_current_codex_20260629.out
reflection: 07-reflections/lrc14-dead-cover-boundary-current-codex-20260629.md
related:
  - HYP-3472
  - HYP-3473
  - HYP-3471
  - HYP-3462
  - HYP-3470
  - HYP-3461
  - HYP-3460
  - HYP-3459
  - HYP-3458
  - HYP-3455
  - HYP-3453
  - HYP-3451
  - HYP-3450
  - HYP-3438
  - HYP-3437
  - HYP-3436
  - HYP-3417
  - HYP-3129
  - THM-523
  - OPEN-Q-108
---

# HYP-3472: LRC14 Dead-Cover Boundary-Current Audit

This is the graph/current sibling of HYP-3471. HYP-3471 showed that every
audited row with a dead component has a rank-`<=2` E/branch survivor gate.
HYP-3451 projected dead components into a branch-coloured blocker graph.
HYP-3472 checks whether those E/branch gates touch the dead-cover projection as
small boundary-current cuts rather than merely existing in the same row.

Known input:

```text
HYP-3451: dead components have a blocker projection graph with low-rank escape data.
HYP-3453: dead_components(row)>0 implies a rank<=2 survivor gate on the 135-row bank.
HYP-3471: dead_components(row)>0 implies a rank<=2 E/branch survivor gate on the same bank.
HYP-3462/HYP-3470: AP84 sidecars are closed enough to serve as named AP packets.
```

## Exact Readout

Script:

```text
04-computation/lrc14_dead_cover_boundary_current_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_dead_cover_boundary_current_codex_20260629.out
```

Finite-bank readout:

```text
rows_audited=135
rows_with_dead_components=130/135
low_rank_e_branch_gates_total=7002
dead_row_low_rank_e_branch_gates=6928
nondead_row_low_rank_e_branch_gates=74
dead_rows_with_e_branch_gate=130/130
dead_rows_without_e_branch_gate=[]

e_branch_gates_touching_dead_projection=6834/6928
e_branch_gates_with_projection_edge_cut=5816/6928
e_branch_gates_with_separating_current=3266/6928

dead_rows_with_touching_gate=130/130
dead_rows_without_touching_gate=[]
dead_rows_with_projection_edge_cut_gate=123/130
dead_rows_without_projection_edge_cut_gate=[
  random_covering_001,
  random_covering_031,
  random_covering_039,
  random_covering_062,
  random_covering_074,
  random_covering_086,
  random_covering_101
]
dead_rows_with_separating_current_gate=121/130
dead_rows_without_separating_current_gate=[
  covering_AP_with_84,
  ap_omit_12_tail_84x01,
  random_covering_001,
  random_covering_031,
  random_covering_039,
  random_covering_062,
  random_covering_074,
  random_covering_086,
  random_covering_101
]
```

AP84/non-AP split:

```text
ap84_dead_rows=13
non_ap_dead_rows=117
ap84_rows_with_projection_edge_cut_gate=13/13
non_ap_rows_with_projection_edge_cut_gate=110/117
ap84_rows_with_separating_current_gate=11/13
non_ap_rows_with_separating_current_gate=110/117
```

The dominant best-per-row current is balanced:

```text
best_current_hist={(1,1): 71, (2,1): 17, (1,2): 15, (2,2): 15, ...}
best_endpoint_kind_hist={B1|E: 35, E|B1: 35, E|B0: 30, B0|E: 30}
```

Tournament Analysis over proof carriers, not runners:

```text
axes=(predicate_retention, dead_projection_payload, cut_strength,
      current_payload, ap_packet_splice, random_gluing_sidecar,
      scalar_firewall)
score_hist={7:1,14:1,52:1,54:1,55:1,63:1,64:1,66:1}
directed_3cycles=0
hamiltonian_path_count=1
hamiltonian_path=
  B00_projection_cut_gate
  -> B01_separating_boundary_current
  -> B02_dead_positive_e_branch_implication
  -> B03_closed_ap84_packet
  -> B04_random031_seven_owner_clause
  -> B05_typed_gate_word
  -> B06_raw_gate_count
  -> B07_raw_dead_fraction
```

## Proof Pull

The universal HYP-3471 implication survives the graph-current test at the
touching level:

```text
dead_components(row)>0
  => rank<=2 E/branch gate touching the dead-cover projection
```

That is stronger than row-local gate existence, and it is exact on the current
bank. HYP-3473 is now the Lean-facing formal interface for the E/branch
producer obligation; HYP-3472 should be added as the dead-cover
projection-touch/cut sidecar after `DeadCoverEBranchSoundness` is instantiated.
The stricter cut claims are not universal:

- Projection edge cut fails only on seven non-AP rows:
  `random_covering_001`, `random_covering_031`, `random_covering_039`,
  `random_covering_062`, `random_covering_074`, `random_covering_086`, and
  `random_covering_101`.
- Separating current fails on those seven rows plus the AP84 base pair
  `covering_AP_with_84` and `ap_omit_12_tail_84x01`; both AP rows still have
  projection edge cuts with `20` removed edges, but no component split.
- `random_covering_031`, the HYP-3455 seven-owner gluing row, is correctly
  placed in the edge-cut exception set. It has an E/branch gate touch but no
  projection edge cut, so it still needs the HYP-3455 gluing clause or another
  owner-current/state-lift sidecar.

Candidate finite lemma ladder:

1. Prove the universal touch lemma:
   `dead row => low-rank E/branch gate touches the dead-cover projection`.
2. Prove projection-edge transfer outside the seven named non-AP exceptions.
3. Prove separating-current transfer outside the seven non-AP exceptions and
   the two AP84 base rows.
4. Discharge the exceptions by HYP-3455, HYP-3462/HYP-3470, owner-current,
   two-adic descent, signed-SPEC, or state-lift debt.

## Computation Definition

1. Join every HYP-3453 row to its HYP-3451 dead-cover blocker projection.
2. For each low-rank E/branch gate, record the adjacent blocker labels it
   exposes on the two bad sides.
3. Remove those labels from the dead-cover incidence graph and test whether the
   gate is a projection cut, a separating cut, or only a local gate touch.
4. Record the current vector of each gate as branch-labelled blocker flow:
   `B0` labels, `B1` labels, balance, removed projection edges, change in
   largest dead-cover component, and row exceptions.
5. Run Tournament Analysis over proof carriers, not runners: projection cut,
   separating current, low-rank E/branch implication, AP84 closed packet,
   random031 seven-owner clause, typed gate word, raw gate count, and raw dead
   fraction.

Assumption challenge:

- Alternate vertices considered: runners, gaps, fixed circle sections, section
  boundaries, wall-crossing events, residues, cover arcs, Fourier modes,
  matroid circuits, survivor gates, dead components, blocker labels,
  owner-current vectors, graph cuts, and proof obligations.
- Chosen vertices for this audit: dead-cover blocker labels, E/branch survivor
  gates, projection-cut witnesses, and proof obligations.
- Preserved LRC predicate: whether a dead-cover obstruction has a graph-current
  route to a low-rank E/branch survivor gate that can feed HYP-3451/HYP-3455.
- Destroyed information: full wall geometry, exact interval lengths, scalar
  color counts, and same-row gate multiplicity beyond the cut/current sidecar.

This is not an LRC14 proof. It is a finite audit meant to determine whether
the HYP-3471 implication can be promoted into a boundary-current transfer
lemma, or whether the remaining non-AP rows must be named as owner-current,
two-adic, signed-SPEC, or state-lift debt.
