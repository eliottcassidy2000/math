---
id: HYP-3528
title: LRC14 random031 proof-contract one-red-clause audit
status: EVIDENCE / executable finite theorem-interface router; not an LRC14 proof
source: codex-2026-06-29 execution of the proof-contract router suggested by HYP-3525, using HYP-3521/HYP-3523/HYP-3524/HYP-3526
tangent: T1528
technique: LTI-528
tournament_technique: LTT-428
script: 04-computation/lrc14_random031_contract_router_one_red_clause_codex_20260629.py
result: 05-knowledge/results/lrc14_random031_contract_router_one_red_clause_codex_20260629.out
reflection: 07-reflections/lrc14-random031-contract-router-one-red-clause-codex-20260629.md
related:
  - HYP-3526
  - HYP-3525
  - HYP-3524
  - HYP-3523
  - HYP-3522
  - HYP-3521
  - HYP-3520
  - HYP-3513
  - HYP-3511
  - HYP-3510
  - HYP-3490
  - HYP-3486
  - OPEN-Q-108
---

# HYP-3528: LRC14 Random031 Proof-Contract One-Red-Clause Audit

## Claim

The random031 proof frontier is now best treated as a finite theorem contract.
The rows of the contract are proof clauses and sidecar obligations, not
runners, arcs, residues, or scalar counts.

The executable router imports HYP-3521, HYP-3523, HYP-3524, HYP-3525, and
HYP-3526 and reports one remaining open clause:

```text
open_clause_names=('residual_pair_close_tail',)
```

All other random031 terminal obligations are reduced to finite evidence plus a
named formalization target:

```text
ordinary_route_emit
free_hole_single_emit
free_hole_doublet_buffer_emit
bypass_transport_emit
bypass_bracket_lift_emit
private_firewall_route_sidecar
vertical_halfturn_guard
```

This does not prove LRC14.  It does make the remaining random031 theorem shape
precise: after transport `(23,93,113)` and branch-boundary lift `(147,169)`,
the only non-streamed owner tail is residual pair `(45,173)`, and HYP-3526
requires route sidecar `R` unless a separate route-reconstruction theorem is
proved.

## Exact Readout

Computed by:

```text
04-computation/lrc14_random031_contract_router_one_red_clause_codex_20260629.py
05-knowledge/results/lrc14_random031_contract_router_one_red_clause_codex_20260629.out
```

Stream imports:

```text
component_tokens=79
terminal_certificate_count=77
emitted_certificate_count=77
component_terminal_hist={
  ordinary_rank2_route_component:64,
  free_hole_single_bracket_packet:10,
  free_hole_doublet_packet:4,
  pure_bypass_owner_boundary_component:1
}
cell_terminal_hist={
  ordinary_rank2_route_component:230,
  free_hole_single_bracket_packet:26,
  free_hole_doublet_packet:14,
  pure_bypass_owner_boundary_component:12
}
max_pending_doublet_components=1
max_owner_carry_width=4
final_owner_carry=(45,173)
doublet_rank_gaps=(1,1)
owner_rank_gap=1
```

Owner carry geometry:

```text
seam_owners=(23,45,93,113,147,169,173)
transport_owners=(23,93,113)
residual_after_transport=(45,147,169,173)
bracket_lift_owners=(147,169)
residual_after_branch_boundary=(45,173)
bypass_rank=45
boundary_component_ranks=(44,46)
```

Route and vertical guardrails:

```text
private_status_by_I=True
private_status_by_Q=True
route_by_IQ=False
route_by_all_existing_plus_IQ=False
route_by_R=True
all_existing_plus_IQ_route_mixed=(1,7)
vertical_halfturn_present_cells=48/282
vertical_pair_class_counts={
  (free_hole,free_hole):2,
  (free_hole,ordinary):14,
  (ordinary,ordinary):8
}
vertical_glued_component_count=69
vertical_mixed_component_count=7
```

## Proof Pull

HYP-3528 proposes this Lean-facing skeleton:

```text
def Random031TerminalContract :=
  ordinary_route_emit
  + free_hole_single_emit
  + free_hole_doublet_buffer_emit
  + bypass_transport_emit
  + bypass_bracket_lift_emit
  + private_firewall_route_sidecar
  + vertical_halfturn_guard
  + residual_pair_close_tail

lemma random031_contract_streams_all_certificates:
  component_events=79 -> terminal_certificates=77

lemma random031_contract_requires_R:
  I/Q proves private cut but terminal route uses R

lemma random031_contract_final_tail:
  after transport and bracket lift, remaining owner tail is exactly (45,173)

remaining_goal=random031_residual_pair_close_tail
```

The next serious proof target is therefore not another scalar scout.  It is the
two-owner residual boundary lemma:

```text
transport + bracket_lift + R + no_hidden_tail_guard
  -> residual (45,173) cannot hide downstream
```

Residue buckets, centered residues, sliced-box volumes, owner counts, and raw
terminal counts are forbidden as replacements for that lemma.

## Tournament Analysis

Vertices are proof-contract clauses and sidecar obligations.

Pairwise observable:

```text
readiness - hidden_tail_risk - scalar_forgetting_risk
```

Switch/gauge: higher contract score; ties use readiness, lower risk, and
clause name.

Fingerprint:

```text
score_hist={-217:1,-7:1,249:1,269:1,274:1,293:1,359:1,370:1,400:1}
directed_3cycles=0
sccs=9 singleton SCCs
hamiltonian_path=
  ordinary_route_emit
  -> free_hole_single_emit
  -> free_hole_doublet_buffer_emit
  -> bypass_transport_emit
  -> bypass_bracket_lift_emit
  -> private_firewall_route_sidecar
  -> vertical_halfturn_guard
  -> residual_pair_close_tail
  -> raw_count_shadow
```

## Assumption Challenge

Candidate vertices considered: runners, arcs, gaps, fixed sections, section
boundaries, wall-crossings, residue buckets, cover arcs, Fourier modes,
matroid circuits, stream states, and proof contracts.

Chosen vertices are proof contracts and sidecar obligations.  This preserves
random031 terminal discharge, owner-tail emission, route `R`, vertical guard,
and forbidden quotient cuts.  It destroys raw row identity, component index
after emission, scalar counts, residue chambers, and hydrotope volumes.

The challenged assumption is that a proof contract should route rows directly.
The better interface routes typed sidecars first and emits terminal tokens only
when the hidden tail can no longer change them.
