---
id: HYP-3523
title: LRC14 random031 spigot-style terminal certificate stream
status: EVIDENCE / finite streaming certificate scheduler; not an LRC14 proof
source: codex-2026-06-29 spigot-inspired continuation of HYP-3522, HYP-3521, HYP-3520, and HYP-3513
tangent: T1523
technique: LTI-523
tournament_technique: LTT-423
script: 04-computation/lrc14_random031_certificate_spigot_codex_20260629.py
result: 05-knowledge/results/lrc14_random031_certificate_spigot_codex_20260629.out
reflection: 07-reflections/lrc14-random031-certificate-spigot-codex-20260629.md
related:
  - HYP-3522
  - HYP-3521
  - HYP-3520
  - HYP-3513
  - HYP-3511
  - HYP-3510
  - HYP-3490
  - HYP-3486
  - HYP-3485
  - HYP-3483
  - HYP-3481
  - THM-523
  - OPEN-Q-108
---

# HYP-3523: LRC14 Random031 Spigot-Style Terminal Certificate Stream

## Claim

The spigot-algorithm prompt suggests a useful proof scheduler for the
`random_covering_031` terminal packet: emit local certificates in witness order,
discard emitted state, and keep only bounded deferred carry when a later local
event can still change the current certificate.

For the current HYP-3521/HYP-3522 random031 packet, this is exact as a finite
audit:

```text
79 component events
77 terminal certificates
max free-hole predigit buffer = 1 component
max owner carry = 4 labels
owner carry reduces after 1 event to residual pair (45,173)
```

Thus the random031 terminal proof can be organized as a streaming invariant:
ordinary route and free-hole single certificates emit immediately; doublet
free-hole packets use a one-component predigit buffer; the bypass emits a
transport word and carries owner debt only until the adjacent branch-boundary
ordinary component reduces it to `(45,173)`.

This does not prove LRC14.  It gives a compact scheduler for the remaining
random031 proof obligations and prevents the proof from retaining the whole
component table after certificates have already emitted.

## Exact Readout

Computed by:

```text
04-computation/lrc14_random031_certificate_spigot_codex_20260629.py
05-knowledge/results/lrc14_random031_certificate_spigot_codex_20260629.out
```

Stream summary:

```text
component_tokens=79
terminal_certificate_count=77
emitted_certificate_count=77
component_type_hist={bypass:1, free_hole:14, ordinary:64}
terminal_class_hist={
  ordinary_rank2_route_component:64,
  free_hole_single_bracket_packet:10,
  free_hole_doublet_packet:4,
  pure_bypass_owner_boundary_component:1
}
action_hist={
  emit_ordinary_route:63,
  emit_ordinary_route+apply_branch_boundary_lift+emit_bypass_owner_boundary_certificate:1,
  emit_free_hole_single:10,
  buffer_free_hole_predigit:2,
  emit_free_hole_doublet:2,
  emit_bypass_transport_predigit:1
}
```

Doublet buffer:

```text
doublet_cluster_count=2
doublet_pending_rank_widths={(8,13):1, (2,3):1}
max_pending_doublet_components=1
```

The two HYP-3511 doublets close at the next component event:

```text
cluster=(8,13) first_rank=30 emit_rank=31 component_indices=(53,5)
cluster=(2,3)  first_rank=75 emit_rank=76 component_indices=(76,77)
```

Owner carry:

```text
seam_owners=(23,45,93,113,147,169,173)
transport_owners=(23,93,113)
branch_boundary_owners=(23,93,147,169)
boundary_component_indices=(27,58)
boundary_component_ranks=(44,46)
bypass_component_rank=45
residual_after_transport=(45,147,169,173)
bracket_lift_owners=(147,169)
residual_after_branch_boundary=(45,173)
owner_carry_rank_width=1
```

The owner events are:

```text
rank=45 open_owner_carry carry=(45,147,169,173)
rank=46 apply_branch_boundary_lift carry=(45,173)
```

## Quotient Result

Target: preserve the stream action at each component event.

```text
component_type:
  mixed_fibers=2, mixed_rows=78
terminal_class:
  mixed_fibers=2, mixed_rows=68
terminal_class_plus_cluster:
  mixed_fibers=3, mixed_rows=68
terminal_class_plus_spigot_state:
  mixed_fibers=0, mixed_rows=0
component_index:
  mixed_fibers=0, mixed_rows=0
```

Terminal class alone forgets whether a free-hole doublet component is a first
predigit or the closing component, and it also forgets whether an ordinary
component applies the bypass branch-boundary lift.  Adding the small spigot
state restores purity without falling back to raw component identity.

## Micro Release Test

A second, smaller quotient test is useful for future proof compression.  The
safe bypass packet is the only carrier in the tested family that both emits
the certified prefix and keeps all typed carry:

```text
full_spigot_packet              score=87
owner_filtration_without_route  score=68  illegal_loss=HYP-3490 route
transport_plus_bracket          score=26  illegal_loss=residual_pair,route
HYP3520_owner_current_only      score=21  illegal_loss=branch_lift,route
terminal_class_only             score=-1  illegal_loss=owner,residual,route
raw_counts                      score=-21 illegal_loss=owner layers,route,free-hole type
raw_seven_owner_shadow          score=-28 illegal_loss=transport,bracket,residual,route
```

This is not a second proof object; it is a compression guardrail.  It says a
candidate quotient may forget raw phase, runner, or count head data only after
the emitted proof packet is certified and the residual/route tail remains
named.

## Proof Pull

HYP-3523 turns the HYP-3521/HYP-3522 terminal interface into five proof tasks:

1. Treat the `77` random031 terminal certificates as a left-to-right stream over
   component events.
2. Prove ordinary route and free-hole single certificates emit immediately and
   require no row-name memory.
3. Prove free-hole doublets need only a one-component predigit buffer, closing
   at the next component event.
4. Prove the bypass owner carry opens as `(45,147,169,173)` and HYP-3522 branch
   bracketing reduces it one event later to `(45,173)`.
5. Attach the HYP-3513 route sidecar `R` or prove route reconstruction; the
   final non-streamed proof object is the residual owner pair `(45,173)`, not
   the full seven-owner seam.

The proof-facing invariant is:

```text
emitted_prefix + predigit_buffer + owner_carry + route_sidecar_R
```

where the predigit buffer has size at most `1` component and the owner carry is
at most `4` labels before reducing to the two-owner residual.

## Tournament Analysis

Vertices are stream states and proof-carrier buffers, not runners or raw arcs.

Pairwise observable: streamability, bounded carry, terminal predicate retention,
and route/firewall sidecar retention.

Switch/gauge: higher proof-facing stream score; ties use the fixed carrier
order.

Fingerprint:

```text
score_hist={10:1, 41:1, 63:1, 80:1, 86:1, 92:1, 98:1, 104:1}
directed_3cycles=0
hamiltonian_path=
  spigot_state_plus_route_R
  -> terminal_class_plus_carry_state
  -> owner_carry_buffer
  -> doublet_predigit_buffer
  -> terminal_certificate_ledger
  -> terminal_class_shadow
  -> component_type_shadow
  -> raw_cell_count_shadow
```

## Assumption Challenge

Candidate vertices considered: runners, gaps, fixed circle sections, section
boundaries, wall-crossing events, residues, cover arcs, Fourier modes, matroid
circuits, proof obligations, and stream states.

Chosen vertices are terminal stream states plus predigit/carry buffers.  This
preserves the random031 terminal-discharge predicate, free-hole doublet
collapse, owner residual, and firewall route sidecar.  It deliberately destroys
raw runner order, raw component index after emission, and scalar cell-count
shadows.

The challenged assumption is that the terminal proof must keep the whole
random031 component table live until the end.  The audit shows all emitted
ordinary/single/doublet certificates can stream away, leaving only bounded
carry plus the residual two-owner lemma and route sidecar.
