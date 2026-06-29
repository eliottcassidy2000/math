---
id: HYP-3523
title: LRC14 random031 spigot-style terminal certificate stream
status: EVIDENCE / finite streaming certificate scheduler plus owner-carry quotient subaudit; not an LRC14 proof
source: codex-2026-06-29 spigot-inspired continuation of HYP-3522, HYP-3521, HYP-3520, HYP-3513, extended by HYP-3525 guarded emission, HYP-3526 route-sidecar dispatch, HYP-3527 proof-contract routing, HYP-3528 one-red-clause execution, and the exact owner-carry companion audit
tangent: T1523
technique: LTI-523
tournament_technique: LTT-423
script: 04-computation/lrc14_random031_certificate_spigot_codex_20260629.py
result: 05-knowledge/results/lrc14_random031_certificate_spigot_codex_20260629.out
companion_script: 04-computation/lrc14_random031_spigot_owner_carry_codex_20260629.py
companion_result: 05-knowledge/results/lrc14_random031_spigot_owner_carry_codex_20260629.out
reflection: 07-reflections/lrc14-random031-certificate-spigot-codex-20260629.md
related:
  - HYP-3528
  - HYP-3527
  - HYP-3526
  - HYP-3525
  - HYP-3524
  - HYP-3522
  - HYP-3521
  - HYP-3520
  - HYP-3513
  - HYP-3512
  - HYP-3511
  - HYP-3510
  - HYP-3494
  - HYP-3493
  - HYP-3490
  - HYP-3486
  - HYP-3485
  - HYP-3483
  - HYP-3482
  - HYP-3481
  - HYP-3402
  - HYP-3034
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
component table after certificates have already emitted.  HYP-3525 supplies the
general guard rule: a proof token may emit only when every hidden tail
compatible with the visible quotient has the same terminal/route/owner target.
Incoming HYP-3526 supplies the concrete route-sidecar guardrail: `I/Q` can be
used as a row-free private-cut interface, but current data does not reconstruct
the HYP-3490 route from `I/Q` or from the tested colored axes plus `I/Q`, so
terminal dispatch still carries route sidecar `R`.
Incoming HYP-3527 then consumes the HYP-3523 stream/carry state as part of the
finite proof-contract router: three interfaces are formal-ready, four clauses
carry required sidecars, and the only open tail lemma is the residual
`(45,173)` no-hidden-tail boundary clause.
Incoming HYP-3528 re-executes that router as a sidecars-to-safe-emission ABI
and makes the remaining red clause unique:
`residual_pair_close_tail = residual (45,173) + R + no_hidden_tail_guard`.
Concurrent THM-578/HYP-3529 is not a random031 dependency, but it reinforces
the same proof discipline on a different LRC14 doublet: exact decomposition
and bounded tail are enough; sharp scalar constants are not substitutes for a
typed residual.

## Exact Readout

Primary stream computation:

```text
04-computation/lrc14_random031_certificate_spigot_codex_20260629.py
05-knowledge/results/lrc14_random031_certificate_spigot_codex_20260629.out
```

Companion owner-carry quotient audit:

```text
04-computation/lrc14_random031_spigot_owner_carry_codex_20260629.py
05-knowledge/results/lrc14_random031_spigot_owner_carry_codex_20260629.out
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

## Owner-Carry Companion

The companion audit checks that the owner-carry emissions are locally stable
before they are treated as proof digits:

```text
mirror_owner_persistent=True
bypass_owner_word_hist={(23,93,113):12}
lower_bypass_matches_transport=True
branch_word_constant=True
branch_intersection_constant=True
carry_widths=(7,4,2,2,2)
bounded_carry_max_after_transport=4
terminal_residual_closed=False
owner_boundary_persistence_word=PDPPOOO
global_flow_minus_bypass=(45,55,147,169,173)
```

The extra owner `55` warns that global flow is not the local bypass emitter.
The safe forgetful relations for the owner target are:

```text
flow_class
owner_boundary_persistence_class
owner_presence_word
seam_debt_word
```

The forbidden shortcuts are:

```text
component_size
endpoint_rank_word
mirror_closed_shadow
raw_owner_count
seam_owner_count
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
safe bypass packet is the only carrier in the tested family that both emits the
certified prefix and keeps all typed carry:

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

HYP-3523 turns the HYP-3521/HYP-3522 terminal interface into six proof tasks:

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
6. Use HYP-3526's route-dispatch check to keep `R` in every terminal emission
   until a genuine route-reconstruction theorem replaces it.
7. Use HYP-3525 guarded emission to state which visible packets may emit route,
   terminal, or owner tokens and which must hold the first missing sidecar.
8. Feed these clauses to HYP-3527's proof-contract router, whose ledger leaves
   exactly one open tail: residual `(45,173)` with route `R`.
9. Use HYP-3528's one-red-clause execution as the closeout target: every
   other contract row is finite evidence plus formalization debt, so the next
   theorem-facing move is only `residual_pair_close_tail`.

The proof-facing invariant is:

```text
emitted_prefix + predigit_buffer + owner_carry + route_sidecar_R
```

where the predigit buffer has size at most `1` component and the owner carry is
at most `4` labels before reducing to the two-owner residual.

## Tournament Analysis

Vertices are stream states and proof-carrier buffers, not runners or raw arcs.

Pairwise observable: streamability, bounded carry, terminal predicate retention,
owner-emission stability, and route/firewall sidecar retention.

Switch/gauge: higher proof-facing stream score; ties use the fixed carrier
order.

Primary stream fingerprint:

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

Owner-carry companion fingerprint:

```text
score_hist={9:2,10:1,27:1,90:1,91:1,92:1,109:1,110:1}
directed_3cycles=0
hamiltonian_path=
  exact_spigot_owner_state
  -> transport_branch_residual_words
  -> owner_boundary_persistence_word
  -> seam_and_bypass_owner_words
  -> flow_class_owner_union_sheet_pgf
  -> bypass_owner_word_only
  -> raw_owner_count_shadow
  -> global_flow_minus_bypass
  -> endpoint_rank_shadow
```

## Assumption Challenge

Candidate vertices considered: runners, gaps, fixed circle sections, section
boundaries, wall-crossing events, residues, cover arcs, Fourier modes, matroid
circuits, proof obligations, output certificates, carry states, quotient
guards, owner rows, free-hole packets, and stream states.

Chosen vertices are stream states with predigit/carry buffers.  They preserve
the random031 terminal predicate because each emitted token is a HYP-3521
terminal certificate, HYP-3511 doublet collapse is handled before emission, and
HYP-3522 owner carry is retained until the branch-boundary lift reduces it.

What this destroys: raw runner order, raw component index after emission, full
row-name memory, and scalar cell-count shadows.  That destruction is legal only
because the stream state keeps the predigit/carry information and because
HYP-3513 says route sidecar `R` is still required unless reconstructed.  The
challenged assumption is that the terminal proof must keep the whole random031
component table live until the end; the audit shows all emitted certificates
can stream away, leaving only bounded carry plus the residual two-owner lemma
and route sidecar.
