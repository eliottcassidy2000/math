---
id: HYP-3526
title: LRC14 spigot route-sidecar dispatch
status: EVIDENCE / finite dispatch guardrail; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3513, HYP-3521, and HYP-3522, inspired by spigot-style bounded-state emission
tangent: T1526
technique: LTI-526
tournament_technique: LTT-426
script: 04-computation/lrc14_spigot_route_sidecar_dispatch_codex_20260629.py
result: 05-knowledge/results/lrc14_spigot_route_sidecar_dispatch_codex_20260629.out
reflection: 07-reflections/lrc14-spigot-route-sidecar-dispatch-codex-20260629.md
related:
  - HYP-3525
  - HYP-3524
  - HYP-3523
  - HYP-3522
  - HYP-3521
  - HYP-3520
  - HYP-3513
  - HYP-3512
  - HYP-3494
  - HYP-3493
  - HYP-3490
  - HYP-3486
  - HYP-3485
  - HYP-3484
  - HYP-3483
  - HYP-3482
  - HYP-3481
  - HYP-3480
  - HYP-3474
  - HYP-3472
  - THM-523
  - OPEN-Q-108
---

# HYP-3526: LRC14 Spigot Route-Sidecar Dispatch

## Claim

HYP-3526 makes the HYP-3513 guardrail operational.  The spigot analogy is:
emit the next terminal certificate from bounded carry state, without returning
to row names.  For the current LRC14 frontier, the emitted object is a terminal
certificate and the carry that cannot currently be removed is route sidecar
`R`.

The row-free private-status route is viable:

```text
I/Q -> multiplicity-one incidence-cut lemma -> private firewall status
```

The terminal route is not currently reconstructable from the tested existing
data:

```text
I/Q does not reconstruct HYP-3490 route
all HYP-3474 axes + I/Q still do not reconstruct HYP-3490 route
R reconstructs HYP-3490 route
```

Therefore terminal dispatch should retain `R` unless a separate route
reconstruction theorem is proved.

Concurrent HYP-3523, HYP-3524, and HYP-3525 are complementary rather than
redundant.
HYP-3523 supplies the random031 certificate stream: `79` component events emit
`77` terminal certificates, with free-hole doublets buffered and the bypass
emitted only with typed owner carry.  HYP-3524 supplies the owner-tail schedule
and chamber canary: transport emits `(23,93,113)`, branch-boundary lift emits
`(147,169)`, and the residual tail is `(45,173)`.  HYP-3525 gives the general
GuardedEmission schema: emit only when visible quotients are constant over
hidden tails, otherwise hold the missing sidecar.  HYP-3526 is the concrete
route-legality canary under that schema: neither `I/Q` nor the existing colored
axes plus `I/Q` replace the HYP-3490 route, so those spigot emissions still
need route sidecar `R` unless a route-reconstruction theorem is proved.

## Exact Readout

Computed by:

```text
04-computation/lrc14_spigot_route_sidecar_dispatch_codex_20260629.py
05-knowledge/results/lrc14_spigot_route_sidecar_dispatch_codex_20260629.out
```

HYP-3513 purity:

```text
private_status_pure_by_I=True
private_status_pure_by_Q=True
route_pure_by_I=False
route_pure_by_Q=False
route_pure_by_IQ=False
route_pure_by_R=True
existing_axes_route_reconstructable=False
all_existing_axes_plus_IQ_route_pure=False
```

The strongest negative route test in the script is:

```text
axes=('K','N','T','S','F','C','M','A','I','Q')
target=h3490_route
fibers=125
max_fiber=7
mixed=1/7
pure=False
```

Random031 terminal facts imported from HYP-3521:

```text
cell_terminal_hist={
  ordinary_rank2_route_component:230,
  free_hole_single_bracket_packet:26,
  free_hole_doublet_packet:14,
  pure_bypass_owner_boundary_component:12
}
component_terminal_hist={
  ordinary_rank2_route_component:64,
  free_hole_single_bracket_packet:10,
  free_hole_doublet_packet:4,
  pure_bypass_owner_boundary_component:1
}
terminal_certificate_count_after_doublet_collapse=77
identity_282=230 ordinary_rank2 + 40 bracketed_free_hole + 12 pure_bypass
identity_77=64 ordinary + 10 free_single + 2 free_doublet + 1 bypass
```

Owner carry imported from HYP-3522:

```text
seam_owners=(23,45,93,113,147,169,173)
transport_owners=(23,93,113)
branch_boundary_owners=(23,93,147,169)
bracket_lift_owners=(147,169)
residual_after_branch_boundary=(45,173)
bypass_size=12
mirror_pair_count=6
```

Dispatch decision:

```text
row_free_dispatch_rule=I/Q may reduce the private cut only before terminal route is supplied
current_route_reconstruction_from_existing_data=False
terminal_dispatch=retain_R
spigot_carry=R
spigot_output=terminal_certificate
```

## Assumption Challenge

Candidate vertices considered: runners, arcs, row names, raw counts, private
bits, quotient axes, incidence sidecars, route sidecars, terminal
certificates, owner-filtration layers, and proof obligations.

The chosen vertices are dispatch compressors.  The quotient preserves the
private-firewall predicate, the random031 terminal route, and the bounded
carry needed to emit terminal certificates.  It destroys row identity, raw
interval geometry, gate order, and scalar counts.

The challenged assumption is that `I/Q` can replace the route.  They cannot on
current data.  `I/Q` are the row-free private-cut interface; `R` is the
terminal carry.

## Proof Pull

This converts the next proof route into a clean interface:

```text
P1. Prove I/Q as a row-free multiplicity-one incidence-cut lemma.
P2. Do not let I/Q stand in for the HYP-3490 five-way route.
P3. For current data, terminal random031 dispatch must retain R.
P4. HYP-3521 emits terminal certificates while HYP-3522 supplies owner carry.
P5. A future route-reconstruction theorem may delete R; until then R is the spigot carry.
```

The immediate LRC14 use is to state terminal lemmas with `R` as an explicit
sidecar rather than burying route selection inside row names or inside the
private-firewall predicate.

## Tournament Analysis

Vertices are dispatch compressors, not runners, arcs, row names, or raw
counts.

```text
pairwise_observable =
  private purity
  + route purity
  + row-free compression
  + bounded carry
  + random031 compatibility

switch =
  higher dispatch score; ties use candidate name

score_hist={-82:1,-53:1,151:1,156:1,159:1,210:1,218:1,226:1}
directed_3cycles=0
sccs=8 singleton SCCs
hamiltonian_path =
  IQ_plus_R_terminal_spigot
  -> terminal_certificate_ledger_plus_R
  -> owner_filtration_plus_R
  -> IQ_without_R_private_only
  -> row_name_exception_list
  -> raw_private_bit
  -> all_colored_plus_IQ_route_reconstruction
  -> raw_count_shadow
```

The winning carrier is `IQ_plus_R_terminal_spigot`: it keeps the row-free
private incidence cut and explicitly carries `R` for terminal route emission.
