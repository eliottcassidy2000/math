---
id: HYP-3137
title: LRC14 generating-function payload atlas
status: EVIDENCE / executable payload atlas and tournament scout; not a proof
source: codex-2026-06-27-S273
tangent: T1202
technique: LTI-263
tournament_technique: LTT-161
script: 04-computation/lrc14_generating_function_payload_atlas_codex_s273.py
result: 05-knowledge/results/lrc14_generating_function_payload_atlas_codex_s273.out
related:
  - HYP-3138
  - HYP-3136
  - HYP-3135
  - HYP-3134
  - HYP-3133
  - HYP-3132
  - HYP-3131
  - HYP-3129
  - HYP-3118
  - HYP-3109
  - HYP-3108
  - HYP-3103
  - HYP-3054
  - HYP-2991
  - THM-077
  - THM-076
  - THM-059
  - OPEN-Q-108
---

# HYP-3137: LRC14 Generating-Function Payload Atlas

## Reservation Claim

This lane claims the S273 generating-function pass requested after the
resolvent packet, A000568 edge-envelope, Lee-Yang root, signed-SPEC, and
coordinate-resurrection frontier.

The intended move is to stop asking whether "a generating function" is useful
in the abstract.  Instead compare which proof payload each generating function
keeps:

```text
miss-count PGF
signed SPEC / resonance-lattice Fourier series
A000568 ordinary count / cycle-index quotient
OCF independence polynomial and Walsh EGF
De Moivre / resolvent elementary-symmetric polynomial
Irving-Omar walk determinant generating function
hard-core / polymer partition function
```

Planned scout: mine repo signals from old OCF/Walsh/independence-polynomial
work and current LRC14 PGF/SPEC/A000568/resolvent work; score candidate
generating functions by retained LRC predicate, preserved coordinates,
quotient legality, root/zero information, coefficient-layer payload,
recursion/descent interface, and terminal proof exit.  Tournament Analysis
will use generating-function carriers as vertices, not runners, roots, raw
polynomials, or scalar evaluations.

Challenged assumption: the whole generating function is not automatically the
right object.  Sometimes the proof payload is a coefficient layer, a root
locus, a log-derivative, a quotient of rooted/unrooted classes, or a signed
tail certificate.

## Executable Readout

`04-computation/lrc14_generating_function_payload_atlas_codex_s273.py`
compares twelve generating-function carriers against ten payload axes and
twelve remaining LRC14 proof obligations.  Output is stored in
`05-knowledge/results/lrc14_generating_function_payload_atlas_codex_s273.out`.

Two maps are the useful product.

**Map 1: carrier payloads.**  The top carriers by retained payload are:

```text
signed_SPEC_resonance_series          total 84
resolvent_elementary_symmetric_GF     total 76
edge_witness_recursion_OGF            total 73
A000568_cycle_index_quotient          total 70
repair_cover_Hilbert_series           total 67
walsh_ocf_support_EGF                 total 67
miss_count_PGF_root_locus             total 64
```

The ranking should not be read as "use only SPEC."  It says signed SPEC is
the load-bearing analytic carrier, while the other two members of the minimal
cover below supply the legal quotient and the zero-motion/core bridge that
SPEC alone forgets.

**Map 2: proof obligations.**  The smallest packet covering all twelve
tracked obligations is:

```text
signed_SPEC_resonance_series
+ A000568_cycle_index_quotient
+ miss_count_PGF_root_locus
```

This triple covers `Q_apex_floor`, `R_safe_positive`, `Rprime_signed_low`,
`SPEC_parseval_tail`, `bounded_core_k8_phi4`, `far_push_zero_motion`,
`edge_tail_tip_children`, `global_consistency_quotient`,
`finite_address_observer_gluing`, `closed_form_uniform_constants`,
`coefficient_middle_layer`, and `destroyed_coordinate_guard`.

Interpretation: the LRC14 endgame should not pick between Lee-Yang/PGF,
SPEC, and A000568/edge quotients.  The proof packet wants exactly their
intersection:

```text
support law + signed tail + global edge quotient + root-motion/core bridge.
```

## Tournament Analysis

Tournament vertices are GF proof carriers, not runners, roots, arcs, or raw
values.  The pairwise observable is majority comparison over retained
proof-payload axes; ties follow the declared Hamiltonian path in the script.

The scout reports:

```text
vertices = 12
edges = 66
score histogram = {0:1, 1:1, ..., 11:1}
directed 3-cycles = 0
Hamiltonian path count = 1
selected path =
  signed_SPEC_resonance_series
  -> edge_witness_recursion_OGF
  -> resolvent_elementary_symmetric_GF
  -> A000568_cycle_index_quotient
  -> repair_cover_Hilbert_series
  -> walsh_ocf_support_EGF
  -> miss_count_PGF_root_locus
  -> trivariate_master_GF
  -> hard_core_independence_partition
  -> IO_walk_determinant_GF
  -> q_pochhammer_modular_tail_GF
  -> raw_scalar_evaluation
```

The acyclic tournament is itself a warning: under these axes, every typed
payload beats raw scalar evaluation.  The interesting non-linearity is not in
the tournament order but in the set-cover map: no single carrier covers the
proof packet, and the minimal full cover needs exactly the SPEC/A000568/PGF
triple.

## New Signals To Measure

- `SPEC_support_sieve`: record the resonance lattice, low modes, L2 tail, and
  sign budget, analogously to the THM-077/THM-076 Walsh support law.
- `edge_recursion_depth_PGF`: count safe edge-witness decks by recursive
  depth and mark tail/tip ownership.
- `global_consistency_class`: the A000568 wedge as quotient legality, between
  raw local sectors and paired endpoint children.
- `PGF_root_trajectory_derivative`: how zeros move under far addition and
  endpoint deletion, not only the nearest-root radius.
- `middle_layer_vector`: pair/triple elementary symmetric layers, cumulants,
  and phi4 sign.
- `k8_reflection_fold_adjoint`: HYP-3138's even fold as a finite lookup that
  preserves the gK8 coordinate while resurrecting odd leakage or naming debt.
- `repair_cover_H(q)`: a Hilbert series for minimal sidecars needed to
  resurrect coordinates destroyed by a quotient.
- `Bravais_resonance_cell`: the low-SPEC lattice cell/chamber; quotient
  shapes without a cell address are unsafe.
- `Savitch_packet_depth`: repeated-squaring style proof reachability for
  joining edge packets in small depth.

## Current Hypothesis

The next LRC14 proof packet should be written as a three-carrier theorem:

1. signed SPEC supplies the support law, exact low resonance, and Parseval
   tail;
2. A000568/edge-witness quotients certify controlled forgetting of local
   tail/tip children into global consistency classes;
3. the miss-count PGF/Lee-Yang root curve supplies the Q/apex factor, far
   push-out monotonicity, and bounded-core k=8 zero-motion bridge.

The De Moivre/resolvent and Walsh/OCF material are not competing routes.  They
are templates for how to package the theorem: keep the middle coefficient
layer before scalar evaluation, and prove a support/cancellation law before
collapsing the generating function.  Post-rebase HYP-3138 makes this more
concrete: the k=8 resolvent carrier now has a candidate
`k8_reflection_fold_adjoint`, where the even fold preserves the gK8 coordinate
and odd leakage is resurrected by a finite lookup or routed to named debt.
