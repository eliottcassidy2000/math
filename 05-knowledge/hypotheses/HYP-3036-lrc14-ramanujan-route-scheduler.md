---
id: HYP-3036
title: LRC14 Ramanujan primitive-period scheduler for coarse-gate route residuals
status: EVIDENCE / residual route scheduler; not a proof
source: codex-2026-06-26-S200
tangent: T1117
script: 04-computation/lrc14_ramanujan_route_scheduler_codex_s200.py
result: 05-knowledge/results/lrc14_ramanujan_route_scheduler_codex_s200.out
related:
  - HYP-3033
  - HYP-3031
  - HYP-3030
  - HYP-3029
  - HYP-3028
  - HYP-3027
  - HYP-3026
  - HYP-3024
  - HYP-3023
  - HYP-3020
  - HYP-2963
  - LTI-184
  - LTT-082
  - T1112
  - T1111
  - OPEN-Q-108
---

# HYP-3036: Ramanujan Primitive-Period Route Scheduler

## Claim

After the HYP-3030 status-topology gate, the remaining route debt can be
scheduled by a small exact-period deck:

```text
coarse ET+Henselian-unit status gate
  + primitive safe-residue counts for 2 <= q <= 13
  -> route-pure scheduler on the S194 residual fibers.
```

The carrier does not use exact `M`, exact safe interval length, full barcode,
or arc-Cech topology.  It tests primitive phases `a/q` and records whether all
runners are at closed distance at least `1/14` from an integer.  This is the
Ramanujan/exact-period pullback of the residual route question: denominator
layers, not runners, carry the route decision.

This is not a proof of LRC14.  It is evidence for a theorem target that
separates post-status open-route scheduling from boundary/open status.

## Computation

Script:

```text
04-computation/lrc14_ramanujan_route_scheduler_codex_s200.py
```

Stored output:

```text
05-knowledge/results/lrc14_ramanujan_route_scheduler_codex_s200.out
```

The script parses the stored S194 result and rebuilds only the `38` packets in
the `15` route-mixed coarse ET+unit residual fibers.  This keeps the dependency
honest: HYP-3030 supplies the full-bank reduction; HYP-3036 tests the residual
scheduler.

Key audit:

```text
coarse_et_unit_gate              fibers=15 mixed_route=15 mixed_status=0 mixed_packets=38
coarse_plus_first_q_2_13         fibers=30 mixed_route=0  mixed_status=0 mixed_packets=0
coarse_plus_first_q_2_14         fibers=30 mixed_route=0  mixed_status=0 mixed_packets=0
coarse_plus_primitive_deck_2_13  fibers=30 mixed_route=0  mixed_status=0 mixed_packets=0
coarse_plus_primitive_deck_2_14  fibers=30 mixed_route=0  mixed_status=0 mixed_packets=0
coarse_plus_ramanujan_trace_2_14 fibers=37 mixed_route=0  mixed_status=0 mixed_packets=0
```

Deck readout:

```text
q_witness_packets=23 unique_deck13=6
covering_packets=15 unique_deck13=1
covering_deck13={(0,0,0,0,0,0,0,0,0,0,0,0): 15}
```

Every residual `Q-WITNESS` row has positive primitive safe mass on some
denominator `q <= 13`.  Every residual `COVERING-MOMENT` row has zero such mass.

## Guardrail

The `q=14` layer is not a theorem-route separator by itself.  Many covering
rows have six primitive safe residues at `q=14`, matching the boundary scale
rather than a direct `q<=13` witness.  The safe-period deck must therefore keep
the cutoff:

```text
q <= 13  direct q-witness scheduler
q = 14   boundary / covering / moment scheduler
q > 14   later covering-moment certificate
```

The integer Ramanujan trace deck `c_q(v)` also separates the tested residuals,
but it is only a diagnostic unless paired with the safe-phase inequality.  Raw
harmonic trace profiles can distinguish packets without proving a lonely-runner
witness.

## Tournament Analysis

Vertices are proof carriers, not runners:

```text
primitive_count_deck_2_13
first_safe_q_2_13
ramanujan_trace_deck_2_14
coarse_et_unit_status_gate
exact_magnitude_cocycle
raw_residue_terminal_word
```

Pairwise observable:

```text
status_purity, route_purity, primitive_period, harmonic_trace,
avoid_exact_M, compression, proof_cost
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1}
directed_3cycles=0
hamiltonian_path_count=1
score_order=
  primitive_count_deck_2_13 >
  first_safe_q_2_13 >
  ramanujan_trace_deck_2_14 >
  coarse_et_unit_status_gate >
  exact_magnitude_cocycle >
  raw_residue_terminal_word
```

## Theorem Target

Prove, familywise and not merely on the S194 residual list:

```text
Inside a coarse ET+Henselian-unit fiber that is already strict-open,
route-mixed Q-WITNESS packets are exactly those with positive primitive
safe-residue mass at some q <= 13.
```

Then covering rows in the same coarse fiber are no longer route ambiguity; they
are routed to the `q=14` boundary/covering-moment layer or to the later
boundary-moment certificate.

## Next Hook

Add `primitive_safe_deck_2_13` and `first_primitive_safe_q_2_13` to the packet
sidecar used by HYP-3027/HYP-3031.  A cached full-bank ledger should verify
that the deck never introduces boundary/open mixing and should report where it
is route-pure beyond the S194 residual set.
