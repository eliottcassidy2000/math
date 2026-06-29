---
id: HYP-3460
title: LRC14 phase-branch color pullback
status: EVIDENCE / exact color-pullback certificate; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-2593/HYP-2595, S359/S363, HYP-3438, HYP-3450, HYP-3455, and HYP-3459
tangent: T1420
technique: LTI-420
tournament_technique: LTT-320
script: 04-computation/lrc14_phase_branch_color_pullback_codex_20260629.py
result: 05-knowledge/results/lrc14_phase_branch_color_pullback_codex_20260629.out
reflection: 07-reflections/lrc14-phase-branch-color-pullback-codex-20260629.md
related:
  - HYP-3459
  - HYP-3458
  - HYP-3457
  - HYP-3456
  - HYP-3455
  - HYP-3454
  - HYP-3453
  - HYP-3452
  - HYP-3451
  - HYP-3450
  - HYP-3439
  - HYP-3438
  - HYP-3436
  - HYP-3422
  - HYP-2991
  - HYP-2989
  - HYP-2595
  - HYP-2594
  - HYP-2593
  - THM-523
  - OPEN-Q-108
---

# HYP-3460: LRC14 Phase-Branch Color Pullback

## Claim

The older coloring route and the newer branch-gate route are compatible
through an exact pullback:

```text
regular circular coloring multiplier t=a/(14V)
phase color b=a mod 14
two-adic branch coordinate u=2t mod 1
branch color = first/second half of t
```

For S3 rows `S=P union {V-e : e in E}`, this lets the HYP-2593/HYP-2595
phase-color CRT reservoir be compared directly with HYP-3438 survivor gates
and HYP-3450 component classes.

## Exact Readout

Script:

```text
04-computation/lrc14_phase_branch_color_pullback_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_phase_branch_color_pullback_codex_20260629.out
```

The main row is the HYP-3455 noncanonical rank-`6` exception:

```text
random_covering_031
S=(12,23,45,55,58,70,84,93,113,120,147,169,173)
V=173
P=(12,)
E=(0,4,26,53,60,80,89,103,115,118,128,150)
```

As a phase-colored CRT row:

```text
Sigma=2184561576461/1192472694400 ~= 1.831959
K=2254
cGP=12
actual_count=282
open_count=282
actual_deficit=41651852906953/1192472694400 ~= 34.928978
8*(k+cGP)+1=193
```

Thus `random031` is not dangerous in the HYP-2595 colored-discrepancy sense:
its exact deficit is far below the candidate `8*(k+cGP)+1` bound, even though
the crude component cutoff `K/Sigma` is large.

The phase counts are mirror-symmetric:

```text
{1:21, 2:16, 3:23, 4:27, 5:22, 6:21, 7:22,
 8:21, 9:22, 10:27, 11:23, 12:16, 13:21}
```

and the phase-branch matrix has no mirror failures under

```text
b -> 14-b, branch0 <-> branch1.
```

Pulling the `282` actual CRT witnesses back through `u=2t mod 1` gives:

```text
component_class_counts={'both_alive':208,'branch0_only':37,'branch1_only':37}
no_component_hits=0
gate_route_counts={
  'edge_singleton_parent_gate':76,
  'edge_survivor_residual':58,
  'mixed_owner_residual':66,
  'owner_current_small_delta':42
}
gate_mask_counts={'both':32,'branch0':105,'branch1':105}
no_gate_hits=40
ambiguous_gate_hits=0
hard_gate_hits={}
hard_component_hits={(43,'branch1',2):6,(54,'branch0',2):6}
```

The key point is `hard_gate_hits={}`.  HYP-3455's two max-delta gates are the
mirror pair at components `43` and `54`, with deltas `(3,4)` and `(4,3)`.
The actual phase-grid witnesses hit neither of those gates.  They only touch
the same hard components through the opposite branch and lower total delta
`2`.

## Proof Pull

This reframes the HYP-3455 obligation.

The seven-owner mirror gate pair is still a real continuous branch-gluing
clause, but it is not a colored-CRT survivor obstruction at `V=173`.  The
colored placement layer routes around it with many exact regular circular
`14`-colorings.

Candidate lemma:

```text
Any branch-colored max-delta survivor-gate obstruction with zero compatible
phase-grid hits must discharge by colored resonance cancellation, a low-rank
component escape, an endpoint-spine/wall lift, owner-current imbalance,
two-adic descent, or signed-SPEC/Rprime debt.
```

For the AP84 comparison rows, the same pullback recovers the canonical side:

```text
ap84_m1: actual=12, boundary_bonus=2, all component hits both_alive,
         hard gates are the finite edge-singleton transient gates.
ap84_m5: actual=68, actual_deficit=-124/77, all component hits both_alive,
         only the rank-one endpoint-phase edge gates are hit.
```

So HYP-3460 separates two kinds of coloring/gate interaction:

1. AP84 side: phase colors hit the named AP endpoint/transient gates.
2. random031 side: phase colors avoid the max-delta noncanonical gates and
   use low-delta branch-compatible routes.

## Tournament Analysis

Vertices are proof carriers, not runners:

```text
regular_circular_coloring_multiplier
phase_color_CRT_reservoir
colored_resonance_discrepancy
phase_branch_pullback
max_delta_gate_avoidance
branch_component_escape_router
raw_phase_counts
```

Stored fingerprint:

```text
score_hist={17:1,50:1,54:1,56:2,58:1,59:1}
directed_3cycles=0
hamiltonian_path=
  phase_branch_pullback
  -> max_delta_gate_avoidance
  -> branch_component_escape_router
  -> colored_resonance_discrepancy
  -> phase_color_CRT_reservoir
  -> regular_circular_coloring_multiplier
  -> raw_phase_counts
```

## Assumption Challenge

Considered vertices:

```text
runners, phase colors, branch colors, endpoint residues, gates, components,
wall crossings, Fourier modes, fixed circle sections, and proof obligations.
```

The chosen pullback preserves both the regular circular-coloring predicate and
the two-adic branch witness predicate.  It destroys raw runner order and raw
gate mass, which is legal only because the lost data is replaced by exact
phase/branch counts, mirror checks, and hard-gate hit counts.

## Status

Relation to HYP-3458 and HYP-3459: HYP-3458 keeps a canonical `35`-state
AP-tail coloring and endpoint-rank subcolor, while HYP-3459 says the AP84
splice is legal only with the full gate/floor/height/endpoint/branch/router/
zipper color packet.  HYP-3460 is the noncanonical random031 sibling: it asks
whether phase-grid colors actually hit the hard branch gates.

HYP-3460 does not prove LRC(14).  It creates a sharper proof interface:
noncanonical gate-gluing debt should be tested against phase-color pullbacks
before being treated as a global obstruction.  For `random031`, that test says
the hard gates are bypassed by exact colored CRT witnesses.
