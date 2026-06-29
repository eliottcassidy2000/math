---
id: HYP-3483
title: LRC14 random031 recursion flow comparator
status: EVIDENCE / exact recursion-sidecar scout; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3482, HYP-3481, HYP-3480, HYP-3479, HYP-3477, HYP-3460, HYP-3455, HYP-3234, HYP-3232, HYP-3230, and HYP-3418
tangent: T1443
technique: LTI-443
tournament_technique: LTT-343
script: 04-computation/lrc14_random031_recursion_flow_comparator_codex_20260629.py
result: 05-knowledge/results/lrc14_random031_recursion_flow_comparator_codex_20260629.out
reflection: 07-reflections/lrc14-random031-recursion-flow-comparator-codex-20260629.md
related:
  - HYP-3481
  - HYP-3480
  - HYP-3479
  - HYP-3478
  - HYP-3477
  - HYP-3476
  - HYP-3460
  - HYP-3455
  - HYP-3453
  - HYP-3451
  - HYP-3450
  - HYP-3438
  - HYP-3234
  - HYP-3232
  - HYP-3230
  - HYP-3418
  - THM-523
  - OPEN-Q-108
---

# HYP-3483: LRC14 Random031 Recursion Flow Comparator

## Claim

The HYP-3481 topology atlas should not be forced into either a pure `n+2`
or pure `n*2` recursion story.

The exact readout says:

```text
phase transport = n*2 / two-adic pullback in u=2t
owner boundary = n+2 / puncture-boundary insertion debt
```

So the useful proof object is a span:

```text
seven-owner forbidden seam + balanced two-adic bypass flow
  -> local random031 terminal packet or named owner/two-adic/SPEC debt
```

## Exact Readout

Model scores:

```text
two_adic_pullback_n_times_2=74
boundary_insertion_n_plus_2=38
controlled_span=130
```

Two-adic bypass:

```text
phase_witnesses=282
hard_gate_hits=0
bypass_hits=12
branch_balance={0:6,1:6}
component_balance={43:6,54:6}
phase_blocks_by_branch={
  0: (7,8,9,10,11,12),
  1: (2,3,4,5,6,7)
}
mirror_pairs=6
```

Boundary insertion:

```text
dead_island_owners=(23,45,93,113)
seam_owners=(23,45,93,113,147,169,173)
bypass_owners=(23,93,113)
seam_only_owners=(45,147,169,173)
rescue_plus_apex=(173,)
```

The seam is the owner-boundary carrier; the bypass is the phase-flow carrier.

## Proof Pull

The random031 terminal theorem should retain the ordered six-hit blocks, not
just the scalar count `12`.  The correct formal shape looks like:

```text
if a seven-owner mirror seam has zero hard phase flux
and the u=2t pullback has two mirror six-hit bypass blocks
then random031 discharges by local island/bypass current
or emits named owner-current / two-adic / signed-SPEC debt.
```

This connects the older `n+2` recursion work to the two-adic descent work
without collapsing either sidecar.

## Tournament Analysis

Vertices are recursion/proof carriers, not runners:

```text
controlled_span_seam_boundary_plus_twoadic_bypass
two_adic_phase_pullback_blocks
mirror_bypass_hit_pairing
seven_owner_boundary_insertion
dead_island_puncture_word
raw_bypass_hit_count
raw_owner_count
```

Fingerprint:

```text
score_hist={16:1,20:1,38:1,47:1,64:1,74:1,130:1}
directed_3cycles=0
hamiltonian_path=
  controlled_span_seam_boundary_plus_twoadic_bypass
  -> two_adic_phase_pullback_blocks
  -> mirror_bypass_hit_pairing
  -> dead_island_puncture_word
  -> seven_owner_boundary_insertion
  -> raw_bypass_hit_count
  -> raw_owner_count
```

## Assumption Challenge

Considered vertices: runners, gaps, fixed circle sections, section boundaries,
wall-crossing events, residues, cover arcs, Fourier modes, dead islands, seam
gates, phase-flow hits, owner labels, recursion operators, and proof
obligations.  The chosen carrier preserves the random031 terminal discharge
predicate by retaining both the owner-boundary seam and the two-adic bypass
flow; it destroys raw runner order and raw gate counts only after replacing
them with these sidecars.
