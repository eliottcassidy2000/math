---
id: HYP-3457
title: LRC14 AP84 finite transient closure
status: EVIDENCE / finite mixed AP-tail transient certificate; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3452, HYP-3454, HYP-3456, and HYP-3439
tangent: T1417
technique: LTI-417
tournament_technique: LTT-317
script: 04-computation/lrc14_ap84_finite_transients_codex_20260629.py
result: 05-knowledge/results/lrc14_ap84_finite_transients_codex_20260629.out
reflection: 07-reflections/lrc14-ap84-finite-transients-codex-20260629.md
related:
  - HYP-3456
  - HYP-3455
  - HYP-3454
  - HYP-3453
  - HYP-3452
  - HYP-3451
  - HYP-3450
  - HYP-3439
  - HYP-3438
  - HYP-3437
  - HYP-3436
  - HYP-3435
  - HYP-3434
  - HYP-3433
  - HYP-3431
  - HYP-3429
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3422
  - HYP-3417
  - HYP-3129
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3457: LRC14 AP84 Finite Transient Closure

## Claim

HYP-3457 closes the finite mixed side of the AP-tail bridge left by
HYP-3452/HYP-3454/HYP-3456.

For

```text
S_m = {1,2,...,11,13,84m},  m=1,2,3,4,
```

the component-cover audit has exactly four low-rank escape windows.  They are

```text
outer_left  = [8/49,(98m-1)/(588m)]         L[B1:7]  R[E:84m]
inner_left  = [(98m+1)/(588m),6/35]         L[E:84m] R[B1:5]
inner_right = [29/35,(490m-1)/(588m)]       L[B0:5]  R[E:84m]
outer_right = [(490m+1)/(588m),41/49]       L[E:84m] R[B0:7].
```

The inner-left window is the best component in each transient row:

```text
m=1: [33/196,6/35]
m=2: [197/1176,6/35]
m=3: [295/1764,6/35]
m=4: [131/784,6/35].
```

Thus the finite side is no longer an unnamed "check m=1..4" instruction.  It is
a four-row exact packet with endpoint labels, lengths, mirror symmetry,
component-cover rank, and dead-projection side data.

## Exact Readout

Script:

```text
04-computation/lrc14_ap84_finite_transients_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_ap84_finite_transients_codex_20260629.out
```

Audit result:

```text
exact_window_failures=[]
closure_failures=[]
rank_drop_failures_for_m_ge_3=[]
inequality_failures=[]
```

Row summary:

```text
m=1: components=26, dead=22, escapes=4, low_rank_escape=4, max_pair=3, projection=1/22
m=2: components=48, dead=44, escapes=4, low_rank_escape=4, max_pair=3, projection=1/44
m=3: components=70, dead=66, escapes=4, low_rank_escape=4, max_pair=2, projection=1/66
m=4: components=88, dead=84, escapes=4, low_rank_escape=4, max_pair=2, projection=1/84
```

The endpoint inequalities explain the phase break.  For `m=1..4`,

```text
(98m+13)/(588m) - 6/35 = (455-98m)/(2940m) > 0,
```

so the high gap still crosses the fixed `B1:5` wall and the best survivor is
mixed `E:84m/B1:5`.  At `m=5` this sign changes, matching HYP-3454's pure
`E:84m/E:84m` endpoint phase.

## Proof Pull

After HYP-3457, the AP-tail bridge has three named pieces:

```text
finite transient packet: HYP-3457
rank-one endpoint interval for m>=5: HYP-3454
mod-35 escape count for all m: HYP-3456
```

The remaining AP-tail work is now the carrier/splice step:

```text
1. import or prove the HYP-3431 fixed-corridor identity as the complete
   low branch-union carrier;
2. splice HYP-3454/HYP-3456/HYP-3457 into the HYP-3439 rank-5 AP-tail descent.
```

This is still not a proof of LRC14.  It only closes the finite AP-tail sidecar
inside the current HYP-3439/HYP-3452 corridor.

## Tournament Analysis

Vertices are finite transient proof obligations, not runners or raw arcs.

```text
pairwise_observable =
  exact window match + endpoint labels + mirror pairing + rank payload
score_hist={21:1,43:1,51:1,52:1,55:1,59:1}
directed_3cycles=0
hamiltonian_path =
  four_window_explicit_packet
  -> component_audit_exact_match
  -> mirror_pairing
  -> mixed_endpoint_clip_inequalities
  -> rank_drop_m3_to_m4
  -> raw_transient_sample_rows
```

Assumption challenge: runners, `m`-values, the four survivor windows, fixed
corridor sections, high-grid wall events, endpoint walls, dead-cover graph
summaries, and proof obligations were considered.  The chosen quotient
preserves the finite AP-tail branch-union survivor predicate and endpoint
labels, but destroys arbitrary non-AP component geometry.  It is only the
finite transient sidecar for the AP-tail bridge.
