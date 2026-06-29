---
id: HYP-3462
title: LRC14 AP84 corridor-splice certificate
status: EVIDENCE / AP-tail bridge splice; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3431, HYP-3439, HYP-3454, HYP-3456, HYP-3457, and HYP-3460
tangent: T1422
technique: LTI-422
tournament_technique: LTT-322
script: 04-computation/lrc14_ap84_corridor_splice_certificate_codex_20260629.py
result: 05-knowledge/results/lrc14_ap84_corridor_splice_certificate_codex_20260629.out
reflection: 07-reflections/lrc14-ap84-corridor-splice-certificate-codex-20260629.md
related:
  - HYP-3460
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

# HYP-3462: LRC14 AP84 Corridor-Splice Certificate

## Claim

HYP-3462 closes the AP84 splice obligation left by HYP-3439:

```text
canonical rank-6 base via HYP-3431,
rank-5 AP-tail descent via HYP-3454/HYP-3456/HYP-3457.
```

It is complementary to HYP-3458, HYP-3459, and HYP-3460: HYP-3458/HYP-3459
name the retained AP-tail color/rank/discrepancy state, HYP-3460 handles the
noncanonical phase-branch color pullback for the random031 sibling, and
HYP-3462 closes the structural AP84 carrier/splice that routes the canonical
tail sidecars into HYP-3439.

The low-core branch-union carrier is exactly the HYP-3431 pair

```text
C1=[8/49,6/35],     C0=[29/35,41/49].
```

The branch-specific low good pieces overlap inside these two corridors:

```text
b0_good=[15/91,13/77] union [29/35,41/49]
b1_good=[8/49,6/35] union [64/77,76/91]
branch_union=C1 union C0.
```

So the correct imported HYP-3431 statement is a complete low branch-union
carrier, not a one-branch-only corridor claim.

On the AP84 tower

```text
S_m={1,2,...,11,13,84m},
```

the one-branch overlap-rescue rank split is:

```text
m=1:  rank 6 core (3,5,7,9,11,13)
m>=2: rank 5 core (5,7,9,11,13)     checked through m=70.
```

The sidecar splice is then:

```text
m=1:   HYP-3431 canonical rank-6 base + HYP-3457 finite windows
m=2..4: rank-5 overlap descent + HYP-3457 finite windows
m>=5: rank-5 overlap descent + HYP-3454 endpoint interval + HYP-3456 floor count
```

This removes the AP84 splice as an unnamed HYP-3439 obligation.  It still does
not prove the global primitive-row theorem; the remaining work is the non-AP
transfer through HYP-3453/HYP-3451/HYP-3455 or named owner/current/state-lift
debt.

## Exact Readout

Script:

```text
04-computation/lrc14_ap84_corridor_splice_certificate_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_ap84_corridor_splice_certificate_codex_20260629.out
```

Carrier:

```text
identity_ok=True
branch_union=('[8/49,6/35]', '[29/35,41/49]')
branch_union_measure=4/245
corridor_lengths=('2/245', '2/245')
```

Rank split:

```text
checked_m=1..70
rescue_rank_hist={5: 69, 6: 1}
route_hist={'canonical_rank6_base': 1, 'rank5_ap_tail_descent': 69}
rank_split_failures=[]
```

AP-tail sidecars:

```text
exact_window_failures=[]
closure_failures=[]
rank_drop_failures_for_m_ge_3=[]
checked_endpoint_failures_m_5_to_70=[]
symbolic_endpoint_containment_failures_m_5_to_420=[]
mirror_failures_m_1_to_70=[]
formula_failures_m_1_to_70=[]
component_audit_failures_m_1_to_70=[]
shift_failures_m_1_to_210=[]
```

## Proof Pull

The AP84 branch is now proof-facing rather than exploratory:

1. Use HYP-3431 as the exact branch-union carrier.
2. Use HYP-3457 for the finite `m=1..4` mixed windows.
3. Use HYP-3454 for the `m>=5` endpoint interval.
4. Use HYP-3456 for the alive-boundary floor count.
5. Use HYP-3437/HYP-3439 only for the one-branch overlap-rank interface:
   rank `6` at `m=1`, rank `5` from `m=2` in the checked AP tail.

The next proof task is therefore not another AP84 endpoint/floor/transient
audit.  It is to move the now-closed AP84 splice into the broader
HYP-3453/HYP-3451 component-gate/conductance route and to discharge the
noncanonical HYP-3455 seven-owner gluing clause or name the residual debt.

## Tournament Analysis

Vertices are splice proof obligations, not runners or raw arcs.

```text
pairwise_observable =
  carrier retention + rank split + endpoint/floor/transient payload
  + scalar penalty
score_hist={19:1,47:1,54:1,56:1,57:2,59:2}
directed_3cycles=0
hamiltonian_path =
  fixed_branch_union_carrier
  -> floor_count_boundary_clock
  -> canonical_rank6_base_split
  -> endpoint_clock_tail
  -> rank5_overlap_descent
  -> finite_transient_sidecar
  -> component_audit_checks
  -> raw_rescue_rank_scalar
```

Assumption challenge: runners, `m`-values, residues, fixed circle sections,
section boundaries, high-grid gaps, wall-crossing events, branch-union
corridors, cover arcs, survivor gates, component graph nodes, odd rescue cores,
endpoint walls, Fourier modes, matroid circuits, and proof obligations were
considered.  The chosen quotient preserves the AP84 branch-union survivor
predicate plus the HYP-3439 rank split.  It destroys arbitrary non-AP row
geometry and is therefore only a bridge splice.
