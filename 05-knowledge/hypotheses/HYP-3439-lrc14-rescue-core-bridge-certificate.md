---
id: HYP-3439
title: LRC14 rescue-core bridge certificate
status: EVIDENCE / AP-tail bridge audit; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3438/HYP-3437/HYP-3436/HYP-3450/HYP-3451
tangent: T1400
technique: LTI-400
tournament_technique: LTT-300
script: 04-computation/lrc14_rescue_core_bridge_certificate_codex_20260629.py
result: 05-knowledge/results/lrc14_rescue_core_bridge_certificate_codex_20260629.out
reflection: 07-reflections/lrc14-rescue-core-bridge-certificate-codex-20260629.md
related:
  - HYP-3457
  - HYP-3456
  - HYP-3454
  - HYP-3453
  - HYP-3452
  - HYP-3451
  - HYP-3450
  - HYP-3438
  - HYP-3437
  - HYP-3436
  - HYP-3435
  - HYP-3434
  - HYP-3431
  - HYP-3429
  - HYP-3427
  - HYP-3425
  - HYP-3422
  - HYP-3418
  - HYP-3415
  - THM-523
  - OPEN-Q-108
---

# HYP-3439: LRC14 Rescue-Core Bridge Certificate

## Claim

HYP-3439 executes the reserved bridge between:

```text
HYP-3437 one-branch overlap rescue cores
HYP-3438 survivor-gate words
HYP-3450/HYP-3451 component-cover escapes
```

The audit is deliberately focused on the proof-facing AP/84m corridor spine.
It is not a full recomputation of the HYP-3450 component bank.  The point is
to decide whether the one rank-`6` overlap rescue core is a broad new
noncanonical obstruction or exactly the canonical corridor-fence case already
handled by HYP-3431.

## Exact Readout

Script:

```text
04-computation/lrc14_rescue_core_bridge_certificate_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_rescue_core_bridge_certificate_codex_20260629.out
```

The HYP-3437 one-branch recap remains:

```text
rows_audited=150
negative_naive_slack_rows=59
negative_rows_with_rescue=59/59
all_row_rescue_rank_hist={0:91, 2:7, 4:2, 5:48, 6:2}
negative_row_rescue_rank_hist={2:7, 4:2, 5:48, 6:2}
max_minimum_rescue_rank=6
max_rank_rows=['covering_AP_with_84', 'canonical_84m_01']
```

The AP/84m bridge joins `14` displayed rows (`13` unique speed rows):

```text
negative_bridge_rows=13/14
bridge_rescue_rank_hist={0:1, 5:11, 6:2}
bridge_route_hist={
  ap_tail_rank5_overlap_core_with_two_color_escape: 11,
  canonical_rank6_corridor_fence_with_two_color_escape: 2,
  nonnegative_one_branch_control_with_two_color_escape: 1
}
rank6_bridge_rows=['covering_AP_with_84', 'ap_omit_12_tail_84x01']
aggregate_survivor_branch_mask_hist={branch0:27, branch1:27}
low_rank_escape_range=[4,22]
max_danger_bridge=ap_omit_12_tail_84x02
```

The canonical row and its duplicate `ap_omit_12_tail_84x01` have the rank-`6`
overlap core:

```text
subset=(3,5,7,9,11,13)
naive_slack=-18586/315315
branch0=563/105105
margin=563/105105
low_rank_escape=4
dead=22/26
max_dead_pair_rank=3
survivor_gates=4
```

Every AP tail `ap_omit_12_tail_84x02` through `x12` still has negative
one-branch naive slack, but the minimum rescue core drops to:

```text
subset=(5,7,9,11,13)
rank=5
```

and each row keeps HYP-3450/HYP-3451 low-rank two-colour escapes.  Incoming
HYP-3452 sharpens this AP-tail clause: for the full canonical tail
`{1,...,11,13,84m}`, `m=1..4` are finite mixed transients, `m>=5` enters the
rank-one `E:84m/E:84m` endpoint phase, paired dead-cover rank is already
`<=2` for `m>=3`, and escape counts follow a mod-`35` Beatty correction.
HYP-3454 closes the endpoint interval, HYP-3456 derives the mod-`35` floor
count, and HYP-3457 closes the finite transient packet.  The nonnegative
control `multi_far_84_154` has `rank=0`, positive naive slack
`14699/504504`, `22` low-rank escapes, and only low danger score `0.078671`.

## Proof Pull

The bridge lemma should be stated as a finite interface, not as a scalar
ranking:

```text
If a primitive covering row has negative one-branch naive slack, then its
overlap-tax rescue core must be paired with either:
  canonical corridor-fence data,
  a bounded rank-<=5 AP-tail/noncanonical rescue core,
  or an explicit two-colour survivor/component escape route.
```

For the AP/84m corridor, HYP-3439 says the only rank-`6` bridge is the
canonical `m=1` fence.  The next theorem target is therefore:

```text
prove the canonical rank-6 base case by HYP-3431/HYP-3450/HYP-3451,
prove the rank-5 AP-tail descent using HYP-3454/HYP-3456/HYP-3457,
then generalize the bridge to arbitrary primitive covering rows.
```

Legal exits remain HYP-3438 survivor gates, HYP-3450/HYP-3451 component
escapes, HYP-3429 endpoint spines, HYP-3427 wall words, HYP-3428 two-adic loss
ledgers, owner-current labels, exact-period/state-lift debt, or signed-SPEC.
Raw rescue rank and raw overlap mass are negative controls unless the endpoint
labels, survivor gates, and component addresses are retained.

## Tournament Analysis

Vertices are bridge proof obligations, not runners or raw arcs.

```text
pairwise_observable =
  predicate retention + one-branch cut exactness + two-colour escape payload
  + survivor-gate route + scalar-firewall safety
score_hist={21:1, 54:1, 55:2, 57:2, 58:1, 59:1}
directed_3cycles=0
hamiltonian_path_count=1
hamiltonian_path =
  canonical_corridor_rank6_bridge
  -> two_color_component_escape
  -> one_branch_overlap_rescue_cut
  -> survivor_gate_word_route
  -> ap_tail_rank5_descent
  -> component_conductance_router
  -> endpoint_spine_wall_handoff
  -> raw_rescue_rank_scalar
```

Assumption challenge: runners, odd blockers, even-half gates, gaps, fixed
circle sections, section boundaries, wall-crossing events, residues, cover
arcs, Fourier modes, matroid circuits, survivor gates, component-cover graph
nodes, endpoint walls, and proof obligations were considered.  The chosen
vertices are bridge obligations between one-branch rescue cuts, survivor-gate
words, and two-colour component escape certificates.  This preserves the
negative-slack relocation predicate and destroys raw runner order and most
interval geometry unless endpoint labels and component addresses are restored.
