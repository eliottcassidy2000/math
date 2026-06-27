---
id: HYP-3098
title: LRC14 two-frontier observer gluing ledger
status: SYNTHESIS / exact scout; not a proof
source: codex-2026-06-27-S258
tangent: T1177
technique: LTI-238
tournament_technique: LTT-136
script: 04-computation/lrc14_two_frontier_gluing_s258.py
result: 05-knowledge/results/lrc14_two_frontier_gluing_s258.out
related:
  - HYP-3097
  - HYP-3096
  - HYP-3095
  - HYP-3094
  - HYP-3093
  - HYP-3092
  - HYP-3091
  - HYP-3090
  - HYP-3089
  - HYP-3088
  - HYP-3085
  - HYP-3083
  - THM-576
  - THM-575
  - THM-573
  - OPEN-Q-108
---

# HYP-3098: LRC14 Two-Frontier Observer Gluing Ledger

## Claim

The current proof frontier should be pushed by alternating two obligations,
not by optimizing either in isolation:

```text
A = polynomial-method witness route
    direct lonely measure / component count / largest arc / finite denominator budget

B = Pascal-equivalence-scissors route
    pair-normalized cap skeleton / cap defect / endpoint and branch scissors
```

The overlap is the hidden object from HYP-3095: an observer-gluing packet.  A
chart may forget a coordinate only when another chart reconstructs it,
annihilates it dually, proves it fiber-constant, descends it, or names the
remaining debt.

The S258 exact scout records both charts on shared rows.  It is not a proof of
LRC14; it is a proof-interface compression test.

## Exact Signals

Task A progress: for fixed non-tight rows, the direct lonely set has exact
component data and therefore a finite denominator-net budget.  Examples from
the scout:

```text
cover tail 12->84:    largest_arc = 3/1960,  D_from_largest = 654
cover tail 12->168:   largest_arc = 23/11760, D_from_largest = 512
near/K33 12->36:      largest_arc = 1/2520,  D_from_largest = 2521
P10+K33:              largest_arc = 1/1960,  D_from_largest = 1961
```

Task A failure feeding Task B: divisor loading already makes the raw direct
time budget explode in the first tested loaded row:

```text
S_B, B=6: largest_arc = 1/5880, D_from_largest = 5881,
          no grid witness for d <= 14 in the scout
```

So raw time cannot be the invariant.  The missing coordinate is the
apex/ruler normalization plus the cap/scissors payload that records why the
loaded apex is a legal finite-principal-part debt rather than a new proof
shape.

Task B progress: the Pascal cap chart is pure pair mass through the first
three sector sizes, and its first failure is exactly one affine pair-mass unit:

```text
j=1,2,3: defect = 0
j=4:     defect = 1/4004, i.e. one unit after multiplying by 4004
j=5:     defect = 1081/76440 = 11891/210 pair units
```

Task B failure feeding Task A: pair mass and positive safe mass do not
distinguish proof routes.  K33 cross-handoff rows are positive-open too, so
the witness route must retain branch chart, active binders, endpoint owners,
and terminal debt before using any scalar cap or measure.

## Reframing

The proof object is neither:

- total lonely measure,
- a raw largest denominator,
- a Pascal count,
- a cap scalar,
- nor a branch label.

It is the overlap map:

```text
arithmetic/CRT chart
  <-> normalized-arc chart
  <-> Pascal-scissors cap chart
  <-> moment/Perron chart
  <-> branch/K33 chart
  <-> formal witness chart
```

The next useful theorem is therefore not just "positive measure" or "pairwise
cap."  It is a gluing statement:

```text
Every THM-573 residual packet either
  (1) has a normalized arc floor whose finite bad-denominator budget is
      compatible with its cap/scissors packet, or
  (2) is rerouted through nested-refinement O2 discharge, or
  (3) is rerouted through cross-handoff O3/K33 state-lift debt, or
  (4) names the first chart overlap that fails.
```

## Tournament Analysis

Tournament vertices are observer charts and proof obligations, not runners.
The pairwise observable is predicate retention, denominator-net survival,
CRT/lift debt, scissors payload, and branch-handoff debt.

The S258 gauge is transitive with Hamiltonian path:

```text
observer_gluing_packet
  > normalized_arc_chart
  > pascal_scissors_chart
  > level7_crt_chart
  > branch_k33_chart
  > safe_mass_scalar
  > raw_denominator_floor
  > raw_pair_count
```

The only edge flip against the chosen input order is
`level7_crt_chart > branch_k33_chart`, reflecting the fact that THM-573 closes
one large arithmetic face before the branch/K33 chart can become terminal.

## Next Ledger Fields

Add these fields to the first `lrc14_observer_gluing_ledger` rows:

```text
source_row_id
crt_c7_lift_status
crt_c2_dyadic_lift_status
direct_lonely_measure
direct_component_count
largest_direct_arc
denominator_net_threshold_D
pascal_pair_mass_unit
triangular_cap_shadow
cap_defect
sector_pair_scissors_signature
grid_class
active_binder_owner_word
endpoint_owner_transition_word
overlap_failure_chart
terminal_exit_or_named_debt
```

The front edge to attack next is either:

1. prove a normalized replacement for `largest_direct_arc` that remains stable
   under divisor loading, or
2. prove the `j=4,5` cap defects are exactly the S3/S4 reflection-Perron
   scissors debt needed by the moment chart.
