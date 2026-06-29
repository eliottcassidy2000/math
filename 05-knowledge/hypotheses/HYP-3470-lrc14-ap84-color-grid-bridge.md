---
id: HYP-3470
title: LRC14 AP84 phase-color grid bridge
status: EVIDENCE / exact AP-tail color-grid discrepancy bridge; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-2593, HYP-2594, HYP-3461, HYP-3460, HYP-3459, HYP-3458, HYP-3456, and HYP-3457
tangent: T1430
technique: LTI-430
tournament_technique: LTT-330
script: 04-computation/lrc14_ap84_color_grid_bridge_codex_20260629.py
result: 05-knowledge/results/lrc14_ap84_color_grid_bridge_codex_20260629.out
reflection: 07-reflections/lrc14-ap84-color-grid-bridge-codex-20260629.md
related:
  - HYP-3461
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
  - HYP-3431
  - HYP-2595
  - HYP-2594
  - HYP-2593
  - THM-523
  - OPEN-Q-108
---

# HYP-3470: LRC14 AP84 Phase-Color Grid Bridge

## Claim

HYP-3470 connects the older phase-color CRT reservoir of HYP-2593/HYP-2594 to
the current AP84 tail sidecars HYP-3454/HYP-3456/HYP-3457.  It is the exact
canonical AP84 `q=14V` placement sidecar under HYP-3459's broader color-packet
legality audit, complementary to HYP-3460's phase-branch pullback, downstream
of HYP-3461's colored-extension gate carrier, and the placement sibling of
HYP-3458's coloring-recursion sidecar.

For

```text
S_m = {1,2,...,11,13,84m},
```

write the HYP-2593 form as

```text
P={1,2,...,11,13},   E={0},   V=84m.
```

The continuous phase-color reservoir collapses:

```text
color 0 has measure 0;
colors 1..13 have the same four live intervals.
```

Those live intervals are

```text
[15/182,13/154]       length 2/1001
[29/70,41/98]         length 1/245
[57/98,41/70]         length 1/245
[141/154,167/182]     length 2/1001
```

with live color measure `426/35035`, total `sigma=426/2695`, and
`K=52`.  Thus the reservoir mass itself does not see the HYP-3456 AP-tail
period-`35` component clock.  The color information lives in the shifted grid
and the closed-boundary correction.

## Exact Readout

Script:

```text
04-computation/lrc14_ap84_color_grid_bridge_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_ap84_color_grid_bridge_codex_20260629.out
```

The script validates the closed interval floor/ceiling count against the
original HYP-2593 CRT predicate on sample rows through `m=70`:

```text
direct_validation_failures=[]
```

The exact total color-grid count `A(m)` has affine period `385`:

```text
A(m+385)-A(m)=5112
total_period385_failures_m_1_to_770=[]
```

The first period-`35` shifts are not constant:

```text
shift_by_35_first35 =
[464,466,466,464,464,466,466,464,464,464,464,464,466,466,464,464,466,466,
 464,464,464,464,464,466,466,464,464,466,466,464,464,464,464,464,466].
```

The residual

```text
R(m)=5112m-385A(m)
```

has checked one-period range:

```text
min=-1542 at m=349
max=882 at m=91
```

Individual live colors need a finer sidecar.  Their uniform affine period is
`5005`:

```text
color_period5005_failures_m_1_to_100=[]
```

while one total period `m=1 -> 386` shifts the colors unevenly:

```text
(0,392,393,393,393,394,394,394,394,394,393,393,393,392).
```

The closed-boundary bonus is a small `7`-clock:

```text
boundary_bonus_hist_m_1_to_385={0:55,2:330}
boundary_bonus_rule_failures_m_1_to_770=[]
bonus=0 iff 7 divides m; otherwise bonus=2.
```

## Proof Pull

HYP-3456 derives the AP-tail escape count from fixed low corridors and moving
high safe gaps.  HYP-3459 says the legal AP84 quotient must keep the color
packet rather than one raw color.  HYP-3460 pulls phase colors back to branch
and gate carriers.  HYP-3461 promotes the AP84 endpoint/floor packet into the
colored-extension gate carrier language.  HYP-3458 records the retained `35`-state
coloring-recursion packet for the escape count.  HYP-3470 is the parallel
CRT-placement sidecar: it supplies the exact `q=14V` colored grid count when
the AP-tail splice wants actual colored grid witnesses rather than branch-union
components.

The useful theorem target is now narrow:

```text
Prove the four fixed live intervals for P={1..11,13}, E={0};
prove the closed color-grid floor/ceiling count;
carry the period-385 total discrepancy and period-5005 color sidecar only when
actual q=14V placement is used.
```

This prevents an illegal scalarization: the live layer mass is the same in
thirteen colors, but the individual colors still carry a real residue sidecar.

## Tournament Analysis

Vertices are color-grid proof carriers, not runners or raw arcs.

```text
pairwise_observable =
  CRT exactness + grid formula + discrepancy period + color sidecar
  + HYP-3456 compatibility
score_hist={17:1,47:1,48:1,52:1,55:1,58:1,60:1}
directed_3cycles=0
hamiltonian_path =
  exact_H2593_CRT_color_predicate
  -> closed_color_grid_floor_formula
  -> total_period385_discrepancy_clock
  -> color_vector_period5005_sidecar
  -> mod7_boundary_bonus_gate
  -> H3456_period35_component_bridge
  -> raw_live_layer_mass_scalar
```

Assumption challenge: runners, phase colors, color residues, fixed live
intervals, interval endpoints, boundary hits, AP-tail `m` residues, HYP-3456
high gaps, Fourier modes, and proof obligations were considered.  The chosen
quotient preserves the exact `q=14V` CRT placement predicate and color-grid
discrepancy, but destroys branch-cover geometry and component adjacency.  The
challenged assumption is that reservoir mass should detect the AP-tail
mod-`35` clock; in this AP84 specialization, mass is blind and the useful color
information is the boundary/grid discrepancy sidecar.
