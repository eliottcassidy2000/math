---
id: HYP-3019
title: LRC14 binding-pair switch carrier
status: SYNTHESIS / exact switch-carrier scout; not a proof
source: codex-2026-06-26-S183
script: 04-computation/lrc14_binding_pair_switch_carrier_codex_20260626.py
result: 05-knowledge/results/lrc14_binding_pair_switch_carrier_codex_20260626.out
related:
  - THM-524
  - HYP-3007
  - HYP-3006
  - HYP-3002
  - HYP-2990
  - HYP-2963
  - LTI-144
  - OPEN-Q-108
---

# HYP-3019: LRC14 Binding-Pair Switch Carrier

## Claim

THM-524 reduces exact LRC gap maximization to pair crossings, but a pair
crossing is theorem-safe only when it keeps the other-runner clearance scan.

The useful proof carrier is therefore not:

```text
binding pair has large pair gap
```

It is the labelled switch packet:

```text
binding pair
sum/difference denominator lane
active blocker deck
other-runner clearance margins
grid/off-grid status
endpoint-owner and packet-family labels
```

A quotient may keep a binding pair only if the active blocker deck and all
clearance inequalities are retained, reconstructed, dual-annihilated, or routed
to an AP/GW boundary, covering branch, unit-petal branch, K33/state-lift branch,
or named F7 residual debt.

## Computation

Script:

```text
04-computation/lrc14_binding_pair_switch_carrier_codex_20260626.py
```

Stored output:

```text
05-knowledge/results/lrc14_binding_pair_switch_carrier_codex_20260626.out
```

The script enumerates THM-524 pair-crossing candidates in `[0,1/2]` for named
LRC14 rows, computes exact `M`, strict safe measure, grid witnesses, witness
switch counts, decoy pair switches, and active optimum labels.

Key exact readouts:

```text
AP13:
  M=1/14, strict_safe_mu=0, grid witnesses at all six units mod 14.
  Optimum active pairs are the literal complement pairs:
    (1,13), (5,9), (3,11).
  Raw pair-gap decoys: 95.

GW 12->24:
  M=1/14, strict_safe_mu=0, same three complement optimum pairs.
  Raw pair-gap decoys: 196.

Petal 10->20:
  M=2/27, strict_safe_mu=1/980, optimum active pair (7,20).

Petal 13->26:
  M=2/27, strict_safe_mu=1/182, optimum active pair (1,26).

K33 12->36:
  M=3/41, strict_safe_mu=1/1260, optimum active pair (5,36).

K33 family drop12,13 add26,36:
  M=3/37, strict_safe_mu=79/8190, optimum active pair (1,36).

Covering 12->84:
  M=7/89, strict_safe_mu=563/105105, no grid witness.
  This is the only off-grid-only witness row in the named audit.

Large-tail decoy stress:
  single far tail 12->200 has 1934 raw-pair decoy times.
  drop6 core add180 has 1788 raw-pair decoy times.
  covering 12->84 has 846 raw-pair decoy times.
```

The decoy counts are the important guardrail.  Many candidate times contain a
pair at distance at least `1/14`, sometimes distance `1/2`, while another
runner collapses the true minimum below threshold.  Thus raw pair gap is a
false quotient; clearance is the LRC predicate.

## Tournament Analysis

Vertices are switch/proof carriers, not runners:

```text
others_clear_binding_switch
active_blocker_deck
binding_denominator_lane
grid_section_unit_witness
complement_sum14_boundary
covering_big_flank_switch
state_lift_nonunit_pair
raw_pair_gap_shadow
snapshot_position_tournament
```

Pairwise observable:

```text
predicate retention
clearance retention
exact denominator scale
off-grid coverage
family separation
decoy resistance
proof maturity
```

Gauge:

```text
carrier A -> B when A wins a majority of observable coordinates;
ties use the declared Hamiltonian path.
```

Fingerprint:

```text
score_hist=[(0,1),(1,1),(2,1),(3,1),(5,3),(7,1),(8,1)]
directed_3cycles=1
scc_sizes=[3,1,1,1,1,1,1]
hamiltonian_path_count=3
tie_path=
  others_clear_binding_switch
  > covering_big_flank_switch
  > grid_section_unit_witness
  > binding_denominator_lane
  > active_blocker_deck
  > complement_sum14_boundary
  > state_lift_nonunit_pair
  > raw_pair_gap_shadow
  > snapshot_position_tournament
```

The nontrivial SCC is healthy: grid, denominator-lane, and active-blocker
information interact and should not be scalar-ranked independently.

## Assumption Challenge

Candidate vertex sets considered:

```text
runners, gaps, fixed circle sections, section boundaries, wall crossings,
residues, cover arcs, Fourier modes, matroid circuits, proof obligations,
and binding-pair switches
```

Chosen vertex set:

```text
binding-pair/proof switches
```

Preserved LRC predicate:

```text
existence of t with min_v ||v t|| >= 1/14, provided the clearance scan is kept
```

Destroyed by raw pair quotients:

```text
active blocker deck, endpoint owners, off-grid witness intervals, exact packet
family labels, and whether the promising pair is blocked by another runner
```

Challenged assumption:

```text
"Since THM-524 reduces the maximum to pair crossings, a large pair crossing
should be enough proof data."
```

The audit refutes that scalar version.  THM-524 gives the right finite support,
but the proof object is a switch with clearance and packet labels attached.

## Next Proof Target

Add the following fields to HYP-2963-style packet records and certificate
manifests:

```text
binding_switch_type
active_pair_shell
denominator_lane
active_blocker_deck
other_clearance_margin
decoy_pair_gap
grid_or_offgrid_status
switch_witness_count
```

Then run a fiber-mixing test: inside each HYP-2963 packet family, check whether
these switch fields are constant, reconstructible, dual-annihilated, or exactly
the coordinate that routes the packet to AP/GW, covering, petal, K33/state-lift,
or F7 debt.
