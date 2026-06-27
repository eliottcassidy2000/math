---
id: HYP-2951
title: LRC14 Haar-Baire taut boundary finite check
status: FINITE SUPPORT / boundary-support lemma scout, not a proof
source: codex-2026-06-24-S146
related:
  - HYP-2948
  - HYP-2949
  - HYP-2950
  - HYP-2947
  - HYP-2942
  - HYP-2940
  - HYP-2937
  - HYP-2932
  - HYP-2920
  - HYP-2908
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_haar_baire_taut_boundary_s146.py
  - 05-knowledge/results/lrc14_haar_baire_taut_boundary_s146.out
---

# HYP-2951: LRC14 Haar-Baire Taut Boundary Finite Check

HYP-2948 suggested the boundary-support lemma:

```text
every reduced threshold-safe strict-Haar-zero row is AP or GW.
```

HYP-2951 adds the first finite Haar-Baire taut-front check around the AP/GW
frontier.  It turns the HBT*/Haar-Baire Wave idea into exact interval fronts:
strict safe components are open intervals, and their endpoints are labelled by
the active boundary owners.

## Computation

Script:

```text
04-computation/lrc14_haar_baire_taut_boundary_s146.py
```

Stored output:

```text
05-knowledge/results/lrc14_haar_baire_taut_boundary_s146.out
```

## Main Finite Results

One-swap AP neighborhood with added value `<=160`:

```text
rows = 1912
boundary_only = 2
positive_open = 1910
covered = 0
```

The only boundary-only rows are:

```text
AP
12->24 = H12:g3 -> D3:g3@24
```

The smallest positive strict Haar mass is:

```text
12->36 = H12:g3 -> D9:g9@36, mass = 1/1260.
```

Two-swap AP neighborhood with both added values in `14..40`:

```text
rows = 27378
boundary_only = 0
positive_open = 27378
covered = 0
```

The two smallest positive rows are exactly the S138 splices:

```text
(10,12)->(20,24), mass = 1/980
(10,12)->(20,36), mass = 4/2205
```

Thus the bounded S138 two-swap neighborhood contributes open Haar-positive
fronts, not new boundary-only packets.

## Boundary-Owner Skeleton

AP and GW have the same six threshold boundary points and the same active owner
pairs:

```text
t=1/14   owners=(1L,13R)
t=3/14   owners=(5L,9R)
t=5/14   owners=(3L,11R)
t=9/14   owners=(3R,11L)
t=11/14  owners=(5R,9L)
t=13/14  owners=(1R,13L)
```

The Goddyn-Wong transfer `H12 -> D3@24` is therefore invisible to the boundary
owner skeleton.  This is exactly why the C27/unital labels from HYP-2937 and
HYP-2942 must remain attached: Haar/Baire boundary owners alone know that GW is
endpoint-tight, but not why the hidden `12->24` transfer is legal.

## Taut Front Readout

The first open rows expose named fronts:

```text
near 12->36:
  (29/70,209/504)   left=5L   right=36R  mid_slack=1/1008
  (295/504,41/70)   left=36L  right=5R   mid_slack=1/1008

petal 10->20:
  (29/98,83/280)    left=7L   right=20R  mid_slack=1/560
  (197/280,69/98)   left=20L  right=7R   mid_slack=1/560

petal 13->26:
  (1/14,27/364)     left=1L   right=26R  mid_slack=1/728
  (337/364,13/14)   left=26L  right=1R   mid_slack=1/728
```

The two-swap rows combine these fronts rather than producing new
boundary-only support.

## HBT* Proof-Order Tournament

Tournament vertices are proof planners, not runners:

```text
HBT*_boundary_rank
Haar-Baire_Wave*
ANYA_interval_taut
CWave_relation_front
BlockA_local_database
Theta_owner_line
FieldD_measure_interp
raw_denominator_A*
```

The pairwise observable is:

```text
endpoint_witness,
haar_interior,
boundary_code,
relation_wave,
lrc_transfer,
dynamic_update,
anti_scalar_guard.
```

The tournament is transitive:

```text
HBT*_boundary_rank
> Haar-Baire_Wave*
> ANYA_interval_taut
> CWave_relation_front
> BlockA_local_database
> Theta_owner_line
> FieldD_measure_interp
> raw_denominator_A*.
```

## New Proof Target

The finite data suggests a sharper lemma:

```text
Boundary-owner skeleton rigidity:
if an AP-neighborhood row has zero strict Haar mass but nonempty threshold
support, then its boundary-owner skeleton is the AP/GW six-pair skeleton,
and the only legal hidden replacement is 12->24 after C27/unital labelling.
```

This does not prove LRC14.  It gives a more concrete finite endpoint for the
HYP-2947/HYP-2948 measurable rank lane:

```text
open front       -> discharge by Haar/Baire strictness
boundary skeleton -> route to C27/unital/state-lift labels
```

## Relation To HYP-2949 And HYP-2950

HYP-2949 is the broader S147 Baire-Haar any-angle carrier.  HYP-2950 reserves
the adversarial counterexample gauntlet.  HYP-2951 is the finite local database
piece suggested by both: it supplies exact AP/GW neighborhood boundary checks
and active-owner skeletons that the gauntlet can use as pass/fail tests.
