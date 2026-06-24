---
id: HYP-2956
title: LRC14 labelled packet classification theorem
status: PROOF-SYNTHESIS / family-sporadic classification target; not a proof
source: codex-2026-06-24-S151
related:
  - HYP-2962
  - HYP-2961
  - HYP-2960
  - HYP-2955
  - HYP-2954
  - HYP-2953
  - HYP-2952
  - HYP-2951
  - HYP-2950
  - HYP-2949
  - HYP-2948
  - HYP-2947
  - HYP-2942
  - HYP-2937
  - HYP-2908
  - THM-523
  - THM-566
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_labelled_packet_theorem_gauntlet_codex_s151.py
  - 05-knowledge/results/lrc14_labelled_packet_theorem_gauntlet_codex_s151.out
external:
  - https://arxiv.org/abs/2606.22636
---

# HYP-2956: LRC14 Labelled Packet Classification Theorem

## Claim

The current LRC14 route should be stated as a labelled packet theorem:

```text
Every primitive 13-speed LRC14 row lies in F0-F5, or is an F6
non-migrating boundary-moment kernel, unless a new F7 Johnson-harmonic
packet sector exists.
```

Consequently, an actual counterexample `M(S)<1/14` must be:

```text
qdiv(S)>=14,
O(S)=empty,
not AP/GW boundary-only,
not a unit-petal/two-block packet,
not a positive K33/open-front packet,
and zero/nonpositive under the boundary-moment bridge.
```

The present gauntlets find no such packet in the audited AP-neighborhood atlas.

Rebase note: HYP-2960 supplies the `qdiv=14` skeleton-gate subclassifier under
this packet theorem, HYP-2961 sharpens the S151 F0-F7 language for strict
counterexamples, and HYP-2962 supplies the executable fixed-margin
packet-signature classifier.  In HYP-2961, AP/GW become equality sporadics,
`qdiv<=14` is a discharge rather than a live family, and the remaining
strict-counterexample work is the L1-L5 decision tree.

## Packet Families

The family split is:

```text
F0 qdiv witness discharge:
  qdiv(S)<14, so t=1/qdiv gives M(S)>1/14.

F1 AP/GW boundary atoms:
  qdiv=14, O(S)=empty, C(S) nonempty, AP/GW owner skeleton.
  These are tight rows, not counterexamples.

F2 unit-petal/two-block discharge:
  positive Haar front with only P10/P13/GW unit labels.

F3 K33 state-lift packet:
  K33/nonunit owner survives after positive or zero-open filtering.
  Positive seeds are already safe but must retain the K33 label for induction.

F4 unlabelled q14 positive front:
  qdiv=14 and O(S) has positive Haar mass without a known petal/K33 key.

F5 covering positive / boundary-moment strictness:
  qdiv>14 or a 14-multiple is present, with O(S) positive.

F6 covering zero-open non-migration kernel:
  qdiv>14, O(S)=empty, C(S)=empty, or a non-AP/GW zero-open kernel.
  This is the only current counterexample-shaped family.

F7 new Johnson-harmonic packet sector:
  a fixed-margin packet-swap audit reveals a non-scalar sector not covered by
  F1-F6.  This is a direct falsifier for the labelled packet theorem.
```

## Exact Sporadic Seeds

`04-computation/lrc14_labelled_packet_theorem_gauntlet_codex_s151.py` audits
the current named seeds with exact S124 `M/q` and S146 Haar-front kernels:

```text
AP, GW 12->24            -> F1
12->36, P10+K33          -> F3
P10, P13, P10+GW         -> F2
12->26                   -> F0
12->96                   -> F4
12->84                   -> F5
```

No named seed lands in F6 or F7.

S150 supplies the larger local gauntlet:

```text
covered qdiv>=14 rows:         0
non-AP/GW boundary-only rows:  0
```

with packet histogram:

```text
AP/GW-boundary-source          2
positive-K33-state-lift        340
positive-covering-off-apex     16488
positive-unit-petal/GW-strip   9632
positive-unlabelled-q14-front  41906
```

## Import From arXiv:2606.22636

Fu, Qin, and Wang prove an inverse-polynomial spectral gap for the lazy swap
chain on binary matrices with prescribed row and column sums.  The relevant
proof pattern is:

```text
binary relation with fixed margins
  -> local 2x2 swaps / two-row heat-bath comparison
  -> reduction to three rows
  -> scalar count sector + Johnson harmonic sectors.
```

For LRC14, the proposed translation is:

```text
binary relation:
  speeds versus packet lenses

packet lenses:
  q clocks, boundary owners, C27 shells, unital pairs, K33 owners,
  exact-period packets, boundary-moment signs

fixed margins:
  source-spectrum labels preserved by packet reductions

scalar count sector:
  qdiv, exact M/Farey node, Haar mass

Johnson harmonic sectors:
  non-scalar owner/carry/state-lift labels
```

The proof obligation is not to import the spectral-gap theorem as a black box.
It is to mimic its architecture: prove that every fixed-margin packet component
reaches F0-F5 by local swaps, unless a three-row Johnson sector is a genuine F7
falsifier.

## Missing Lemmas

The classification reduces the remaining proof to four lemmas:

```text
1. Fixed-Margin Packet Encoding.
   Every reduced LRC14 residual has a binary packet matrix whose margins retain
   the source-spectrum predicate and the labels used by HYP-2953/HYP-2954.

2. Packet Swap Connectivity.
   Within fixed packet margins, local row/column swaps connect all rows without
   changing the target family except through declared boundary faces.

3. Three-Row Johnson Audit.
   Every three-row packet component splits into scalar q/Farey/Haar count
   sectors and non-scalar Johnson sectors already represented by C27/unital/K33
   labels.  A missing sector is exactly F7.

4. Boundary-Moment Kernel Positivity.
   Every F6 covering zero-open non-migrating packet has positive
   boundary-moment/gK8 slack or constructs the HYP-2908/THM-572 state lift.
```

If these four lemmas hold, LRC14 follows from the packet classification.

## Tournament Analysis

Vertices are packet families and proof obligations, not runners or arcs.

Assumption challenge:

```text
alternate vertex sets considered:
  runners, gap intervals, fixed denominator-14 circle sections,
  section-boundary owner events, wall-crossing events, residues,
  cover arcs, Fourier modes, matroid/circuit-like packet dependencies,
  and proof obligations.

chosen vertex set:
  proof obligations and packet families.

preserved LRC predicate:
  the observer-source condition "some phase makes every speed at distance
  >=1/14", together with exact scale/Haar status and the labels needed to
  route a row to discharge, boundary atom, K33/state lift, or covering moment.

destroyed information:
  raw runner order, many magnitudes inside an already-safe open packet,
  and unlabelled apex tournament shadows.

challenged assumption:
  the hard object is not a tournament on runners; it is a fixed-margin
  labelled packet component whose non-scalar sectors must be classified.
```

Pairwise observable:

```text
which object preserves the observer-source LRC predicate and the labels needed
to route a non-migrating packet.
```

Switch/gauge:

```text
majority over source retention, boundary retention, non-scalar sector
visibility, endpoint strength, finite auditability, and anti-scalarization.
```

S151 gives a transitive tournament:

```text
F6_boundary_moment_kernel
> three_row_Johnson_sectors
> fixed_margin_packet_swaps
> K33_state_lift_endpoint
> boundary_owner_rigidity
> unit_petal_discharge
> qdiv_scalar_gate
> raw_row_or_scalar
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3cycles=0
```

## Status

This is a rigorous classification target, not a completed proof.  Its value is
that it names exactly what a counterexample must look like: an F6
boundary-moment kernel, unless an F7 Johnson sector exposes a missing label.
