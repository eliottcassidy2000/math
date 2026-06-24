---
id: HYP-2997
title: LRC14 cocycle normal-form atlas
status: PROOF-INTERFACE / cohomological synthesis and tournament carrier atlas; not a proof
source: codex-2026-06-24-S167
related:
  - HYP-2996
  - HYP-2995
  - HYP-2994
  - HYP-2993
  - HYP-2992
  - HYP-2991
  - HYP-2990
  - HYP-2989
  - HYP-2987
  - HYP-2986
  - HYP-2985
  - HYP-2978
  - HYP-2974
  - HYP-2969
  - HYP-2963
  - HYP-2942
  - HYP-2937
  - HYP-2930
  - HYP-2171
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_cocycle_normal_form_atlas_codex_s167.py
  - 05-knowledge/results/lrc14_cocycle_normal_form_atlas_codex_s167.out
---

# HYP-2997: LRC14 Cocycle Normal Form Atlas

HYP-2997 rewrites the current LRC14 proof stack in cocycle normal-form
language.  It sits after HYP-2992's cocycle-sheaf exactness target, HYP-2995's
packet-cocycle `omega_Q` atlas, HYP-2994's obstruction ledger, and the S168
HYP-2996 residual-section packet-grid claim: those name the complex, the broad
carrier family, the legal cochain exits, and the finite residual-section bank;
this pass asks for the first nonzero residual class once those exits fail.
HYP-2990 gave the controlled-kernel zipper rule, and HYP-2991 identified the
local Haar mixed coordinate

```text
zeta(T)=T00-T01-T10+T11.
```

The rebased POKE coordination layer supplies the terminal reading: F7 is now
the harmonic/state-lift residual sector, and the HYP-2963 Haar-grid sweep is
the finite packet bank where the first-nonzero cocycle ledger should be tested.

The next abstraction is that all load-bearing proof data should be treated as
cochains on a labelled packet complex.  A quotient is sound only when each
forgotten cochain descends, is a coboundary, is annihilated by a dual
certificate, or is routed to a named residual group.

## Definitions

The proof packet complex `P` is not the original runner set.  It is the finite
and semi-finite complex whose cells are proof states:

```text
0-cells: labelled rows, endpoint cells, exact M/Farey nodes, proof packets
1-cells: quotient moves, certificate handoffs, wall crossings, comparisons
2-cells: commuting proof squares, Haar squares, tournament triples, chart overlaps
```

A `k`-cochain is any observable assigned to `k`-cells of `P`: endpoint owner
currents, exact-period characters, Haar curls, Farey defects, C27 carries,
Toeplitz moment functionals, tournament comparison weights, or state-lift
residual labels.

The cocycle condition is:

```text
delta a = 0
```

on every declared cell whose boundary preserves the target LRC predicate.

The coboundary/certificate condition is:

```text
a = delta phi
```

on a quotient fiber, or else a dual certificate pairs to zero with `a`.

The no-free-slider rule becomes a cohomological quotient rule:

```text
If pi:P->Q forgets a direction, then that direction must be
fiber-constant, reconstructed from retained labels, a coboundary on the
fiber, killed by a dual certificate, or mapped to a named residual class.
```

## Computation

Script:

```text
04-computation/lrc14_cocycle_normal_form_atlas_codex_s167.py
```

Stored output:

```text
05-knowledge/results/lrc14_cocycle_normal_form_atlas_codex_s167.out
```

The script defines `13` cocycle carriers:

```text
labelled_packet_total_cocycle
haar_zipper_2cocycle
endpoint_owner_boundary_cocycle
farey_excess_mediant_1cocycle
c27_carry_lift_1cocycle
fejer_toeplitz_dual_coboundary
ramanujan_exact_period_character_cocycle
tope_cocircuit_wall_cocycle
tournament_path_h1_cocycle
boundary_moment_multichart_cocycle
state_lift_obstruction_class
curried_section_derivative_cocycle
raw_scalar_shadow
```

For each carrier it records:

```text
degree
base space
cochain
cocycle equation
coboundary / certificate / residual exit
LRC pull
first nonzero counterexample family
anchors
retention vector
```

## Tournament Analysis

Tournament vertices are cocycle channels / proof obligations, not runners.
The pairwise observable is coordinate-wise comparison of the retention vector:

```text
predicate
base_fiber
gauge_invariance
coboundary_test
dual_annihilator
local_to_global
formalizable
residual_named
anti_scalar_guard
```

The switch/gauge orients `A -> B` when `A` has more strictly larger retention
coordinates than `B`; ties use the declared Hamiltonian path.

The S167 fingerprint is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:1}
directed_3cycles=0
SCC_sizes=[1,1,1,1,1,1,1,1,1,1,1,1,1]
Hamiltonian_path_count=1
```

Top-to-bottom carrier order:

```text
labelled_packet_total_cocycle
haar_zipper_2cocycle
farey_excess_mediant_1cocycle
endpoint_owner_boundary_cocycle
tope_cocircuit_wall_cocycle
c27_carry_lift_1cocycle
state_lift_obstruction_class
fejer_toeplitz_dual_coboundary
ramanujan_exact_period_character_cocycle
boundary_moment_multichart_cocycle
curried_section_derivative_cocycle
tournament_path_h1_cocycle
raw_scalar_shadow
```

This transitivity is intentional and conservative.  It says that when the
vertex set is changed from "proof methods" to "obstruction classes", the
packet cochain itself is the controlling object, raw scalar shadows are the
negative control, and all other lanes are typed first-nonzero channels.

## Repository Scan

The script scanned `16` recent and historical anchors, including HYP-2991,
HYP-2990, HYP-2989, HYP-2987, HYP-2986, HYP-2985, HYP-2978, HYP-2974,
HYP-2969, HYP-2963, HYP-2942, HYP-2937, HYP-2930, `beta1_upper_bound.out`, and
`heisenberg_frustration_bridge.out`.

Largest marker-hit channels:

```text
tournament_path_h1_cocycle                 428
tope_cocircuit_wall_cocycle                348
endpoint_owner_boundary_cocycle            343
fejer_toeplitz_dual_coboundary             281
c27_carry_lift_1cocycle                    218
curried_section_derivative_cocycle         212
ramanujan_exact_period_character_cocycle   206
haar_zipper_2cocycle                       192
state_lift_obstruction_class               163
farey_excess_mediant_1cocycle              109
```

The scan is not evidence for correctness, but it shows that the repo already
uses cocycle-like objects everywhere: endpoint currents, tournament `H^1`,
Haar curls, C27 carries, Fourier dual coboundaries, exact-period characters,
and named state-lift residuals.

## AP/Goddyn-Wong Closure Profile

In this language, AP and Goddyn-Wong are not merely two rows with the same
small scalar invariants.  They are boundary packets where all current
cocycle channels close before the state-lift residual:

```text
qdiv/Farey      e=14p-q is zero at the tight 1/14 boundary
endpoint/tope   regular-open safe mass is zero; boundary cocircuit debt remains
Haar zipper     mixed zeta stops at a boundary atom
C27/Jacobsthal  AP is base tiling; GW is the unique gated D3 acceleration
Fejer/Toeplitz  PSD-failure certificates vanish because there is no open safe interval
tournament      raw iso class is insufficient; labelled optimum/Farey packet is useful
state lift      no F7/THM-572 residual is invoked
```

The common DNA is therefore not "same residues" or "same tournament class".
It is closed cocycle accounting.

## Counterexample Family Normal Form

A primitive LRC14 counterexample would need a first nonzero class:

```text
F_zeta       unpaired mixed Haar curl after all zipper teeth fail
F_wall       no open tope and no AP/GW boundary cocircuit
F_farey      first nonzero excess class beyond the 1/14 boundary
F_carry      nonunit C27 monodromy, usually K33/D9
F_psd        harmonic moment cone fails to certify or fails to descend
F_period     exact-period character survives divisor quotient
F_cocircuit  wall packet not discharged by open-tope or equality atom
F_tournament non-potential tournament pressure class
F_chart      boundary-moment transition class fails to close
F_lift       irreducible named state-lift obstruction
F_section    frozen-coordinate derivative refuses to commute
```

If every listed class is a coboundary, dual-annihilated, or exits to AP/GW,
the packet has no remaining typed place to hide.  The only legal last bucket is
then the named `state_lift_obstruction_class`, not an anonymous exception.

## Assumption Challenge

This pass considered runners, arcs, fixed circle sections, section boundaries,
endpoint walls, residues, cover arcs, Fourier modes, Haar tiles, divisor
packets, matroid cocircuits, chart overlaps, proof obligations, and cocycle
channels as possible Tournament Analysis vertices.

Chosen vertices:

```text
cocycle channels / proof obligations
```

Preserved predicate:

```text
strict witness, AP/GW equality, positive certificate, or named state-lift residual
```

Destroyed information:

```text
raw runner identity and continuous scalar order
```

Retained information:

```text
exact M/Farey, endpoint owners, carry, Haar curl, period labels, residual class
```

The challenged assumption is that the useful tournament must live on the
objects of the original problem.  Here the useful tournament lives on the
obstruction classes.

## Theorem Target

Candidate cocycle normal-form theorem:

```text
For every reduced LRC14 packet, the retained proof data form a cocycle in the
packet complex.  If the packet has no strict safe interval, then each
non-boundary cocycle is either a coboundary on the quotient fiber, killed by a
dual certificate, transferred through a typed zipper tooth, or mapped to the
named F7/THM-572 obstruction class.  AP and GW are exactly the zero-open
boundary packets where all channels close before F7.
```

Next finite test:

```text
Build HYP-2963 packet-level cocycle ledgers and record the first nonzero class
for every low-frontier row.  No row should land in raw_scalar_shadow.
```
