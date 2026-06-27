---
id: HYP-3025
title: LRC14 closed arc-Cech nerve carrier
status: SYNTHESIS / exact topology proof-interface carrier; not a proof
source: codex-2026-06-26-S174
script: 04-computation/lrc14_arc_cech_nerve_carrier_codex_s174.py
result: 05-knowledge/results/lrc14_arc_cech_nerve_carrier_codex_s174.out
related:
  - HYP-3024
  - HYP-3023
  - HYP-3022
  - HYP-3021
  - HYP-3020
  - HYP-3018
  - HYP-3016
  - HYP-3015
  - HYP-3014
  - HYP-3013
  - HYP-3012
  - HYP-3011
  - HYP-3010
  - HYP-3009
  - HYP-3008
  - HYP-2997
  - HYP-2986
  - HYP-2975
  - HYP-2974
  - HYP-2970
  - HYP-2963
  - HYP-2990
  - THM-572
  - OPEN-Q-108
---

# HYP-3025: LRC14 Closed Arc-Cech Nerve Carrier

This supersedes the S174 reservation after S188 claimed `HYP-3024` for the
fiber-zipper convergence audit.  The arc-Cech carrier treats that convergence
audit as a predecessor side channel, not as the same proof object.

## Claim

The exact topology carrier for an LRC14 row is not the runner set, the sequence
name, or the runner-level nerve.  It is the closed Cech nerve of the individual
threshold danger arcs on the time circle, together with the open-cell tope
split and the boundary cocircuit facets.

At threshold `1/14`, each speed contributes several danger arcs.  Collapsing
all arcs owned by the same runner is a quotient.  That quotient is useful only
if it carries an explicit Betti-defect side channel:

```text
runner_quotient_betti_defect =
  |beta0(closed_arc_nerve) - beta0(closed_runner_nerve)|
+ |beta1(closed_arc_nerve) - beta1(closed_runner_nerve)|.
```

The closed arc nerve preserves the actual circular cover predicate.  For a
good cover of the circle by proper arcs, `beta1=1` is the full-cover signal;
`beta1=0` leaves positive safe components.  The open arc nerve preserves a
different predicate: it records which components are glued only by endpoint
equality walls.

## Computation

Script:

```text
04-computation/lrc14_arc_cech_nerve_carrier_codex_s174.py
```

Output:

```text
05-knowledge/results/lrc14_arc_cech_nerve_carrier_codex_s174.out
```

Named-row audit:

```text
AP                 closed_arc_betti=(1,1), open_arc_betti=(6,0), safe_mu=0
GW_12_to_24        closed_arc_betti=(1,1), open_arc_betti=(6,0), safe_mu=0
K33_12_to_36       closed_arc_betti=(2,0), open_arc_betti=(8,0), safe_mu=1/1260
petal_10_to_20     closed_arc_betti=(2,0), open_arc_betti=(8,0), safe_mu=1/980
petal_13_to_26     closed_arc_betti=(2,0), open_arc_betti=(8,0), safe_mu=1/182
P10_plus_K33       closed_arc_betti=(4,0), open_arc_betti=(10,0), safe_mu=4/2205
covering_12_to_84  closed_arc_betti=(8,0), open_arc_betti=(8,0), safe_mu=563/105105
fibbinary_first13  closed_arc_betti=(38,0), open_arc_betti=(40,0), safe_mu=66077/399840
moser_first13      closed_arc_betti=(64,0), open_arc_betti=(72,0), safe_mu=4264747/40348854
```

AP and Goddyn-Wong are the only named full covers.  They also share the same
boundary cocircuit current:

```text
boundary_safe_owner_sums_mod_14 = (0, 0, 0, 0, 0, 0).
```

The one-swap AP neighborhood through `add<=160` has:

```text
primitive_one_swaps=1911
zero_open_one_swaps=1
zero_open = drop 12, add 24
smallest_positive_safe_mu = 1/1260 at drop 12, add 36
```

Thus the closed arc nerve separates AP/GW from K33 in the same place as the
old exact safe-component audit, but it gives a more structural reason: K33 has
two closed danger components and no cover cycle.

## Tournament Analysis

Vertices are proof carriers, not runners:

```text
arc_cech_good_cover_nerve
endpoint_tope_cocircuit_wall
taut_owner_current
safe_component_interval_measure
runner_quotient_nerve
fejer_toeplitz_dual_certificate
automaton_gap_sidecar
raw_speed_or_sequence_scalar
```

Pairwise observable:

```text
predicate equivalence
topology exactness
quotient-defect visibility
endpoint locality
dual handoff
computability
maturity
```

The carrier tournament is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
score_path=endpoint_tope_cocircuit_wall
  > arc_cech_good_cover_nerve
  > taut_owner_current
  > safe_component_interval_measure
  > runner_quotient_nerve
  > fejer_toeplitz_dual_certificate
  > automaton_gap_sidecar
  > raw_speed_or_sequence_scalar
```

The score path puts endpoint topes slightly before the Cech nerve because the
endpoint wall language sees cocircuit ownership.  The Cech nerve supplies the
global cover topology that the wall language alone can fragment.

## Necessary Conditions

1. **Closed-cover condition.**  A zero-open LRC14 packet must have closed
   danger-arc Cech `beta1=1`, unless it is being discharged by a stronger
   dual certificate or a named residual sector.

2. **Boundary-gluing condition.**  If the open nerve has multiple components
   glued only at endpoints, the packet must retain boundary cocircuit facets
   and owner-current sums.  AP/GW have six such boundary points with owner
   sums `0 mod 14`.

3. **Runner-quotient condition.**  Any argument that collapses individual
   danger arcs to runners must retain `runner_quotient_betti_defect` or prove
   that the cover predicate is constant on the quotient fiber.

4. **Positive-row exit condition.**  Named positive rows have closed arc
   `beta1=0`, so they must exit by safe-component interval measure, Fejer /
   Toeplitz, Ramanujan / twist, Haar / discrepancy, K33 state-lift, or another
   labelled certificate.  They are not closed-cover equality atoms.

5. **F7 good-cover quotient-defect condition.**  A genuine future F7 residual
   should be forced to exhibit a runner-level homology shadow that cannot be
   lifted to the individual-arc good-cover nerve, or else supply a new
   state-lift / harmonic obstruction.

## Packet Fields

Add these to HYP-2963 packet records or a sidecar manifest:

```text
closed_arc_cech_beta
open_arc_component_count
boundary_cocircuit_facet_word
boundary_owner_sum_word_mod_14
runner_quotient_betti_defect
private_arc_count
private_runner_count
safe_tope_count
arc_cech_exit_route
```

These fields are meant to sit below the recent automaton and sequence-shadow
fields from HYP-3011/HYP-3012.  A finite-state word is allowed to help route a
packet only after the exact danger-cover topology is still attached.

## Proof Target

The next theorem shape is:

```text
Closed Arc-Cech Dichotomy.
Every primitive zero-open LRC14 packet either:
  (1) has the AP/GW closed arc-Cech cover cycle and boundary owner-current law,
  (2) emits a named K33/state-lift or Fejer/Ramanujan/Haar dual exit, or
  (3) is the first genuine F7 good-cover quotient defect.
```

This is not a proof of LRC14.  It is a better carrier for the proof: it makes
the actual full-cover topology exact before any runner, automaton, or scalar
quotient is allowed to forget it.
