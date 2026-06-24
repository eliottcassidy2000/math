---
id: HYP-2986
title: LRC14 certificate handoff atlas and zipper theorem target
status: PROOF-INTERFACE / labelled packet handoff atlas; not a proof
source: codex-2026-06-24-S164b
artifacts:
  - 04-computation/lrc14_certificate_handoff_atlas_codex_s164.py
  - 05-knowledge/results/lrc14_certificate_handoff_atlas_codex_s164.out
related:
  - HYP-2985
  - HYP-2984
  - HYP-2983
  - HYP-2982
  - HYP-2981
  - HYP-2980
  - HYP-2979
  - HYP-2978
  - HYP-2977
  - HYP-2976
  - HYP-2975
  - HYP-2974
  - HYP-2973
  - HYP-2972
  - HYP-2969
  - HYP-2968
  - HYP-2965
  - HYP-2963
  - HYP-2956
  - HYP-2908
  - THM-572
  - OPEN-Q-108
---

# HYP-2986: LRC14 Certificate Handoff Atlas

The current LRC14 frontier is no longer a search for one scalar obstruction.
It is a gluing problem.  The repo has strong local certificates:

```text
q-witnesses and twist ladders
Haar/Baire positive-open fronts
endpoint bridge and taut-current certificates
danger-count moment duals
Ramanujan exact-period projectors
Fejer/Toeplitz interval certificates
analytic-sieve/Kaczynski smoothing packets
HYP-2908/THM-572 state-lift endpoint
```

HYP-2986 says the next proof object should be a handoff atlas: a finite list
of proof carriers, packet families, and allowed quotient arrows.  An arrow is
valid only when it preserves the LRC predicate or explicitly retains the label
that would otherwise be destroyed.

## Computed Artifact

Script:

```text
04-computation/lrc14_certificate_handoff_atlas_codex_s164.py
```

Stored output:

```text
05-knowledge/results/lrc14_certificate_handoff_atlas_codex_s164.out
```

The script defines certificate carriers with retention vectors for:

```text
LRC predicate
exact scale
phase/period
topology
endpoint owners
packet family
dual certificate
formal checkability
residual routing
```

It then runs Tournament Analysis on proof carriers rather than runners.  The
carrier tournament has:

```text
score_hist={0:1,1:1,2:1,4:3,6:1,7:1,8:1}
directed_3cycles=1
SCC_sizes=[1,1,3,1,1,1,1]
Hamiltonian_path_count=3
```

The tie/retention Hamiltonian path is:

```text
labelled_source_packet
> exact_interval_fejer
> ramanujan_exact_period
> endpoint_bridge_graph
> twist_ladder_blocker
> danger_count_moment_dual
> analytic_sieve_kaczynski
> tournament_state_lift
> raw_scalar_shadow
```

The nontrivial SCC is useful rather than alarming: it records that Ramanujan
exact-period data, endpoint bridges, and twist ladders can overtake each other
depending on the gauge.  That is exactly why a proof must state which
predicate is preserved at each handoff.

## Handoff Table

The script records seven packet handoffs:

```text
qdiv<14 / direct q-witness
AP/Goddyn-Wong boundary atom
K33 / 12->36 and K33 splices
unit petal and P10+GW strip
covering boundary-moment packets
lcm-tail / finite-denominator wall
SOURCE-SPECTRUM-UNKNOWN / F7
```

Each handoff records:

```text
current carrier -> exit carrier
status
retained predicate
fields destroyed without a guardrail
open theorem arrow
```

This is the concrete quotient-discipline object requested by the user's
Robbins/Robin guardrail theme.  It turns "do not forget labels" into a finite
bridge contract.

## Zipper Theorem Candidate

The atlas isolates six open arrows:

```text
O1 source-kernel exclusion:
   every primitive row enters this atlas.

O2 formal interval backend:
   replace prototype pi/trig intervals by checkable balls.

O3 family compression:
   selected Fejer/twist certificates lift to packet templates.

O4 admissible smoothing:
   analytic-sieve/Kaczynski quotients retain approach labels.

O5 state-lift construction:
   any zero-open non-AP/GW residual builds HYP-2908/THM-572.

O6 F7 definition:
   if a Johnson-harmonic residual sector exists, name its preserved predicate.
```

The fixed-margin swap-chain paper `arXiv:2606.22636` is a useful imported
proof pattern here.  Its abstract proof shape is: keep row/column margins as
fibers, compare the full chain with a low-row heat-bath core, reduce to three
rows, then split the three-row function space into a count sector and Johnson
harmonic sectors.  The LRC14 translation is direct enough to be a guardrail:
source packets are fixed margins, the desired reduction is to a finite packet
core, and any surviving F7 sector must be a named harmonic residual with a
preserved LRC predicate rather than an anonymous failure of Fejer/twist/moment
certificates.

The proposed theorem is:

```text
If O1-O6 hold, then every primitive LRC14 row either has a strict
witness/dual certificate in the atlas, is the AP/Goddyn-Wong equality atom,
or constructs the forbidden HYP-2908/THM-572 state lift.
```

This would prove LRC14.  HYP-2986 does not prove O1-O6.  Its contribution is
that it localizes the remaining proof debt and forbids future scalar shortcuts
from hiding one of those six arrows.

## Creative Proof Reading

The new picture is a zipper rather than a ladder.  One side is primal:
twists, Haar-open intervals, endpoint bridges, and lift packets.  The other
side is dual: danger-count moments, Toeplitz/Fejer PSD failure, Ramanujan
primitive projectors, and analytic smoothing.  The zipper closes only when the
teeth match packet-by-packet.

AP/Goddyn-Wong are the zipper stops: they are equality atoms, not strict
counterexamples.  K33, petals, covering rows, and lcm tails are not separate
mysteries; they are places where a different tooth is load-bearing:

```text
K33       -> incidence/state-lift label or Fejer certificate
petal     -> C27 unit-hole/two-block label or Fejer certificate
covering  -> endpoint/moment/lift chart or Fejer certificate
lcm-tail  -> exact-period/twist ladder or Kaczynski approach label
F7        -> new sector, unless it constructs the forbidden state lift
```

The useful next proof pass is therefore not "find a better scalar."  It is to
choose one open arrow and make it formal.  The highest-leverage candidates are
O3 family compression for the Fejer certificates and O5 state-lift
construction for zero-open non-AP/GW residuals.

## Assumption Challenge

This pass explicitly rejects runner vertices as the default.  The relevant
vertices are:

```text
proof carriers
packet fibers
endpoint-owner bridges
exact-period modes
dual certificates
smoothing kernels
state-lift obligations
```

The quotient preserves the LRC strict-witness predicate only if it retains the
packet labels named in the handoff table.  It destroys raw runner order,
floating shadows, endpoint ownership, or denominator prime-power labels only
after those fields are reconstructed, annihilated by a dual certificate, or
parked in a residual bucket.
