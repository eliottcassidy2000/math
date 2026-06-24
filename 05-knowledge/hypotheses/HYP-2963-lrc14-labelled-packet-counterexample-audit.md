---
id: HYP-2963
title: LRC14 labelled-packet counterexample audit
status: PROOF-CLASSIFIER / bounded labelled-packet audit refining HYP-2961 and supporting HYP-2962/HYP-2956; not a proof
source: codex-2026-06-24-labelled-packet-audit
related:
  - HYP-2962
  - HYP-2961
  - HYP-2960
  - HYP-2956
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
  - HYP-2940
  - HYP-2937
  - HYP-2908
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_labelled_packet_counterexample_classifier_codex_20260624.py
  - 05-knowledge/results/lrc14_labelled_packet_counterexample_classifier_codex_20260624.out
external:
  - https://arxiv.org/abs/2606.22636
---

# HYP-2963: LRC14 Labelled-Packet Counterexample Audit

## Claim

The breakthrough theorem should be a labelled packet theorem:

```text
primitive 13-speed residual S
  -> labelled packet P(S)
  -> one of finitely many theorem buckets:

     q-witness family
     AP/Goddyn-Wong boundary family
     unit-petal sporadic family
     K33/state-lift sporadic family
     covering boundary-moment family
     source-spectrum unknown sporadic
```

The packet must retain enough data that no false quotient can hide a
counterexample:

```text
P(S) =
  exact M(S), binding denominators, q_threshold,
  Farey branch and excess e=14p-q,
  strict Haar/Baire safe mass and component count,
  finite boundary debt,
  C27 shell transfer,
  S145 packet route/rank,
  K33/state-lift flag,
  source family and covering class.
```

The target theorem is:

```text
Every primitive LRC14 counterexample candidate emits P(S), and P(S)
lands in one of the safe labelled families or in an explicitly named
sporadic state-lift obligation.
```

In the bounded classifier run here, the unknown bucket is empty.

This is an audit-level refinement of the incoming HYP-2961 classifier and a
bounded companion to the fixed-margin HYP-2962 packet theorem and the HYP-2956
labelled packet classification theorem.  HYP-2961 names the conceptual
strict-counterexample families; HYP-2962 makes "family" mean a fixed-margin
labelled packet class; HYP-2956 names the F0-F7 theorem buckets; HYP-2963
supplies an exact labelled packet emitter over a bounded adversarial bank and
confirms that the bank routes through the same family/sporadic partition.  It
also complements the HYP-2955 packet-migration gauntlet by checking exact `M`,
Haar/Baire mass, C27 transfer, S145 route/rank, and K33/state-lift flags in the
same packet.

## Computation

Script:

```text
04-computation/lrc14_labelled_packet_counterexample_classifier_codex_20260624.py
```

Stored output:

```text
05-knowledge/results/lrc14_labelled_packet_counterexample_classifier_codex_20260624.out
```

Default run:

```text
single_limit=180
two_swap_limit=36
alias_depth=4
lcm_tail_max=5
workers=8
```

The audited bank contains named adversaries, shell aliases, magnitude liars,
covering rows, AP single-swaps through `180`, and AP two-swaps with added
values through `36`.

Main census:

```text
audited rows              21913
below threshold           0
tight rows                2: AP, GW 12->24
M<=2/27 low packets       7
unknown packets           0

Q-WITNESS                 14676
BOUNDARY-AP-GW            2
BOUNDARY-PETAL-SPORADIC   4
K33-STATE-LIFT            3
COVERING-MOMENT           7228
```

The `M<=2/27` layer remains the familiar S145/S148 packet:

```text
AP
GW 12->24
near/K33 12->36
petal 10->20
petal 13->26
P10+GW
P10+K33
```

The classifier also records a larger K33-family marker outside that low layer:

```text
drop(12,13)->add(26,36), M=3/37, q_threshold=14, safe_mu=79/8190.
```

This is useful because it says K33 is not only one near-miss point.  It is a
family label: once a nonunit K33 packet appears, the proof should route it to
the HYP-2908/THM-572 state-lift lane rather than try to discharge it as a raw
open interval.

## Family Definitions

### A. q-Witness Family

If `q_threshold(S)<=13`, the direct denominator witness gives:

```text
M(S) >= 1/q_threshold(S) > 1/14.
```

This is an infinite safe family.  It includes many residue-liars and shell
aliases that would look dangerous after a magnitude-blind quotient.

### B. AP/Goddyn-Wong Boundary Family

The boundary family is:

```text
q_threshold=14,
M=1/14,
strict Haar/Baire safe mass=0,
labelled packet route in {tight AP floor, tight GW floor}.
```

The bounded classifier finds exactly:

```text
AP
GW 12->24
```

These are tight atoms, not counterexamples.  They are the zero-open boundary
core that the global theorem must classify.

### C. Unit-Petal Sporadic Family

The unit-petal family carries C27 unit-visible `p=2` structure and rank-0
packet discharge.  In the audit this includes:

```text
petal 10->20
petal 13->26
P10+GW
drop(10,13)->add(20,26)
```

The first three are in the `M<=2/27` frontier.  The last is a larger positive
open petal splice, so it is not a counterexample but it shows the petal label
also forms a family.

### D. K33 / State-Lift Sporadic Family

The K33 family is the nonunit packet branch:

```text
near/K33 12->36
P10+K33
drop(12,13)->add(26,36)
```

The theorem role is not "K33 rows are counterexamples."  They are positive-open
or loose rows whose labelled nonunit incidence should construct the
TournamentStateLift endpoint.  If that lift has packet value `7`, THM-572 blocks
it.

### E. Covering Boundary-Moment Family

Rows with `q_threshold>14`, multiple-of-14 pressure, or positive strict Haar
interior route to the covering boundary-moment branch.  In the audit this is
the largest non-q bucket:

```text
COVERING-MOMENT rows = 7228
```

The intended global theorem here is the HYP-2954/HYP-2953 bridge:

```text
adaptive exact-period boundary packet
  -> boundary/missed-depth quotient
  -> gK8 / L_y moment positivity
  -> positive source interval or K33 state-lift.
```

### F. Source-Spectrum Unknown Sporadic

This is the only bucket that would represent a new obstruction:

```text
q_threshold>=14,
strict Haar/Baire safe mass=0,
not AP/GW,
not unit-petal,
not K33/state-lift,
not covering-moment positive.
```

The bounded classifier found:

```text
SOURCE-SPECTRUM-UNKNOWN = 0.
```

## Relation To The Swap-Chain Paper

Fu, Qin, and Wang's June 2026 paper
[Spectral Gap for the Binary Fixed-Margin Swap Chain](https://arxiv.org/abs/2606.22636)
proves an inverse-polynomial spectral gap by comparing the swap chain to
two-row heat-bath moves, reducing to a three-row inequality, and then splitting
the proof into a scalar count sector and Johnson-harmonic non-scalar sectors.

The LRC14 import is methodological, not a direct theorem transfer:

```text
binary matrix count sector      -> q_threshold, exact M, Farey excess
three-row reduction             -> AP/GW/petal/K33 local packet atlas
Johnson harmonic sectors        -> owner-labelled C27/K33/boundary modes
heat-bath comparison            -> packet-preserving local replacement moves
```

The lesson is that scalar counts become honest only after conditioning on the
labels that define the fiber.  That is exactly the labelled-packet discipline:
do not compare rows after forgetting C27 ownership, boundary status, or K33
incidence.

## Tournament Analysis

Vertices are counterexample-classification routes, not runners.

Pair observable:

```text
exact scale retention,
boundary retention,
packet-owner retention,
state-lift visibility,
covering visibility,
anti-scalar guard.
```

Switch:

```text
A -> B iff A's retention vector is lexicographically >= B's.
```

Tie Hamiltonian path:

```text
COUNTEREXAMPLE
-> Q-WITNESS
-> BOUNDARY-AP-GW
-> BOUNDARY-PETAL-SPORADIC
-> K33-STATE-LIFT
-> COVERING-MOMENT
-> SHELL-ALIAS-LOOSE
-> MAGNITUDE-LIAR-LOOSE
-> SOURCE-SPECTRUM-UNKNOWN
```

Fingerprint in the stored run:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
c3=0
scc=(1,1,1,1,1,1,1,1,1)
Hamiltonian paths=1
```

The classifier tournament is deliberately transitive.  Its job is not to model
chaos; its job is to state which packet labels are allowed to be forgotten
later in the proof.

## Missing Global Step

HYP-2963 does not prove LRC14.  It sharpens the global proof obligation:

1. Prove every primitive residual can be reduced to a packet with the fields
   listed above.
2. Prove the zero-open `q=14` packet core is exactly AP/GW after C27/unital
   labels are retained.
3. Prove unit-petal packets discharge by finite C27/two-block rigidity.
4. Prove nonunit K33 packets construct the HYP-2908/THM-572 state-lift or are
   already positive-open.
5. Prove covering packets have positive boundary-moment/gK8/L_y image, or else
   collapse into the same state-lift endpoint.

That is the labelled packet theorem.

## Assumption Challenge

This session explicitly rejects the assumption that tournament vertices must be
runners or arcs.  The considered vertex sets included runners, gaps, fixed
circle sections, section boundaries, wall-crossing events, residues, cover arcs,
Fourier modes, C27 shell transfers, K33 incidence packets, and proof
obligations.

The chosen quotient uses proof-route packet vertices.  It preserves the LRC14
predicate only by retaining exact scale, open-vs-boundary status, and owner
labels.  It destroys raw speed identity and therefore cannot be applied before
the q-threshold, Farey, Haar/Baire, C27, and K33 fields are attached.
