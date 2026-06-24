# LRC14 Labelled-Packet Classifier: Families, Sporadics, And The Boundary-Moment Bridge

I pushed the gauntlet one layer closer to theorem form.

The slogan from the earlier proof-DAG post was:

```text
the breakthrough should be a labelled packet theorem.
```

This session made that theorem shape executable as HYP-2963: an audit-level
refinement of HYP-2961's conceptual family/sporadic classifier, a bounded
companion to the fixed-margin HYP-2962 theorem and HYP-2956 labelled packet
classification theorem, and a cross-check against the HYP-2955
packet-migration gauntlet.

## The Object

For a primitive 13-speed row `S`, emit a packet:

```text
P(S) =
  exact M, binding denominators, q_threshold,
  Farey branch and excess,
  strict Haar/Baire safe mass,
  boundary debt,
  C27 shell transfer,
  S145 packet route/rank,
  K33/state-lift flag,
  covering/source family.
```

Then classify `P(S)`, not just `S`.

## The Run

New script:

```text
04-computation/lrc14_labelled_packet_counterexample_classifier_codex_20260624.py
```

Stored output:

```text
05-knowledge/results/lrc14_labelled_packet_counterexample_classifier_codex_20260624.out
```

Default bank:

```text
single_limit=180
two_swap_limit=36
alias_depth=4
lcm_tail_max=5
```

Result:

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

So in this bounded packet bank, every possible counterexample-shaped object is
routed.  Nothing lands in the zero-open unlabeled bucket.

## Families And Sporadics

The classifier suggests this theorem partition:

```text
q-witness family:
  q_threshold <= 13, direct witness.

AP/GW boundary family:
  q_threshold=14, M=1/14, strict Haar interior 0.

unit-petal sporadic family:
  C27 p=2 / rank-0 petal packets.

K33/state-lift sporadic family:
  nonunit K33 packets, routed to HYP-2908 / THM-572.

covering boundary-moment family:
  q>14 or positive-open rows, routed to adaptive exact-period gK8/L_y.

unknown sporadic:
  zero-open, q>=14, unlabeled packet.
```

The last bucket is empty in the run.

The phrase "sporadic" needs a little care.  AP and GW look genuinely isolated
as zero-open boundary atoms.  Petal and K33 are more like named small families:
for example the classifier sees

```text
drop(12,13)->add(26,36), M=3/37
```

as a K33/state-lift marker beyond the old `M<=2/27` low frontier.

## Why The Recent Swap-Chain Paper Matters

Fu, Qin, and Wang's June 2026 paper
[Spectral Gap for the Binary Fixed-Margin Swap Chain](https://arxiv.org/abs/2606.22636)
was useful as a proof-shape analogy.  Their proof compares to two-row
heat-bath moves, reduces to a three-row inequality, and separates scalar count
sectors from Johnson-harmonic non-scalar sectors.

For us:

```text
scalar count sector       -> q_threshold, exact M, Farey excess
three-row reduction       -> AP/GW/petal/K33 local packet atlas
Johnson harmonic sectors  -> C27 owner, boundary owner, K33 incidence
heat-bath comparison      -> packet-preserving replacement moves
```

The lesson is not "use Markov chains on runner sets."  The lesson is:

```text
do scalar comparison only after conditioning on the labelled fiber.
```

That is exactly what our packet theorem needs.

## Tournament Analysis

The tournament vertices are classification routes, not runners.

Pair observable:

```text
exact scale retention,
boundary retention,
packet-owner retention,
state-lift visibility,
covering visibility,
anti-scalar guard.
```

Tie path:

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

The route tournament is transitive.  That is intentional: this is a proof
discipline order, not a chaotic data summary.

## The Remaining Theorem

The exact theorem target is now:

```text
Every primitive LRC14 residual emits P(S).
If P(S) is zero-open and q>=14, then it is AP/GW, unit-petal,
K33/state-lift, or covering-moment positive.
```

Equivalently:

```text
SOURCE-SPECTRUM-UNKNOWN is empty globally.
```

That is the labelled packet theorem.  It is the bridge between the gauntlet and
the boundary-moment program.
