---
id: HYP-2941
title: LRC14 triangular affine-operator carrier
status: PROOF-INTERFACE / scalar equality refuted, labelled carrier retained
source: codex-2026-06-24-S139
related:
  - HYP-2940
  - HYP-2939
  - HYP-2938
  - HYP-2937
  - HYP-2936
  - HYP-2935
  - HYP-2934
  - HYP-2933
  - HYP-2932
  - HYP-2931
  - HYP-2930
  - HYP-2908
  - THM-572
---

# HYP-2941: LRC14 Triangular Affine-Operator Carrier

The prompt's equality

```text
x*(2x-1) = 2*log2(x) + 1 - x
```

is exactly equivalent to

```text
x^2 = 1/2 + log2(x).
```

This equation has no positive real solution.  The script
`04-computation/lrc14_triangular_affine_operator_codex_s139.py` stores output at
`05-knowledge/results/lrc14_triangular_affine_operator_codex_s139.out` and
computes the gap

```text
g(x)=x^2 - 1/2 - log2(x).
```

Its minimum occurs at

```text
x = 1/sqrt(2 ln 2) = 0.849321800288019...
g(x) = 0.456964333972033... > 0.
```

So the scalar equality is false.  The useful proof object is the labelled
carrier behind it: quadratic/product growth dominates logarithmic address
depth, and the composition order must not be erased.

## Cubic Factor Packet

The proposed cubic

```text
C(x)=((x-1)/2)*(x/2)*(x-1/2)
```

has roots `0`, `1/2`, and `1`.  It should be read as a signed three-wall packet
analogous to Bonferroni/depth labels, not as a substitute for `log2(x)`.
For example:

```text
x=3/4: C(x)=-3/256,  x^2-1/2-C(x)=19/256.
x=2:   C(x)=3/4,     x^2-1/2-C(x)=11/4.
```

The cubic packet can be useful only after its sign region and forgotten
geometry are named.

## Triangular and Perfect-Number Lane

The expression

```text
x*(2x-1)
```

is exactly the triangular number `T_{2x-1}`.  If `x=2^(r-1)` and `2^r-1` is a
Mersenne prime, then this triangular number is the even perfect number

```text
2^(r-1)*(2^r-1).
```

For LRC14 this belongs to the product side channel.  On the unit-excess chain

```text
p/q = p/(14p-1),
```

the product lane is

```text
p*q = p*(14p-1),
```

which is the same family as `p*(2p-1)` but at apex `14` rather than apex `2`.
It is not the binding denominator `q`; it is the area/coimage/Kpq side channel.

## Affine Composition Readout

Let

```text
a(x)=x/2,    b(x)=x+1.
```

For left-to-right words in `a,b`, each `b` contributes `2^{-h}`, where `h` is
the number of future halvings.  Thus the word records a depth profile, not just
letter counts.  The noncommutator is

```text
ab(x)-ba(x)=1/2.
```

For the staircase word `(ba)^n`, the `b` depths are

```text
n, n-1, ..., 1,
```

so the depth sum is the triangular number `T_n`.  For the block word `b^n a^n`,
the depths are all `n`, so the depth sum is `n^2`.  Same letter counts can
produce triangular, square, or zero depth totals depending on composition
order.

This is the LRC14 lesson: a scalar carrier destroys the order label that the
proof needs.  The affine-depth profile is a conservative analogue of the
current `p+q`/`p*q`/C27/K33 packet split.

## Tournament Analysis

Tournament vertices are proof carriers, not runners.  The pairwise observable is

```text
(branch retention, composition-order retention, unit preservation,
 finite checkability, state-lift fit, scalar-decoy resistance).
```

The conservative tournament is transitive:

```text
exact M/Farey branch
> C27 two-swap splice
> Kpq/K33 owner packet
> affine depth profile
> product triangular/perfect lane
> cubic barrier carrier
> log/quadratic curvature guardrail
> raw scalar equality numerology.
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
c3=0
```

## Proof Target

The refined interface is:

```text
After exact M/Farey branch, C27 transfer, and Kpq/K33 owner data are attached,
try to assign an affine-depth packet to any low-gap non-AP/GW residual.  Unit-
visible depth entries should force the C27 petal/two-swap splice; nonunit depth
entries should expose the K33/octahedral/Clebsch state-lift packet.
```

This is not a proof of LRC14.  It is a guardrail and a packet-construction
proposal: triangular/perfect-number structure is useful only as a labelled
product-depth carrier, while the false log/quadratic equality must not enter as
a theorem.
