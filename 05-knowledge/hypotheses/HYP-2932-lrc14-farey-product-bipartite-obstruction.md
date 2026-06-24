---
id: HYP-2932
title: LRC14 Farey-product bipartite obstruction ledger
status: PROOF-INTERFACE / new local split; not a proof of LRC14
source: codex-2026-06-23-S131
related:
  - HYP-2934
  - HYP-2933
  - HYP-2931
  - HYP-2930
  - HYP-2928
  - HYP-2926
  - HYP-2920
  - HYP-2917
  - HYP-2908
---

# HYP-2932: Farey products are complete-bipartite obstruction ledgers

The prompt's product observation gives the S130 payload `p*q` a precise
structural role:

```text
p/q  ->  K_{p,q},      p*q = |E(K_{p,q})|.
```

This does **not** replace the LRC14 denominator.  For an exact optimum
`M(S)=p/q`, the theorem-level gap remains

```text
M(S)-1/14 = (14p-q)/(14q).
```

Thus `q` is still the binding scale.  The product ledger instead records the
incidence depth of the Farey escape.  Since `K_{p,q}` is nonplanar exactly when
`min(p,q)>=3`, the first product-minor obstruction is the Kuratowski
`K_{3,3}` wall.

## Computation

The script
`04-computation/lrc14_farey_bipartite_obstruction_codex_s131.py` stores output
at
`05-knowledge/results/lrc14_farey_bipartite_obstruction_codex_s131.out`.

It verifies the ordinary Farey fact:

```text
F_3: nonzero new complete-bipartite ledgers stay planar.
F_4: 3/4 -> K_{3,4}, which contains K_{3,3}.
```

For the LRC14 unit-excess chain

```text
e = 14p-q = 1,      q = 14p-1,
```

the product ledger is:

```text
1/13 -> K_{1,13}   star / q-threshold parent
2/27 -> K_{2,27}   planar two-block strip
3/41 -> K_{3,41}   first K_{3,3} wall
```

The `3/41` wall is exactly the S128 near-miss row
`{1,...,11,13,36}`.

## Row-bank result

On the S130 `749`-row AP/GW/petal/single-replacement bank, the exact bucket
counts are:

```text
tight-floor     star       2
unit-excess     star      54
unit-excess     two-block  2
unit-excess     K33-wall   1
nonunit-excess  star     319
nonunit-excess  two-block 315
nonunit-excess  K33-wall  56
```

The unit-excess rows split cleanly:

```text
1/13  -> K_{1,13}: many coarse q-threshold parent rows
2/27  -> K_{2,27}: swap 10->20 and swap 13->26
3/41  -> K_{3,41}: near-miss 12->36 only
```

So the product ledger does not characterize tightness, but it does isolate the
first genuinely three-owner escape from the planar one- and two-owner escapes.

## Tournament Analysis

S131 deliberately challenges the assumption that tournament vertices must be
runners or arcs.  It considers fractions, bipartition sides, incidence edges,
`K_{3,3}` minors, LRC rows, residues, wall crossings, and proof ledgers.

The chosen proof-ledger tournament has vertices:

```text
q-binding, sum-recursion, Kpq-product, K33-wall, power-stress, raw-iso
```

The pairwise observable is a lexicographic role score:

```text
theorem safety,
additive locality,
product side-channel strength,
obstruction signal,
magnitude stress,
false-positive resistance.
```

The resulting tournament is transitive:

```text
q-binding > sum-recursion > Kpq-product > K33-wall > power-stress > raw-iso
```

This preserves which quotient can serve which part of the proof and destroys
exact runner geometry.  That is intentional: it is a proof-interface hierarchy.

## Proof target

HYP-2930 said every non-AP/GW `q=14` survivor should either have nonunit
excess or enter a non-tight unit-excess Farey child.  HYP-2932 refines the
unit-excess branch:

```text
e=0:       AP/GW floor candidates.
e=1,p=1:   coarse parent 1/13, already q-threshold loose.
e=1,p=2:   planar two-block strip 2/27, to be killed by petal/two-block rigidity.
e=1,p>=3:  K_{3,3} product-minor wall, a finite three-owner obstruction target.
```

The new theorem should **not** claim nonplanarity is itself a contradiction.
Instead, it should prove that any remaining `q=14` non-AP/GW atom either falls
into the `p=2` two-block strip handled by petal rigidity, or reaches `p>=3`
where the `K_{3,3}` ledger supplies the finite three-owner obstruction packet
that can feed the HYP-2908 state-lift / forbidden-H endpoint.

Follow-up HYP-2934 sharpens the `p=2` side: `2/27` is exactly the
`C=27=2*14-1` summand-unit shell branch, while `3/41` is the first branch that
leaves that shell and crosses the `K_{3,3}` incidence wall.
