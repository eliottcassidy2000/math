# LRC14 Summand/Multiplicand Farey Bridge

**Session:** codex-2026-06-23-S133
**Artifacts:** `04-computation/lrc14_summand_multiplicand_bridge_codex_s133.py`,
`05-knowledge/results/lrc14_summand_multiplicand_bridge_codex_s133.out`,
HYP-2934.

S131 made the prompt's graph reading precise:

```text
p/q -> K_{p,q},  p*q = |E(K_{p,q})|.
```

S133's job was to reconnect that with the older summand graph and multiplicand
graph notes.  The important correction is that there are two additive-looking
objects in play.

The formal Farey pair `(p,q)` is a summand pair at node `p+q`.  This is useful
as a recursion ledger, but it is not the LRC second-gap shell.  The LRC shell is
the denominator itself when `q=2n-1`.  For `n=14`, that means `q=27`.

So the unit-excess chain is not just:

```text
1/13, 2/27, 3/41, ...
```

It is:

```text
1/13: q-threshold parent.
2/27: the exact C=27 summand-unit quotient.
3/41: the first product-incidence K33 wall.
```

That is the session's main split.

## What Changed

Before this pass, S131 had already separated `2/27` as planar and `3/41` as the
first `K_{3,3}` wall.  S133 sharpens the planar side: `2/27` is not merely
"not yet nonplanar."  It is exactly the HYP-2083/HYP-2161 object:

```text
C=27 = 2*14-1,
summand shells {a,27-a},
unit-visible shells plus gcd 3 and gcd 9 nonunit holes.
```

That means the two S131 petal rows `10->20` and `13->26` should not be routed
to the same proof machine as the `12->36` near-miss.  They belong to the
second-gap/petal/lift problem.

The near-miss `12->36` has `M=3/41`.  It has left the `C=27` shell and crossed
the first place where the complete-bipartite incidence expansion has three
left owners:

```text
3/41 -> K_{3,41}.
```

That is why it is the correct first candidate for a finite three-owner
obstruction packet.

## Multiplicand Graph Versus Kpq

The multiplicand graph and `K_{p,q}` also needed to be separated.

For `2/27`, the product node is:

```text
54 = 2 * 3^3,
factor fibers: (2,27), (3,18), (6,9).
```

So the ordinary product graph sees a rich `3`-adic tower.  But the complete
bipartite blow-up attached to the specific Farey fiber is only `K_{2,27}`.
There is no `K_{3,3}` wall.

For `3/41`, the product node is just:

```text
123 = 3 * 41,
factor fibers: (3,41).
```

The product node is simpler, but the incidence blow-up `K_{3,41}` has a
`K_{3,3}` minor.  This is a nice warning: divisor richness and incidence rank
are different coordinates.

## Proof Shape

The better LRC14 split is now:

```text
e=0:        AP/GW floor candidates.
e=1,p=1:    coarse 1/13 parent.
e=1,p=2:    C=27 second-gap summand-unit/petal branch.
e=1,p>=3:   K33 multiplicand-incidence branch.
```

The `p=2` branch wants shell machinery: unit/nonunit gcd strata, petal rigidity,
lift denominators, endpoint owners, and CRT conservativity.

The `p>=3` branch wants packet machinery: three independent owners, K33
incidence, and eventually a HYP-2908-style state lift into the forbidden
tournament-conflict category.

The practical theorem target is:

```text
Every remaining q=14 non-AP/GW atom either reduces to the C=27 p=2
shell branch, or crosses the p>=3 K33 packet wall.
```

This is still a proof interface, not a proof.  But it tells us where to put the
next force: do not try to make product nonplanarity do the whole job, and do
not send the `2/27` petals to the same machine as the `3/41` near-miss.

## Tournament Analysis

The tournament vertices were proof carriers, not runners:

```text
q-binding,
C27-shell,
Farey p+q,
product fiber,
Kpq blowup,
K33 packet,
raw apex iso.
```

The tournament is transitive:

```text
q-binding > C27-shell > Farey p+q > product fiber
          > Kpq blowup > K33 packet > raw apex iso.
```

The product-only gauge flips nine carrier pairs.  That is the anti-scalarization
lesson in miniature: product data is useful, but only after the denominator
shell and incidence labels have had their say.
