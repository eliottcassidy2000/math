---
id: HYP-2946
title: LRC14 Farey-perfect Kuratowski carrier
status: PROOF-INTERFACE / product-perfect and graph-minor guardrail, not a proof
source: codex-2026-06-24-S144
related:
  - HYP-2945
  - HYP-2944
  - HYP-2943
  - HYP-2942
  - HYP-2941
  - HYP-2940
  - HYP-2939
  - HYP-2938
  - HYP-2937
  - HYP-2935
  - HYP-2934
  - HYP-2933
  - HYP-2932
  - HYP-2931
  - HYP-2930
  - HYP-2908
  - THM-572
---

# HYP-2946: LRC14 Farey-Perfect Kuratowski Carrier

The prompt's graph-theoretic correction is the right guardrail:

```text
p/q  ->  K_{p,q},       p*q = |E(K_{p,q})|,
```

but forbidden-minor structure is not mediant averaging.  The disjoint union
`K5 + K3,3` is nonplanar because it already contains the two known obstruction
cores.  It is not a third minimal obstruction: deleting one component leaves
the other nonplanar component.

The graph iteration that matters is minor/subdivision transitivity:

```text
J <= H <= G  implies  J <= G.
```

So the LRC14 proof interface should keep three layers separate:

```text
exact M/Farey branch
K_{p,q} product incidence
K5/K3,3 minor-transitive obstruction core
```

## Computation

The script
`04-computation/lrc14_farey_perfect_kuratowski_codex_s144.py` stores output at
`05-knowledge/results/lrc14_farey_perfect_kuratowski_codex_s144.out`.

It verifies the key `F_3`/`F_4` split:

```text
F_3:  2/3 -> K_{2,3},  p*q = 6,   planar, first even-perfect product
F_4:  3/4 -> K_{3,4},  p*q = 12,  nonplanar by inherited K_{3,3}
```

Thus the first perfect-product Farey fraction is not the first
complete-bipartite Kuratowski wall, and the first complete-bipartite
Kuratowski wall is not a perfect product.

## Perfect-Number Edge Loads

The even perfect number formula

```text
2^(r-1) * (2^r - 1)
```

is exactly an edge count:

```text
|E(K_{2^(r-1), 2^r-1})|.
```

The first rows are:

```text
2/3       -> 6            planar K_{2,3}
4/7       -> 28           K_{3,3} inherited
16/31     -> 496          K_{3,3} inherited
64/127    -> 8128         K_{3,3} inherited
4096/8191 -> 33550336     K_{3,3} inherited
```

There is also a coarse product fiber warning: `1/N` has product `N`, so
product-perfect alone is too coarse.  The meaningful perfect-number lane is
the nontrivial Mersenne factor fraction

```text
2^(r-1)/(2^r-1),
```

not the star fraction `1/N`.

## Mediant Cross Terms

If `a/b` and `c/d` have mediant `(a+c)/(b+d)`, then the product graph edge
load is

```text
(a+c)(b+d) = ab + cd + ad + bc.
```

The terms `ad` and `bc` are typed cross-incidences.  They are not graph
averages.  For the LRC14 unit-excess route,

```text
mediant(1/14, 2/27) = 3/41,
```

and the edge-load decomposition is

```text
14 + 54 + 27 + 28 = 123 = |E(K_{3,41})|.
```

This is exactly the first unit-excess `K_{3,3}` wall from HYP-2932, the
near-miss `12->36`.

## New POKE Hypotheses

1. **F3/F4 decoupling.**  `2/3` in `F_3` is the perfect-product planar seed;
   `3/4` in `F_4` is the first complete-bipartite K33 wall.  Perfect product
   and product-minor obstruction are distinct signals.
2. **Mersenne Kpq ladder.**  Even perfect numbers are complete-bipartite edge
   loads.  After `r=2`, their nonplanarity is inherited from `K_{3,3}`, not a
   new obstruction species.
3. **Mediant cross-term ledger.**  A mediant product is parent edge load plus
   typed cross-incidence.  Any proof using mediants should keep those cross
   terms visible.
4. **Minor-transitive proof iteration.**  The graph iteration layer should be
   minors of minors or subdivisions of subdivisions, not weighted graph
   mixtures.
5. **Apex-N wall rule.**  For an apex `N` unit-excess chain `p/(Np-1)`, the
   first complete-bipartite wall is `3/(3N-1)`.  At `N=14` this recovers
   `3/41`.
6. **Carrier order.**  Exact `M`/Farey branch and C27/unital labels precede
   Kpq, perfect-number, mediant, polyhedral, tiling, pi, and flower analogies.

## Tournament Analysis

Tournament vertices are proof carriers, not graphs:

```text
exact M/Farey branch
C27/unital branch-local pair grammar
Kuratowski minor-transitive core
Kpq product incidence ledger
Mersenne/perfect edge-load lane
mediant cross-term ledger
polyhedral/tiling recursion labels
pi/flower/unital analogy labels
raw graph average or scalar product
```

The pairwise observable is:

```text
branch retention,
unit preservation,
obstruction minimality,
state-lift fit,
scalar-decoy resistance.
```

The resulting conservative tournament is transitive:

```text
exact M/Farey branch
> C27/unital branch-local pair grammar
> Kuratowski minor-transitive core
> Kpq product incidence ledger
> Mersenne/perfect edge-load lane
> mediant cross-term ledger
> polyhedral/tiling recursion labels
> pi/flower/unital analogy labels
> raw graph average or scalar product.
```

## Proof Target

Refine the HYP-2932/HYP-2934 unit-excess split as follows:

```text
p=1:  star parent, q-threshold loose
p=2:  planar two-block strip / C27 petal branch
p>=3: K_{3,3} incidence packet
```

The new addition is the guardrail that perfect-number products are only edge
loads.  They can stress the product side channel, but they do not supply a new
forbidden-minor theorem.  The actual graph-minor proof step must construct a
minor-transitive `K_{3,3}` or `K5` packet, or route the typed packet into the
HYP-2908/THM-572 state-lift endpoint.

This is not a proof of LRC14.  It is a POKE carrier ordering and a set of
new finite hypotheses to test against the current low-gap frontier.
