---
id: HYP-2945
title: LRC14 Farey mediants and graph minors do not commute
status: COMPUTATIONAL SCOUT / graph-operation guardrail; not a proof
source: codex-2026-06-24
tags: [lrc14, farey, kuratowski, graph-minors, mediants, bipartite]
related:
  - HYP-2944
  - HYP-2934
  - HYP-2932
  - HYP-2908
results:
  - 04-computation/lrc14_farey_product_perfect_kuratowski_codex.py
  - 05-knowledge/results/lrc14_farey_product_perfect_kuratowski_codex.out
---

# HYP-2945: Farey mediants and graph minors do not commute

The product reading:

```text
a/b -> K_{a,b},  a*b = |E(K_{a,b})|
```

is useful only if graph operations are kept distinct from Farey operations.

The decisive example is:

```text
2/3 -> K_{2,3}, planar
1/1 -> K_{1,1}, planar
mediant(2/3,1/1) = 3/4 -> K_{3,4}, contains K33.
```

So two planar Farey parents can mediant into the first complete-bipartite
nonplanar child.

This is not a paradox.  It says:

```text
mediant-taking is not minor-taking,
mediant-taking is not subdivision-taking,
edge-count averaging is not obstruction averaging.
```

Kuratowski/Wagner works because minor and subdivision containment are
transitive:

```text
J <= H <= G  implies  J <= G.
```

Farey mediants do not supply that transitive graph-obstruction order.

## Edge-Count Alias Guardrail

The script checks first Farey product appearances of several graph edge
counts:

```text
K33 edge count 9:              F9,  1/9 -> K_{1,9}, planar
K5 edge count 10:              F5,  2/5 -> K_{2,5}, planar
K5 disjoint K33 edge count 19: F19, 1/19 -> K_{1,19}, planar
```

The first actual `K33` minor in the Farey-product `K_{a,b}` reading is:

```text
F4: 3/4 -> K_{3,4}.
```

Therefore `K5`, `K33`, and `K5` disjoint `K33` are not points on a numerical
edge-count continuum.  `K5` and `K33` are the two irreducible nonplanarity
cores; the disjoint union is nonplanar only because it already contains them.

## LRC14 Use

The `K_{p,q}` packet in the LRC14 Farey branch should be used as a typed
incidence carrier:

```text
p=1: star/threshold parent
p=2: planar two-block/petal strip
p>=3: K33 three-owner packet
```

It should not be used as raw product numerology.  Product values such as `9`,
`10`, `19`, or `28` can appear as edge-count aliases without carrying the
corresponding obstruction core.

The live proof handoff remains:

```text
exact M/Farey branch
> C27 shell labels
> Kpq incidence side
> K33 state-lift packet.
```
