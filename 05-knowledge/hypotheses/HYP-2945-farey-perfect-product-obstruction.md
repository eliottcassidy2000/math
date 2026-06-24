---
id: HYP-2945
title: Farey perfect-product obstruction guardrail
status: PROOF-INTERFACE / arithmetic-graph synthesis; not a proof
source: codex-2026-06-24-S143
related:
  - HYP-2940
  - HYP-2932
  - HYP-2931
  - HYP-2930
  - HYP-2221
  - HYP-2220
  - HYP-2908
  - THM-572
---

# HYP-2945: Farey Perfect-Product Obstruction Guardrail

The prompt's graph reading

```text
a/b  ->  K_{a,b},      a*b = |E(K_{a,b})|
```

has a second arithmetic control family: the even perfect numbers occur as
Farey-product rows on the `n=2` unit-excess chain

```text
a/(2a-1).
```

More precisely, the Euclid-Euler perfect number

```text
N = 2^(p-1)(2^p-1)
```

is the product of the reduced Farey term

```text
2^(p-1)/(2^p-1),
```

when `2^p-1` is prime.  These fractions are Farey neighbors of `1/2` because

```text
(2^p-1) - 2*2^(p-1) = -1.
```

This makes the perfect-number rows an exact `n=2` unit-excess product chain,
parallel to the LRC14 unit-excess chain

```text
a/(14a-1).
```

For `a=2^k` and prime `q=n*a-1`, the two-prime-power product satisfies

```text
sigma(aq)/(aq) = n(2a-1)/(na-1)
                = 2 + (2-n)/(na-1).
```

Thus `n=2` is exactly perfect, while the LRC14 prime-`q` shadow is deficient by

```text
2 - sigma(a(14a-1))/(a(14a-1)) = 12/(14a-1).
```

This is the disciplined connection to HYP-2220/HYP-2221: perfect numbers are
not a scalar proof of LRC14, but they are the exact fixed-point control case
for the same unit-excess product shape.

## Computation

The script
`04-computation/farey_perfect_product_obstruction_codex_s143.py` stores output
at
`05-knowledge/results/farey_perfect_product_obstruction_codex_s143.out`.

It verifies the first ordinary Farey seam:

```text
F_3: 2/3 -> product 6 -> K_{2,3}, planar, perfect product.
F_4: 3/4 -> product 12 -> K_{3,4}, first reduced K_{3,3} wall.
```

The first Euclid-Euler product rows are:

```text
2/3       product 6          first_F 3     K_{2,3} planar
4/7       product 28         first_F 7     K33-wall
16/31     product 496        first_F 31    K33-wall
64/127    product 8128       first_F 127   K33-wall
4096/8191 product 33550336   first_F 8191  K33-wall
```

Every displayed row has abundancy `sigma(N)/N=2`.  The first product `6`
appears in `F_3` as `2/3`; the later perfect products also have star-factor
duplicates such as `1/28`, but the Euclid-Euler reduced factor pair is the
earliest high-incidence product address.

## Kuratowski/Wagner Guardrail

The graph product ledger must keep minor labels attached.

```text
K5 has 10 edges, but 2/5 -> K_{2,5} also has 10 edges and is planar.
K3,3 has 9 edges, but 3/3 is not a reduced Farey term; 1/9 is a star.
```

The first reduced complete-bipartite Farey term containing `K_{3,3}` is

```text
3/4 -> K_{3,4}.
```

The disjoint union `K5 disjoint-union K3,3` has `(edges,vertices)=(19,11)`,
the unreduced density mixture of `(10,5)` and `(9,6)`, but it is nonplanar
only because it already contains the two irreducible obstruction cores.
Mediants and weighted edge-density averages do not create new forbidden
minors.  The true graph iteration is minor/subdivision closure, which is
transitive.

## Tournament Analysis

S143 explicitly challenges the assumption that tournament vertices should be
runners or arcs.  Alternate vertex sets considered:

```text
Farey fractions, products, graph isomorphism types, minor cores,
factor pairs, divisor atoms, density ratios, wall crossings,
and proof obligations.
```

The chosen vertices are carrier roles:

```text
minor_closure, kuratowski_core, K33_incidence, perfect_abundancy,
unit_excess_chain, farey_level, product_edges, raw_edge_density.
```

The pairwise observable is: for each criterion, which role better preserves
the target predicate.  The gauge is majority vote over:

```text
planarity_obstruction,
perfect_fixed_point,
lrc_farey_transfer,
anti_scalar_guard,
density_bookkeeping.
```

The tie Hamiltonian path is:

```text
minor_closure > kuratowski_core > K33_incidence > perfect_abundancy
> unit_excess_chain > farey_level > product_edges > raw_edge_density.
```

The majority tournament is not transitive.  Its fingerprint is:

```text
score histogram = {0:1, 1:1, 2:1, 3:1, 5:2, 6:2}
directed 3-cycles = 2
nontrivial SCC = {K33_incidence, farey_level, product_edges, unit_excess_chain}
Hamiltonian paths = 5
```

The directed cycles are:

```text
K33_incidence -> product_edges -> farey_level -> K33_incidence
unit_excess_chain -> product_edges -> farey_level -> unit_excess_chain
```

This is the actionable warning.  In this thread, Farey level, product edges,
unit-excess address, and K33 incidence form a real middle packet.  None of
them should be scalarized alone before the exact `M`/Farey branch and graph
minor labels are retained.

## Proof-Route Use

For the LRC14 program, HYP-2945 says:

1. Keep exact `M=p/q` and the Farey excess `e=14p-q` first.
2. Use the product `pq` as a labelled complete-bipartite incidence side
   channel, not as the theorem denominator.
3. Treat even perfect products as the exact `n=2` fixed-point control family
   for the same unit-excess shape.
4. Treat the LRC14 chain `p/(14p-1)` as a deficient shadow of that perfect
   control, with defect `12/(14p-1)` in the prime-`q` power-of-two model.
5. Use `K_{p,q}` only after attaching the minor predicate: `p=2` remains a
   planar two-block branch, while `p>=3` crosses the `K_{3,3}` incidence wall.

The live proof target remains the HYP-2932/HYP-2940 route:

```text
exact Farey branch
-> C27 shell for p=2
-> K33 incidence for p>=3
-> possible HYP-2908 / THM-572 forbidden-H state lift.
```

HYP-2945 adds an arithmetic control family and a scalarization guardrail; it
does not prove the terminal state lift.

## Next Tests

1. Add a row-bank annotation for `p/(14p-1)` unit-excess rows that records
   whether `p` is a power of two and whether `14p-1` is prime.
2. Compare divisor-lattice abundancy shadows with the `K_{p,q}` rank on the
   S130/S136 low-frontier rows.
3. Search whether any LRC14 packet construction can use the `n=2` perfect
   chain as a calibration family for a `p=2` versus `p>=3` split.
4. Build a graph-minor-preserving version of the product ledger that stores
   edge count, bipartite rank, and obstruction core separately.

## See Also

`04-computation/farey_perfect_product_obstruction_codex_s143.py`;
`05-knowledge/results/farey_perfect_product_obstruction_codex_s143.out`;
`07-reflections/farey-perfect-product-obstruction-codex-s138.md`;
HYP-2940; HYP-2932; HYP-2931; HYP-2930; HYP-2221; HYP-2220; HYP-2908; THM-572.
