---
id: HYP-2998
title: Farey-Fibonacci additive-basis carrier and representation-economy split
status: SYNTHESIS / finite scout and proof-interface carrier; not a proof
source: codex-2026-06-24-S169
artifacts:
  - 04-computation/farey_fibonacci_additive_basis_s169.py
  - 05-knowledge/results/farey_fibonacci_additive_basis_s169.out
  - 07-reflections/farey-fibonacci-additive-basis-carrier-s169.md
related:
  - HYP-2934
  - HYP-2933
  - HYP-2932
  - HYP-2931
  - HYP-2930
  - HYP-2219
  - HYP-2218
  - HYP-2091
  - HYP-1920
  - HYP-1902
  - HYP-2997
  - HYP-2995
  - HYP-2990
  - HYP-2982
  - OPEN-Q-108
---

# HYP-2998: Farey-Fibonacci Additive-Basis Carrier

The old Goldbach, Fermat polygonal, Zeckendorf, and Farey-mutation threads
become one useful proof interface when they are read as representation
hypergraphs with different proof economies.

```text
Goldbach / ternary Goldbach: many prime branches and local residue factors.
Fermat polygonal numbers: bounded arity pays a finite local-obstruction bill.
Zeckendorf / Fibonacci: sparse recurrence plus a no-adjacent normal form.
Farey mutations: one fraction address split into q, p+q, p*q, and power stress.
```

The practical claim is a controlled-kernel rule for sequence shadows: before a
quotient forgets representation data, it must declare which economy it is using.
Smoothing, bounded arity, normal form, and fraction-address ledgers preserve
different coordinates.

## Computation

The scout

```text
04-computation/farey_fibonacci_additive_basis_s169.py
```

stores output at

```text
05-knowledge/results/farey_fibonacci_additive_basis_s169.out
```

It records the user's Fibonacci arrangement as the sparse-carry diagonal

```text
F_{n+1} = sum_k binom(n-k,k)
```

with rows

```text
1
1
1+1
1+2
1+3+1
1+4+3
1+5+6+1
...
```

Equivalently, it counts legal supports with no adjacent carries.  This is the
finite combinatorial face of Zeckendorf: the recurrence does not merely cover
integers, it supplies a canonical sparse support.

The same script isolates the golden Farey/Stern-Brocot spine

```text
p/q = F_i/F_{i+1}.
```

On this all-ones continued-fraction spine, the additive Farey payload is exactly
the next Fibonacci number:

```text
p+q = F_{i+2}.
```

The product payload is the complete-bipartite area

```text
p*q = |E(K_{p,q})|,
```

matching HYP-2932/HYP-2934.  The power payloads are ordered magnitude-stress
channels, not proof denominators.  In the finite table, `q^p` wins for `1/2`
and `2/3`, while `p^q` wins starting at `3/5`.  That flip is useful: powers
detect the ordered address that sum and product can partially hide.

The finite additive-basis audit also recovers the old S501 readout through
`N=500`:

```text
binary Goldbach missing evens: none in the range
ternary Goldbach missing odds: none in the range
3..8-gonal max summands: exactly 3..8 in the range
Zeckendorf atoms <= 500: 13 atoms, max digit count 6
```

## Tournament Analysis

The vertices are proof carriers, not runners:

```text
zeckendorf_sparse_normal_form
farey_address_vector
fermat_polygonal_bounded_arity
farey_product_Kpq_area
ternary_goldbach_smoothing
binary_goldbach_pair_graph
farey_power_stress_channel
raw_scalar_rep_count
```

The pairwise observable is the tuple

```text
(retained coordinates, proof power, scalar safety, LRC transfer)
```

with the listed carrier order as the tie Hamiltonian path.  The resulting
tournament is transitive:

```text
score_hist=[(0,1),(1,1),(2,1),(3,1),(4,1),(5,1),(6,1),(7,1)]
directed_3cycles=0
Hamiltonian path:
  zeckendorf_sparse_normal_form
  > farey_address_vector
  > fermat_polygonal_bounded_arity
  > farey_product_Kpq_area
  > ternary_goldbach_smoothing
  > binary_goldbach_pair_graph
  > farey_power_stress_channel
  > raw_scalar_rep_count
```

This ranking is not a claim that normal forms are always "better" than
smoothing.  It says that for LRC-style quotient guardrails, a canonical support
or explicit fraction-address vector is safer than a raw representation count.

## Relation To Earlier Threads

S501 already put Goldbach, Fermat polygonal, and Zeckendorf under one
representation-hypergraph grammar.  HYP-2998 adds the missing Farey operation
axis:

```text
q       = theorem-level binding scale,
p+q     = additive Stern-Brocot recursion ledger,
p*q     = K_{p,q} incidence / product area,
p^q,q^p = ordered magnitude-stress tests.
```

This refines HYP-2931/HYP-2932/HYP-2934 rather than replacing them.  The LRC14
unit-excess chain still keeps `q` as the binding denominator, uses `p+q` as the
linear recursion ledger, and uses `p*q` only after the `K_{p,q}` incidence label
is retained.  The new Fibonacci observation says that on the all-ones
Stern-Brocot spine, the additive ledger is not arbitrary: it is exactly the next
Fibonacci carry level.

## Proof Target

For future LRC/OCF sequence-shadow uses, attach a representation-economy label
to every quotient:

```text
SMOOTHING: keep local residue/singular-series side data.
BOUNDED_ARITY: keep the finite summand budget and residue invoice.
NORMAL_FORM: keep the carry automaton and forbidden adjacency/collision rule.
FAREY_ADDRESS: keep q, p+q, p*q, and power-stress as separate clocks.
```

The candidate theorem schema is:

```text
A representation quotient is admissible only when the target predicate is
constant on the selected economy's fibers, or the forgotten coordinate is
reconstructed, annihilated by a dual certificate, exact as a cocycle, descended
to a smaller family, or emitted as named residual debt.
```

For LRC14, this should become a new field in packet classifiers beside the
HYP-2995 carrier atlas and the HYP-2997 cocycle normal-form atlas: every
sequence shadow must state whether it is using Goldbach smoothing, Fermat
bounded arity, Zeckendorf normal form, or Farey address retention.  Otherwise
it is just another raw scalar representation count, and should be treated as
unsafe.
