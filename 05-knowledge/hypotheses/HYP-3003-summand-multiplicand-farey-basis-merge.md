---
id: HYP-3003
title: Summand/multiplicand Farey additive-basis merge
status: SYNTHESIS / proof-interface carrier; not a proof
source: codex-2026-06-24-S171
script: 04-computation/summand_multiplicand_farey_basis_merge_codex_s171.py
result: 05-knowledge/results/summand_multiplicand_farey_basis_merge_codex_s171.out
related:
  - HYP-3000
  - HYP-2999
  - HYP-2998
  - HYP-3001
  - HYP-2934
  - HYP-2935
  - HYP-2083
  - HYP-2067
  - HYP-1902
  - HYP-2218
  - HYP-2219
  - LTI-153
---

# HYP-3003: Summand/Multiplicand Farey Additive-Basis Merge

## Claim

The Goldbach / ternary Goldbach / Fermat polygonal / Zeckendorf / Fibonacci /
Farey-mutation threads should be read through two operation graphs over the
same labelled packet:

```text
summand graph:       x+y=N          dense antidiagonal fibers
multiplicand graph:  x*y=N          sparse hyperbola / factor fibers
Farey packet:        p/q            root packet whose shadows are p+q and p*q
```

HYP-2998, HYP-2999, and HYP-3000 already identify the representation-economy
axis.  HYP-3003 adds the missing two-graph normal form: `p+q` is not just a
number, it is a summand-graph antidiagonal; `p*q` is not just a scalar product,
it is a factor/Kpq incidence fiber.  The power lanes `q^p,p^q` remain stress
tests unless the root packet and forgotten fibers are retained.

## Fibonacci Row

The user-specified row

```text
1, 1, 1+1, 1+2, 1+3+1, 1+4+3, 1+5+6+1, ...
```

is the `d=2` Pascal/path row

```text
row(n,k) = C(n-1-k,k),      sum_k row(n,k)=F_n.
```

The scalar Fibonacci number is therefore a shadow.  The proof packet is the
whole row vector: rank counts of independent sets in a path, with Zeckendorf
normal form obtained by adding the no-adjacent carry rule.

## Computation

Script:

```text
04-computation/summand_multiplicand_farey_basis_merge_codex_s171.py
```

Stored output:

```text
05-knowledge/results/summand_multiplicand_farey_basis_merge_codex_s171.out
```

The scout checks:

- Fibonacci/Pascal rows through `n=10`;
- binary Goldbach pair counts and ternary Goldbach triple counts on sample
  targets;
- Fermat polygonal bounded arity through `300` for sides `3..8`;
- Zeckendorf normal forms for selected Fibonacci/tournament-side numbers;
- summand pair counts versus multiplicand factor-pair counts on shared nodes;
- golden Farey packets `F_i/F_{i+1}` where `p+q=F_{i+2}`;
- LRC14 unit-excess packets `p/(14p-1)`.

Representative LRC14 unit-excess readout:

```text
p=1, q=13:  p+q=14,  p*q=13   -> q-parent/AP-GW side
p=2, q=27:  p+q=29,  p*q=54   -> C27 petal/two-block side
p=3, q=41:  p+q=44,  p*q=123  -> first Kpq/K33 product wall
```

This refines the older HYP-2934/HYP-2935 lesson: `2/27` belongs to the
summand/C27 second-gap branch, while `3/41` is the first product-incidence wall.

## Tournament Analysis

Vertices are proof currencies / representation fibers, not runners:

```text
exact_farey_packet
summand_antidiagonal_fiber
multiplicand_hyperbola_fiber
zeckendorf_path_carry
fibonacci_pascal_row
fermat_polygonal_invoice
ternary_goldbach_smoothing
binary_goldbach_pair_sieve
farey_power_stress
```

Pairwise observable:

```text
(exact_packet, additive_fiber, product_fiber, carry_normal_form,
 bounded_arity, smoothing_entropy, local_residue, lrc_transfer,
 anti_scalar_guard)
```

Gauge: coordinate majority; tie Hamiltonian path is the listed order.  The
scout gives a transitive fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
```

## Guardrail

A representation/Farey quotient is theorem-safe only after declaring:

```text
root Farey packet: p/q, q, excess, route labels
summand shadow:    p+q as antidiagonal fiber / additive-pinch row
multiplicand:      p*q as factor/Kpq incidence fiber
carry law:         Pascal row and Zeckendorf no-adjacent normal form when used
proof economy:     smoothing, bounded arity, normal form, product incidence, or stress test
```

If one of these coordinates is forgotten, it must be reconstructed, killed by a
dual certificate, exact/coboundary, descended to a smaller family, or routed to
a named residual sector.  This is the summand/multiplicand version of the
current no-free-slider rule.
