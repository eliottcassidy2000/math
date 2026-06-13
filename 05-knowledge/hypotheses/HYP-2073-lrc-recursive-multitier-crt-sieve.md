---
id: HYP-2073
status: PROGRESS - exact residual-packet invariant found; proof norm still open
source: codex-2026-06-02-S562
related:
  - HYP-1844
  - HYP-1866
  - HYP-1868
  - HYP-1942
  - HYP-1952
  - HYP-2036
  - HYP-2064
  - HYP-2069
  - HYP-2072
  - THM-369
  - THM-391
---

# HYP-2073: Multi-tier CRT sieves recurse on residual packets with conserved frontier product

## Statement

A multi-tier CRT sieve should not be modeled as one flat list of divisibility
tests.  It is a recursion on product-tree residual packets:

```text
node = (base denominator n, skipped quotient label, valuation vector of scale).
```

At each node, a local tier either gives a THM-369/S561 witness or exports a
residual packet to child valuation nodes.  The exported packet can shrink the
visible Archimedean gap, but exact hard rows suggest it does so by creating
matching endpoint/frontier debt.

The first invariant is the raw product

```text
Gap_A(packet) * boundary_count(packet),
```

where `Gap_A = max_gap / (1/n)` and `boundary_count` is the number of exact
endpoint witnesses in the S356 interval audit.  On the audited dyadic residual
lifts, multiplying the scale by `2` divides `Gap_A` by `2`, doubles the boundary
count, and preserves the product.

## Evidence

`lrc_multitier_crt_sieve_s562.py` audits the gate-packet family

```text
{1} union {scale*q : 1 <= q < n, q != skip}
```

for the known `n=14`, `n=17`, and `n=18` packets.

```text
n=14 skip=6:
  scale 7  -> 14 -> 28
  gap/th 5/924 -> 5/1848 -> 5/3696
  boundary 84 -> 168 -> 336
  product 5/11 throughout

n=17 skip=8:
  scale 17 -> 34 -> 68
  gap/th 1/272 -> 1/544 -> 1/1088
  boundary 450 -> 900 -> 1800
  product 225/136 throughout

n=18 skip=8:
  scale 9 -> 18 -> 36
  gap/th 1/176 -> 1/352 -> 1/704
  boundary 176 -> 352 -> 704
  product 1 throughout
```

The `n=17` dyadic control is the surprise: even without a squarefree CRT base
denominator, the skip-8 residual packet translates in a new dyadic endpoint
tier with the same product conservation.  Thus the phenomenon is not only
"two primes interact"; it is "a generated residual packet has a local norm that
survives valuation translation."

All audited packets are sieve-complete in the THM-369 zero-branch sense:
`missing_zero_branches = -`.  They are not open-cover candidates; they are
positive-gap residual packets whose visible gap has been traded for boundary
debt.

## Recursive Sieve Reading

S561/HYP-2072 gives the first two-tier algorithmic split for `k+1=2q`:

```text
mod-q oracle -> generated residual -> apex/ratio-cover subtiers.
```

HYP-2073 adds the next recursion layer:

```text
generated residual -> product-tree valuation node
                   -> child node after local lift
                   -> conserved/weighted frontier product.
```

For `n=14`, the odd payload is `7` and the lift moves in the `2`-tree:

```text
v(scale) = (2:0,7:1), (2:1,7:1), (2:2,7:1).
```

For `n=18`, the odd payload is `3^2` and the lift again moves in the `2`-tree:

```text
v(scale) = (2:0,3:2), (2:1,3:2), (2:2,3:2).
```

For `n=17`, the prime payload is `17`, but the residual endpoint packet still
admits dyadic translation:

```text
v(scale) = (2:0,17:1), (2:1,17:1), (2:2,17:1).
```

This suggests an abstract proof strategy: a counterexample cannot win by
repeatedly lifting residual packets if each lift preserves a positive product
mass in at least one local tree or product-building direction.

## Tournament Analysis

Vertices were recursive CRT packet states, not runners:

```text
n14:scale7:skip6, n14:scale14:skip6, ...
```

Pairwise observable:

```text
(gap/th, gap*boundary, scale valuation depth, boundary)
```

Switch/gauge: smaller visible gap wins; ties prefer conserved lower debt and
deeper valuation translation.  The nine-state audit is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3_cycles=0
sccs=[1,1,1,1,1,1,1,1,1]
hamiltonian_paths=1
```

This transitivity is itself diagnostic.  Bare valuation translation is a
ledger, not yet a cyclic obstruction.  Cycles should reappear only after
endpoint owners, pressure leaves, wall-crossing events, or cross-prime coupling
are retained.

## Assumption Challenge

The session considered runners, quotient labels, skipped gate labels, product
p-adic nodes, endpoint owners, boundary witnesses, fixed circle sections,
wall-crossing events, residues, cover arcs, Fourier/Gabor modes, matroid
circuits, and proof obligations.

Chosen quotient:

```text
vertices = whole residual gate packets / proof-obligation states.
```

Predicate preserved:

```text
THM-369 sieve completeness, skipped quotient label, valuation translation,
exact gap ratio, exposed boundary count, and raw frontier product.
```

Information destroyed:

```text
individual runner identity, endpoint ownership, exact interval adjacency,
pressure-leaf structure, and cross-prime cancellation signs.
```

Challenged assumption:

```text
CRT sieve tiers are just prime factors of n, or tournament vertices must be
runners.
```

The `n=17` dyadic control refutes the narrow version: a recursive tier can be
introduced by the residual endpoint denominator, not only by the base `n`.

## Open Work

1. Replace raw `boundary_count` by the HYP-1868 product-building frontier mass,
   retaining owner labels and cross-prime signs.
2. Extend S561's generated residual set so every residual packet carries a
   canonical valuation node and a weighted frontier ledger.
3. Test non-dyadic local lifts `*p`; the current exact rows only test the
   dyadic translation edge.
4. Prove a descent/no-escape lemma: every infinite residual recursion path has
   positive product-building mass, so it cannot converge to both zero visible
   gap and zero endpoint debt.

## Files

- `04-computation/lrc_multitier_crt_sieve_s562.py`
- `05-knowledge/results/lrc_multitier_crt_sieve_s562.out`
- `07-reflections/lrc-recursive-multitier-crt-sieve-s562.md`
