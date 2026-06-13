---
id: THM-361
name: product-sum-defect-normal-form
status: PROVED
date: 2026-05-31
session: codex-2026-05-31-S365
lean: TournamentH7.ProductSum
depends_on: []
---

# THM-361: Product-Sum Defect Normal Form

## Statement

Let `M` be a finite multiset of positive integers and let

```text
P(M) = product of entries of M,
S(M) = sum of entries of M.
```

Remove all `1`s from `M` and call the remaining multiset `C`.  If `m` is the
number of removed `1`s, then

```text
P(M) = S(M)
```

if and only if

```text
P(C) - S(C) = m.
```

Conversely, every finite multiset `C` of integers at least `2` with
`d = P(C)-S(C) >= 0` gives a product-sum multiset by adjoining exactly `d`
copies of `1`.

In particular:

1. The only two-entry positive product-sum multiset is `{2,2}`.
2. The only product-sum multiset with all entries distinct and at least two
   entries is `{1,2,3}`.

## Proof

Removing `1`s does not change the product and subtracts exactly `m` from the
sum:

```text
P(M) = P(C),
S(M) = S(C) + m.
```

Therefore `P(M)=S(M)` is equivalent to `P(C)=S(C)+m`, i.e.
`P(C)-S(C)=m`.  The converse is the same identity read backwards.

For two entries, `a+b=ab` gives

```text
(a-1)(b-1)=1,
```

so `a=b=2`.

Now assume all entries are distinct.  If no entry is `1`, then for two entries
the preceding equation forces the repeated pair `{2,2}`, and for three or more
entries the smallest possible distinct core is `{2,3,4}`, where

```text
2*3*4 - (2+3+4) = 15 > 0.
```

Increasing any core entry by `e>0` increases the product by `e` times the
product of the other entries, while increasing the sum only by `e`; with at
least three entries at least `2,3,4`, the defect can only grow.  Hence a
distinct no-`1` product-sum tuple is impossible.

If a distinct product-sum multiset contains `1`, remove it.  The remaining
distinct core `C` has defect `1`.  A one-entry core has defect `0`.  A two-entry
core `{a,b}` has

```text
P(C)-S(C)=ab-a-b=(a-1)(b-1)-1,
```

so defect `1` is equivalent to `(a-1)(b-1)=2`, giving `{a,b}={2,3}`.  A core
with at least three distinct entries has defect at least `15`, by the previous
paragraph.  Thus the only distinct positive solution is `{1,2,3}`.

## Interpretation

Product-sum equations are not isolated coincidences.  They are the zero-defect
interface between additive and multiplicative modes:

```text
defect(C) = product(C) - sum(C).
```

The entry `1` is exactly additive slack: each copy of `1` repairs one unit of
multiplicative excess without changing the product.

This explains why the first low-dimensional resonances are:

```text
2+2 = 2*2 = 4,
1+2+3 = 1*2*3 = 6.
```

The first is a diagonal binary resonance; the second is the unique distinct
positive resonance.

## Lean formalization

The list-level core of this theorem is formalized in
`04-computation/lean/TournamentH7/TournamentH7/ProductSum.lean` and audited
from `TournamentH7.Verify`.

Formalized statements include:

```lean
Tournament.ProductSum.product_sum_iff_core
Tournament.ProductSum.pad_core_product_sum
Tournament.ProductSum.core_defect_eq_ones_of_product_sum
Tournament.ProductSum.two_entry_product_sum
Tournament.ProductSum.no_three_ge_two_product_sum
Tournament.ProductSum.one_cons_two_entry_product_sum
Tournament.ProductSum.three_entry_distinct_product_sum
```

This covers the core defect normal form, the padding converse, the natural
subtraction corollary, the positive ordered two-entry uniqueness `(2,2)`, the
three-entry obstruction when all entries are at least `2`, the explicit
`1,a,b` classification, and the ordered distinct positive ternary uniqueness
up to `List.Perm [1,2,3]`. The multiset-level quotient and arbitrary-length
distinct positive uniqueness theorem remain future formalization targets.

## Related

- `04-computation/lean/TournamentH7/TournamentH7/ProductSum.lean`.
- `05-knowledge/results/lean_product_sum_s366.out`.
- `04-computation/natural_mode_graph_s365.py`.
- `05-knowledge/results/natural_mode_graph_s365.out`.
- `07-reflections/summand-graph-fermat-zeckendorf.md`.
- `07-reflections/product-graph-sc-spine-fractal-dimensions.md`.
