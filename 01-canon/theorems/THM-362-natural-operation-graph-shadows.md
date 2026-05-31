---
id: THM-362
name: natural-operation-graph-shadows
status: PROVED
date: 2026-05-31
session: codex-2026-05-31-S365
depends_on:
  - THM-361
---

# THM-362: Natural Operation Graph Shadows

## Statement

Let `N+ = {1,2,3,...}`.

Define the additive operation graph `A` by putting a loopless simple directed
edge `x -> z` when there exists `y in N+` with `x+y=z`.  Then

```text
x -> z in A  iff  x < z.
```

Consequently, the induced graph on `{1,...,N}` is the transitive tournament,
and its transitive reduction is the successor chain.

Define the multiplicative operation graph `M` by putting a loopless simple
directed edge `x -> z` when there exists `y in N+` with `xy=z`.  Then

```text
x -> z in M  iff  x divides z and x < z.
```

Thus `M` is the proper-divisor DAG.  On `{1,...,N}`, its transitive reduction
has exactly the edges

```text
x -> xp
```

with `p` prime and `xp <= N`.

A positive `k`-tuple is a collision between the `k`-ary additive and
multiplicative operation hyperedges exactly when

```text
x_1 + ... + x_k = x_1 * ... * x_k.
```

If all entries equal to `1` are deleted and the remaining nonunit core is
`c=(c_1,...,c_m)`, with product `P(c)` and sum `S(c)`, then the deleted number
of ones is forced to be `P(c)-S(c)`.  Hence product-sum tuples are exactly
nonunit cores with nonnegative defect, expanded by that many copies of `1`.

For every product-sum tuple with product `P`, division by `P` gives the
Egyptian-fraction identity

```text
sum_i 1/(P/x_i) = 1.
```

For fixed arity `k`, there are only finitely many positive sorted solutions.

## Proof

For addition, if `x+y=z` with `y>=1`, then `x<z`.  Conversely, if `x<z`, then
`y=z-x` is a positive integer and witnesses `x->z`.  Thus the finite induced
simple graph on `{1,...,N}` has exactly one directed edge from the smaller
number to the larger number for every unordered pair.  This is the transitive
tournament, and its cover graph is the successor chain.

For multiplication, if `xy=z`, then `x` divides `z`; after loops are removed,
`x<z`.  Conversely, if `x` is a proper divisor of `z`, then `y=z/x` is a
positive integer larger than `1`, so `x->z`.

The transitive reduction of divisibility is the prime-multiplication relation.
If `z/x` is prime, then there is no divisor strictly between `x` and `z`.  If
`z/x` is composite, write `z/x=ab` with `a,b>1`; then

```text
x -> xa -> xab=z
```

is a nontrivial directed path, so the edge `x->z` is not a cover.

The collision statement is the definition of the two operation hyperedges
landing at the same target: the additive target is `sum_i x_i`, and the
multiplicative target is `prod_i x_i`.

For the core-defect statement, suppose a product-sum tuple has `q` entries
equal to `1` and nonunit core `c`.  The product is unchanged by those ones,
while the sum increases by `q`; hence

```text
P(c) = S(c) + q,
```

so `q=P(c)-S(c)`.  Conversely, if `P(c)-S(c)=q>=0`, then adding exactly `q`
ones to `c` gives `sum = S(c)+q = P(c) = product`.

Dividing `sum_i x_i=P` by `P` gives

```text
sum_i x_i/P = sum_i 1/(P/x_i) = 1.
```

For finiteness at fixed arity, sort a collision as

```text
1 <= a_1 <= ... <= a_k.
```

Let `P_{k-1}=a_1...a_{k-1}`.  Since each `a_i <= a_k`,

```text
P_{k-1} a_k = a_1+...+a_k <= k a_k.
```

Thus `P_{k-1} <= k`, so the first `k-1` sorted entries have bounded product.
There are only finitely many such prefixes, and once a prefix has product `P`
and sum `S`, the last entry is forced by `x(P-1)=S`.

## Verification Record

`04-computation/natural_operation_graphs_s365.py` verifies the finite graph
statistics and enumerates sorted product-sum collisions through arity `20`.
The stored output is
`05-knowledge/results/natural_operation_graphs_s365.out`.

## Notes

The theorem separates two objects that are easy to conflate: the unlabeled
simple graph shadow of an operation, and the labeled operation cospan/hypergraph
with its fibers.

The additive simple graph forgets the complement `y=z-x` and collapses to
ordinary order.  Multiplication survives this forgetting because divisibility
is already a sparse partial order.  The equations `x+y=xy`, `x+y+z=xyz`, and
their higher-arity analogues therefore belong most naturally to the two-colored
operation hypergraph, not merely to the simple edge shadow.

## Related

- THM-361: product-sum defect normal form.
- HYP-1820: operation-shadow critical pairs.
- HYP-1821: natural-mode defect interface.
- HYP-1822: natural-number operation graphs are labeled cospans.
- `07-reflections/natural-operation-graphs-s365.md`.
