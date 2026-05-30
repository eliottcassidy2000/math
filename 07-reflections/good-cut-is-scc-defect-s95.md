# Good-Cut Is SCC Defect S95

The good-cut coordinate looked path-geometric: choose a Hamiltonian path, mark
which path cuts are crossed by upward tiles, and count them.

The hidden invariant is simpler:

```text
good_cut_count = n - number_of_strong_components.
```

A non-good cut is exactly a boundary between strong components. Along any
directed Hamiltonian path, tournament SCCs appear as contiguous blocks in
their forced total order. Internal cuts of an SCC must be crossed by some
backward edge; otherwise the right side could not reach the left side.

This closes the older merged-coordinate question. Good-cut count descends to
ordinary tournament classes first, then to merged classes because complement
preserves SCC count while reversing the component order.

It also explains why endpoint private pivots were good-cut pure. There was no
special endpoint miracle in the purity itself: every tournament class at the
child level is already pure. The endpoint-specific mystery that remains is
selection:

```text
why do private pivots concentrate so strongly near g = n-1?
```

That question is now about endpoint transfer preferring strongly connected or
nearly strongly connected children, not about whether good-cut height is a
valid quotient coordinate.

## Useful Reframing

The good-cut bucket distribution is the distribution of SCC defects:

```text
g = 0      transitive tournaments
g = n-1    strongly connected tournaments
g = n-c    tournaments with c strong components
```

So the old good-cut triangle is a strong-component-count triangle for
fixed-path tilings. The missing bucket `g=1` is no longer mysterious: it would
require exactly `n-1` SCCs, hence one SCC of size `2`, which cannot happen in a
tournament.

For endpoint transfer, this suggests a sharper rank story:

```text
endpoint insertion tends to create private pivots in children with few SCCs;
merged parity then keeps only the self-complementary odd buckets.
```

The next refinement is HYP-1789: in the merged quotient, private-pivot failure
appears only in the strongly connected SC layer. Decomposable SC boundary
nodes still behave triangularly.
