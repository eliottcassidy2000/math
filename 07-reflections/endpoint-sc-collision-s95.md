# Endpoint SC Collision S95

The merged endpoint boundary is the SC child set. Private pivots account for
most of that boundary, but not all of it.

The collision probe splits SC children into two blocks:

```text
decomposable SC children:  private through 6->7
strongly connected SC:     private plus collision columns
```

The non-private SC sequence is:

```text
[0, 0, 0, 1, 14]
```

All of it lives in the top good-cut bucket, i.e. the strongly connected layer.
Even more sharply, every non-private SC column found has support exactly `3`:

```text
5->6: {3: 1}
6->7: {3: 14}
```

So the first failure of triangularity is ternary. It is not a broad many-parent
collapse. The collision columns are also independent in the tested range:

```text
nonprivate_sc:   [0, 0, 0, 1, 14]
nonprivate_rank: [0, 0, 0, 1, 14]
```

## Meta Connection

This looks like the endpoint-transfer version of the recurring cubic
obstruction:

```text
directed 3-cycles
beta_3 generators
cubic Krawtchouk residuals
support-3 SC endpoint collisions
```

The same pattern may be appearing in different coordinates: when a binary
quotient cannot stay triangular, its first irreducible obstruction is often a
three-way relation.

## New Working Problem

Reduce merged endpoint full rank to two pieces:

1. decomposable SC private pivots give a triangular boundary block;
2. strongly connected SC support-3 collisions form a ternary hypergraph whose
   incidence must have the remaining rank.

If the support-3 hypergraph has a recursive description, HYP-1773's merged
rank conjecture might reduce to a cubic-rank statement on the strongly
connected SC spine.

## Geometry Check

The follow-up script `endpoint_collision_geometry_s95.py` checks whether these
support-3 columns are literal triangles in the parent merged arc-flip
metagraph.  They are not.  At `6->7`, the 14 owner triples have induced edge
counts

```text
{0: 1, 1: 6, 2: 5, 3: 2}
```

So the collision block should be treated as endpoint-transfer incidence
geometry.  Parent-metagraph adjacency is a useful feature, but it is not the
relation itself.

The same script found a replacement triangular mechanism: the support-3
collision hypergraph leaf-peels completely in the tested range.

```text
peel_removed: [0, 0, 0, 1, 14]
peel_core:    [[], [], [], [], []]
```

This is stronger than simply observing full rank.  It suggests the collision
block may have a recursive pivot order even when no column is initially
private in the original merged transfer matrix.
