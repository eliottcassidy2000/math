---
id: HYP-1793
status: OPEN
source: endpoint_collision_geometry_s95.py
session: codex-2026-05-30
---

# HYP-1793: SC Collision Hypergraph Is Leaf-Peelable

## Statement

The non-private self-complementary collision columns in merged endpoint
transfer form a 3-uniform incidence hypergraph that is leaf-peelable.

Equivalently, after restricting to the collision block, every nonempty
subcollection of collision columns has some parent row appearing in exactly one
remaining column.

If true, HYP-1791 follows immediately: ordering the leaf-peeling removals
gives a triangular incidence submatrix over `F_2`, so the collision columns
are linearly independent.

## Evidence

`endpoint_collision_geometry_s95.py` computes the support-3 hypergraph on
non-private SC columns and greedily removes columns with a degree-one owner.

For transitions `2->3` through `6->7`:

```text
nonprivate_sc: [0, 0, 0, 1, 14]
peel_removed:  [0, 0, 0, 1, 14]
peel_core:     [[], [], [], [], []]
```

At `6->7`, the owner-degree spectrum is:

```text
{1: 5, 2: 7, 3: 4, 5: 1, 6: 1}
```

The hypergraph is not merely linear: pair intersections among the 14 triples
include six pairs sharing two owners:

```text
pair_intersections = {0: 53, 1: 32, 2: 6}
```

So the observed independence is not explained by disjointness or ordinary
linearity.  Leaf-peelability is the stronger visible mechanism.

## Interpretation

HYP-1792 says the collision relation is not parent-metagraph clique geometry.
HYP-1793 says it may still be triangular, but in the correct incidence
category.  The "private pivot" mechanism fails column-by-column, yet a weaker
private owner appears after previous collision columns are peeled away.

This refines the merged endpoint-rank problem:

```text
decomposable SC boundary   -> private pivots
strong SC collision block  -> leaf-peelable 3-uniform incidence hypergraph
```

The next proof target is to construct the leaf order without exhaustive
enumeration, perhaps by sorting collision columns by endpoint-exposure height,
owner `H`, or SCC-defect tuple.
