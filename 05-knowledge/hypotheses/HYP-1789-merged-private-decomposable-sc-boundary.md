---
id: HYP-1789
status: OPEN
source: endpoint_private_goodcut_s95.py; THM-354
session: codex-2026-05-30
---

# HYP-1789: Merged Private Pivots Cover the Decomposable SC Boundary

## Statement

For complement-merged endpoint transfer `n -> n+1`, every self-complementary
child node with more than one strongly connected component is a private odd
child column.

Equivalently, using THM-354:

```text
if child node U is SC and good_cut(U) < (n+1)-1,
then U is a private odd column of the merged endpoint-transfer matrix.
```

Thus any self-complementary child node that is not private must be strongly
connected.

## Evidence

Exact transitions `2->3` through `6->7`:

```text
child n=3:
  all SC by good-cut:         {0: 1, 2: 1}
  private SC by good-cut:     {0: 1, 2: 1}
  non-private SC by good-cut: {}

child n=4:
  all SC by good-cut:         {0: 1, 3: 1}
  private SC by good-cut:     {0: 1, 3: 1}
  non-private SC by good-cut: {}

child n=5:
  all SC by good-cut:         {0: 1, 2: 1, 4: 6}
  private SC by good-cut:     {0: 1, 2: 1, 4: 6}
  non-private SC by good-cut: {}

child n=6:
  all SC by good-cut:         {0: 1, 3: 1, 4: 1, 5: 9}
  private SC by good-cut:     {0: 1, 3: 1, 4: 1, 5: 8}
  non-private SC by good-cut: {5: 1}

child n=7:
  all SC by good-cut:         {0: 1, 2: 1, 4: 7, 6: 79}
  private SC by good-cut:     {0: 1, 2: 1, 4: 7, 6: 65}
  non-private SC by good-cut: {6: 14}
```

The only missing private pivots are in the top good-cut bucket, i.e. the
strongly connected child layer.

The same probe also found that every non-private SC column has support `3` in
the tested range:

```text
5->6: {3: 1}
6->7: {3: 14}
```

This support-3 refinement is split out as HYP-1790.

## Interpretation

The merged endpoint boundary is exactly the SC child set, but not every SC
child must be private. HYP-1789 says the non-private part of the boundary is
concentrated entirely inside the strongly connected layer.

So the merged rank problem splits:

1. decomposable SC child nodes give triangular private pivots;
2. strongly connected SC child nodes form the only collision zone;
3. the odd-Smith/full-rank phenomenon must resolve those top-layer collisions.

This sharpens HYP-1788 from "private pivots are good-cut pure" to "private
pivot failure only begins after the SCC-defect has vanished."

## Possible Proof Route

A decomposable SC child has a nontrivial ordered SCC decomposition. Endpoint
deletion may interact asymmetrically with the first and last component blocks,
leaving exactly one merged parent with odd deletion parity. Strongly connected
SC children have no component boundary to anchor this deletion, so multiple
parents can collide in the same odd column.

## Questions

1. Can the missing strongly connected SC nodes be characterized by endpoint
   deletion profiles of support size `>1`?
2. Are the missing counts `0,0,0,1,14,...` controlled by a known SC strongly
   connected tournament sequence?
3. Does the merged odd-Smith full-rank phenomenon reduce to proving full rank
   on the strongly connected SC collision block?
