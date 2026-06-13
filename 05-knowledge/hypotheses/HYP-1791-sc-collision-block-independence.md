---
id: HYP-1791
status: OPEN
source: endpoint_sc_collision_s95.py
session: codex-2026-05-30
---

# HYP-1791: SC Collision Block Independence

## Statement

In complement-merged endpoint transfer, the non-private self-complementary
collision columns are linearly independent over `F_2`.

Equivalently, the support-3 collision hypergraph from HYP-1790 has full column
rank when viewed as a parent-by-collision incidence matrix.

## Evidence

Exact transitions `2->3` through `6->7`:

```text
nonprivate_sc:   [0, 0, 0, 1, 14]
nonprivate_rank: [0, 0, 0, 1, 14]
```

The whole merged SC block has rank:

```text
sc_rank: [1, 2, 3, 9, 30]
```

At `6->7`, this means the 14 ternary collision columns add 14 independent
directions inside a top-SC block of rank 27.

The follow-up geometry probe `endpoint_collision_geometry_s95.py` gives a
candidate mechanism stronger than raw rank: the support-3 collision hypergraph
leaf-peels completely through `6->7`.

```text
peel_removed: [0, 0, 0, 1, 14]
peel_core:    [[], [], [], [], []]
```

If this persists, HYP-1791 can be proved by triangularizing the collision
incidence matrix in the leaf-peeling order.

## Interpretation

Private columns give immediate triangular pivots. Non-private SC columns lose
privacy but may still behave like independent ternary constraints. This is
precisely the kind of structure needed for full merged rank to survive after
private coverage drops.

If HYP-1791 holds, the merged endpoint-rank problem can be split into:

1. private decomposable SC pivots;
2. private strongly connected SC pivots;
3. independent support-3 SC collision columns;
4. remaining NSC even-boundary columns that can help row rank but do not
   change the SC boundary.

## Questions

1. Is there a canonical orientation on each support-3 collision triple that
   makes the collision block triangular after ordering by `H` or SCC defect?
2. Does the collision block remain independent integrally, i.e. with odd Smith
   invariant factors?
3. Is full merged endpoint rank equivalent to independence of this collision
   block plus private coverage of the decomposable SC boundary?
