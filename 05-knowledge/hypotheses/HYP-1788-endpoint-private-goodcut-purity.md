---
id: HYP-1788
status: PROVED
source: endpoint_private_goodcut_s95.py; goodcut_class_purity_s95.py; THM-354
session: codex-2026-05-30
---

# HYP-1788: Endpoint Private Good-Cut Purity

## Statement

Every private odd endpoint-transfer child column in the tournament fixed-path
quotients is good-cut pure.

More explicitly, if a child class or merged child node `D` is reached oddly
from exactly one parent row in endpoint transfer, then every tiling inside `D`
has the same good-cut count.

For the complement-merged quotient, the private child is also necessarily
self-complementary by the boundary corollary of THM-356.

This purity statement is now a corollary of THM-354: good-cut count equals
`n - scc(T)`, so it descends to every tournament isomorphism class, not only
to private endpoint pivots.

## Evidence

Exact endpoint-transfer/private-column computation for transitions `2->3`
through `6->7` gives:

```text
unmerged private columns:      [2, 4, 10, 29, 133]
unmerged owners covered:       [1, 2, 4, 12, 56]
unmerged good-cut width zero:  [2, 4, 10, 29, 133]

merged private columns:        [2, 2, 8, 11, 74]
merged owners covered:         [1, 2, 3, 9, 28]
merged good-cut width zero:    [2, 2, 8, 11, 74]
```

The merged private-child kind spectra are:

```text
2->3: {'SC': 2}
3->4: {'SC': 2}
4->5: {'SC': 8}
5->6: {'SC': 11}
6->7: {'SC': 74}
```

The SC membership is forced by parity. The good-cut purity was first observed
in the private-pivot scan and is now explained by THM-354.

## Good-Cut Bucket Distribution

Private columns are not confined to one height, but they skew toward high
good-cut buckets:

```text
unmerged 6->7 private columns by good-cut: {0: 1, 2: 5, 3: 4, 4: 19, 5: 47, 6: 57}
merged   6->7 private columns by good-cut: {0: 1, 2: 1, 4: 7, 6: 65}
```

The merged quotient is especially sharp: by `6->7`, almost every private pivot
lies in the maximal good-cut bucket, with a small residue on the transitive
and low-height self-complementary spine.

## Interpretation

Endpoint transfer is a deletion/insertion recursion. Good-cut count is a Morse
height on the fixed-path tiling cube, and THM-354 identifies that height as
strong-component defect. HYP-1788 says that the visible triangular pivots of
endpoint transfer land in definite SCC-defect fibers.

This would connect three previously separate facts:

1. tournament automorphism groups have odd order;
2. complement-merged odd buckets are exactly self-complementary buckets;
3. endpoint-transfer rank witnesses appear to choose definite SCC-defect child
   fibers.

## Questions

1. Is the high-good-cut concentration a theorem after removing the transitive
   spine and its first endpoint descendants?
2. Can private pivots be constructed by a canonical "maximal barrier crossing"
   endpoint signature?
3. Do merged rows without private pivots fail because several SC pivots with
   the same SCC-defect collide after complement identification?
