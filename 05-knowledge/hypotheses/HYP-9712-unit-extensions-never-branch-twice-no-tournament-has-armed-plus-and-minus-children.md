---
id: HYP-9712
title: "Unit extensions never branch twice: no tournament T has both an armed +1 child and an armed -1 child, where a child adds one vertex with an arbitrary out-set and has H = 2H(T) +- 1, and armed means having unit extensions of both signs; so no stable plus/minus fan pair, fixed or re-chosen, generates the binary tree of odd numbers"
status: >
  OPEN (CONJECTURE), mac-mini-2026-11-15-slides, wave 26 of chanpak_waves_20261113.
  FINITE-EXACT (THM-5182 (D)): true for every tournament of order <= 9 (all 191,536 order-9 classes as seeds,
  children to order 10), with arbitrary out-sets re-chosen at every node.
  Two-vertex analogue (H -> 4H + eps, eps in {+-1, +-3}): seeds of orders 3-6 have depth <= 1 except H = 3
  (depth 2).
  Pre-registered P77, P78 REFUTED (arming is rarer and shallower than registered), P79, P80, P86 CONFIRMED.
session: mac-mini-2026-11-15-slides
date: 2026-11-15
source: 05-knowledge/results/chanpak_waves_20261113/README.md#wave-26-fans-slides-and-reservoirs--unit-extensions-never-branch-twice-a-mod-3-law-for-passive-channels-and-marked-two-block-completeness
scripts:
  - 05-knowledge/results/chanpak_waves_20261113/scripts/armed.c (cut formula, census), tree2.c (adaptive trees),
    childpair.c (children by sign), xsize.c (out-set sizes), base4.c (two-vertex steps)
related:
  - THM-5182 (D), chanpak_leaf_fans_20261115 (fan law, cross-port law, movable-port census: the leaf case),
    HYP-9702
---

# HYP-9712 - unit extensions never branch twice

**Setting.** A unit extension of T adds one vertex x with any out-set X ⊆ V(T) such that H(T + x) = 2H(T) ± 1.
By the odd-cycle expansion,

    H(T + x) = H(T) + 2W(X),   W(X) = Σ_{u∈X, v∉X} Π(u,v),

where Π(u,v) sums H(T − P) over the even-order paths P from u to v. So a unit extension is a cut of value
(H ± 1)/2 in the weighted digraph Π. T is **armed** when cuts of both values exist.

**Conjecture.** No tournament has both an armed +1 child and an armed −1 child.

A tree of unit moves would generate every odd number in [(h − 1)2^k + 1, (h + 1)2^k − 1] from a seed of count h:
the owner's "full binary tree" from a stable plus/minus fan pair. The conjecture says no such tree has depth 2,
whatever the out-sets and however they are re-chosen at each node.

**Evidence (exhaustive, orders ≤ 9).**

| seed order | eligible seeds (H ≠ 3, 11) | armed | with an armed child of each sign |
|---:|---:|---:|---:|
| 5 | 8 | 5 | 0 |
| 6 | 50 | 18 | 0 |
| 7 | 448 | 50 | 0 |
| 8 | 6,880 (all classes) | 134 | 0 (8 with an armed +1 child, 18 with an armed −1 child) |
| 9 | 191,536 (all classes) | 308 | 0 (10 with an armed +1 child, 27 with an armed −1 child) |

- About half of all unit out-sets are leaves or reverse leaves. The rest have intermediate sizes, so this goes
  beyond the leaf model of chanpak_leaf_fans_20261115.
- The children freeze their sign. Typically a +1 child has only +1 unit extensions and a −1 child only −1 ones.
  Examples at order 6:
  - H = 19: the +1 children (39) admit only +1, the −1 children (37) only −1;
  - H = 33: the +1 children (67) are not armed, the −1 children (65) are;
  - H = 9: the +1 children (19) are armed, the −1 children (17) are not.
- With the branches whose target is 7 or 21 skipped, the deepest trees are 3 → 5 → {9, 11} → {17, 19, 23}
  (depth 3, from the 3-cycle).
- Arming gets rarer: the armed fraction among eligible classes is 0.625, 0.36 and 0.11 at orders 5, 6 and 7. Among
  all classes it is 1.9% at order 8 and 0.16% at order 9.

**Consequence for proof architecture.** A binary or base-b tree of exact affine moves H → bH + ε cannot be run
with one-vertex (or, at small order, two-vertex) attachments. Any digit tree has to take multi-vertex steps, and
their intermediate defects are not units. This generalizes the leaf-model obstructions of
chanpak_leaf_fans_20261115 (the cross-port law and the repeat-stable classification) to arbitrary out-sets and
adaptive choices.

**What a proof might use.** The leaf case has a mechanism: after a unit leaf at q, a second port picks up a 3^k
channel unless q's out-neighbourhood is tiny (CROSS_PORTS). For general cuts there is a hint in the small cases.
Many unit cuts of a child of the same sign are X itself, a subset of X, or X plus the newborn. In the order-5
example with X = {2, 3} the +1 cuts of the child are {2, 3}, {3} and {2, 3, 5}. A proof could try to show that
same-sign cuts are inherited this way, while every opposite-sign cut of the child forces its own children to
inherit that sign.

**Predictions.**
- *P86 (CONFIRMED).* The conjecture holds for all seeds of order 9 (191,536 classes).
- *P87.* Among armed tournaments of order n, the fraction whose +1 children are armed and the fraction whose −1
  children are armed both tend to 0. Their sum stays at most 1 for every n.
