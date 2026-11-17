# Toothpick diagrams for A000568: the extension graph of tournament classes, its nonuniform contacts and unequal growth

Owner: *"for each term in A000568 a matching toothpick diagram can be constructed to represent how the isomorphism classes
relate to each other … unequal growth or nonuniform contacts."* The toothpick sequence (A139250) grows uniformly: every
exposed end sprouts one perpendicular toothpick per step. Tournament classes grow through one-vertex extension, where a
class of order `n` has between 1 and `n` distinct parents (its vertex deletions) and a class of order `n − 1` has up to
`2^{n−1}` children; the diagram below makes that growth visible with exactly `A000568(n)` toothpicks in generation `n`.

## 1. The extension graph, orders 3–9 (FINITE-EXACT; [receipts/extgraph.out](receipts/extgraph.out), `extgraph.py`, nauty `labelg`)

Every class and every one-vertex deletion canonically labelled; a *contact* is a pair (class, distinct deletion class).

| n | classes | contacts | distinct parents per class (min, mean, max) | children per class of order n−1 (min, mean, max) | children of the transitive `M_{n−1}` | tree branching (max, childless, mean) |
|---|---|---|---|---|---|---|
| 3 | 2 | 2 | 1, 1.00, 1 | 2, 2.0, 2 | 2 (max) | 2, 0, 2.0 |
| 4 | 4 | 7 | 1, 1.75, 2 | 3, 3.5, 4 | 4 (max) | 4, 1, 2.0 |
| 5 | 12 | 29 | 1, 2.42, 4 | 6, 7.25, 9 | 9 (max) | 9, 1, 3.0 |
| 6 | 56 | 215 | 1, 3.84, 6 | 7, 17.9, 24 | 21 | 21, 4, 4.67 |
| 7 | 456 | 2591 | 1, 5.68, 7 | 14, 46.3, 58 | 48 | 48, 22, 8.14 |
| 8 | 6880 | 50357 | 1, 7.32, 8 | 11, 110.4, 121 | 106 | 106, 215, 15.1 |
| 9 | 191536 | 1656619 | 1, 8.65, 9 | 22, 240.8, 248 | 227 | 227, 2977, 27.8 |

- **Nonuniform contacts.** At order 9 a class has on average 8.65 distinct parents of a possible 9 (96 %: deletions are
  almost never isomorphic), while the number of children of an order-8 class runs from 22 to 248 of a possible 256 — the
  most fertile class realises 97 % of its extension patterns as distinct classes. Every class has at least one child (a
  sink padding), so there are no dead ends in the full graph. R-3 confirmed.
- **The transitive trunk.** The transitive tournament has the most children only through order 4; from order 5 it is
  beaten (21 < 24, 48 < 58, 106 < 121, 227 < 248) by classes with trivial automorphism group. In the canonical-parent
  tree (each class attached to its lexicographically least canonical deletion) the transitive tournament is nevertheless
  the largest branch at every order, because it is the least string: 227 of the 191,536 order-9 classes hang on `M_8`.
- **Unequal growth.** The tree branching at order 8 → 9 has maximum 227, mean 27.8, and 2,977 of the 6,880 order-8 classes
  (43 %) have no tree child at all (R-4 confirmed). First creators have fewer tree children than other classes at every
  order (4.3 vs 13.6 at order 6, 22.0 vs 50.4 at order 8; [receipts/tree_creators.out](receipts/tree_creators.out)).

## 2. The diagram

Construction (`toothpick.py`): one toothpick per class; generation `n` perpendicular to generation `n − 1`, as in the
toothpick sequence; each class touches its tree parent at an interior point of the parent toothpick, the children of one
parent at distinct points along it (the nonuniform contacts), each child leaving to its own side of the plane in a slice
proportional to the size of its subtree; first creators in red. Figures: [figures/toothpick_7.svg](figures/toothpick_7.svg)
(532 toothpicks, orders 1–7, 370 first creators) and [figures/toothpick_8.svg](figures/toothpick_8.svg) (7,412 toothpicks,
5,827 first creators). The matching with A000568 is exact by construction: generation `n` has `A000568(n)` toothpicks and
the contacts are the extension relation restricted to least deletions; the full relation (every deletion) is the graph of
§1 and is not planar in any useful sense at order 9 (1.66 million contacts).

## 3. Reading the growth

The toothpick sequence's growth is dictated by geometry (ends get blocked); here it is dictated by symmetry and by the
count. Children per class are many when the automorphism group is trivial, few when it is large (the regular tournaments
of orders 5 and 7, the Paley tournament); the counts those children carry are the values of the parent's quadratic form
([GRAMMAR.md](GRAMMAR.md) §1), 4 to 97 distinct values among the 128 corners at order 7. The spectrum grows faster than the
number of classes: 85 % of the order-9 classes carry a count that no smaller tournament has.
