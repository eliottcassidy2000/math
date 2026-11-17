# First creators as the grammar's new words: most classes are new, the genealogy is degree-blind, and 63

Owner: *"for every odd number other than 7 and 21 there is an isomorphism class or set of isomorphism classes that relate to
the first appearance of those counts at the smallest size, and there must be some relationship between these
constructions and themselves."* THM-4922 defines the first creators (classes of least order with a given count), proves
their fibres converse-closed and strong (≤ 11), and finds the genealogy (every first creator at orders 4–10 except the pair
creating 63 has a vertex deletion that is a first creator). Here the genealogy is typed by the deleted vertex, and the
first creators are placed in the grammar of [GRAMMAR.md](GRAMMAR.md) and the tree of [TOOTHPICK.md](TOOTHPICK.md).

## 1. Census, orders 3–9 ([receipts/creators9.out](receipts/creators9.out), `creators9.c`, 6 s)

| order | new counts | first creators | share of all classes | with a first-creator deletion | with such a deletion of out-degree 1 or n−2 | with a digit-0 vertex |
|---|---|---|---|---|---|---|
| 3 | 1 | 1 | 50 % | 0 (the 3-cycle, orphan) | 0 | 1 |
| 4 | 1 | 1 | 25 % | 1 | 1 | 1 |
| 5 | 4 | 6 | 50 % | 6 | 4 | 3 |
| 6 | 12 | 33 | 59 % | 33 | 27 (82 %) | 7 (21 %) |
| 7 | 58 | 328 | 72 % | 328 | 235 (72 %) | 24 (7.3 %) |
| 8 | 243 | 5457 | 79 % | 5455 (the 63 pair orphan) | 3017 (55 %) | 92 (1.7 %) |
| 9 | 1200 | 163676 | 85 % | 163676 | 63060 (38.5 %) | 610 (0.37 %) |

The 1,520 values with least order ≤ 9 and the orphans (3 at order 3; the converse pair creating 63 at order 8) agree
with THM-4922. **First creators are the rule, not the exception**: 85 % of the order-9 classes carry a count no smaller
tournament has. The spectrum grows faster than the class count because the quadratic forms of GRAMMAR.md §1 take new
values on most corners of most cubes.

## 2. The genealogy is nearly blind to the type of the deleted vertex (R-5 refuted)

Out-degree of the deleted vertex `v` over all genealogy edges `(F, v)` (`F − v` a first creator):

| order | d = 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| 7 | 165 | 449 | 621 | 449 | 165 | | |
| 8 | 1846 | 6331 | 10814 | 10814 | 6331 | 1846 | |
| 9 | 35176 | 144462 | 296548 | 372259 | 296548 | 144462 | 35176 |

Symmetric (first creators are converse-closed) and close to binomial: at order 9 the counts divided by `C(8, d)` are
4397, 5159, 5295, 5318 for `d = 1, 2, 3, 4` — a vertex of a first creator yields a first creator on deletion at nearly
the same rate whatever its out-degree, only slightly less at the extremes. The mode is the middle degree, not a sole
exit or sole entry; the share of first creators with a sole-exit or sole-entry genealogy edge falls from 82 % (order 6)
to 38.5 % (order 9). **So the recursive grammar of first creators is the full extension grammar, not its additive
(near-join, tower, dead-end) part**: the additive productions, which are the ones that encode the holes by their 2-adic
structure (GRAMMAR.md §3 (b)), account for a vanishing share of new words — the digit-0 share of first creators is 21 %,
7.3 %, 1.7 %, 0.37 % at orders 6–9, in line with the 95 % sterile share of all strong classes (THM-5264 (F)).

## 3. The 63 orphans ([receipts/sixty_three.out](receipts/sixty_three.out))

The two first creators of 63 (order 8, a converse pair, both strong) have deletion counts `1, 11, 25, 25, 33, 33, 37, 37`
— every deletion's count was born earlier (1, 11, 25, 33, 37 at orders 1, 5, 6, 6, 6), so neither is grown from a first
creator (THM-4922's orphans). But they are not orphans of the grammar: one of the pair (out-degrees `1, 1, 2, 4, 4, 5, 5, 6`)
has a digit-0 vertex `z`, so `63 = h(S − z) + B_t(S − z)` is a dead-end append over an order-7 base of count 25 or 11; its
converse (out-degrees `1, 2, 2, 3, 3, 5, 6, 6`) has none — a digit-0 vertex is a sole exit, its converse a sole entry.
63 = 3·21 = 9·7 is the atom through both holes (waves 48, 55), the Zsigmondy exception and the period-7 Mandelbrot count
(THM-5265 (E)); as a word of the grammar it is a single dead-end append.

## 4. In the toothpick tree ([receipts/tree_creators.out](receipts/tree_creators.out))

The tree parent (least canonical deletion) of a first creator is itself a first creator for 17 %, 27 %, 44 %, 61 %, 73 % of
the first creators of orders 5–9 — below the base rate of first creators among all classes at the parent order (25 %,
50 %, 59 %, 72 %, 79 %), because the transitive trunk, never a first creator after order 1, is the least string and
collects the most children. First creators have fewer tree children than other classes (22.0 vs 50.4 at order 8). R-7
refuted below order 8.
