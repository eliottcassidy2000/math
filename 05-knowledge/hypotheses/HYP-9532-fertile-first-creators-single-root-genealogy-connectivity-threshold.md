---
id: HYP-9532
title: "First creators beyond order 11: (A) every first creator has an extension that is a first creator, and at least half of its 2^n extension patterns are fertile except for the two creators of 63; (B) every first creator other than the two creators of 63 contains a chain of first creators, one at each order, starting from the 3-cycle; (C) the share of new values reached only from first-creator parents tends to 1; (D) the generation threshold k*(n) = max arc-strong connectivity of a non-first-creator of order n is unbounded; (E) no odd value above 63 is born at its Busch bound; (F) whether every tournament is a first creator, a first creator plus a vertex, an ordinal sum, a substitution, or an ordinal sum re-connected by at most two arcs"
status: "(A) FINITE-EXACT at orders 3-10 (pre-registered R3 and R4 CONFIRMED at order 10). (B) FINITE-EXACT through order 10. (C) EXPLORATORY: 0.350, 0.425, 0.470, 0.515 at orders 8-11. (D) k*(n) = 1 for n <= 8, 2 for n = 9, 10, and at least 3 at n = 11 (pre-registered R11 CONFIRMED); OPEN beyond. (E) FINITE-EXACT through 95,095, conditional beyond on HYP-9481 (a). (F) FINITE-EXACT through order 11; OPEN question beyond."
session: kind-pasteur-2026-11-08-S17
date: 2026-11-08
related: [THM-4942, THM-4922, HYP-9512, HYP-9481, THM-1370, THM-1900, THM-4744]
script: 05-knowledge/results/first_creator_generation_20261108/ext.c; 05-knowledge/results/first_creator_generation_20261108/lam11.c
---

# HYP-9532 - fertile first creators, a single genealogical root, and the connectivity threshold

Packet: [README](../results/first_creator_generation_20261108/README.md),
[PROOF](../results/first_creator_generation_20261108/PROOF.md). The proved and finite-exact facts are THM-4942.

## (A) Fertility

An extension pattern of a first creator T is **fertile** if the extension is again a first creator.
- At orders 3-10 every first creator has a fertile pattern.
- All but the two creators of 63 have at least half their patterns fertile. The two creators of 63 have 77 of 256.
- The minimum at order 10 is 564 of 1,024.

Conjecture: both statements hold at every order. Together with HYP-9512 (A), the minimal set is then
self-propagating in both directions: it grows from itself, and every member has offspring in it.

## (B) A single root

At order 9, the 152 first creators that extend a creator of 63 all have another first-creator parent. Through order
10, every first creator except the creators of 63 therefore lies on a chain of first creators that starts at the
3-cycle. These are Moon chains inside the minimal set.

Conjecture: this holds at every order. Then the 63 pair is a side branch, and the 3-cycle is the only genealogical
root.

## (C) Parents of new values

Consider the share of new values whose every realization from one order down has a first-creator parent. It rises:
0.350, 0.425, 0.470 and 0.515 at orders 8-11.

Conjecture: it tends to 1, in line with the first-creator share of all classes (72%, 79%, 85.5% and 88.8% at orders
7-10).

## (D) The connectivity threshold

Let k*(n) be the largest arc-strong connectivity of a non-first-creator of order n. Then:
- k*(n) = 1 for n <= 8;
- k*(n) = 2 at orders 9 and 10;
- k*(11) >= 3, by exhaustive search. f_3(11) = 11,475 < max H(10). Four of the seven 3-arc-strong values below
  15,745 (11,475, 15,213, 15,273 and 15,561) lie in S(10).

Conjecture: k*(n) is unbounded.

Heuristic. For fixed k the floor f_k(n) is at most exponential in n, through chain constructions. Meanwhile max
H(n - 1) grows faster than any exponential. So the k-th floor eventually falls below max H(n - 1). Low-H
k-arc-strong tournaments then meet the attained range, which needs a completeness input such as HYP-9481 (a).

## (E) Values born at their Busch bound

Through 95,095, n(h) = B(h) holds exactly for 3, 5, 9, 11, 13, 17, 19, 23, 35, 39 and 63.

Conjecture: there are no others.

For every h > 95,095, m(h) <= B(h) - 2:
- B(h) >= 22 there, because f(22) = 78,125.
- f grows exponentially, while max H, and hence m, outruns it. For example, B(h) >= 25 once h > 545,835 = max H(12).

So n(h) = B(h) would give n(h) >= m(h) + 2, and HYP-9481 (a) (n(h) <= m(h) + 1) excludes it. The conjecture
therefore follows from HYP-9481 (a).

## (F) The five-operation grammar

Through order 11, every tournament is one of:
- a first creator;
- a first creator plus a vertex;
- an ordinal sum;
- a substitution;
- an ordinal sum re-connected by at most two arcs.

At order 11 the seven non-first-creator classes with lambda >= 3 are all substitutions. Their floor is C3[R5, R5, 1].

Question: does this hold at every order? If not, how many re-connecting arcs are needed, as a function of n? A first
failure would be a prime, 3-arc-strong non-first-creator with no first-creator deletion.
