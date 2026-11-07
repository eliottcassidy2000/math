---
id: HYP-9512
title: "First creators beyond order 10: (A) every odd h other than 7, 21 and 63 has a first creator obtained from a first creator of order n(h) - 1 by adding one vertex, and 63 is the only orphan value; (B) no non-first-creator ever has every vertex deletion a first creator (no hypo analogue); (C) the fraction of fully critical first creators (every deletion a first creator) tends to 1; (D) at odd orders the SC classes are equidistributed between H = 1 and H = 3 (mod 4), exactly (1/1, 4/4, 44/44, 1376/1376 at n = 3, 5, 7, 9); (E) the blue share among values not forced black by THM-791 stays near 0.83"
status: "(A)-(C) FINITE-EXACT through order 10 (pre-registered Q1-Q3 CONFIRMED at order 10); OPEN beyond. (D) FINITE-EXACT at n = 3, 5, 7, 9; OPEN in general (not found in the repository; THM-791 covers even orders). (E) EXPLORATORY: 84% (odd orders <= 9), 82.5% (order 10, h = 1 mod 4). Pre-registered Q5 (mixed fibers the majority at order 10) FAILED, because THM-791 forces every new value = 3 (mod 4) at even order to be black. AUDIT: pending."
session: kind-pasteur-2026-11-07-S16
date: 2026-11-07
related: [THM-4922, THM-791, THM-344, THM-4870, HYP-9481, HYP-9359, THM-4748]
script: 05-knowledge/results/first_creators_genealogy_20261107/classgen.c; 05-knowledge/results/first_creators_genealogy_20261107/order10.c
---

# HYP-9512 - first-creator genealogy, the orphan 63, and the blue/black split

Packet: [README](../results/first_creators_genealogy_20261107/README.md),
[PROOF](../results/first_creators_genealogy_20261107/PROOF.md). The proved facts are THM-4922.

## (A) Genealogy

First creators grow from first creators by adding one vertex, with one exception through order 10: the converse pair
of first creators of 63 at order 8.
- 63 is order 7's isolated hole (klein).
- It is a Mersenne number 2^6 - 1, and such values need at least three back arcs (THM-4748).
- It is an outlier in HYP-9359.

Conjecture: 63 is the only orphan value. THM-4922 (2) proves the genealogy above n max H(n - 2); the open range is
below that threshold, where holes like 63 live. If (A) holds, the new values of each order are exactly the values of
one-vertex extensions of the previous order's first creators, plus the single orphan 63. That is a recursive generation
of the spectrum, in the direction the owner proposed.

## (B), (C) The hypohamiltonian analogue

A vertex is critical when its deletion stays in the class of first creators, in THM-4870's sense of non-critical
reductions. The fully critical first creators are 24%, 39%, 50% and 62.8% of all first creators at orders 7-10.

Conjecture: the share tends to 1. Also, "not a first creator, yet every deletion is one" never happens. This is the
analogue of the non-existence of hypohamiltonian and hypotraceable tournaments.

## (D), (E) The two kinds

A value is **blue** if some first creator is self-converse, and **black** if not.
- At even least order, every value 3 (mod 4) is black (THM-791).
- Elsewhere about 83% of values are blue.
- At odd orders the SC classes split exactly evenly between H = 1 and 3 (mod 4).

The owner's proposed dichotomy (a single class or a single converse pair) is false: it holds only for small values and
near each order's maximum (THM-4922 (3)).

**Lesson from Q5.** A prediction about the share of a kind must first remove the cells a known theorem forces. Here
THM-791 forces the black cell at even orders.
