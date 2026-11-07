---
id: HYP-9512
title: "First creators beyond order 10: (A) every odd h >= 5 other than 7, 21 and 63 has a first creator obtained from a first creator of order n(h) - 1 by adding one vertex, and 63 is the only orphan value with h >= 5; (B) no non-first-creator of order n >= 3 has every vertex deletion a first creator (no hypo analogue); (C) the share of fully non-critical first creators (every deletion a first creator) tends to 1; (E) the blue share among values not forced black by THM-791 has a limit near 0.83"
status: "(A)-(C) FINITE-EXACT through order 10 (pre-registered Q1-Q3 CONFIRMED blind at order 10); OPEN beyond. The former (D), the odd-order self-converse balance mod 4, is now PROVED as THM-4922 (4) (audit proof). (E) EXPLORATORY: the blue share is 0.72, 0.66, 0.84 and 0.825 at orders 7-10 (82% in aggregate). Pre-registered Q5 FAILED, because THM-791 forces every new value = 3 (mod 4) at even order to be black. AUDIT: independent, PASS WITH CORRECTIONS (MISTAKE-735)."
session: kind-pasteur-2026-11-07-S16
date: 2026-11-07
related: [THM-4922, THM-791, THM-344, THM-4748, THM-4870, THM-4889, HYP-9481, HYP-9359, 05-knowledge/results/session224_first_occurrence_20261107/README.md]
script: 05-knowledge/results/first_creators_genealogy_20261107/classgen.c; 05-knowledge/results/first_creators_genealogy_20261107/order10.c
---

# HYP-9512 - first-creator genealogy, the orphan 63, and the value colours

Packet: [README](../results/first_creators_genealogy_20261107/README.md),
[PROOF](../results/first_creators_genealogy_20261107/PROOF.md),
[AUDIT](../results/first_creators_genealogy_20261107/audit/AUDIT.md). The proved facts are THM-4922.

## (A) Genealogy

First creators grow from first creators by adding one vertex. The only exception at orders 4-10 is the converse pair
of first creators of 63 at order 8.
- These are the Fibonacci fan F_6 and its converse (THM-4748, THM-344, the Session24 multiplicities).
- Each is a transitive 7-tournament plus an apex, which is why one deletion has H = 1.
- 63 is order 7's isolated hole (klein), and HYP-9359 asks why.

Conjecture: 63 is the only orphan value with h >= 5. THM-4922 (2) proves the genealogy above n max H(n - 2); the open
range lies below that threshold, where holes like 63 live.

If (A) holds, the new values of each order n are exactly the new values among the one-vertex extensions of the first
creators of order n - 1, plus 63 at n = 8. That is a recursive generation of the spectrum, in the direction the owner
proposed.

## (B), (C) The hypohamiltonian analogue

A vertex is non-critical, in THM-4870's sense, when its deletion stays in the class of first creators; k(T) counts the
non-critical vertices. The fully non-critical first creators are 24%, 39%, 50% and 62.8% of all first creators at
orders 7-10.

Conjecture: the share tends to 1. Also, "not a first creator, yet every deletion is one" never happens at order >= 3.
This is the analogue of the non-existence of hypohamiltonian and hypotraceable tournaments (Rédei, Camion).

## (E) The value colours

A value is blue if some first creator is self-converse, and black if not. These value colours are distinct from the
metagraph's line colours.
- At even least order, every value 3 (mod 4) is black (THM-791).
- Elsewhere the blue share is 0.72, 0.66, 0.84 and 0.825 at orders 7-10.
- At odd orders the self-converse classes are exactly balanced mod 4 (THM-4922 (4), proved).

**Lesson from Q5.** A prediction about the share of a kind must first remove the cells a known theorem forces. Here
THM-791 forces the black cell at even orders.
