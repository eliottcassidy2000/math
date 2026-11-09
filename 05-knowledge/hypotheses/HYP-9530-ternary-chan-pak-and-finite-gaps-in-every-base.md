---
id: HYP-9530
title: "Packing counts in every base: (A) ternary Chan-Pak - the base-3 packing counts I(Omega(T), 3) take every value congruent to 1 mod 3 except 10 and 31, the same alpha-sums 3 and 10 that give 7 and 21 in base 2; (B) every base has finitely many gaps; (C) the carry-free layer has infinitely many singletons, while binary carries leave only 3, 5, 11, 13"
status: >
  OPEN.
  2026-11-09 update (kind-pasteur-2026-11-09-S18, THM-5070):
  - (A) is FINITE-EXACT to values <= 10^6 (alpha-sums <= 333,333): port chains with blocks of order <= 5 attain every
    v = 1 (mod 3) there except 10 and 31. The audit extended the 15-generator family to 3*10^6. HYP-9640 (A) reduces
    (A) to the cofiniteness of an explicit 15-matrix orbit.
  - (B) is FINITE-EXACT to alpha-sum 30,000 for b = 4..8: the gaps there are exactly those below 81 (THM-5070 (6),
    HYP-9640 (E)).
  Earlier status, kept for the record:
  (A): FINITE-EXACT and complete for alpha-sums <= 9678 (values <= 29,035): the gaps there are exactly 10 and 31.
  THM-4940 (3), (4), (9) gave this to alpha-sum 1799 (values <= 5398). THM-4950 (opus-2026-11-09-S87) extends it with
  the complete strong order-10 census; its pre-registered Q1-Q3 are CONFIRMED. THM-4950 also reduces (A) to the
  strong realisation of the irreducibles of M3 = {v = 1 (mod 3)} minus {10, 31}: the lanes P, Q, 10P, 31P, 31Q and
  five sporadic values. It is OPEN beyond. Pre-registered P1 (exact tables to alpha-sum 80) is CONFIRMED; the
  extension to 1799 is constructive and was not pre-registered.
  (B): consistent with the exact tables to alpha-sum 80 for b <= 90.
  (C): FINITE-EXACT below 8192: 64 carry-free singletons, all self-converse, and every one with a >= 8 has a strong
  non-intercyclic binary twin (session to 1024, audit to 8192). The binary half beyond that is HYP-9477 (A).
  Audited with THM-4940 (PASS WITH CORRECTIONS; MISTAKE-739).
session: opus-2026-11-08-S86
source: 05-knowledge/results/packing_bases_repunits_20261108/README.md
related:
  - THM-4940
  - THM-4950
  - HYP-9477
  - THM-1370 (h-spectrum-omits-7-21-all-n)
  - THM-4884
---

# HYP-9530 - ternary Chan-Pak and finite gaps in every base

**Setting.** As in THM-4940:
- H_b(T) = I(Omega(T), b) is the number of b-coloured odd-cycle packings;
- Gap_b is the set of alpha-sums m with 1 + b m not a value;
- A is the set of odd-cycle counts of intercyclic tournaments.

## (A) Ternary Chan-Pak

**Statement.** Gap_3 = {3, 10}. Every integer congruent to 1 (mod 3) except 10 and 31 is the number of 3-coloured
odd-cycle packings of some tournament.

**What is proved.**
- 10 and 31 are permanent gaps. 10 is THM-4940 (3), at all orders. 31 is THM-4940 (4) with b = 3: its candidate
  vectors (10, 0), (7, 1) and (4, 2) are excluded at all orders by the complete inventory.
- Below 242 these are the only gaps: THM-4940 (6), pre-registered P1, CONFIRMED.
- Below 5399 they are still the only gaps: THM-4940 (9), constructive, with every other value attained by a tournament
  of order <= 9.
- Below 29,036 they are still the only gaps (THM-4950 (3), opus-2026-11-09-S87). This is FINITE-EXACT by explicit
  witnesses: every value is an order join of strong tournaments of order <= 10. The complete strong order-10 census
  has 10,070 values, independently re-derived by the audit. S86's other pieces add nothing below 29,038.
- **Reduction (THM-4950 (1), (2)).** (A) is equivalent to: every irreducible of the monoid
  M3 = {v = 1 (mod 3)} minus {10, 31} is H_3 of a strong tournament. The irreducibles are:
  - the primes = 1 (mod 3) other than 31;
  - the products of two primes = 2 (mod 3) other than 10;
  - 10p, 31p and 31pq;
  - 40, 250, 310, 9610, 29791.
  The first irreducible not reached by the current pieces is 29,038 = 2 * 14,519.

**Comparison with base 2.**
- Base 2 and base 3 have the same known gap alpha-sums {3, 10}: base 2 up to alpha-sum 2^34 (all-odd coverage),
  base 3 up to alpha-sum 9678 (THM-4950; 1799 in THM-4940 (9)). Chan-Pak and (A) assert this at every alpha-sum.
- The ternary pair 10, 31 was first observed at n <= 6 by `04-computation/ip_hierarchy.py` (opus-S71e, March 2026),
  which asked whether they are "gaps analogous to H = 7 and H = 21". (A) is that question with an all-order proof of
  both gaps and exact data to alpha-sum 80.
- By THM-4940 (7) they are the only bases that can share them.
- In both, every hole of A other than 3 and 10 is rescued by carries, at least up to alpha-sum 80.

**Numerology, typed and not used.**
- The second gap 1 + 10b is the repunit 111_(b+2) exactly for b in {2, 3}: 21 = 111_4 and 31 = 111_5.
- The first gap is 7 = 111_2 in base 2.
- MISTAKE-559 applies: no mechanism is claimed.

**Falsifier.** A value 1 + 3m (m not 3 or 10) that is not a 3-coloured packing count of any tournament. Every
m <= 9678 is already settled (THM-4950).

## (B) Finite gap sets

**Statement.** For every b >= 2, Gap_b is finite.

**What is known.**
- **Lower bound.** Gap_b contains every hole of A below b, so |Gap_b| grows at least like a positive proportion of b
  (Session 16's density bound).
- **Exact data to alpha-sum 80.** Bases 4 to 8 have at most five gaps: {3, 14, 27}, {3, 10, 27}, {3, 10, 14, 21, 46},
  {3, 10, 14, 29}, {3, 35, 52}. Each gained at most one new gap on (44, 80] (pre-registered P3, CONFIRMED).
- **Mechanism suggested by the rescue sets.** As m grows, the number of carried decompositions sum alpha_k b^(k-1) = m
  grows, and the realisable vectors with alpha_2 >= 1 cover more of them.

**Falsifier.** A base with gaps at arbitrarily large alpha-sums.

## (C) Carry-free singletons are infinite; binary carries kill all but four

**Statement.**
- Infinitely many a in A have exactly one strong intercyclic realiser:
  - 26 below 1024 (largest 933);
  - 64 below 8192 (audit), all self-converse.
  - Above order 9 they occur only at odd orders.
- Every such a >= 8 has another strong realiser of H = 1 + 2a in base 2, coming from carries.

**What is known.**
- The second half is HYP-9477 (A) restricted to H = 1 + 2a with a in A.
- It is FINITE-EXACT below 625 (THM-4884), where the binary strong singletons are exactly 3, 5, 11, 13, i.e.
  a = 1, 2, 5, 6.
- For the carry-free singletons it is FINITE-EXACT below 8192: each of the 60 with a >= 8 has an explicit strong
  twin with two disjoint cyclic triangles. These come from THM-4940 (8) and the audit's `aud_twins_hi.out`.

**Reading.** In base 2, carries both fill holes and destroy uniqueness:
- holes: Chan-Pak, OPEN, says only 3 and 10 survive;
- uniqueness: HYP-9477 (A), OPEN, says only a <= 7 survive. This is verified below 625, and for the carry-free
  singletons below 8192 by the audit.

The known binary exceptions of both kinds lie in the carry-free range a <= 10.

**Falsifier.** For the first half: a finite list, i.e. no carry-free singletons above some bound. For the second half:
a carry-free singleton a >= 8 whose binary value has a single strong class.
