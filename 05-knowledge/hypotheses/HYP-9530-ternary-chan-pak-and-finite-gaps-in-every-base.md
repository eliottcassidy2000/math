---
id: HYP-9530
title: "Packing counts in every base: (A) ternary Chan-Pak - the base-3 packing counts I(Omega(T), 3) take every value congruent to 1 mod 3 except 10 and 31, the same alpha-sums 3 and 10 that give 7 and 21 in base 2; (B) every base has finitely many gaps; (C) the carry-free layer has infinitely many singletons, while binary carries leave only 3, 5, 11, 13"
status: >
  OPEN.
  (A): pre-registered P1 CONFIRMED. It is FINITE-EXACT and complete at all orders for alpha-sums <= 80, i.e. values
  <= 241.
  (B): consistent with the exact tables to alpha-sum 80 for b <= 90.
  (C): FINITE-EXACT below 1024 (26 carry-free singletons); the binary half is HYP-9477 (A).
session: opus-2026-11-08-S86
source: 05-knowledge/results/packing_bases_repunits_20261108/README.md
related:
  - THM-4940
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
- Below 242 these are the only gaps (THM-4940 (6), pre-registered P1, CONFIRMED).

**Comparison with base 2.**
- Base 2 (Chan-Pak) and base 3 have the same gap alpha-sums {3, 10}.
- By THM-4940 (7) they are the only bases that can share them.
- In both, every hole of A other than 3 and 10 is rescued by carries, at least up to alpha-sum 80.

**Numerology, typed and not used.**
- The second gap 1 + 10b is the repunit 111_(b+2) exactly for b in {2, 3}: 21 = 111_4 and 31 = 111_5.
- The first gap is 7 = 111_2 in base 2.
- MISTAKE-559 applies: no mechanism is claimed.

**Falsifier.** A value 1 + 3m (m not 3 or 10) that is not a 3-coloured packing count of any tournament. Every m <= 80
is already settled.

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
- Infinitely many a in A have exactly one strong intercyclic realiser. Below 1024 there are 26, all self-converse,
  the largest being 933.
- Every such a >= 8 has another strong realiser of H = 1 + 2a in base 2, coming from carries.

**What is known.**
- The second half is HYP-9477 (A) restricted to H = 1 + 2a with a in A. It is FINITE-EXACT below 625 (THM-4884), where
  the binary strong singletons are exactly 3, 5, 11, 13, i.e. a = 1, 2, 5, 6.

**Reading.** In base 2, carries both fill holes (Chan-Pak keeps only 3, 10) and destroy uniqueness (Carmichael-type
singletons keep only a <= 7). Both families of binary exceptions lie in the carry-free range a <= 10.

**Falsifier.** For the first half: a finite list, i.e. no carry-free singletons above some bound. For the second half:
a carry-free singleton a >= 8 whose binary value has a single strong class.
