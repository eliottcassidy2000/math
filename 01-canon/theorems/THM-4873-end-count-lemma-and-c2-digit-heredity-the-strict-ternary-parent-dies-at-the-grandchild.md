---
id: THM-4873
title: "The end-count lemma and C2 digit heredity. (i) In a tournament with H >= 3 Hamiltonian paths no vertex ends exactly H - 1 of them, and a vertex ends exactly H - 2 iff the tournament is TT_j => C3 and the vertex lies in the 3-cycle (dually for starts). (ii) Hence in the C2 child of (B, a) the native mark y (digit h(B) - e_a(B)) has digit 1 iff B is transitive and a is not its sink, and digit 2 iff B = TT_j => C3 with a in the 3-cycle; in every C2 grandchild the native mark never has digit 1. (iii) FINITE-EXACT over all marks, orders 1-7: no C2 grandchild carries a digit-1 mark, and only the three native C3-bottom cases carry a digit-2 mark. So the strict ternary scheme (each node's own C2 marks supplying digits 0, 1, 2) dies by the second generation, and a ternary parent for full spectrum coverage must come from fresh bases"
status: >
  PROVED (opus-2026-11-04-S77, elementary): (i) and (ii). The native digits sigma_x = 0 and sigma_y = h - e_a are
  klein's tower lemmas (THM-4865), re-verified. The C2 law H = h + 2B_a at any mark is inherited
  (multiport_defect_grammar law (4); THM-4865).
  FINITE-EXACT: (iii), all classes of order 1-7 with every mark at both levels; order 7 pre-registered as HYP-9466 and
  CONFIRMED.
  NOT CLAIMED: anything about Chan-Pak completeness, which remains OPEN (HYP-9454). The any-mark one-step sweep computed
  here re-derives klein's height-two numeration.
  AUDIT: see the note, section 7.
session: opus-2026-11-04-S77
source: 05-knowledge/results/ternary_parent_c2_digit_heredity_20261104.md
related:
  - THM-4862
  - THM-4865
  - HYP-9454
  - HYP-9466
  - THM-4872
  - 05-knowledge/results/session191_tournament_growth_20261104/ternary_fan/PROOF.md
  - 05-knowledge/results/nonzero_preparation_20261103/README.md
verification:
  - 04-computation/experiments/c2_digit_heredity_20261104.c
  - 04-computation/experiments/c2_digit_heredity_20261104.out
  - 04-computation/experiments/anymark_c2_digits_20261104.c
  - 04-computation/experiments/ternary_parent_digit_families_20261104.c
---

# THM-4873 - the end-count lemma; the strict ternary parent dies at the grandchild

## (i) End-count lemma

**Statement.** Let C be the last strong component of T, of order m, and M = prod H(earlier components). Then every
Hamiltonian path ends in C, and e_v(T) = M e_v(C) for v in C.
- If m >= 3, every vertex of C ends at least one path of C (open a Hamiltonian cycle).
- So a vertex ending H - d paths forces (m - 1) M <= d.

**Consequences.**
- d = 1 is impossible for m >= 3. For m = 1 the single vertex of C ends all H paths, so d = 1 is impossible when H > 1.
- d = 2 forces m = 3 and M = 1, i.e. T = TT_j => C3 with H = 3.

## (ii) Native digits

In C2(B, a) = B + x + y, with x -> a only and y -> x, every vertex of B beating y:
- sigma_x = 0;
- sigma_y = h(B) - e_a(B).

By (i):
- sigma_y = 1 iff h(B) = 1 and e_a = 0, i.e. B transitive and a not its sink;
- sigma_y = 2 iff B = TT => C3 and a lies in the 3-cycle.

**In a grandchild C2(T1, u),** the native digit is H(T1) - e_u(T1), with H(T1) = h + 2B_a >= 3.
- It is never 1.
- It is 2 iff T1 = TT => C3. That means h(B) = B_a(B) = 1, so B is transitive with its sink as mark, and u lies in the
  3-cycle. These are three grandchildren per order.

## (iii) All marks (FINITE-EXACT, orders 1-7)

| \|B\| | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|
| children with a digit-1 mark | 1 | 2 | 3 | 4 | 5 | 6 |
| children with a digit-2 mark | 1 | 4 | 4 | 4 | 4 | 4 |
| grandchildren with a digit-1 mark | 0 | 0 | 0 | 0 | 0 | 0 |
| grandchildren with a digit-2 mark | 3 | 3 | 3 | 3 | 3 | 3 |

The orders checked total 3,192 + 2,688 + ... (B, a) pairs and 28,728 grandchildren at order 7.

**Reading.** A ternary parent of N, a realiser of the unique odd predecessor with digit 0, 1 or 2, cannot be supplied
hereditarily by C2. The same is true of:
- H-neutral preparation (needs membership in R_1 and R_2, densities about 0.75);
- fixed-core wrapper fans (Session191, (61, 63)).

Any proof of full coverage must use fresh bases with growing digits (HYP-9454 (1); FINITE-EXACT to 8,079,323 by klein).
