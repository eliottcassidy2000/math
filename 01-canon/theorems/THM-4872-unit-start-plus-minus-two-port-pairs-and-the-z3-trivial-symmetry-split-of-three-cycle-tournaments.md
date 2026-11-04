---
id: THM-4872
title: "(A) The unit-start +-2 port pair of order 7. Among strong rooted 7-vertex tournaments with one Hamiltonian path from the root, exactly one pair has equal two-slot count F = 157 and path counts 29, 31. Its classes are 4 strong arc reversals apart (38 shortest routes, none keeping the unit root start, which forces distance 7; 2 keeping every intermediate wrapper below the start). So Session189's wrappers H(D_m) = h + 2^m F differ by exactly 2 at every clock. Such pairs occur sporadically at every order 7-14 (1, 2, 5, 1, 6, 10, 4, 8). (B) The sole-exit-free tournaments with exactly three Hamiltonian cycles split by symmetry as F(2n - 9) = F(2n - 10) (Aut = Z_3, C3-block cores) + F(2n - 11) (asymmetric); those with four or five cycles are asymmetric"
status: >
  (A) FINITE-EXACT (opus-2026-11-04-S76): exhaustive over all 2^21 labelled tournaments of order 7 and over all classes
  to order 9; orders 10-14 by the FUNNELS F1/F2 grammar (PROVED complete there), which reproduces F(2n - 4). The wrapper
  formula is Session189's (PROVED there), re-verified on 220 cases.
  (B) PROVED from THM-4869 (B1, B4, B5) and THM-4868 (E); FINITE-EXACT at orders 6-14.
  PRE-REGISTERED: HYP-9465 P1-P3 (orders 13-14), all CONFIRMED.
  AUDIT: see the note, section 7.
session: opus-2026-11-04-S76
source: 05-knowledge/results/unit_start_pm2_port_pairs_symmetry_split_mrb_20261104.md
related:
  - THM-4869
  - THM-4868
  - THM-4870
  - THM-4862
  - HYP-9465
  - HYP-9464
verification:
  - 04-computation/experiments/plus_minus_two_port_routes_20261104.c
  - 04-computation/experiments/plus_minus_two_port_routes_20261104.out
  - 04-computation/experiments/unit_start_pm2_pairs_20261104.c
  - 04-computation/experiments/unit_start_pm2_pairs_20261104.out
  - 04-computation/experiments/unit_start_pm2_pairs_grammar_20261104.c
  - 04-computation/experiments/unit_start_pm2_pairs_grammar_20261104_14.out
  - 04-computation/experiments/fewcycle_symmetry_and_unit_pairs_20261104.c
  - 04-computation/experiments/fewcycle_symmetry_and_unit_pairs_20261104_14.out
---

# THM-4872 - unit-start +-2 port pairs; the Z_3 / trivial split of the three-cycle counts

## (A) The order-7 pair

**Profile.** A rooted tournament (T, v) has profile (h, s, q, F): H(T); Hamiltonian paths from v; H(T - v); spanning
two-slot covers with the first path ending at v (ADAPTIVE.md section 2). The two classes, with root 0:

| | out-rows | (h, s, q, F) |
|---|---|---|
| A | 8 73 3 4 15 31 61 | (29, 1, 15, 157) |
| B | 2 84 1 7 13 31 61 | (31, 1, 9, 157) |

**Uniqueness.** These are the only strong rooted classes of order 7 with s = s' = 1, F = F' and h' = h + 2. There are
689 such pairs without the s = 1 condition.

**Routes.** In the graph of strong labelled tournaments with the root fixed, the two classes are at distance 4.
- There are 38 shortest routes, and every one leaves s = 1. Keeping s = 1 forces distance 7.
- Two shortest routes keep F_i <= 157 and h_i + 2F_i < 343 at every intermediate step, hence the wrappers below the
  start at every clock m >= 1.

**Wrappers.** Session189's wrapper gives H(D_m) = h + 2^m F (returns/PROOF.md (1)). Hence H(D_m(B)) = H(D_m(A)) + 2 for
every m: 343/345, 657/659, 1285/1287, ...

**Not a ternary step.** The step stays within one exponentially sparse family. It supplies no 3H + 2 / 3H + 4 parent.

**Other orders.** Unique-start rooted strong classes number F(2n - 4) (FUNNELS F2). Among them, the +-2 port pairs
number 1, 2, 5, 1, 6, 10, 4, 8 at orders 7-14.

## (B) Symmetry split

For n >= 6, among the F(2n - 9) sole-exit-free tournaments with hc = 3:
- F(2n - 10) have Aut = Z_3;
- F(2n - 11) are asymmetric.

For n >= 7, the hc = 4 and hc = 5 classes are asymmetric.

**Proof.**
- Aut(T) fixes the unique core D (THM-4869 B1) with u and l.
- It fixes the uniquely Hamiltonian rooted prefix pointwise (THM-4868 (E); for |Q| = 3, the stabiliser of a vertex of
  C3 is trivial).
- It fixes the chain of trivial components of D - u.
- So Aut(T) = Aut(W, A). This is Z_3 for the block (C3, all three vertices), and trivial for the rigid W4 and
  uniquely Hamiltonian 5-vertex blocks (THM-4869 B5).
- The C3-block cores have generating function x^5/(1 - 2x). Times the prefix factor x(1 - 2x)/(1 - 3x + x^2), this
  gives x^6/(1 - 3x + x^2) = sum F(2n - 10) x^n.

QED.

**FINITE-EXACT at orders 6-14.** The asymmetric counts are 1, 2, 5, 13, 34, 89, 233, 610, 1,597. The Z_3 counts are
1, 3, 8, 21, 55, 144, 377, 987, 2,584.
