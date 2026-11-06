---
id: HYP-9480
title: "Finite O/Z tail banks and root changes: (A) residue-resolved Poisson law - in every residue class c, a seed's tail image misses the fraction exp(-kappa lambda_c) of the odd integers, lambda_c = parameters per odd integer of the class, kappa = 1/2 for the root (converse involution) and 1 otherwise; it contains the local factors of confined seeds, the root's density 1/2 and the mod-7 bias of the one-apex image; (B) bank members are positively correlated (nesting JZ = OJ, shared halves, dilated roots), so a bank misses more than the product of its members; (C) every finite bank misses a positive proportion of odd integers, so no finite bank covers all large odd counts through its O/Z extensions; (D) root-change trichotomy - mean multiplicity bounded / ~(pi^2/6) ln N / ~ln^2 N for 0 / 1 / 2 root changes from the rooted edge; one root change leaves infinitely many misses, two leave finitely many"
status: "(A) CONFIRMED out of sample (P1: 109/109 unseen seeds within 0.004; P2: root 0.498538; P6a: root per class mod 7/9/21/63 within 0.0022 at 2^25 and 2^26; P6b: 5908 classes of 102 unseen seeds within 0.00084); (B) EXPLORATORY, not pre-registered (bank A misses 6x the product of its members' miss rates); (C) OPEN, P3 CONFIRMED (0.021894, 0.021721); (D) OPEN, P4 CONFIRMED (F_(<=2) has no odd miss in (3,780,967, 10^8]); audit PENDING"
session: kind-pasteur-2026-11-06-S13
date: 2026-11-06
related: [THM-4888, THM-4878, THM-4885, THM-4886, HYP-9474, HYP-9478, session205/finite-oz-tail-bank/ARITHMETIC_REDUCTION.md, session206/parameter-capacity/PARAMETER_CAPACITY.md, session206/weighted-parameter-capacity/WEIGHTED_CAPACITY.md, session206/root-oz-noncofiniteness/ROOT_OBSTRUCTION.md, session199/hidden-structure/PROOF.md]
script: 05-knowledge/results/finite_oz_bank_poisson_root_changes_20261106/verify.py; 05-knowledge/results/finite_oz_bank_poisson_root_changes_20261106/density.c; 05-knowledge/results/finite_oz_bank_poisson_root_changes_20261106/residue_law.c; 05-knowledge/results/finite_oz_bank_poisson_root_changes_20261106/root_changes.c; 05-knowledge/results/finite_oz_bank_poisson_root_changes_20261106/PREREGISTRATION.md; 05-knowledge/results/finite_oz_bank_poisson_root_changes_20261106/PREREGISTRATION_P6.md
---

# HYP-9480 - finite O/Z tail banks: residue-resolved Poisson law, positive correlation, root-change trichotomy

Owner's question: can a finite bank of unit-start tournaments cover all sufficiently large odd counts through its O/Z
extensions? This hypothesis records the quantitative answer (no, with a positive miss proportion) and the structure
that would be needed instead (two root changes). Proved ingredients: THM-4888; Sessions 205/206. Details, tables and
outcomes: `05-knowledge/results/finite_oz_bank_poisson_root_changes_20261106/PROOF.md`.

## (A) Residue-resolved Poisson law (single tail images)

For a seed s, a window of odd integers and an odd modulus M, let lambda_c be the number of parameters m with H_s(m) in
the window and in the class c, per odd integer of the class. Then the image misses the fraction
miss_c = exp(-kappa lambda_c) of the class, with kappa = 1/2 for the root (its converse involution pairs the
parameters) and kappa = 1 for every other seed. Summed over the classes, the intensity is 2c_s per odd integer
(Session 206's capacity), so:
- with M = 1 the law reads 1 - exp(-2c_s);
- for the root it gives 1 - exp(-ln 2) = 1/2;
- THM-4888 (C) confinement is lambda_c = 0 off the allowed classes;
- as X grows, lambda_c tends to the limit given by Session 206's weighted capacity: 2g (c_0 [c = h mod g] + c_1
  [c = h - e mod g]) with g = gcd(R, d, p), uniform on the allowed coset.

The law also accounts for the mod-7 bias of the one-apex image recorded in THM-4878 (5). At finite scales the root's
parameter intensity is non-uniform mod 7, and the law converts it, class by class, into miss rates ranging from 0.46
to 0.56. This works in every window tested, although the pattern changes from window to window. By Session 206's
weighted capacity the intensities become uniform as X grows, so the bias is a finite-size effect. It does shrink: the
mod-7 span is 0.133 at 2^20 and 0.080 at 2^26.

Exceptions are exactly the extra symmetries of THM-4888 (B) and (D):
- the root (kappa = 1/2);
- d = 1 seeds, which lie inside the root image;
- R = 2d seeds, dilated roots with a partial involution;
- balanced seeds, which are nested.

Evidence, all pre-registered: P1 and P2 (commit d36f57c903); P6a and P6b (commit e3fdff4ceb).
- P1: 109 unseen seeds, all within 0.004 of the local-factor prediction, 108 within 2 x 10^-4. The maximum is 0.0020
  at (81, 99, 27, 45), where 9 | gcd(d, R). Post hoc, the mod-9 refinement fits it, and P6b includes that refinement.
- P2: the root's density at [2^27, 2^28) is 0.498538.
- P6a: the root at [2^25, 2^26) and [2^26, 2^27); every class mod 7, 9, 21 and 63 is within 0.0022, and the mod-7
  spans are 0.090 and 0.080.
- P6b: 102 unseen generic seeds at [2^24, 2^25); all 5908 live classes mod 63 are within 0.00084, with mean error
  4.2 x 10^-5 against 0.0025 for the class-blind law.
- Systematic part: the root sits about 0.001 above the law in every class and creeps toward it (density 0.4974 in
  [2^20, 2^21), 0.4985 in [2^27, 2^28)).

## (B) Bank members are positively correlated (exploratory)

The three redundancy identities of THM-4888 (D) are:
- JZ = OJ (KNOWN, Session 205's audit), so seeds W O J are nested in W J;
- seeds that agree in (D - E, U, d), resp. (D, U, d), share their odd-m, resp. even-m, images;
- R = 2d seeds are dilated root images.

Under these identities a bank's union misses more than the product of its members' residue-resolved miss rates. In
the window [2^22, 2^23), bank A misses 0.0219 of the odd integers, while the product of its 85 members' miss rates
per class mod 3 is 0.0036. Without its 35 balanced (nested) seeds the union is unchanged and the product is
0.0157-0.0169 at every residue or position resolution. The 34 generic seeds still miss 24% more than their product.
Not pre-registered.

## (C) Finite banks miss a positive proportion

Conjecture: every finite physical bank misses a positive proportion of the odd integers. The heuristic lower bound is
the product of the members' residue-resolved miss rates; correlations (B) only add misses.

Evidence: P3, pre-registered. Bank A (the root and the 84 physical J-born seeds with d > 1 and U <= 64) misses
0.021894 at [2^26, 2^27) and 0.021721 at [2^27, 2^28). Its intensity is 2 sum c_s = 8.29, 7.60 with the root's
involution and 5.15 without its 35 nested seeds, while -ln(miss) is only 3.82. By class mod 3 the miss rates are
0.006 / 0.040 / 0.019; the class-1 deficit is THM-4888 (C).

Proved special cases (Sessions 205/206):
- capacity-deficient banks;
- d = 1-only banks;
- residue-deficient banks;
- banks failing the modularly weighted capacity test (for example the root with JJ and OZJ).

Scale: a bank whose effective intensity is Lambda should miss a proportion of order exp(-Lambda), so a large bank
can cover every odd N up to an astronomically large bound and still miss a positive proportion. S82's one-apex x
order-<=7 bank family (one free block, mean multiplicity 48.2; smallest class mean mod 3 41.4) covers every odd
N <= 2^35 except 7 and 21 (HYP-9478 P5-P6). That is consistent with (C): its first misses would be expected only far
beyond 10^15. Finite censuses cannot settle the question, only structure can.

A proof for large banks needs anti-concentration for Q along the parameter tree. Neither capacity test can suffice:
for bank A the weighted capacity exceeds every class modulo every odd p <= 105 by a factor of at least 4.7. Collision
structure cannot suffice either.

## (D) Root-change trichotomy and the two-root-change conjecture

F_k is the set of Hamiltonian-path counts of t . w_1 (J w_2) ... (J w_(k+1)): at most k root changes from the rooted
edge t.

Masses (words with H <= X):
- k = 0: (ln 2) X + O(sqrt X) (Session 206);
- k = 1: (pi^2/12) X ln X - 1.17 X + o(X), numerically for X = 10^5..10^7; the constant matches S82's THM-4885 (B);
- k = 2: about 0.50 X ln^2 X (numerical fit).

Misses:
- k = 0: positive density;
- k = 1: growing about 1.58x per dyadic window (density about N^(-1/3)), like S82's F1 (ratios 1.5876-1.5879 near
  2^(2/3) = 1.5874 at 2^32-2^34, same mass constant pi^2/12); 224,380 up to 10^8;
- k = 2: sporadic, almost all 1 (mod 3) (THM-4888 (C)). P4 (pre-registered): exactly 140 odd misses up to 10^8,
  none above 3,780,967; 14 of the 15 above 10^5 are 1 (mod 3), the exception being 136,607.

Conjecture (D):
1. One root change leaves infinitely many misses.
2. Two root changes leave finitely many: every odd N > N_0 is H(t . w_1 J w_2 J w_3).

(D)(2) plus a finite check would prove Chan-Pak. One free block can never suffice because its multiplicity is bounded
(that is (C)).
