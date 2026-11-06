---
id: HYP-9480
title: "Finite O/Z tail banks and root changes: (A) residue-resolved Poisson law - in every residue class c, a seed's tail image misses the fraction exp(-kappa lambda_c) of the odd integers, lambda_c = parameters per odd integer of the class, kappa = 1/2 for the root (converse involution) and kappa = 1 for generic seeds (J-born, d > 1, R != 2d, not nested in a seed with d = 1 or R = 2d); it contains the local factors of confined seeds, the root's density 1/2 and the slowly decaying mod-7 bias of the one-apex image; (B) bank members are correlated exactly through the redundancy identities of THM-4888 (D) (nesting JZ = OJ, shared halves, dilated roots), so a bank misses more than the product of its members; (C) every finite bank misses a positive proportion of odd integers, so no finite bank covers all large odd counts through its O/Z extensions; (D) root-change trichotomy - mean multiplicity bounded / ~(pi^2/6) ln N / ~ln^2 N for 0 / 1 / 2 root changes from the rooted edge; one root change leaves infinitely many misses, two leave finitely many"
status: "(A) CONFIRMED out of sample: for the root (kappa = 1/2) P2 (0.498538) and P6a (every class mod 7/9/21/63 within 0.0022 at 2^25 and 2^26) are decisive, a no-collision model missing by more than 0.2; for generic seeds (kappa = 1) P1 (109/109 within 0.004) and P6b (all 5908 classes mod 63 of 102 generic P1 seeds within 0.00084) are CONFIRMED but discriminate weakly - the no-collision model miss_c = 1 - lambda_c also passes P1 and fails P6b only on 16 high-intensity classes, and generic seeds never exceed per-class intensity 0.2425; the Poisson form's mean error 4.2e-5 against 4.3e-4 is reported, not pre-registered. (B) EXPLORATORY, not pre-registered: no correlation beyond the identities of THM-4888 (D) (audit). (C) OPEN, P3 CONFIRMED (0.021894, 0.021721). (D) OPEN, P4 CONFIRMED (F_(<=2) has no odd miss in (3,780,967, 10^8]). Audit PASS WITH CORRECTIONS (4 majors, 12 minors; MISTAKE-727)"
session: kind-pasteur-2026-11-06-S13
date: 2026-11-06
related: [THM-4888, THM-4878, THM-4885, THM-4886, HYP-9474, HYP-9478, session205/finite-oz-tail-bank/ARITHMETIC_REDUCTION.md, session206/parameter-capacity/PARAMETER_CAPACITY.md, session206/weighted-parameter-capacity/WEIGHTED_CAPACITY.md, session206/root-residue-capacity/ROOT_RESIDUE_CAPACITY.md, session206/root-oz-noncofiniteness/ROOT_OBSTRUCTION.md, session199/hidden-structure/PROOF.md]
script: 05-knowledge/results/finite_oz_bank_poisson_root_changes_20261106/verify.py; 05-knowledge/results/finite_oz_bank_poisson_root_changes_20261106/density.c; 05-knowledge/results/finite_oz_bank_poisson_root_changes_20261106/residue_law.c; 05-knowledge/results/finite_oz_bank_poisson_root_changes_20261106/root_changes.c; 05-knowledge/results/finite_oz_bank_poisson_root_changes_20261106/PREREGISTRATION.md; 05-knowledge/results/finite_oz_bank_poisson_root_changes_20261106/PREREGISTRATION_P6.md
---

# HYP-9480 - finite O/Z tail banks: residue-resolved Poisson law, correlation through identities, root-change trichotomy

Owner's question: can a finite bank of unit-start tournaments cover all sufficiently large odd counts through its O/Z
extensions? This hypothesis records the conjectural answer (no: a positive miss proportion) and the structure that
would be needed instead (two root changes). Proved ingredients: THM-4888; Sessions 205/206. Details, tables and
outcomes: `05-knowledge/results/finite_oz_bank_poisson_root_changes_20261106/PROOF.md`. Audit:
`.../audit/AUDIT.md` (PASS WITH CORRECTIONS).

## (A) Residue-resolved Poisson law (single tail images)

For a seed s, a window of odd integers and an odd modulus M, let lambda_c be the number of parameters m with H_s(m) in
the window and in the class c, per odd integer of the class. Then the image misses the fraction
miss_c = exp(-kappa lambda_c) of the class:
- kappa = 1/2 for the root, whose converse involution pairs the parameters;
- kappa = 1 for generic seeds. These are the J-born seeds with d > 1 and R != 2d that are not balanced seeds W O J
  whose parent W J has d = 1 or R = 2d. All P1 and P6b seeds are generic.

Summed over the classes, the intensity is 2c_s per odd integer (Session 206's capacity), so:
- with M = 1 the law reads 1 - exp(-2c_s);
- for the root it gives 1 - exp(-ln 2) = 1/2;
- THM-4888 (C) confinement is lambda_c = 0 off the allowed classes;
- as X grows, lambda_c tends to the limit given by Session 206's weighted capacity: 2g (c_0 [c = h mod g] + c_1
  [c = h - e mod g]) with g = gcd(R, d, p), uniform on the allowed coset.

The law also accounts for the mod-7 bias of the one-apex image recorded in THM-4878 (5).
- At the tested scales the root's parameter intensity is non-uniform mod 7. The law converts it, class by class,
  into miss rates from 0.46 to 0.56.
- This works in every window tested, although the pattern changes from window to window.
- By Session 206's root-residue-capacity theorem the intensities become uniform as X grows, so the bias is a
  finite-size effect. It does shrink: the mod-7 span is 0.133 at 2^20, 0.080 at 2^26 and 0.074 at 2^27.

Known exceptions to kappa = 1, each from THM-4888 (B) or (D); the list is not shown to be complete. At [2^22, 2^23)
mod 63, these seeds miss more than exp(-lambda_c), by 0.02 to 0.14 per class:
- d = 1 seeds, which lie inside the root image;
- R = 2d seeds, dilated roots with a partial involution;
- O/Z descendants of R = 2d seeds, for example (5, 12, 0, 3) = O(JJ);
- balanced seeds W O J whose parent W J is exceptional, which inherit its collisions.

A balanced seed nested in a generic seed follows kappa = 1: (53, 53, 17, 17), nested in (27, 35, 9, 17), deviates by
1e-4.

Evidence, all pre-registered: P1 and P2 (commit d36f57c903); P6a and P6b (commit e3fdff4ceb).
- P2 and P6a (the root) are decisive.
  - P2: the root's density at [2^27, 2^28) is 0.498538.
  - P6a: at [2^25, 2^26) and [2^26, 2^27), every class mod 7, 9, 21 and 63 is within 0.0022, and the mod-7 spans are
    0.090 and 0.080.
  - A no-collision model is off by more than 0.2.
- P1: all 109 unseen generic seeds lie within 0.004 of the local-factor prediction, and 108 within 2 x 10^-4.
  - The maximum, 0.0020, is at (81, 99, 27, 45), where 9 | gcd(d, R). Post hoc, the mod-9 refinement fits it, and P6b
    includes that refinement.
  - P1 tests the intensities and the confinement, not the Poisson form: the no-collision model also passes (maximum
    0.0031).
- P6b: 102 generic P1 seeds at [2^24, 2^25). They were not used in the P6 exploration, but their overall densities in
  this window were known from P1.
  - All 5908 live classes mod 63 are within 0.00084.
  - Mean error 4.2 x 10^-5; the no-collision model gives 4.3 x 10^-4 and the class-blind law 0.0025.
  - The no-collision model fails the registered criterion only on 16 high-intensity classes.
- No sharper generic test exists: among the 18,080 generic J-born seeds with U <= 4096, the largest limiting per-class
  intensity mod 63 is 0.2425, reached at (9, 11, 3, 5) and its multiples. There the Poisson and no-collision
  predictions differ by 0.029.
- Systematic part: in every class mod 7 and mod 9 the root sits 0.0005 to 0.002 above the law and creeps toward it.
  Its density is 0.4974 in [2^20, 2^21) and 0.4985 in [2^27, 2^28).

## (B) Bank members are correlated exactly through the identities of THM-4888 (D) (exploratory)

The three redundancy identities of THM-4888 (D) are:
- JZ = OJ (KNOWN, Session 205), so seeds W O J are nested in W J;
- seeds that agree in (D - E, U, d), resp. (D, U, d), share their odd-m, resp. even-m, images;
- R = 2d seeds are dilated root images.

Under these identities a bank's union misses more than the product of its members' residue-resolved miss rates. In
the window [2^22, 2^23), bank A misses 0.0219 of the odd integers, while the product of its 85 members' miss rates
per class mod 3 is 0.0036. Without its 35 balanced (nested) seeds the union is unchanged and the product is
0.0157-0.0169 at every residue or position resolution.

The rest is accounted for (audit):
- Once their six shared-half pairs (D)2 are merged, the 34 generic seeds behave as independent half-images: 0.0908
  against their union miss of 0.0909.
- The generic and dilated-root blocks are independent of each other, class by class.
- The remaining gap lies inside the dilated-root block (D)3.

Not pre-registered.

## (C) Finite banks miss a positive proportion

Conjecture: every finite physical bank misses a positive proportion of the odd integers. The heuristic lower bound is
the product of the members' residue-resolved miss rates, and the identities (B) only add misses.

Evidence: P3, pre-registered. Bank A (the root and the 84 physical J-born seeds with d > 1 and U <= 64) misses
0.021894 at [2^26, 2^27) and 0.021721 at [2^27, 2^28).
- Its intensity 2 sum c_s is 8.29; it is 7.60 with the root's involution and 5.15 without its 35 nested seeds.
  -ln(miss) is only 3.82.
- By class mod 3 the miss rates are 0.006 / 0.040 / 0.019. The class-1 deficit comes mostly from the confinement of
  THM-4888 (C) (heuristic).

Proved special cases (Sessions 205/206):
- capacity-deficient banks;
- d = 1-only banks;
- residue-deficient banks;
- banks failing the modularly weighted capacity test (for example the root with JJ and OZJ).

Scale: under a Poisson heuristic, a bank of effective intensity Lambda misses a proportion of order exp(-Lambda).
- S82's one-apex x order-<=7 bank family has one free block and mean multiplicity 48.2; its smallest class mean mod 3
  is 41.4. It covers every odd N <= 2^35 except 7 and 21 (HYP-9478 P5-P6).
- A Poisson model would place its first misses anywhere from about 6 x 10^9 (with bank A's halving of the intensity)
  to about 6 x 10^18 (no loss).
- S82's F2 is far from Poisson, so this census does not test (C).
- It shows only that a positive miss proportion can lie beyond any feasible census.

A proof for large banks needs anti-concentration for Q along the parameter tree. Neither capacity test can suffice:
for bank A the weighted capacity exceeds every class modulo every odd p <= 105 by a factor of at least 4.7 (3.39
without the nested seeds). Collision arguments are limited too:
- (B) of THM-4888 excludes the full root involution for every other seed;
- but R = 2d seeds keep a partial involution, d = 1 seeds lie inside the root image, and other collision structures
  are not excluded.

## (D) Root-change trichotomy and the two-root-change conjecture

F_k is the set of Hamiltonian-path counts of t . w_1 (J w_2) ... (J w_(k+1)): at most k root changes from the rooted
edge t.

Masses N_k(X) count the words with exactly k root changes and H <= X. Trailing O letters of the final block are
deleted; intermediate blocks keep all their letters.
- k = 0: (ln 2) X + O(sqrt X) (Session 206).
- k = 1: (pi^2/12) X ln X - 1.17 X + o(X), numerically for X = 10^5..10^8.
  - The constant pi^2/12 is THM-4885 (B)'s. Its unordered pairs have mass (pi^2/24) X ln X, so F_1 has twice that.
  - For "at most one" root change, the second constant is -0.48.
- k = 2: about 0.500 X ln^2 X - 2.03 X ln X + 4.13 X (numerical fit). The mean multiplicity is about 204 at
  N = 4 x 10^6.

Misses:
- k = 0: positive density.
- k = 1: growing about 1.58x per dyadic window (density about N^(-1/3)); 224,380 up to 10^8. This is like S82's F1,
  whose ordered pairs have the same mass constant: ratios 1.5876-1.5879, near 2^(2/3) = 1.5874, at 2^32-2^34.
- k = 2: sporadic, almost all 1 (mod 3), where the confinement of THM-4888 (C) removes mass (heuristic).
  - P4 (pre-registered): exactly 140 odd misses up to 10^8, none above 3,780,967.
  - 14 of the 15 above 10^5 are 1 (mod 3); the exception is 136,607.

Conjecture (D):
1. One root change leaves infinitely many misses.
2. Two root changes leave finitely many: every odd N > N_0 is H(t . w_1 J w_2 J w_3).

An effective form of (D)(2), with an explicit N_0, plus a finite check below N_0 would prove Chan-Pak. Conjecturally
(C), one free block never suffices: its mean multiplicity is bounded (Session 206).
