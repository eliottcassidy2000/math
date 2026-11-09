---
id: THM-4940
title: "Odd-cycle packing counts in every base. For a tournament T and an integer b >= 2 put H_b(T) = I(Omega(T), b) = sum_k alpha_k(T) b^k, the number of b-coloured odd-cycle packings; H_2 = H by OCF. (1) H_b is multiplicative over order joins, invariant under converse, congruent to 1 mod b, and monotone under vertex deletion. (2) Block lemma: for 0 <= a < b, 1 + ab is a value iff a lies in A, the set of odd-cycle counts of tournaments with no two disjoint cycles (Session 16's intercyclic spectrum). (3) The 7-hole is base-free: 1 + 3b is never a value, for every b >= 2 (indeed 1 + 3x is never I(Omega(T), x) for real x > 0, x != 1). (4) The 21-hole is a carry obstruction: 1 + 10b is a value iff b is 4 or 8. (5) The repunit 111_b = b^2 + b + 1 is a value in base b iff b + 1 lies in A; this fails exactly for b = 2, 9, 13, 16, 19, 20, 26, ... (the base-2 failure is 7 = 111_2, and the alpha-sum-10 failure is 91 = 111_9, not 21 = 111_4, which is a base-4 value). (6) Exact gap tables at all orders for alpha-sums m <= 80: bases 2 and 3 miss only m = 3 and 10 (the values 7, 21 and 10, 31), bases 4 to 8 miss at most five, and every base b >= 81 misses exactly the 31 holes of A below 81. (7) Only the bases 2 and 3 can have gap set {3, 10}. (8) Carry-free singletons: below 1024 exactly 26 values of A have a unique strong intercyclic realiser; in base 2 each of the 22 with a >= 8 acquires a second, non-intercyclic strong realiser of H = 1 + 2a (THM-4884 below 625; explicit order-9 witnesses for 323 <= a <= 933, audit), so among them only a <= 7, i.e. H = 3, 5, 11, 13, stay binary singletons. (9) Constructively, base 3 attains every value congruent to 1 mod 3 up to 5398 except 10 and 31, already at order <= 9, so Spec_3 is exact on [1, 5398]"
status: >
  PROVED:
  - (1), (2), (3) (all orders), (5);
  - (4): b >= 11 by argument, 2 <= b <= 10 by complete finite exhaustion;
  - (7): b >= 15 by (2), 4 <= b <= 14 by the complete finite tables of (6).
  FINITE-EXACT, complete at all orders in the stated range:
  - (6) for alpha-sums m <= 80 and bases b <= 90. Moon's sharp cycle count bounds the order of every strong component.
    Independently, the audit's enumerator finds the next level, order 19, empty, so Moon's deletion lemma alone gives
    completeness.
  - (8) below 1024. By Session 16's classification, intercyclic strong tournaments with alpha_1 < 1024 have order
    <= 13. The base-2 twins come from THM-4884's census below H = 625 and from explicit order-9 witnesses for
    625 < H < 2048 (audit/aud_compare.out).
  - The audit extends (8) to 8192: 64 carry-free singletons, all self-converse, each with a >= 8 having a strong
    non-intercyclic binary twin (audit/aud_twins_hi.out).
  - (9): FINITE-EXACT, by explicit witnesses. It is exploratory, not pre-registered: it was found after
    PREREGISTRATION-2's P6-P8 failed. The audit's follow-up reproduced it independently and found the base-3 range
    already attained at order <= 9.
  KNOWN and credited:
  - OCF H = I(Omega(T), 2) (THM-002; Grinberg-Stanley);
  - the chain constraint alpha_k >= 1 => alpha_1 >= k (THM-338, THM-343);
  - three pairwise-intersecting odd cycles force a fourth (THM-343, THM-029);
  - THM-1370 (h-spectrum-omits-7-21-all-n): 21 is not an H-value;
  - Session 16's classification of intercyclic tournaments (05-knowledge/results/feedback_word_hamilton_spectrum_20261020/intercyclic_rays.md,
    audited ACCEPT), resting on Bang-Jensen-Bessy, "Small cycle transversals in tournaments", Theorem 1.8;
  - Moon's theorems: a strong tournament of order n >= 4 has a strong subtournament of order n - 1 (an (n-1)-cycle,
    by pancyclicity), and at least n - l + 1 cycles of each length l (Moon 1966; the sharp count of MISTAKE-560);
  - THM-4884 / THM-4810: the binary strong singletons below 625 are exactly 3, 5, 11, 13.
  PRIOR OBSERVATIONS credited:
  - HYP-1147 = HYP-1109 (opus-S71e, March 2026, `04-computation/ip_hierarchy.py`): the "moving gap at 1 + 3x" in
    I(Omega, x), checked exhaustively at n <= 6 and stated "for all x >= 0". (3) proves it at every order for every
    integer x >= 2 and every real x > 0 except x = 1. The stated range was too wide: at x = 0 the value 1 is attained,
    and at x = 1 the value 4 is attained (C3 => C3).
  - The same script also found I(Omega, 3) = 31 absent at n <= 6, and asked whether 10 and 31 are "gaps analogous to
    H = 7 and H = 21". That is the ternary pair of HYP-9530 (A) and of (4) at b = 3, observed at n <= 6.
  - 07-reflections/cubic-family-x6-coloring.md (oracle-2026-05-16-S2): I(Omega, 6) = 19 = 1 + 6 * 3 is never attained,
    by the same alpha = (3, 0) obstruction. This is the b = 6 case of (3).
  - 07-reflections/fugacity-axis-five-points.md (oracle-2026-05-17-S1): the degree-1 alpha_1 gaps
    {3, 10, 14, 17, 20, 21, ...}, sampled at n = 7. This is an early view of the holes of A.
  - HYP-886 (kind-pasteur-S62, x = 3) and HYP-892 (opus-S72, every x): I(Omega, x) = 1 (mod x).
  - HYP-1144 = HYP-1106 (opus-S71e): the base-2 row of (4), sampled at n = 7 and proved by THM-1370.
  - HYP-1015 (kind-pasteur-S65): x^2 + x + 1 is never an independence polynomial Omega(T), i.e. alpha = (1, 1) is
    impossible.
  - THM-467: I(Omega(T), t) = eps([det(I + XA)/det(I - XA)]^(t/2)) for even t. Its (d)1 is the strong-component
    factorisation at every x, which is the first and third items of (1).
  - THM-4754: H = 3 sum J_k alpha_k + I(Omega, -1), with base-4 repunits J.
  - 07-reflections/the-tournament-alphabet.md defines the same object H_q(T) = I(Omega(T), q) (the "q-tournament")
    and notes that it is odd iff q is even. Its readings are ANALOGY.
  - The base-2 join law a, a' -> a + a' + 2aa' on (H - 1)/2 is older than Session 210: it appears in
    paired_prime_collatz_20261003.md, tiling_modular_atoms_20261003.md and Session142.
  NEW: all-order proofs of (3)-(5) and (7); the exact second-hole set {4, 8}; the block lemma; the exact all-order
  gap tables and rescue sets to alpha-sum 80; the repunit criterion; the carry-free singleton list. Also an exhaustive
  confirmation of Session 16's intercyclic spectrum below 1024 by two generators that share one canonical-form
  routine; the audit's independent enumerator confirms it below 8192. A brief literature search found nothing earlier
  on H_b for b != 2.
  NOT CLAIMED: completeness of any base-b spectrum beyond the tables. Binary completeness is Chan-Pak (OPEN); ternary
  completeness is HYP-9530 (A). No mechanism is claimed for the projective-plane or repunit numerology recorded in the
  source note (MISTAKE-559).
  AUDIT: independent audit PASS WITH CORRECTIONS (2 majors, both on wording: the base-2 half of (8) was stated beyond
  its proof, now closed by the audit's witnesses; a README row typed "alpha_1 = 3 impossible" PROVED). There were 15
  minors. The audit reproduced every number with an independent enumerator. All corrections are applied; MISTAKE-739.
  Report: 05-knowledge/results/packing_bases_repunits_20261108/audit/AUDIT.md.
session: opus-2026-11-08-S86
source: 05-knowledge/results/packing_bases_repunits_20261108/README.md
related:
  - THM-002
  - THM-029
  - THM-338
  - THM-343
  - THM-1370 (h-spectrum-omits-7-21-all-n)
  - THM-1862
  - THM-4094
  - THM-4754
  - THM-4810
  - THM-4884
  - THM-4892
  - HYP-9530
  - 05-knowledge/results/feedback_word_hamilton_spectrum_20261020/intercyclic_rays.md
verification:
  - 05-knowledge/results/packing_bases_repunits_20261108/few_cycle_strong.c
  - 05-knowledge/results/packing_bases_repunits_20261108/base_spectra.py
  - 05-knowledge/results/packing_bases_repunits_20261108/intercyclic_gen.c
  - 05-knowledge/results/packing_bases_repunits_20261108/carry_free_layer.py
  - 05-knowledge/results/packing_bases_repunits_20261108/alpha_census.c
  - 05-knowledge/results/packing_bases_repunits_20261108/single_apex_spectrum.py
  - 05-knowledge/results/packing_bases_repunits_20261108/audit/aud_fewcycle.c
  - 05-knowledge/results/packing_bases_repunits_20261108/audit/aud_spectra.py
  - 05-knowledge/results/packing_bases_repunits_20261108/audit/aud_compare.py
  - 05-knowledge/results/packing_bases_repunits_20261108/audit/aud_twins_search.c
  - 05-knowledge/results/packing_bases_repunits_20261108/audit/run_audit.sh
  - 05-knowledge/results/packing_bases_repunits_20261108/ternary_attainment.py
  - 05-knowledge/results/packing_bases_repunits_20261108/audit/aud_followup.py
---

# THM-4940 - odd-cycle packing counts in every base

## Setting

- **Packings.** Omega(T) is the conflict graph of the directed odd cycles of T. alpha_k(T) counts the sets of k
  pairwise vertex-disjoint odd cycles, with alpha_0 = 1.
- **Packing count in base b.** For an integer b >= 1, H_b(T) = I(Omega(T), b) = sum_k alpha_k(T) b^k counts the odd-cycle
  packings with each cycle coloured in one of b colours. H_2 = H (OCF).
- **Spectrum and alpha-sums.** Spec_b = {H_b(T)}. Writing H_b(T) = 1 + b m, the integer m = sum_{k>=1} alpha_k b^(k-1)
  is the alpha-sum of T in base b. Gap_b is the set of m >= 0 with 1 + b m not in Spec_b.
- **The carry-free layer.** A = {alpha_1(T) : alpha_2(T) = 0}, the odd-cycle counts of intercyclic tournaments (no
  two disjoint cycles; equivalently no two disjoint cyclic triangles).
  - By Session 16, A = {0} u {alpha_1(w)} u {3 * 2^s : s >= 1} u {9 * 2^s - 2 : s >= 0}.
  - Here alpha_1(w) runs over the one-apex words w (a transitive tournament plus one apex). It equals
    sum over i < j with w_i = 1, w_j = 0 of c(j - i), where c(1) = 1 and c(d) = 2^(d-2).
  - The holes of A begin 3, 10, 14, 17, 20, 21, 27, 28, 29, 33, 35, ...

## Statement

1. **Basic properties.**
   - H_b(A => B) = H_b(A) H_b(B), so H_b is the product over the strong components.
   - H_b(T^op) = H_b(T).
   - H_b(T) = 1 (mod b).
   - H_b(T - v) <= H_b(T).
2. **Block lemma.** For 0 <= a < b, 1 + a b lies in Spec_b iff a lies in A. So
   Spec_b ∩ [1, b^2) = {1 + a b : a ∈ A, a < b}.
3. **The base-free 7-hole.** For every integer b >= 2, 1 + 3b is not in Spec_b.
   - More generally, for every real x > 0 with x != 1, no tournament has I(Omega(T), x) = 1 + 3x.
   - At x = 1 the value 4 = I(Omega(C3 => C3), 1) is attained, and at x = 0 the value is 1.
4. **The 21-hole is a carry obstruction.** 1 + 10b lies in Spec_b iff b is 4 or 8.
   - Witnesses:
     - b = 4: a strong 6-tournament with alpha = (6, 1), giving H_4 = 41.
     - b = 8: C3 => C3, with alpha = (2, 1), giving H_8 = 81 = 9^2.
   - In every other base, every decomposition of the alpha-sum 10 that satisfies the chain constraint is
     unrealisable. The list below omits the vectors already excluded by alpha_2 <= C(alpha_1, 2) or by
     alpha_3 >= 1 => alpha_2 >= 3: (2, 4), (6, 0, 1), (4, 1, 1) at b = 2, and (2, 2) at b = 4.
     - b = 2: (10, 0), (8, 1), (6, 2), (4, 3);
     - b = 3: (10, 0), (7, 1), (4, 2);
     - b = 5, 6, 7: (10, 0) and (10 - b, 1);
     - b >= 9: (10, 0).
5. **Repunits in all bases.** For b >= 2, 111_b = b^2 + b + 1 lies in Spec_b iff b + 1 lies in A.
   - **Failing bases.** It fails exactly for b + 1 in the holes of A: b = 2, 9, 13, 16, 19, 20, 26, 27, 28, 32, 34, 36,
     37, 40, ...
   - **Density.** By Session 16 (4), at least (1/8 - o(1)) 2^(L-1) of the integers below 2^(L-1) are holes of A, so
     the failing bases have positive upper density.
   - **Repunit primes.** Among the bases with 111_b prime, the failures are b = 2 (7), 20 (421), 27 (757), 54 (2971),
     66 (4423), 80 (6481), 90 (8191), ...
   - **Where 21 sits.** The failure tied to 21's alpha-sum 10 is 111_9 = 91. By contrast 21 = 111_4 IS a base-4 value
     (alpha = (5, 0)).
6. **Exact gap tables** (all orders; alpha-sums m <= 80):

   | b | Gap_b, m <= 80 | missing values 1 + b m |
   |---|---|---|
   | 2 | 3, 10 | 7, 21 |
   | 3 | 3, 10 | 10, 31 |
   | 4 | 3, 14, 27 | 13, 57, 109 |
   | 5 | 3, 10, 27 | 16, 51, 136 |
   | 6 | 3, 10, 14, 21, 46 | 19, 61, 85, 127, 277 |
   | 7 | 3, 10, 14, 29 | 22, 71, 99, 204 |
   | 8 | 3, 35, 52 | 25, 281, 417 |
   | >= 81 | the 31 holes of A below 81 | |

   **Rescue sets.** For each hole a <= 80 of A, the bases b in which 1 + a b is nevertheless attained:
   - a = 3: none;
   - a = 10: exactly {4, 8};
   - every other hole: always contains 2 and 3. For example a = 14 is rescued exactly in {2, 3, 5, 8, 12}.
7. **Only bases 2 and 3 can have gap set {3, 10}.** For every b >= 4, Gap_b contains an element outside {3, 10}, or
   misses 10:
   - b = 4 and b = 8 miss 10;
   - b = 5, 6, 7 contain 27, 14, 14;
   - 9 <= b <= 14 contain 14 or 17;
   - b >= 15 contains 14, by (2).
8. **Carry-free singletons.**
   - Below 1024, exactly 26 values a in A have a unique strong intercyclic realiser, all self-converse: a = 1, 2, 5, 6,
     8, 15, 32, 45, 89, 106, 176, 185, 229, 310, 323, 369, 414, 429, 544, 604, 680, 699, 746, 838, 859, 933.
   - In every base b > 933 the first-block values 1 + a b therefore have a unique strong realiser.
   - In base 2, every strong tournament with H <= 15 is intercyclic, and the binary strong singletons below 625 are
     H = 3, 5, 11, 13 (THM-4884). These are a = 1, 2, 5, 6.
   - Carries give every carry-free singleton with 8 <= a < 1024 a second strong realiser in base 2:
     - for a < 312 by THM-4884's census;
     - for the twelve values 323 <= a <= 933 by explicit strong order-9 tournaments with two disjoint cyclic
       triangles (audit).
   - The audit extends both halves to 8192: there are 64 carry-free singletons, and each of the 60 with a >= 8 has an
     explicit carried binary twin at orders 6-11.
9. **Constructive attainment** (FINITE-EXACT, explicit witnesses).
   - **Base 3.** Every value v = 1 (mod 3) with v <= 5398, other than 10 and 31, is H_3(T) for some tournament T of
     order <= 9, taken from the complete census of orders <= 9.
   - **Base 2.** Every odd v <= 2925 other than 7 and 21 is H(T) of an order join of explicit strong pieces: the strong
     classes of order <= 9 or with alpha_1 <= 80, the strong one-apex tournaments, and Session 16's rays.
     - Every odd v <= 2881 except 7 and 21 is already attained at order <= 9.
     - This half is a control, inside the known all-odd coverage to 2^35.
   - **Consequence.** With (3) and (4): Spec_3 ∩ [1, 5398] = {v <= 5398 : v = 1 (mod 3)} minus {10, 31}, i.e.
     Gap_3 ∩ [0, 1799] = {3, 10}.
   - **Where the construction stops.** It first fails at 5401 = 11 * 491. A tournament with H_3 = 5401 would have a
     single nontrivial strong component, of order >= 10, with alpha_1 > 80 and two disjoint odd cycles. Whether one
     exists is OPEN.
   - No gap is claimed beyond (3), (4) and (6).

## Proof

**(1).**
- **Join.** No odd cycle crosses an order join, so Omega(A => B) is the disjoint union of Omega(A) and Omega(B), and
  independence polynomials multiply over disjoint unions.
- **Converse.** Converse reverses every cycle, preserving Omega.
- **Mod b.** alpha_0 = 1.
- **Deletion.** Omega(T - v) is an induced subgraph of Omega(T), and every independent set of it is one of Omega(T).

**(2).**
- If H_b(T) = 1 + a b with a < b, then sum_{k>=1} alpha_k b^(k-1) = a < b. This forces alpha_k = 0 for k >= 2 and
  alpha_1 = a, so T is intercyclic and a is in A.
- Conversely, an intercyclic T with alpha_1 = a has H_b = 1 + a b.

**(3).**
- **Reduction to (3, 0).** Suppose alpha_1 + alpha_2 x + alpha_3 x^2 + ... = 3 with x > 0.
  - alpha_k >= 1 for some k >= 2 forces alpha_1 >= 2, since two disjoint odd cycles give two vertices of Omega.
    Then alpha_1 >= 2 and alpha_2 x >= x.
  - If alpha_1 = 2, then alpha_2 <= 1 and alpha_3 = 0, so 2 + alpha_2 x = 3, which needs alpha_2 = 1 and x = 1.
  - If alpha_1 >= 3, then alpha_1 = 3 and every higher alpha_k vanishes. If alpha_1 <= 1, the sum is at most 1.
  - So for x != 1 the only candidate is alpha = (3, 0).
- **(3, 0) is impossible.** It would be an intercyclic tournament with exactly three odd cycles. Three pairwise
  intersecting odd cycles lie in one strong component. On 3, 4 or at least 5 vertices that component carries 1, 2 or
  at least 4 odd cycles (THM-343).

**(4).**
- **b >= 11.** An alpha-sum of 10 < b is alpha = (10, 0), i.e. an intercyclic tournament with H = 1 + 2 * 10 = 21,
  impossible by THM-1370.
- **2 <= b <= 10.** Every realisation has alpha_1 <= 10.
  - Each strong component therefore has alpha_1 <= 10, and by Moon's sharp count (m(8) = 6 + 4 + 2 = 12 > 10) has
    order <= 7.
  - The complete inventory of strong classes with alpha_1 <= 10 (`few_cycle_strong.c`; it reproduces the order <= 7
    census) and multiplicativity (1) list every packing polynomial with alpha_1 <= 10.
  - The decomposition table in `base_spectra.out` checks every b and every decomposition. Only (6, 1) at b = 4 and
    (2, 1) at b = 8 are realised.

**(5).**
- 111_b = 1 + b(b + 1), so we need sum_{k>=1} alpha_k b^(k-1) = b + 1.
- alpha_3 = 0, since b^2 > b + 1. alpha_2 <= 1, and alpha_2 = 1 would force alpha_1 = 1, impossible by the chain
  constraint (HYP-1015).
- So alpha = (b + 1, 0), which exists iff b + 1 is in A.
- The list of failing bases is A's hole list shifted by one (Session 16's formula, confirmed below 1024 by
  `carry_free_layer.py`).

**(6).**
- **Completeness.** Let H_b(T) = 1 + b m with m <= 80. Then alpha_1(T) <= 80, so every strong component has
  alpha_1 <= 80, and so has order <= 18, because m(19) = 81.
- **The inventory.** `few_cycle_strong.c 80` lists every strong class with alpha_1 <= 80, with its full packing vector.
  - Completeness: alpha_1 does not increase under deletion, and Moon's deletion lemma gives every strong class of
    order >= 4 a strong one-vertex deletion.
  - Size: 5,324 classes of orders 3-18.
  - Checks:
    - it reproduces the full census at orders <= 9;
    - its alpha_1 <= 44 part reproduces the separate budget-44 run;
    - its intercyclic part reproduces A below 81;
    - the audit's independent enumerator (`audit/aud_fewcycle.c`, a different canonical form and cycle counter)
      finds the identical multiset of (order, alpha, H).
  - Independently, the extension step finds no strong 19-class with alpha_1 <= 80. So completeness also follows from
    Moon's deletion lemma alone, without the cycle-count floor.
- **The tables.** `base_spectra.py` forms, base by base, the multiplicative closure of the strong values H_b(C) up to
  1 + 80b, which is valid by (1). The rescue sets come from the same tables.

**(7).** This is read off (4), (6) and (2).

**(8).**
- **Completeness.** By Session 16's classification, a strong intercyclic tournament of order n has
  alpha_1 >= 2^(n - 4):
  - one-apex reduced words start with 1 and end with 0, and the first-last pairs alone give 2^(n-4);
  - the two rays are larger.
  So alpha_1 < 1024 forces n <= 13, and `intercyclic_gen.c 13` is exhaustive there by Moon's deletion lemma and
  heredity.
- **Checks.** Its counts at orders <= 9 match the census, and the set of alpha_1 values below 1024 equals Session 16's
  formula.
- **The binary statement.**
  - A strong tournament with two disjoint cycles has order >= 6. It then has alpha_1 >= 6 (Moon) and alpha_2 >= 1, so
    its binary alpha-sum is >= 8.
  - Hence for H = 1 + 2a with a <= 7, every strong realiser is intercyclic.
  - For a < 312, the rest is THM-4884's census.
  - For the twelve singletons 323 <= a <= 933, the order-9 census contains a strong tournament with H = 1 + 2a and two
    disjoint cyclic triangles (`audit/aud_compare.out` section 4).
  - This completes (8) below 1024. Beyond it, the audit's extension to 8192 finds a twin for each of the 38 singletons
    in [1024, 8192), by randomised search at order 10 (order 11 for a = 7845). Each twin is re-verified in C and in
    Python: H by subset DP, strong, two disjoint cyclic triangles (`audit/aud_twins_hi.out`).

**(9).**
- **Attainment.**
  - `ternary_attainment.py` forms the multiplicative closure of the piece values, which is valid by (1). It prints
    factorisations as spot checks.
  - The audit's follow-up (`audit/aud_followup.py`) evaluates H_3 directly on all 198,948 classes of orders <= 9 and
    finds every value = 1 (mod 3) up to 5398 except 10 and 31. In base 2 the direct evaluation reaches 2881, and the
    products reach 2925.
- **The two gaps** 10 and 31 are (3) and (4) at b = 3.

QED.
