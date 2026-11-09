---
id: HYP-9640
title: "Chan-Pak in both bases as cofiniteness of an explicit matrix-semigroup orbit, and one-defect digraph spectra. (A) The 15 explicit base-3 port-chain generators (nonnegative 3 x 3 integer matrices) carry (1,0,1) to states with every first coordinate v = 1 (mod 3), v >= 866; with the finite censuses this implies ternary Chan-Pak. (B) The 56 base-2 port-chain matrices of blocks of order <= 4 reach every odd v >= 493; this implies binary Chan-Pak. (C) With exactly one reciprocal pair and no missing pair, the digraph path totals are exactly the integers >= 2 other than 4 and 7. (D) With exactly one missing pair and no reciprocal pair, every nonnegative integer is a total. (E) For bases 4-8 the gap set is exactly its part below 81 (verified below alpha-sum 30,000)."
status: >
  OPEN.
  (A): FINITE-EXACT to 10^6 (THM-5070 (4)).
  (B): FINITE-EXACT to 2*10^5.
  (C): the exclusions are PROVED at every order (THM-5071 (4b)); attainment FINITE-EXACT below 198 (orders <= 7).
  (D): FINITE-EXACT below 130 (orders <= 7).
  (E): FINITE-EXACT below alpha-sum 30,000 (THM-5070 (6)).
session: kind-pasteur-2026-11-09-S18
source: 05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/README.md
related:
  - THM-5070
  - THM-5071
  - HYP-9530
  - THM-4950
  - THM-4952
  - THM-4875
  - THM-1370 (h-spectrum-omits-7-21-all-n)
---

# HYP-9640 - Chan-Pak as a matrix-semigroup orbit, and one-defect spectra

## (A) Ternary

**Statement.** Let G_3 be the 15 matrices listed in
05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/out/generators_base3_15.txt. Each is a port-chain
block of order <= 4 (THM-5070 (1)) with x = 3. Every v = 1 (mod 3) with v >= 866 is the first coordinate of
M_r ... M_1 (1, 0, 1) for some M_i in G_3.

**Why it implies ternary Chan-Pak.**
- Below 866, the values other than 10 and 31 are attained by known tournaments: already by order-5 port chains
  (THM-5070 (4)), and earlier by THM-4940 (9) and THM-4950 (3).
- 10 and 31 are never values (THM-4940 (3), (4)).

**Evidence.** The statement holds on [866, 10^6], with 64,420,757 reachable states.

**Falsifier.** A value v = 1 (mod 3), v >= 866, outside the orbit.

**Route.** THM-5070 (5) excludes an affine digit system along the tower. A proof needs a second mechanism, for
instance a set of reachable states, stable under a few blocks, whose first coordinates fill intervals.

## (B) Binary

**Statement.** The 56 base-2 block matrices of order <= 4 reach every odd v >= 493. The 13-matrix subset
`generators_base2_13.txt` reaches every odd v >= 961.

**Why it implies Chan-Pak.** The odd values below 493 other than 7 and 21 are attained (exhaustive censuses;
THM-4875 (D)).

**Evidence.** Holds to 2*10^5. THM-4875 (D) is a different family that reaches 2^32.

## (C) One reciprocal pair

**Statement.** The totals of loopless semicomplete digraphs with exactly one reciprocal pair are exactly {2, 3, 5, 6}
together with every integer >= 8.

**Status.** Exclusion of 1, 4 and 7 is PROVED at all orders (THM-5071 (4b)). Attainment is checked below 198.

**Falsifier.** An integer >= 8 never attained.

## (E) The small gap sets are all the gaps (bases 4-8)

**Statement.** For b = 4, 5, 6, 7, 8, Gap_b is exactly the set of its elements below 81:
- b = 4: {3, 14, 27};
- b = 5: {3, 10, 27};
- b = 6: {3, 10, 14, 21, 46};
- b = 7: {3, 10, 14, 29};
- b = 8: {3, 35, 52}.

This sharpens HYP-9530 (B) for these bases.

**Evidence.** It holds below alpha-sum 30,000 (THM-5070 (6)).

**Falsifier.** Any alpha-sum m > 30,000 with 1 + b m attained by no tournament.

## (D) One missing pair

**Statement.** Every integer h >= 0 is the total of a tournament with one arc deleted.

**Evidence.** Checked below 130 at orders <= 7.

**Falsifier.** A missing value.
