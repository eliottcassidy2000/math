---
id: THM-5312
title: "Gluing targets and one shared vertex. (A) If induced witnesses of odd primes p, q (not both 3) share at most one vertex, every tournament containing them has H >= p + q + 7. (B) Overlap inclusion-exclusion: H >= p + q - r + (a-1)(q-r) + (b-1)(p-r) - r(a-1)(b-1), correcting the intersection term of CONDITIONAL.md. (C) Target map: for odd 23 <= N <= 123 and rows p + q - c, c in {1, -1, -3}, the census of orders <= 9 decides which prime pairs sit as induced subtournaments in a count-N tournament; the (J) row admits a pair at only 18 targets, and 27, 57, 67, 71, 91, 103, 105, 115 admit none in any row. (D) The conflict graphs of value 7 and 21: K_3 only for 7; six for 21 (shield of K_3, P_4, four joins, three of them 11 v 11)"
status: >
  PROVED: (A) for all orders (case analysis on the strong cores, Moon's vertex-pancyclicity, and four explicit mixed
  odd cycles); (B) (pack inclusion-exclusion; the intersection of the two classes carries every pack of W); (D) (finite
  enumeration of complements with at most four edges on at most ten vertices, omega_forms.py).
  FINITE-EXACT + CITED: (C), all tournaments of orders 3-9 (gentourng). Completeness uses Moon-Busch f(10) = 125 (Busch,
  EJC 2006, CITED) and THM-1370 monotonicity. Independent OCF replay for N < 45.
  FINITE-EXACT: replay of (A) and (B) on every cover at orders <= 9. The printed CONDITIONAL.md bound fails 135 times at
  order 9 (MISTAKE-775).
  NOT CLAIMED: any physical gluing; HYP-9862 (one-vertex supermultiplicativity, OPEN); Chan-Pak (OPEN).
session: klein-2026-11-16-join
date: 2026-11-16
source: 05-knowledge/results/chanpak_gluing_targets_one_vertex_20261116/README.md
scripts:
  - 05-knowledge/results/chanpak_gluing_targets_one_vertex_20261116/jglue.c, jcontain.c, onevx.c, supmult.c
  - 05-knowledge/results/chanpak_gluing_targets_one_vertex_20261116/verify_rows_ocf.py, omega_forms.py, make_tables.py
related:
  - CONDITIONAL.md (conditional_determinant_vitali_20261116; one-vertex bound p+q+3 sharpened, quantitative bound corrected)
  - fresh_unit_prime_gluing_20261116/PRIME_GLUE.md (no overlap <= 3 for primes >= 11; deletion-form families), SPINE_AND_HALF_INDEX.md
  - chanpak_waves_20261113 wave 43 (two-apex gluing only (3,3); {7,21} = {h(3), h(3 o 1)})
  - THM-1370 (7, 21 excluded; strong floor), THM-4940 (alpha-sums, 21-hole carry obstruction), THM-4114 (Mobius positivity), THM-202 (P_4), THM-343
---

# THM-5312 - gluing targets and one shared vertex

**(A) One shared vertex.** Let p, q be odd primes, not both 3. Let A and B be induced subtournaments
with H(A) = p, H(B) = q and |A ∩ B| ≤ 1. Then every tournament containing them has H ≥ p + q + 7.
- Equality is attained by (3, 5).
- For (3, 3) the bound is 9 = p + q + 3.
- Proof: [ONE_VERTEX.md](../../05-knowledge/results/chanpak_gluing_targets_one_vertex_20261116/ONE_VERTEX.md).
  - Disjoint cores give pq ≥ p + q + 7.
  - A core cycle avoiding the shared vertex adds a mixed pair to the bowtie.
  - Otherwise both cores are one-apex at the shared vertex, and four distinct mixed odd cycles exist.

**(B) Overlap inclusion-exclusion.** Write W = A ∩ B, r = H(W), a = H(A − W) and b = H(B − W). Then

    H(A ∪ B) ≥ p + q − r + (a−1)(q−r) + (b−1)(p−r) − r(a−1)(b−1).

- The printed version in CONDITIONAL.md has (a−1)(b−1) in the last term. It is false: at order 9 there
  is H = 99 < 107 with p = 33, q = 9, r = a = b = 3.

**(C) Target map.** Let N be odd with 23 ≤ N ≤ 123, c ∈ {1, −1, −3}, and p ≤ q odd primes ≠ 7 with
p + q = N + c. Whether some tournament of count N contains induced subtournaments of counts p and q is
decided by the census of orders ≤ 9.
- In the (J) row (c = 1) exactly 18 targets admit a pair:
  25, 33, 35, 37, 39, 45, 51, 53, 59, 61, 69, 73, 79, 81, 95, 117, 119, 121.
- Across all rows, 27, 57, 67, 71, 91, 103, 105 and 115 admit none.
- 11 is never a (J) leg in this range.
- Full table:
  [TARGET_MAP.md](../../05-knowledge/results/chanpak_gluing_targets_one_vertex_20261116/TARGET_MAP.md).

**(D) Conflict graphs of value 7 and 21.** I(G, 2) = 7 only for K_3. I(G, 2) = 21 for exactly six
graphs:
- K_10, K_8 − e, K_6 minus a 2-matching, and K_3 ∨ (K_2 ⊔ K_1) are joins; every one splits as
  3 ∨ 19 and 5 ∨ 17, and the first three also as 11 ∨ 11;
- K_3 ⊔ K_1 is the triangle shield of 7;
- P_4 is neither.

K_3 is never a connected component of Ω(T).
