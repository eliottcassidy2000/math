# Which prime pairs can sit inside a tournament of count N = p + q − c

**Question.** CONDITIONAL.md's overlapping prime-gluing hypothesis asks for a tournament with count
p + q − c, built from witnesses of the odd primes p and q (both ≠ 7), for each row c ∈ {1, −1, −3}.
- c = 1 is the gluing (J): count p + q − 1, a conflict-graph join.
- c = −1 and c = −3 give p + q + 1 and p + q + 3.

The weakest version that keeps the inputs is **induced retention**. The output must contain an
induced subtournament with count p and one with count q. Nothing is assumed about how the two overlap
or whether they cover the output.

This note decides induced retention completely for every target N ≤ 123.

## Completeness lemma (PROVED, using Busch's theorem)

Let p ≤ q be odd primes with p + q = N + c, N ≥ 23 and c ∈ {1, −1, −3}. Suppose a tournament T with
H(T) = N contains induced X and Y with H(X) = p and H(Y) = q. Then some tournament of order at most 9
with count N contains such a pair, provided N < 125.

*Proof.*
- **Cores.** H is the product of H over the strong components (THM-1370 (1)). Since p is prime, X has
  exactly one nontrivial strong component K_X, and H(K_X) = p. Likewise K_Y with H(K_Y) = q.
- **One component.** K_X is strongly connected, so it lies inside one strong component S of T. Likewise
  K_Y lies inside a strong component S′.
  - If S ≠ S′, then by monotonicity under vertex deletion H(T) ≥ H(S)H(S′) ≥ pq.
  - But pq > p + q + 3 ≥ N unless p = q = 3, which needs N ≤ 9.
- **The component carries all of N.** H(S) divides N and H(S) ≥ q ≥ (N − 3)/2 > N/3, so H(S) = N.
- **The component is small.** The strong floor f(n), the least count of a strong n-tournament, is
  nondecreasing for n ≥ 5 (THM-1370 (3)). Moon–Busch gives f(n) = min{3^a 5^b : 2a + 3b = n − 1}
  (Busch, EJC 2006, CITED), so f(10) = 125. Hence N < 125 forces |S| ≤ 9, and S itself is the required
  tournament. ∎

For N < 75 the cited theorem is not needed. The order-9 census gives f(9) = 75 directly, and THM-1370's
monotonicity extends it to all larger orders.

## Census (FINITE-EXACT)

`jcontain.c` computes H of every induced subtournament of every tournament of orders 3–9 (all 191,536
classes at order 9, from `gentourng`). In `ROWS=1` mode it lists every admissible prime pair of every
row that sits in the tournament. Output: `rows_census.out`, 338 distinct (N, p, q, c) records.

Without `ROWS`, it classifies (J)-row pairs by the strongest form present (`jrow_forms.out`):
- C: contained anywhere;
- K: the two subtournaments cover the host;
- D: vertex exchange, A = V − a and B = V − b.

`verify_rows_ocf.py` is an independent check for every N < 45. It does not use the path DP: it
enumerates the directed odd cycles and evaluates H = I(Ω, 2) (THM-002). Its table agrees with the C
census exactly (`verify_rows_ocf.out`).

## Result

1. **(J) row.** Only 18 of the 51 odd targets in [23, 123] admit any induced pair:
   25, 33, 35, 37, 39, 45, 51, 53, 59, 61, 69, 73, 79, 81, 95, 117, 119, 121.
   - The other 33 admit none: 23, 27, 29, 31, 41, 43, 47, 49, 55, 57, 63, 65, 67, 71, 75, 77, 83, 85,
     87, 89, 91, 93, 97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 123.
   - Allowing the small legs 3 and 5 changes nothing in this row.
2. **The leg 11.** The prime 11 never sits as a (J) leg below 125, with any partner. 13 sits only in
   13 + 13 − 1 = 25.
   - (11, 11) → 21 is the forced failure noted in fresh_unit_prime_gluing_20261116/SPINE_AND_HALF_INDEX.md.
   - 11 and 13 are the binary strong singletons (THM-4940 (8)).
3. **All three rows.** Even with all three rows and all primes ≠ 7, eight targets admit no induced pair
   at all: **27, 57, 67, 71, 91, 103, 105, 115.**
   - So induced retention fails for CONDITIONAL.md's hypothesis at these counts, whichever row the
     reduction picks. At N = 67 its rule chooses c = −3, and none of the pairs of 64 (3+61, 5+59, 11+53,
     17+47, 23+41) sits in 67.
4. **Shape of the successes.**
   - The (J) successes are near-diagonal: q/p ≤ 1.6 for the covers, at most 1.83 overall.
   - Every pair that sits as a cover also sits as a vertex exchange (A = V − a, B = V − b), at the same
     least order.
   - This agrees with the exact deletion-form hosts of fresh_unit_prime_gluing_20261116/PRIME_GLUE.md,
     credited: 13,13→25 at order 6; 17,19→35, 23,17→39, 29,23→51 and 31,23→53 at order 7.
   - Two pairs sit without forming a cover: 17,29→45 (order 8) and 29,53→81 (order 9).
5. **63 = 7·9.** 63 is a count (THM-1370, order 8), yet no Goldbach pair of 64 sits in it.

**What this does and does not say.**
- It does not forbid the counts: every target above is a Hamilton count (THM-1370 covers [1, 609] except
  7 and 21).
- It forbids keeping prime witnesses as induced subtournaments. A gluing rule that works at these targets
  must rebuild ("prepare") at least one witness so that the input is no longer present in the output.
- PRIME_GLUE.md proves that an overlap of at most three vertices fails for primes ≥ 11. This map is the
  target-side complement: below 125, no overlap of any size helps at the 33 (resp. 8) listed targets.

## Table (generated by `make_tables.py` → `target_map_table.out`)

The bracket gives the forms present (C, K, D) and the least order at which the pair appears.

| N | Goldbach pairs p+q=N+1 (p,q>=11) | induced (J) pairs [forms, least order] | rows -1 / -3 realised (any primes != 7) |
|--:|:--|:--|:--|
| 23 | (11,13) | **none** | 11+11+1 |
| 25 | (13,13) | (13,13) [CDK 6] | 11+13+1, 11+11+3 |
| 27 | (11,17) | **none** | **none** |
| 29 | (11,19) (13,17) | **none** | 13+13+3 |
| 31 | (13,19) | **none** | 11+19+1 |
| 33 | (11,23) (17,17) | (17,17) [CDK 7] | 11+19+3 |
| 35 | (13,23) (17,19) | (17,19) [CDK 7] | 17+17+1 |
| 37 | (19,19) | (19,19) [CDK 7] | 11+23+3 |
| 39 | (11,29) (17,23) | (17,23) [CDK 7] | 19+19+1, 17+19+3 |
| 41 | (11,31) (13,29) (19,23) | **none** | 17+23+1 |
| 43 | (13,31) | **none** | 19+23+1 |
| 45 | (17,29) (23,23) | (17,29) [C 8], (23,23) [CDK 7] | 19+23+3 |
| 47 | (11,37) (17,31) (19,29) | **none** | 23+23+1 |
| 49 | (13,37) (19,31) | **none** | 19+29+1, 17+29+3, 23+23+3 |
| 51 | (11,41) (23,29) | (23,29) [CDK 7] | 17+31+3, 19+29+3 |
| 53 | (11,43) (13,41) (17,37) (23,31) | (23,31) [CDK 7] | 19+31+3 |
| 55 | (13,43) (19,37) | **none** | 23+31+1 |
| 57 | (11,47) (17,41) (29,29) | **none** | **none** |
| 59 | (13,47) (17,43) (19,41) (23,37) (29,31) | (29,31) [CDK 8] | 29+29+1 |
| 61 | (19,43) (31,31) | (31,31) [CDK 8] | 23+37+1 |
| 63 | (11,53) (17,47) (23,41) | **none** | 23+37+3 |
| 65 | (13,53) (19,47) (23,43) (29,37) | **none** | 31+31+3 |
| 67 | (31,37) | **none** | **none** |
| 69 | (11,59) (17,53) (23,47) (29,41) | (29,41) [CDK 8] | **none** |
| 71 | (11,61) (13,59) (19,53) (29,43) (31,41) | **none** | **none** |
| 73 | (13,61) (31,43) (37,37) | (37,37) [CDK 7] | 31+41+1, 23+47+3, 29+41+3 |
| 75 | (17,59) (23,53) (29,47) | **none** | 37+37+1 |
| 77 | (11,67) (17,61) (19,59) (31,47) (37,41) | **none** | 37+37+3 |
| 79 | (13,67) (19,61) (37,43) | (37,43) [CDK 8] | **none** |
| 81 | (11,71) (23,59) (29,53) (41,41) | (29,53) [C 9], (41,41) [CDK 8] | **none** |
| 83 | (11,73) (13,71) (17,67) (23,61) (31,53) (37,47) (41,43) | **none** | 41+41+1 |
| 85 | (13,73) (19,67) (43,43) | **none** | 31+53+1, 29+53+3, 41+41+3 |
| 87 | (17,71) (29,59) (41,47) | **none** | 31+53+3, 37+47+3 |
| 89 | (11,79) (17,73) (19,71) (23,67) (29,61) (31,59) (37,53) (43,47) | **none** | 43+43+3 |
| 91 | (13,79) (19,73) (31,61) | **none** | **none** |
| 93 | (11,83) (23,71) (41,53) (47,47) | **none** | 31+59+3 |
| 95 | (13,83) (17,79) (23,73) (29,67) (37,59) (43,53) | (37,59) [CDK 8] | 41+53+1, 47+47+1 |
| 97 | (19,79) (31,67) (37,61) | **none** | 43+53+1, 47+47+3 |
| 99 | (11,89) (17,83) (29,71) (41,59) (47,53) | **none** | 37+59+3 |
| 101 | (13,89) (19,83) (23,79) (29,73) (31,71) (41,61) (43,59) | **none** | 41+59+1 |
| 103 | (31,73) (37,67) (43,61) | **none** | **none** |
| 105 | (17,89) (23,83) (47,59) (53,53) | **none** | **none** |
| 107 | (11,97) (19,89) (29,79) (37,71) (41,67) (47,61) | **none** | 53+53+1 |
| 109 | (13,97) (31,79) (37,73) (43,67) | **none** | 53+53+3 |
| 111 | (11,101) (23,89) (29,83) (41,71) (53,59) | **none** | 37+71+3 |
| 113 | (11,103) (13,101) (17,97) (31,83) (41,73) (43,71) (47,67) (53,61) | **none** | 43+67+3 |
| 115 | (13,103) (19,97) (37,79) (43,73) | **none** | **none** |
| 117 | (11,107) (17,101) (29,89) (47,71) (59,59) | (59,59) [CDK 8] | 53+61+3 |
| 119 | (11,109) (13,107) (17,103) (19,101) (23,97) (31,89) (37,83) (41,79) (47,73) (53,67) (59,61) | (59,61) [CDK 9] | 59+59+1 |
| 121 | (13,109) (19,103) (43,79) (61,61) | (61,61) [CDK 9] | 47+73+1, 59+59+3 |
| 123 | (11,113) (17,107) (23,101) (41,83) (53,71) | **none** | 61+61+1, 47+73+3 |
