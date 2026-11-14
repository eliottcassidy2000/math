# Chan–Pak waves: prediction, test, verdict (session mac-mini-2026-11-13-cuts, sixth prompt)

**Chan–Pak OPEN.** Each wave states predictions, runs the cheapest decisive test, and records the verdict. Data:
complete strong spectra Sp(n) of orders ≤ 9 (`gentourng` + `hl`, see the busch_staircase note). Notation: in base
b, I(Ω(T), b) = 1 + b·m with m = c1 + b·c2 + b²·c3 + …, where c_k counts sets of k pairwise vertex-disjoint odd
cycles (THM-5070). Gap_b is the set of unattained m.

## Wave 1: pairwise-meeting odd cycles

- **Prediction P11.** Every attained H is attained by a tournament whose odd cycles pairwise meet (Ω complete).
- **Verdict: REFUTED.** Below 125, the odd counts with no Ω-complete realizer are exactly:

  7, 21, 29, 35, 41, 43, 55, 57, 59, 67, 71, 75, 77, 83, 89, 93, 101, 105, 109, 111, 113, 119 (22 values).

  The list is FINITE-EXACT. Ω is complete iff the tournament has no two vertex-disjoint 3-cycles, since any
  non-transitive set contains a triangle. Such a realizer has a single nontrivial strong component, and Busch's floor
  b(10) = 125 caps its order at 9.

- **Theorem A (PROVED, elementary).** For m < b, m ∈ Gap_b ⟺ 2m + 1 is not the count of an Ω-complete tournament.
  *Proof.* If c2 ≥ 1 then m ≥ b. So m < b forces c2 = c3 = … = 0, i.e. Ω ≅ K_m, and then H = 1 + 2m. ∎
- **Corollary.** For every base b ≥ 60, Gap_b contains the 22 m-values of the list above:
  {3, 10, 14, 17, 20, 21, 27, 28, 29, 33, 35, 37, 38, 41, 44, 46, 50, 52, 54, 55, 56, 59}.
  This is consistent with THM-5070's tables for b ≤ 8, including m = 14 (N = 29) in bases 4, 6, 7, m = 27 (55) in
  bases 4, 5, and m = 21 (43) in base 6. The number of base-b gaps grows with b; HYP-9530 (B) asks only for
  finiteness at each b.
- **Prediction P13 (OPEN).** The Ω-complete spectrum misses a positive proportion of odd numbers (about 35% below
  125), so |Gap_b| ≥ c·b. The route is a structure theorem for tournaments with no two disjoint triangles. They
  include the one-centre (hub) family, but not only it: for example the regular 5-tournament, H = 15, has no hub.

## Wave 2: the universal gaps

- **Theorem B (PROVED).** m = 3 ∈ Gap_b for every b ≥ 2.
  *Proof.* m = 3 forces c1 = 3 and c2 = 0; in base 2 the alternative c1 = 1, c2 = 1 is impossible. So there are three
  pairwise-meeting odd cycles and no others. They lie in one strong component C, and every other component is
  trivial. But C has 1 odd cycle if |C| = 3 and 2 if |C| = 4 (the unique strong 4-tournament R_2). If |C| ≥ 5 it has
  at least |C| − 2 ≥ 3 triangles (Moon) plus a 5-cycle (pancyclicity), so at least 4 odd cycles. ∎
- **Theorem C (PROVED; uses 21 ∉ S).** m = 10 ∈ Gap_b for every b ≥ 9.
  *Proof.* For b ≥ 11, Theorem A applies and N = 21 is not attained at all. For b = 9 or 10, the decompositions
  10 = c1 + b·c2 force c2 = 0, since a disjoint pair needs c1 ≥ 2. So Ω ≅ K_10 again. ∎
  - With THM-4940 / THM-5070 for b ≤ 8: m = 10 is a gap for every b except **4 and 8**. In base 8 the exception is
    C3 ⇒ C3, with m = 2 + b = 10. In base 4 it is a profile c1 = 6, c2 = 1.
  - *Correction.* The brick_moves note listed base 3 as an exception; base 3 does miss N = 31 (m = 10).

## Wave 3: the brick law

- **Prediction P1** (from the brick_moves note). A count shared by two-component compositions (a, b) and (c, d) of n,
  with a ≤ b and c ≤ d, is an S-unit re-association with S = primes(Sp(a)) ∪ primes(Sp(c)).
  - **Verdict: CONFIRMED at n = 13.** (4,9)/(5,8): 373 shared counts; (4,9)/(6,7): 249; (5,8)/(6,7): 320. 0
    violations. Together with n = 11, 12 that is 1,481 shared counts with 0 violations.
- **Stronger prediction:** the fixed set {3, 5, 11, 13} suffices for all n.
  - **Verdict: REFUTED.** At n = 13, 425 = 5·85 = 17·25 moves a 17, a prime of Sp(6). At n = 14 there are 360 and
    363 failures. The movable bricks grow with the smallest core size, exactly as P1 says.
- **Prediction P7:** 273 = Φ3(16) is a Fermat-block cut 17 + 2^j·L.
  - **Verdict: REFUTED for sources of order ≤ 9.** No realizer of 17 up to order 9 has a power-of-two port. The
    ports of R_k are 2^(k+1) − 1, then a descending interior sequence, then 2^k + 2 and 2^k + 1; for example R_4 has
    31, 69, 44, 28, 18, 17. Powers of two occur only as q(R_2) = 32 (a left-padding root) and the port 16 of R_3. 273
    is attained by other means.

## Next waves (queued predictions)

- **P13.** Classify tournaments with no two disjoint triangles; prove the Ω-complete spectrum has density < 1.
- **P14.** P1 at n = 14 with the full prime set: (5,9)/(6,8), (5,9)/(7,7), (6,8)/(7,7).
- **P15.** The exponent pattern of Gap_b. For b ≥ 60, every m < b with 2m + 1 in the 22-list is a gap. Predict which
  of them survive for m ≥ b, extending THM-5070's 30,000 bound to b = 9–12.
