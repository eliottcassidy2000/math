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

## Wave 4: structure of pairwise-meeting tournaments (FINITE-EXACT, orders 5–9)

Strong tournaments with no two vertex-disjoint triangles (Ω complete):

| order n | total | with a hub (one vertex on every triangle) | hub-free | path counts of hub-free ones |
|---|---:|---:|---:|---|
| 5 | 6 | 3 | 3 | 13, 15 |
| 6 | 11 | 7 | 4 | 25, 33 |
| 7 | 20 | 15 | 5 | 49, 69 |
| 8 | 37 | 31 | 6 | 97, 141 |
| 9 | 70 | 63 | 7 | 193, 285 |

**Exact fits.** The total is 2^(n−3) + n − 3, of which 2^(n−3) − 1 have a hub and **n − 2 are hub-free**. Every
hub-free one has a triangle transversal of size 2, and its count is **1 + 3·2^(n−3)** or **9·2^(n−4) − 3**. The first
sequence is THM-5126's line 1 + 3·2^s for s ≥ 2; the line's only missing point is 7.

- **Prediction P16.** The pattern holds for every n ≥ 5. Then the pairwise-meeting spectrum is the hub (one-centre)
  spectrum W together with {1 + 3·2^s} and {9·2^s − 3}. By Theorem A, for large b the base-b gaps below b are the
  m with 2m + 1 outside that set: about half of all m < b, if HYP-9469's one-centre density 1/2 holds. Base 2 has
  only two gaps; large bases have about b/2 below b.
  - *Test:* order 10, where the prediction is 8 hub-free classes with counts 385 and 573.
  - The τ ≤ 2 part should be a directed-triangle packing/covering (Erdős–Pósa type) fact; this has not been checked
    against the literature.

## Wave 5: the Two-Hub Conjecture (headline)

- **Prediction.** Two hubs suffice for Chan–Pak, while one hub misses about half.
- **Verdict: CONFIRMED far beyond expectation.** Order ≤ 9 tournaments with a triangle transversal of size ≤ 2
  attain every odd N < 700 except 7 and 21. The exact counter `scripts/twohub.c` (below) then shows that every odd
  N ≤ **22,156,805** except 7 and 21 is the path count of a two-hub tournament on at most 18 vertices (2,327,673 at
  16 vertices). The first miss
  at word length ≤ 14 is 2,327,675, a length-cap effect. Length ≤ 12 already covers every odd N ≤ 200,000 except 7
  and 21.

**The counter.** The spine is TT_m. Hubs u1 and u2 have words x, y ∈ {0,1}^m (x_i = [u1 → i]) and ε = [u1 → u2]. A
Hamiltonian path has the form S0 h1 S1 h2 S2, with increasing spine segments S0, S1, S2 and junction arcs at the
hubs. So H is a **rational series in the four-letter alphabet (x_i, y_i)**: a 36-state transfer (2 hub orders × 3 × 3
× 2 segment states), with H = u_ε·M_w·v0. It agrees with brute force on 300 random configurations (0 mismatches).

**Two-Hub Conjecture (TWO-HUB CHAN–PAK).** Every odd N ∉ {7, 21} is H(T) for a tournament T with two vertices
whose deletion leaves a transitive tournament. This is a strengthening of Chan–Pak and a reduction of it to one
explicit 4-letter numeration.

**Proof route.** By OCF, H = 1 + 2(c1 + c2 + c12) + 4b. Here c1, c2 count odd cycles through u1 only or u2 only,
c12 those through both, and b the disjoint pairs. When every u1-cycle is disjoint from every u2-cycle, b = c1·c2 and

    H = (1 + 2c1)(1 + 2c2) + 2·c12,

a product of two one-hub counts (each from the half-dense W) plus an even cross-cycle correction. The aim is a
sub-automaton in which c12 moves in steps of 1 over an interval while c1 and c2 stay fixed; combined with W·W this
gives a digit argument for every large odd N.

**Predictions.**
- **P17. CONFIRMED.** The horizon grows about tenfold per two extra spine vertices. At word length ≤ 16 (18
  vertices), every odd N ≤ **22,156,805** except 7 and 21 is a two-hub count; the first miss is 22,156,807. The run
  took 4 minutes.
- **P18.** A sub-automaton exists with H′ = 3H + d for d ∈ {0, 2, 4}, closed under itself: the clone-clock ternary
  digits realized inside two hubs. *Test:* search 1–2 letter suffixes with affine action on H for a fixed prefix
  class.
- **P19.** No odd N other than 7 and 21 needs more than two hubs. *Test:* the two-hub counter at word length 16
  combined with the 2^32 verification of THM-4875 (D).

## Next waves (queued predictions)

- **P13.** Classify tournaments with no two disjoint triangles; prove the Ω-complete spectrum has density < 1.
- **P14.** P1 at n = 14 with the full prime set: (5,9)/(6,8), (5,9)/(7,7), (6,8)/(7,7).
- **P15.** The exponent pattern of Gap_b. For b ≥ 60, every m < b with 2m + 1 in the 22-list is a gap. Predict which
  of them survive for m ≥ b, extending THM-5070's 30,000 bound to b = 9–12.
