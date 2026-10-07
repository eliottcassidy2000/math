---
id: THM-4591
title: "Integer cycles of the Terras map T(x) = x/2, (3x+1)/2 on Z with period L <= 301,993 are exactly the five known ones, {0}, {-1}, {1,2}, {-5,-7,-10} and the 11-cycle through -17, so their periods are exactly {1, 2, 3, 11}; as intervals 3^k : 2^L they are the octave, fifth, fourth and whole tone (Gersonides' four, gap 1; THM-4484) and the apotome 2187:2048 (gap 139); Ellison's five exceptional shapes x >= 12 are exactly the five shapes with the largest expected cycle count Lyn(L,k)/|2^L - 3^k| and carry no integer cycle, but carry the known 3x+5 cycles (27,17) and eight 3x+23 cycles on the negative integers (19,12)"
status: >
  FINITE-EXACT (two independent methods). (a) Lyndon-word enumeration (k >= L/2, sufficient by the lemma) with exact divisibility by 2^L - 3^k for every
  L <= 40 (3.06e10 words, count-checked against the necklace formula). (b) Minimal-element scan: every |x| below the
  lemma's bound (7.22e9 positive, 1.03e10 negative, the extremes over L <= 301,993, at 125743/79335 and 176251/111202)
  iterated to return or descent, with 0 capped and 0 overflowed runs. PROVED: the lemma (product identity) and the
  Gersonides classification (THM-4484). DICTIONARY: the interval names. For positive cycles much stronger results are
  KNOWN (Simons-de Weger, Hercher); the new content is the self-contained two-sign census to L <= 301,993 (THM-4484 had
  L <= 90) and the Ellison-shape table. Platonic-solid readings of these numbers are NUMEROLOGY (base rate p = 0.78).
session: mac-mini-2026-10-07-golden (platonic reader; census re-read and checked against THM-4484 here)
source: 05-knowledge/results/golden_collatz_resonance_20261007.md
scripts:
  - 04-computation/experiments/golden_20261007_readers/platonic/cycles/ (necklace.c + necklace_*.txt; minscan_mt.c + minscan_301993.txt; bounds3.py; heuristic*.py; nearcyc*.py)
related:
  - THM-4484 (free and sporadic cycles: Gersonides 2-1, 3-2, 4-3, 9-8; -17 sporadic, 1 of 30 necklaces; census p <= 90)
  - THM-4512 (Ellison's bound and exceptions), HYP-9230 (the same commas as resonance frequencies of class densities)
---

# THM-4591 — the five Pythagorean cycles

## Statements

1. **Lemma (PROVED).** Around a cycle, `Π_(odd x_i) (3 + 1/x_i) = 2^L`.
   * A positive cycle therefore has `3^k < 2^L` and `min x ≤ 1/(2^(L/k) − 3)`.
   * A negative cycle has `3^k > 2^L` and `min |x| ≤ 1/(3 − 2^(L/k))`.
2. **Census (FINITE-EXACT).** For every period `L ≤ 301,993` the integer cycles of `T` on `Z` are exactly:

   | cycle | shape `(L, k)` | `2^L − 3^k` | interval `3^k : 2^L` |
   |---|---|---|---|
   | `{0}` | (1, 0) | 1 | octave 2:1 (read `2^L : 3^k`) |
   | `{−1}` | (1, 1) | −1 | fifth 3:2 |
   | `{1, 2}` | (2, 1) | 1 | fourth 4:3 (`2^2 : 3`) |
   | `{−5, −7, −10}` | (3, 2) | −1 | whole tone 9:8 |
   | `{−17, …, −136}` | (11, 7) | −139 | apotome 2187:2048 |

   The −17 cycle is a single necklace out of 30 for its shape.
3. **Ellison's shapes (FINITE-EXACT for `L ≤ 1500`, ranking over all `k`; DICTIONARY).** Put `E(L,k) = Lyn(L,k)/|2^L − 3^k|`, the naive expected number of integer cycles. The five shapes with `L ≥ 12` and largest `E` are exactly Ellison's exceptions:

   | shape | `E` | gap | rational cycles realised (reduced denominators `< |gap|`) |
   |---|---|---|---|
   | (19, 12) | 0.371 | `−7153 = −23·311` | 23 (8), 311 (117): eight `3x+23` cycles on the negative integers (positive `3x−23` cycles), since `2^19 < 3^12` |
   | (16, 10) | 0.077 | `6487 = 13·499` | 499 (41); never 13 |
   | (27, 17) | 0.062 | `5077565 = 5·71·14303` | 5 (2), 71 (5), 355 (16), …: the `3x+5` cycles with minima 187 and 347 (KNOWN, Lagarias 1990) |
   | (13, 8) | 0.061 | `1631 = 7·233` | 233 (15); never 7 |
   | (14, 9) | 0.043 | `−3299` | none |

   * The next shape is (12, 7), with `E = 0.035`.
   * The expected total over the five is 0.61, and none carries an integer cycle.
   * This is the cycle-count reading of Ellison's inequality: both measure `|2^x − 3^y|/2^x`.

## Proof of 1

* Each odd step multiplies by `(3 + 1/x_i)/2`, and each even step by `1/2`. Going once around the cycle gives `Π (3 + 1/x_i) = 2^L`.
* All factors are positive, and `3 + 1/x` is monotone in `x`. So the factor at the minimal `|x|` (an odd element) is the largest factor of a positive cycle and the smallest factor of a negative cycle.
* Comparing it with the geometric mean `2^(L/k)` gives `3 + 1/m ≥ 2^(L/k)` (positive) and `3 − 1/m ≤ 2^(L/k)` (negative). (Corrected after audit B: the earlier text said "largest" for both signs.)
* The extremes over `L ≤ 301,993` are at the best approximations `125743/79335` (positive side) and `176251/111202` (negative side) of `log2 3`. ∎

## Remarks

* **Pythagorean reading (DICTIONARY).** The odd Terras step is a fifth up (`×3/2`) and the even step an octave down (`×1/2`). A cycle is a closed walk on the spiral of fifths whose net interval `3^k/2^L` is close enough to 1 for the cycle equation to have an integer solution.
  * The four cycles with gap 1 (Gersonides, Catalan, Størmer's 3-smooth superparticular ratios) are integral for every word.
  * The apotome cycle is a Diophantine accident (THM-4484).
* **{2,3,11}.** The cycle periods are 1, 2, 3 and 11, and 11/7 is a semiconvergent (mediant) of `log2 3`. Every link from the period 11 to the other 11s (Golay length, `PSL(2,11)`, eleven squares) is NUMEROLOGY. The eleven-squares 11 is `L_5`, through THM-4592.
* **Platonic solids.** Polyhedral numbers land among the numerators and denominators `≤ 120` of the convergents and semiconvergents of `log2 3` at the base rate (`p = 0.78`).
  * 12-TET ↔ the icosahedron's 12 vertices has no equivariant form, because `A_5 × C_2` has no element of order 12 (PROVED).
  * The comma `3^12/2^19` is Ellison's (19, 12), and it is also the 12-mode of HYP-9230. That 12 is structural; the icosahedral 12 is not.

**Audit (2026-10-07, independent audit B).**
* CONFIRMED:
  * the scan bounds, which are exact maxima over `L ≤ 301,993` (7,216,102,492.69 and 10,295,871,816.1);
  * the scan logic;
  * an independent full re-run with signed 128-bit code, giving exactly −1, −5, −17 and 1 (0 capped, 0 overflowed);
  * the Lyndon tally (30,594,555,229 words);
  * the Ellison ranking over all `k`, and the rational-cycle table.
* Corrected above (MISTAKE-584): the negative-cycle proof sentence and the sign of the 3x+23 cycles.
* The novelty of the negative side has not been checked against published 3x−1 verifications.
