---
id: THM-4580
title: "Zeroless powers of two: the trailing digits of 2^n form a 5-adic tree in which every zeroless class has 4 or 5 zeroless lifts according to one parity bit; Z_k = #(zeroless k-digit multiples of 2^k) satisfies Z_(k+1) = (9 Z_k + Delta_k)/2 with Delta_k a cyclotomic-unit trace; the growth rate lies in [4.47848, 4.52386] and the exponent set has 5-adic dimension in [0.93156, 0.93782]; every 2^n with 87 <= n < 1.1e11 contains a 0; no argument using finitely many leading or trailing digits can settle the conjecture"
status: "PROVED: lift lemma, bijection, recursion, unit formula (checked for m <= 12), doubling criterion, finite-digit obstruction. FINITE-EXACT: Z_k for k <= 40 (Z_1..Z_26 = OEIS A181610; Z_1..Z_9 re-enumerated independently here); verification of 87 <= n < 1.1e11 (session reader's C verifier, validated against Python big integers, end states checked against exact 2^N mod 10^288; independently re-verified here for n < 2e9 with separate code, reproducing OEIS A031142 records 24-38). PROVED (computer-assisted): the growth and dimension brackets. The conjecture itself (no zeroless 2^n with n > 86) is OPEN. KNOWN context: OEIS A007377 records a check to 1e10 (Radcliffe 2022); A031142's record table (Griffiths 2012), if complete, implies the conjecture for n < 7.88e12."
session: mac-mini-2026-10-07-oaimath3 (owner prompt "investigate whether every power of two above 2^86 contains a zero"; work by the session's zeroless reader, re-checked here)
source: 05-knowledge/results/oai3_two_orbits_twos_and_threes_20261007.md
scripts:
  - 04-computation/experiments/oai3_20261007_readers/zeroless/ (verify.c, verify2.c, zk_mitm.c, fourier_levels.c, checks.py (ALL CHECKS PASSED), leading_measure.py, heuristic.py, run logs)
  - 04-computation/experiments/oai3_20261007_zeroless_verify.c (+ .out: independent verifier, n < 2e9)
related:
  - HYP-9219 (growth exactly 9/2, dimension log_5(9/2))
  - THM-4556 (iv) (the 2-adic clock of 3: the same lifting-the-exponent step); THM-4554 (Moran-type sieve); THM-4072 / THM-3848 (Mahler 3/2 finite-prefix obstruction)
  - J. Lagarias, Ternary expansions of powers of 2, J. LMS 79 (2009); R. Saye, J. Integer Seq. 25 (2022) (ternary lift lemma); W. Narkiewicz (1980)
  - OEIS A007377, A181610, A031142
---

# THM-4580 — the 5-adic tree of zeroless powers of two

## Statements

1. **Lift lemma (PROVED).**
   * Let `T_k = 4·5^(k−1)`. Lifting the exponent gives `2^(T_k) = 1 + 5^k u` with `5 ∤ u`.
   * For `n ≥ k + 1` and `j = 0, …, 4`, the numbers `2^(n + jT_k)` share their last `k` digits, and their digit `k + 1` runs through the five digits of one parity.
   * So a residue class whose last `k` digits are zeroless has 5 zeroless lifts if that parity is odd, and 4 if it is even.
   * Base 3 (Saye) keeps exactly 2 of 3 lifts at every node. Base 10 branches on a parity bit.
2. **Bijection (PROVED).** Because 2 is a primitive root mod `5^k`, the `k`-digit tails of `2^n` (`n ≥ k`) are exactly the `k`-digit strings divisible by `2^k` and prime to 5. Hence `Z_k`, the number of zeroless tail classes, equals the number of zeroless `k`-digit multiples of `2^k`. This is OEIS A181610.
3. **Recursion and unit formula (PROVED).**
   * Let `q = x/2^k` and `Δ_k = #(odd q) − #(even q)`. Then `Z_(k+1) = (9Z_k + Δ_k)/2`, and `(2/9)^k Z_k = 1 + Σ_(j<k) Δ_j 2^j/9^(j+1)`.
   * `Δ_(m−1) = 2^(1−m) Tr P_m(ζ_(2^m))`, where `P_m(ζ) = ζ^((10^m−1)/9) Π_(i ≤ m−4) (1 − ζ^(9·10^i))/(1 − ζ^(10^i))` is a cyclotomic unit. Checked for `m ≤ 12`.
4. **Data (FINITE-EXACT).**
   * `Z_k` for `k ≤ 40`, for example `Z_40 = 119333906141890097435122400`. `Δ_k` for `k ≤ 39`.
   * `|Δ_k|^(1/k) ≈ 1.6`, far below the `≈ 2.12` of a random parity.
   * `c = lim (2/9)^k Z_k = 0.8876940431151482645` (NUMERICAL, `±3·10^(−19)`).
5. **Brackets (PROVED, computer-assisted).**
   * The growth rate lies in `[4.47848, 4.52386]`.
   * `dim_H` of the 5-adic exponent set lies in `[0.93156, 0.93782]`.
   * `#{n ≤ x : 2^n zeroless} ≪ x^0.93783`, the decimal analogue of Narkiewicz's ternary bound.
6. **Doubling criterion (PROVED).** For zeroless `x`, `2x` is zeroless iff `x` contains none of the blocks `51, 52, 53, 54` and does not end in 5. A zero arises exactly at a 5 that receives no carry. This explains runs such as `n = 31, …, 37`.
7. **Verification (FINITE-EXACT).** Every `2^n` with `87 ≤ n < 1.1·10^11` contains a 0. For `n ≥ 957` the 0 lies among the last 251 digits.
   * The rightmost-zero records reproduce OEIS A031142 entries 24–41, including `n = 103233492954` (249 zeroless trailing digits) and `n = 109171987836` (250).
   * An independent verifier written in this session confirms `n < 2·10^9` (records `1757, …, 781717865`).
8. **Finite-digit obstruction (PROVED, elementary).** For every `k` the zeroless tail classes are nonempty and the exponent set is perfect. Every leading-digit set has positive measure. So no argument using finitely many trailing digits, or finitely many leading digits, can settle the conjecture.
   * A proof must couple 5-adic and archimedean information through the middle digits.
   * Collatz has the same shape: 2-adic parity information must be coupled with the size condition `2^K` versus `3^d` (THM-4564 (6): the debt's archimedean size against its 2-adic valuation).

## Heuristic (HEURISTIC)

* The model is `P(2^n zeroless) ≈ τλ·0.9^(D(n))`, where `D(n)` is the digit count, `τ = 5c/4 = 1.10962` is the 5-adic end factor, and `λ = 1.08453` is the Benford end factor (leading digits).
* It predicts 33.9 cases with `n ≤ 86`; there are 36.
* It predicts 2.31 cases beyond 86, so an empty tail was a 10–15% event a priori. Beyond `1.1·10^11` it predicts about `10^(−1.5·10^9)`.
* The joint frequency of zeroless leading and trailing digit blocks factorizes. For fixed block lengths this is PROVED (Weyl). Up to `n = 10^10` it holds NUMERICALLY deep into the tail, answering Lagarias's (2009) question empirically in base 10.

## Not claimed

* The conjecture itself.
* That the growth rate is exactly 9/2 (HYP-9219).
* That openai/math helps: none of its 722 manuscripts concerns decimal digits.
