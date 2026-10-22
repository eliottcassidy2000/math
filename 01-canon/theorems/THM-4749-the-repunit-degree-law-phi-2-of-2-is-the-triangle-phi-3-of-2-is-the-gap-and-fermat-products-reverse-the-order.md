---
id: THM-4749
title: "The repunit degree law: for fixed h, the number N(n,h) of n-vertex tournament classes with h Hamiltonian paths grows like n^w(h), where w(h) is the longest factorization of h into strong values (all odd numbers >= 3 except 7 and 21; KNOWN, Session24). For the binary repunits this is read off the cyclotomic factorization 2^k -+ 1 = prod Phi_d(2): Phi_2(2) = 3 = H(C_3) divides 2^k + 1 exactly for odd k, Phi_3(2) = 7 is a gap and divides 2^k - 1 exactly when 3 | k, and 2^(2^j) - 1 is a product of Fermat numbers. Hence N(n, 2^k - 1)/N(n, 2^k + 1) -> 0 for k = 5, 6, 7, 9, 13, 15 and -> infinity for k = 4, 8, 12, 16; the degrees tie at k = 10, 11, 14 (FINITE-EXACT degrees for k <= 16, exact polynomials for h < 2100)"
status: >
  KNOWN: the polynomial law (Session24, 05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md:
  degree = maximum number of nontrivial strong components, onset n >= floor(3 log_5 h)). PROVED: the cyclotomic corollaries
  (statements 2-5), given the strong value set V below 80,265 (FINITE-EXACT, THM-4097/4102/4104 and the S50/S51 inventories).
  FINITE-EXACT: the degree table for k <= 16 (every prime factor of 2^k -+ 1 is below 80,265), the exact polynomials for
  h = 2^k -+ 1 < 2100, and the counts N(n, 2^k -+ 1) for 8 <= n <= 14, k <= 11. Conditional for k >= 17 on V-membership of
  the large prime factors (the Chan-Pak completeness conjecture, THM-1370 OPEN part). Session opus-2026-10-21-S51.
session: opus-2026-10-21-S51
source: 05-knowledge/results/repunit_rarity_hyp9357_20261021.md
scripts:
  - 04-computation/experiments/repunits_hyp9357_20261021/repunit_degrees.py (+ .out)
  - 04-computation/experiments/repunits_hyp9357_20261021/repunit_rarity.py (+ .out)
related:
  - THM-1370 (7 and 21 are never H)
  - THM-4744 (strong inventories; every odd d <= 623 except 7, 21 is a strong count)
  - THM-4748 (the Mersenne side needs three reversed arcs; the Fibonacci fan)
  - HYP-9359 (the arithmetic-support law at fixed n and fixed order)
---

# THM-4749 — the repunit degree law

## Setting

* `N(n, h)` is the number of isomorphism classes of `n`-vertex tournaments with exactly `h` Hamiltonian paths.
* `V` is the set of strong values: odd `d ≥ 3` other than 7 and 21.
  * FINITE-EXACT below 80,265: THM-4097/4102/4104. Below 2100 it is reproved by the complete strong inventory in `repunit_rarity.py`.
  * For all `d` it is the Chan–Pak completeness conjecture.
* `w(h)` is the largest `r` with `h = d_1 ⋯ d_r`, every `d_i ∈ V`; and `w(1) = 0`.

## Statements

1. **The degree law (KNOWN, Session24).**
   * `N(n, h) = Σ_r Σ C(n − M + r, r)`. The inner sum runs over ordered sequences of `r` strong classes of orders `m_i ≥ 3` with `∏ H = h`, and `M = Σ m_i`.
   * It is a polynomial in `n` of degree `w(h)` from `n = ⌊3 log_5 h⌋` on.
   * The leading coefficient is `1/w(h)!` times the number of ordered sequences of `w(h)` strong classes with product `h`.
2. **The triangle on the Fermat side (PROVED).**
   * For odd `k`, `3 = Φ_2(2) = H(C_3)` divides `2^k + 1` and does not divide `2^k − 1`.
   * For odd `3 ≤ k ≤ 17`, `(2^k + 1)/3 ∈ V`, so `w(2^k + 1) ≥ 2`.
   * Hence, whenever `2^k − 1` is prime or blocked, `deg N(·, 2^k − 1) ≤ 1 < 2 ≤ deg N(·, 2^k + 1)` (k ≤ 17: `k = 3, 5, 7, 9, 13, 17`). Blocked means every factorization uses 7 or 21.
3. **The gap on the Mersenne side (PROVED).** `7 = Φ_3(2)` divides `2^k − 1` iff `3 | k`, and never divides `2^k + 1`.
   * `63 = 3^2·7`: its factorizations `3·21`, `9·7`, `3·3·7` all use a gap, so `w(63) = 1` and `N(n, 63) = 2(n − 7)` for `n ≥ 7`. The order-8 atoms are THM-344.
   * `511 = 7·73` is blocked: `w(511) = 1` against `w(513) = 4`.
   * `4095 = 3^2·5·7·13` absorbs the 7 into `35`: `w(4095) = 4`.
4. **Fermat products reverse the order (PROVED).**
   * `2^(2^j) − 1 = F_0 F_1 ⋯ F_(j−1)`, with `F_i = 2^(2^i) + 1 = H(NT_(2^i+2))`, so `w(2^(2^j) − 1) ≥ j`.
   * `2^(2^j) + 1 = F_j` is prime for `j ≤ 4`.
   * So `N(n, 15)/N(n, 17)`, `N(n, 255)/N(n, 257)` and `N(n, 65535)/N(n, 65537)` tend to infinity.
5. **The table (FINITE-EXACT for k ≤ 16).** Here `w_−` stands for `w(2^k − 1)` and `w_+` for `w(2^k + 1)`.

   | k | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
   |---|---|---|---|---|---|---|---|---|---|---|---|---|---|
   | `w_−` | 2 | 1 | 1 | 1 | 3 | 1 | 3 | 2 | 4 | 1 | 3 | 2 | 4 |
   | `w_+` | 1 | 2 | 2 | 2 | 1 | 4 | 3 | 2 | 2 | 2 | 3 | 4 | 1 |

   * In the large-`n` limit, `2^k − 1` is the rarer value for `k = 5, 6, 7, 9, 13, 15`, and the commoner for `k = 4, 8, 12, 16`.
   * For `k = 10, 11, 14` the degrees tie and the leading coefficients decide. The exact polynomials give:
     * `k = 10`: third differences 24 against 18, so 1023 is commoner;
     * `k = 11`: second differences `184 = 2·S(23)·S(89)` against `516 = 2·S(683)`, so 2049 is commoner.
   * For odd `5 ≤ k ≤ 40` the Mersenne side never has the larger degree; the degrees tie at `k = 11, 23, 29`. This is conditional for `k ≥ 17`.
6. **The finite record (FINITE-EXACT).**

   | `n` | 15 / 17 | 31 / 33 | 63 / 65 | 127 / 129 | 255 / 257 | 511 / 513 | 1023 / 1025 | 2047 / 2049 |
   |---|---|---|---|---|---|---|---|---|
   | 8 | 20 / 6 | 10 / 23 | 2 / 25 | 10 / 40 | 28 / 30 | 32 / 35 | 0 / 0 | 0 / 0 |
   | 10 | 42 / 10 | 18 / 49 | 6 / 59 | 34 / 136 | 188 / 163 | 290 / 363 | 602 / 543 | 990 / 1029 |
   | 14 | 110 / 18 | 34 / 125 | 14 / 151 | 82 / 424 | 1176 / 455 | 946 / 2427 | 4770 / 3133 | 7426 / 10168 |

   Every `N(·, 2^k ± 1)` with `h < 2100` is the polynomial of degree `w(h)` from the stated onset (`repunit_degrees.py`).

## Proofs

**1.** Session24, §2. The condensation of a tournament is a transitive order of strong components, `H` is multiplicative over them, and placing the singleton components gives the binomials. Busch's bound gives the onset. ∎

**2.** For odd `k`, `2^k ≡ 2 ≡ −1 (mod 3)`. So `3 | 2^k + 1` and `2^k − 1 ≡ 1 (mod 3)`.
* For odd `3 ≤ k ≤ 17`, `(2^k + 1)/3` is one of 3, 11, 43, 171, 683, 2731, 10923, 43691. Each is odd, at least 3, not 7 or 21, and below 80,265, hence in `V`. So `2^k + 1 = 3·(2^k + 1)/3` is a two-factor factorization, and `w ≥ 2`.
* When `2^k − 1` is prime, its only factorization is itself, so `w = 1`. When it is blocked, the same holds. Primes in this range: 31, 127, 8191, 131071.
* Statement 1 then gives the degrees. ∎

**3.** `ord_7(2) = 3`, so `7 | 2^k − 1` iff `3 | k`. Since `ord_7(2)` is odd, −1 is not a power of 2 modulo 7, so `7 ∤ 2^k + 1`.
* For 63: the divisors of 63 lying in `V` are 3, 9 and 63, and each cofactor `21, 7, 1` either is a gap or ends the factorization. So 63 has no factorization with two or more factors.
* For 511: the divisors in `V` are 73 and 511, and the cofactor of 73 is 7.
* `N(n, 63) = 2(n − 7)` from statement 1 and the two order-8 atoms (THM-344; S50). ∎

**4.** `2^(2^j) − 1 = (2^(2^(j−1)) − 1)(2^(2^(j−1)) + 1)`, and by induction it equals `F_0 ⋯ F_(j−1)`. Every `F_i` is in `V`, since `F_i ≠ 7, 21` and `F_i ≤ 65537 < 80,265` for `i ≤ 4`. So `w ≥ j`. For `j ≤ 4`, `F_j` is prime, so `w(F_j) = 1`. ∎

**5–6.** `repunit_degrees.py` computes `w` from the factorizations. It checks every polynomial on the exact counts from the complete strong inventory below 2100: the `w`-th difference is a positive constant, and the next difference vanishes from the onset on. `repunit_rarity.py` prints the class counts. ∎

## What this does and does not explain

* **The large-`n` answer.** This is the exact large-`n` reason that `31, 63, 127, 511` are rarer than `33, 65, 129, 513`, and that `15, 255` are commoner than `17, 257`. The leading term counts orderings of strong components, and the triangle `C_3` is the cheapest component (`H = 3`, three vertices).
* **The small-`n` part.** At `n = 8, 9`, where the rarity was first seen (S50), the counts are dominated by single strong components of orders 6–8. There the strong spectra themselves are thin at 31, 63 and 127 (`S_7(63) = 0`, `S_7(127) = 1`). That small-order part is measured in HYP-9359 and the note, §5.
