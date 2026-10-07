---
id: THM-4561
title: "The Sierpinski skew-Hadamard tower in the Littlewood problem: its rows are Walsh rows plus sparse dyadic spike trains, every doubling step is a Morse step x -> x(1 +- z^m) up to one flipped entry, so every row has merit factor < 2 and the rows are Thue-Morse-like (max F -> 0, sup/sqrtN -> infinity, up to a perturbation sketch); the Rudin-Shapiro switching D H D of the same tower has every row sup <= (3 + 2 sqrt2) sqrt(2N), i.e. merit factor >= 1/(33 + 24 sqrt2)"
status: "PROVED (elementary: closed form by induction, a two-line autocorrelation transfer identity, Cauchy-Schwarz, Golay/Davis-Jedwab pieces); 'max F -> 0' PROVED up to a perturbation sketch (one flipped entry per step is an O(N^(-1/2)) relative change of (R, s)). NUMERICAL: the rates. FINITE-EXACT: H_16 has no bent switching. Found by the session's reader of openai/math #076 (Littlewood polynomials) and re-checked independently. INDEPENDENTLY AUDITED 2026-10-06 (audit A: doubling, closed form, transfer identity, Walsh bound, Golay claim and the (3 + 2 sqrt2) sqrt(2N) bound re-derived; PASS WITH CORRECTIONS, applied: the F-range, wording, typing of the column remark)."
session: mac-mini-2026-10-06-oaimath2
source: 05-knowledge/results/oai2_openai_math_second_reading_20261006.md
scripts:
  - 04-computation/experiments/oai2_20261006_littlewood_tower.py (+ .out, ALL CHECKS PASSED)
related:
  - HYP-9162 / THM-4557 (the tower T_k; Aut = F_21)
  - openai/math #076 (ultraflat Littlewood polynomials, unrefereed) and its two companions (merit factor -> infinity; two-sided flat)
  - J. A. Davis, J. Jedwab, Peak-to-mean power control in OFDM, Golay complementary sequences, and Reed-Muller codes (1999)
---

# THM-4561 — the tower's rows are Thue–Morse; Rudin–Shapiro switching flattens them

## Setting

* The tower: `H_2 = [[1, 1], [−1, 1]]`, `H_(2N) = [[H, H], [−H^T, H^T]]`. Each `H_(2^k)` is skew Hadamard; normalized, it is the doubly regular tournament `T_k` of THM-4557.
* Rows and columns as polynomials: `r_j(z) = Σ_y H[j, y] z^y` and `c_j` similarly; `w = z^N`.
* Merit factor `F = N² / (2 Σ_(u ≥ 1) C_u²)` and `R = ‖x‖₄⁴ / N² = 1 + 1/F`.
* For a sequence `x` of length `m`, `s = Σ_(u=1)^(m−1) C_u C_(m−u) / ‖x‖₄⁴`.

## Statements (PROVED)

1. **Doubling.**
   * Rows of `H_(2N)` are `r_j (1 + w)` and `(r_j − 2z^j)(1 − w)`.
   * Columns are `c_j ∓ w r_j`, with `c_j = 2z^j − r_j`.
   * So rows take Morse steps. On columns the step has the Rudin–Shapiro/Golay form, but `(c_j, r_j)` is an anti-aligned pair, and numerically it does not flatten (NUMERICAL: max column `F` = 4, 4, 1.10, 0.68, 0.46, 0.30, 0.23 for `k = 2..8`).
2. **Closed form.** `r_i(y) = w_i(y) − 2 Σ_(l : i_l = 1) [y ≡ i (mod 2^l)] (−1)^(Σ_(m ≥ l) i_m y_m)`, with `w_i(y) = (−1)^(popcount(i & y))`. That is, a Walsh row plus sparse dyadic spike trains. (Induction on `k` using 1.)
3. **Morse transfer.** For `y = (x, εx)`: `‖y‖₄⁴ = 6‖x‖₄⁴ + 8εS` and `S_y = ‖x‖₄⁴ + 4εS`. Hence `R′ = (3/2 + 2εs)R` and `s′ = f(εs)`, with `f(t) = (1 + 4t)/(6 + 8t)`.
4. **Merit factor below 2.**
   * Cauchy–Schwarz gives `|s| ≤ (1 − 1/R)/2`, hence `F(x(1 ± z^m)) ≤ 2F/(F + 1) < 2` for `len x ≥ 2`.
   * **Every row of `H_(2^k)`, `k ≥ 2`, has `F < 2`.** Columns can reach `F = 4` at `k = 2, 3`.
   * **No two rows are Golay complementary, for `k ≥ 2`.** Write the rows of `H_(2N)` as in 1.
     * Upper rows vanish at the `N` roots of `z^N = −1`, and lower rows at the `N` roots of `z^N = 1`. So two rows from the same half share zeros.
     * For a mixed pair, `z = 1` gives `|u(1)|² + |v(1)|² = 4 r_j(1)² ∈ {0, 4N²}`. A Golay pair of length `2N` would need `4N`, so this fails for `N ≥ 2`.
5. **Walsh rows.** `f([−1/4, 1/4]) = [0, 1/4]`, and two steps multiply `R` by at least `7/4 + t ≥ 3/2`. So `F(Walsh row of length 2^k) ≤ 1/((3/2)^⌊k/2⌋ − 1)`.
6. **Tower rows are Thue–Morse-like (PROVED up to a sketch).** The same transfer holds with one flipped entry per lower-half step, an `O(m^(−1/2))` relative change of `(R, s)`. So `max_row F = O((2/3)^(k/2)) → 0` and `sup/√N ≥ √R → ∞`.
   * NUMERICAL: max row `F` is 0.236, 0.160, 0.127, 0.093, 0.073, 0.055 for `k = 8..13`.
   * Best Walsh row: `F ≈ 0.905·(4/(1+√17))^k`. The tower's best row is about 1.47 times that.
7. **Rudin–Shapiro switching flattens.**
   * Let `D = diag((−1)^(Σ y_m y_(m+1)))`. Then `D H D` is skew Hadamard and renormalizes to the same `H`, since row 0 is all ones.
   * **Every row of `D H_(2^k) D` has `sup_(|z|=1) |r| ≤ (3 + 2√2)·√(2N) ≈ 8.24 √N`, hence `F ≥ 1/(33 + 24√2) ≈ 0.0149`.**
   * *Proof.* By 2, `RS ⊙ r_i` is `±RS ⊙ w_i` (a Golay sequence, `|·|² ≤ 2N`) plus, for each `l ≥ 1` with `i_l = 1`, a term `2z^(i mod 2^l) G_l(z^(2^l))`.
   * On the progression `y ≡ i (mod 2^l)` the Rudin–Shapiro form restricts to a path form plus a linear form. So `G_l` is a Davis–Jedwab Golay sequence of length `N/2^l`.
   * Sum: `√(2N)(1 + 2Σ_(l ≥ 1) 2^(−l/2)) = (3 + 2√2)√(2N)`. ∎
   * NUMERICAL: the worst row has `sup ≤ 3.0√N` (`k ≤ 12`), and `F ∈ [0.70, 3.20]` for `3 ≤ k ≤ 12` (at `k = 12`: `[0.753, 3.001]`, mean 1.30).

## FINITE-EXACT and conditional remarks

* `H_16` (tower) has no `±1` vector `d` with `|H_16 d| ≡ 4` entrywise (exhaustive; Sylvester's `H_16` has 448 with `d_0 = 1`).
  * So no row/column negation of the tower's `H_16` is regular. We have not matched its class against Hall's five classes of order 16.
  * Every switching of `T_4` has some row sum of absolute value at most 2.
* **CONDITIONAL on #076 (trivial).** Row 0 of `D H D` is `d` itself (up to sign), so a switching by an ultraflat `d` makes one row ultraflat.
* **Open.** Is there, for every `k`, a switching of `T_k` all of whose rows are two-sided flat (`c√N ≤ |r| ≤ C√N`)? Ultraflat?
  * #076's sequences are not a dyadically closed family. Davis–Jedwab Golay is such a family, and it gives only the upper bound.

## Not claimed

Nothing about the merit-factor or ultraflatness problems themselves. The tower is a specific dyadic family, and its rows sit in the Riesz-product (singular) regime. The RS-switched rows are sums of `√2`-flat Golay pieces, with proven `sup ≤ 8.24√N` (observed about `3.0√N`).
