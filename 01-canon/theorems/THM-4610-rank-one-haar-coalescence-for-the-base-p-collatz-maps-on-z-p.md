---
id: THM-4610
title: "Rank-one Haar coalescence on Z_p: for every odd prime p and the base-p Collatz map C_p(x) = x/p on pZ_p, ((p+1)x + p - i)/p on x = i mod p, i != 0 (on positive integers x/p or ceil((p+1)x/p); Carnielli's T_p, in Hasse's class and Moller's family with multiplier p+1; p = 2 gives the Terras map and THM-4581), two Haar p-adic orbits related by u = (p+1)^k v + e with (p+1)^max(0,-k) e in Z merge at equal time almost surely. The flip skeleton of the debt k is a fair simple random walk (k itself is lazy); at every level h >= 1 the normalized offset is multiplied by the branch factor (1/p or (p+1)/p) of a fresh uniform reference digit, up to a bounded additive term; a level weight s^|k| with s < 1 makes the accumulated drift summable; an explicit descent path gives accessibility from every integer offset; P(no merge by T) >= c T^(-1/2). This is the almost-sure part of HYP-9244's requirement 1 (rank one, d >= 3) for the maps C_p, including Z_7, Z_11 and Z_13"
status: >
  PROVED: statements 1-5 (elementary probability: the exact integer pair chain; V_n minus its accumulated additive terms
  is a supermartingale, and Fatou is applied to V >= 0; optional stopping, Levy's extension of Borel-Cantelli, Levy's 0-1
  law). The route is HYP-9244's requirement-1 sketch (THM-4581's architecture with a level weight s^|k|), whose per-step
  drift audit G checked in rank1_weight_check.py (mac-mini). New here: the odd-p valuation-run lemma 1(e), the
  Borel-Cantelli plus injectivity step that replaces THM-4581's parity argument, the explicit descent of statement 4, and
  the visit-length bookkeeping L_h. FINITE-EXACT: the local lemmas were also checked exhaustively on boxes of states x
  digits for p = 2, 3, 5, 7, 11, 13; the descent path for 2 <= |e| <= 20000. KNOWN: the remark after statement 5 (the
  3x - 1 barrier; Akin 2004's "the method is too general") and the cycles of C_p (Carnielli's Table 1, computed by
  K. R. Matthews), re-verified here. NUMERICAL: statement 6 (merge-time tails, cycles and repunit lines on Z_7, Z_11,
  Z_13). AUDITED 2026-10-08 by the session's independent audit: proofs sound and complete; corrections applied
  (MISTAKE-649). Session opus-2026-10-08-S22.
session: opus-2026-10-08-S22 (owner: "keep aiming at remaining collatz steps, think about Z_7, Z_11 and Z_13")
source: 05-knowledge/results/zp_rank_one_coalescence_20261008.md
scripts:
  - 04-computation/experiments/zp_rank_one_20261008/zp_chain_checks.py (+ .out) (A)-(F): table, local lemmas, descent, level weight, cycles, dictionary
  - 04-computation/experiments/zp_rank_one_20261008/zp_basins_survey.py (+ .out): cycles and basins, both signs (remark after statement 5), boundaries
  - 04-computation/experiments/zp_rank_one_20261008/zp_tails.py (+ .out): merge-time tails (statement 6)
  - 04-computation/experiments/zp_rank_one_20261008/zp_repunit_lines.py (+ .out): repunit lines (statement 6)
related:
  - THM-4581 (p = 2: the Terras map; this proof follows its architecture)
  - HYP-9244 (Polya trichotomy; requirement 1, rank one for d >= 3); THM-4606, THM-4607, THM-4608, THM-4609, THM-4611 (mac-mini, the other rows)
  - 05-knowledge/results/debt_rank_least_requirements_20261008.md (section 6 item 1; the words "Not yet written" were in section 4 at commit 97b640656); 04-computation/experiments/rank_20261008/audit_G/rank1_weight_check.py (the level-weight drift); rank_20261008/rank1_two_valued.out (row Z_3 (1, 4, 4) is C_3)
  - THM-4590 (p = 2: o(x) cuts; the negative-integer basins of 3x + 1, i.e. the cycles of 3x - 1)
  - Matthews and Watts (1984); Hasse; Moller (1978); W. Carnielli, Some natural generalizations of the Collatz problem, Applied Mathematics E-Notes 15 (2015) 207-215 (T_d; Table 1 by K. R. Matthews); Bruschi arXiv:0810.5169 (the undivided form p C_p); Tremblay arXiv:1907.13302 (cycles of T_d); Goncalves-Greenfeld-Madrid arXiv:2111.06170 (Remarks 1.6, 1.7); Kontorovich-Lagarias arXiv:0910.1944 (section 1.1); Akin, Why is the 3x + 1 problem hard?, Contemp. Math. 356 (2004); Durrett, Probability (Thms 4.3.1, 4.3.4)
---

# THM-4610 — rank-one Haar coalescence for the base-p Collatz maps on Z_p

## Setting

* `p` is an odd prime and `P = p + 1`. For `x ∈ Z_p` let `δ(x) ∈ {0, …, p − 1}` be its last digit, `x mod p`.
* **The map.**

      C(x) = x/p                    if δ(x) = 0,
      C(x) = (P x + p − δ(x))/p     if δ(x) ≠ 0.

  * On positive integers `C(x) = x/p` or `⌈P x/p⌉`. For `p = 2` this is the Terras map `x/2`, `(3x+1)/2`.
  * In the Matthews–Watts form `(m_i x + r_i)/p` the multipliers are `m_0 = 1` and `m_i = P` for `i ≠ 0`, all `≡ 1 mod p` (translation-only), with two values (rank one). It contracts: `Λ = (1/p) ln(1/p) + ((p−1)/p) ln(P/p) < 0`.
* **Digits.** Let `y` be Haar on `Z_p`, `v_n = C^n(y)`, `β_n = δ(v_n)` and `F_n = σ(y mod p^n)`.
  * **Fact D.** `β_n` is uniform on `{0, …, p − 1}` and independent of `F_n`, and `v_n` is Haar given `F_n`.
  * **Injectivity.** Two points of `Z_p` with the same digit sequence are equal.
  * Proof: on each class `i + pZ_p`, `C` is an affine bijection onto `Z_p` with a `p`-adic unit multiplier, and it scales Haar measure by `p`. If `x ≠ x'` share `n` digits, then `C^n x − C^n x'` has valuation `v_p(x − x') − n`; so the digits eventually differ. ∎
* **The pair chain.** Let `u_0 = P^(k_0) y + e_0`, where `k_0 ∈ Z`, `A_0 := P^max(0,−k_0) e_0 ∈ Z` and `(k_0, A_0) ≠ (0, 0)` (an *admissible* start). Put `u_n = C^n(u_0)`.
  * For every `n`, `u_n = P^(k_n) v_n + e_n` with `A_n = P^max(0,−k_n) e_n ∈ Z`.
  * The digit of `u_n` is `j_n = (β_n + A_n) mod p`, because `P ≡ 1 mod p`.
  * With `i = β_n`, `j = j_n` and `h = −k` when `k < 0`:

    | `(i, j)` | `k ≥ 1` | `k = 0` | `k ≤ −1` |
    |---|---|---|---|
    | `(0, 0)` | `(k, A/p)` | `(0, A/p)` | `(k, A/p)` |
    | `(≠0, ≠0)` | `(k, [PA + (p−j) − P^k (p−i)]/p)` | `(0, [PA + i − j]/p)` | `(k, [PA + P^h (p−j) − (p−i)]/p)` |
    | `(0, ≠0)` | `(k+1, [PA + p − j]/p)` | `(1, [PA + p − j]/p)` | `(k+1, [A + P^(h−1) (p−j)]/p)` |
    | `(≠0, 0)` | `(k−1, [A − P^(k−1) (p−i)]/p)` | `(−1, [PA − (p−i)]/p)` | `(k−1, [PA − (p−i)]/p)` |

  * Every numerator is divisible by `p`, since `j ≡ i + A` and `P ≡ 1`.
  * `(0, 0)` is absorbing, and absorption means `u_n = v_n`: a merge at equal time.
  * `k_n − k_0` counts the steps where `v` is divided and `u` multiplied, minus the reverse.
* **Notation.**
  * Level `h = |k|`.
  * Normalized offset `f = e P^(−max(k, 0))`, so `f = e` when `k ≤ 0`.
  * Reference digit `ρ_n`: `β_n` if `k_n ≥ 1`, `j_n` if `k_n ≤ −1`.
  * Branch factor `F(ρ)`: `1/p` if `ρ = 0`, `P/p` otherwise.
  * A *flip* is a step that changes `k`, i.e. exactly one of `i`, `j` is 0.
  * `μ_h = v_p(P^h − 1) = 1 + v_p(h)` for `h ≥ 1`.
  * `c = (p − 1)/p`.
  * For `θ > 0`:
    * `κ(θ) = (1/p) p^(−θ) + ((p − 1)/p)(P/p)^θ`. It is strictly convex with `κ(0) = κ(1) = 1`; `κ(1) = 1` because the multipliers sum to `p²` (`1 + (p − 1)(p + 1) = p²`). So `κ(θ) < 1` for every `θ ∈ (0, 1)`. Fix such a `θ`.
    * `g(s) = (1/p) p^(−θ) s + (1/p)(P/p)^θ s^(−1) + ((p − 2)/p)(P/p)^θ`. Since `g(1) = κ(θ) < 1`, there is `s ∈ (0, 1)` with `g(s) ≤ 1`. Fix one.
    * At `θ = 1/2` the least such `s` is 0.788, 0.762, 0.754 for `p = 7, 11, 13`. The same equation at `p = 2` gives THM-4581's `s = 0.8966`.

## Statements

1. **Local lemmas (PROVED; also exhaustive on boxes of states × digits for p = 2, 3, 5, 7, 11, 13).**
   * (a) **Factor law.** At every step `|f_(n+1)| ≤ F |f_n| + c`, where:
     * at levels `h ≥ 1`, `F = F(ρ_n)`;
     * at `k = 0`, a non-flip step has `i = j = 0` (`F = 1/p`) or both digits nonzero (`F = P/p`), and a departure (flip at `k = 0`) has `F = 1/p`.
   * (b) **Direction.** At a level `h ≥ 1`, a flip moves toward 0 iff `ρ_n ≠ 0`, at cost `P/p`. It moves away iff `ρ_n = 0`, and then pays `1/p`.
   * (c) **Flip law.** Conditionally on `F_n`, `ρ_n` is uniform. If `p ∤ A_n`, a step is an away flip, a toward flip or a non-flip with probabilities `1/p`, `1/p` and `(p − 2)/p`; a non-flip then has both digits nonzero. If `p | A_n`, then `i = j` and no flip occurs.
   * (d) **Runs at `k = 0`.** If `p | A_n` at `k_n = 0`, then `e_(n+1) = e_n/p` or `(P/p) e_n` exactly, and `v_p` drops by one. So the chain stays at `k = 0` for exactly `v_p(e)` steps of exact multiplication.
   * (e) **Valuation runs at `h ≥ 1`.** Suppose `p | A_n` and `w = v_p(A_n)`.
     * If `w < μ_h`, then `v_p(A_(n+1)) = w − 1` for every digit.
     * If `w ≥ μ_h`, then `v_p(A_(n+1)) < μ_h` for exactly `p − 1` of the `p` digits.
     * So a maximal stretch with `p | A` at level `h` lasts at most `μ_h − 1 + G` steps, with `G` conditionally `Geom((p−1)/p)`.
2. **Drift and return drift (PROVED).** Let `V_n = |f_n|^θ s^(|k_n|)`.
   * At a level `h ≥ 1`: `E[V_(n+1) | F_n] ≤ V_n + c^θ s^(h−1)`.
   * At `k = 0`: `E[V_(n+1) | F_n] ≤ λ_0 V_n + c^θ`, with `λ_0 = max(κ(θ), (2/p) p^(−θ) s + ((p − 2)/p)(P/p)^θ) < 1`. When `p | A_n` the bound is `κ(θ) V_n` exactly.
   * Let `τ_0 < τ_1 < …` be the successive arrival times at `k = 0` (with `τ_(j+1) = τ_j + 1` after absorption). Then

         E[ |e_(τ_(j+1))|^θ | F_(τ_j) ]  ≤  λ_0 |e_(τ_j)|^θ + C_p(θ),     C_p(θ) < ∞.

3. **Theorem: almost-sure coalescence (PROVED).** From every admissible start the chain is absorbed almost surely. Equivalently, for Haar-almost every `y ∈ Z_p`, `C^n(P^(k_0) y + e_0) = C^n(y)` for some `n`.
4. **Accessibility by explicit descent (PROVED).** Every state `(0, e)` with `e ∈ Z` reaches `(0, 0)` along a finite digit string:
   * `p | e`: digit 0 gives `(0, e/p)`.
   * `e ≥ 2`, `p ∤ e`, with `ē = e mod p`: digit `p − ē` gives `(−1, A′)` with `A′ = (Pe − ē)/p`. Then digit 0 while `p | A`, ending at some `A*`. Then digit 0 once more gives `(0, ⌈A*/p⌉)`, and `1 ≤ ⌈A*/p⌉ ≤ (Pe − 1 + p² − p)/p² < e`.
   * `e ≤ −2`, `p ∤ e`: digit 0 gives `(1, (Pe + p − ē)/p)`. Then digit 0 while `p | A`, ending at some `A*`. Then digit `(−A*) mod p` gives `(0, ⌈A*/p⌉ − 1)`, with `e < ⌈A*/p⌉ − 1 ≤ −1`.
   * `e = 1`: digits `0, p − 2` give `(0, 1) → (1, 2) → (0, 0)`.
   * `e = −1`: digits `1, 0` give `(0, −1) → (−1, −2) → (0, 0)`.
5. **Consequences (PROVED).**
   * (a) **Translates.** For every nonzero integer `e`, Haar `y` and `y + e` merge at equal time almost surely. Absorption by time `K` depends only on `y mod p^K`, so for each `e` the integers `n` with `C^t(n) = C^t(n + e)` for some `t ≤ K` have density tending to 1 as `K → ∞`.
   * (b) **Repunit lines.**
     * For an integer `E ≥ 1`, `R_E = p^E − 1` rises for `E − 1` steps to `x_E = p P^(E−1) − 1`. For `E ∈ Z_p` put `x_E = p P^(E−1) − 1`, using the `p`-adic exponential (`R_E` itself has no meaning for `p`-adic `E`).
     * `x_(E−D) = P^(−D) x_E + P^(−D) − 1`, so with reference orbit `y = x_E` the pair `(x_E, x_(E−D))` is the admissible start `(−D, 1 − P^D)`.
     * `E ↦ P^(E−1)` is a measure isomorphism from `Z_p` onto `1 + pZ_p` (p odd), so `x_E` is Haar on `p − 1 + p²Z_p`.
     * Hence for each `D ≥ 1` the orbits of `x_E` and `x_(E−D)` merge at equal time for Haar-almost every `E ∈ Z_p`. Absorption by time `n ≥ 2` depends only on `E mod p^(n−2)`, so the integer exponents `E` with `C^n(x_E) = C^n(x_(E−D))` for some `n ≤ N` (an exit `R_E ⇝ R_(E−D)` with shift `D`, certified by time `N + E − 1`) have density tending to 1 as `N → ∞`.
   * (c) **Rate, lower bound.** For every non-absorbed admissible start, `P(no merge by T) ≥ c T^(−1/2)`. The upper bound `C T^(−1/2)(log T)^2` of THM-4581 (4) should follow from the same sketch; it is not claimed here.
   * (d) **p = 2.** The argument specializes to `p = 2` and reproves THM-4581 3 there, with two differences: `μ_h = v_2(3^h − 1)` (1 for odd `h`, `2 + v_2(h)` for even `h`) replaces `1 + v_p(h)`, and the `λ_0` obtained is weaker than THM-4581's. There are no lazy steps, and `g(s) = 1` is THM-4581's equation for `s`. Statement 4 is replaced by THM-4581 3(iii).
   * (e) **Basin boundaries have density zero.** Let `B` be any set of positive integers that is a union of grand orbits of `C_p` restricted to `Z_(>0)` (for example the basin of one cycle). Then for every `e ≥ 1` the set `{n : n ∈ B xor n + e ∈ B}` has natural density 0. This is the base-`p` analogue of THM-4590's `o(x)` cuts.
   * **Remark (the barrier; KNOWN).** Statements 1–4 and 5(a)–(e) are statements about Haar measure on `Z_p`. They are invariant under the measure-preserving conjugation `x ↦ −x`, so they hold verbatim for `F_p(x) = −C_p(−x)`, which is `x/p` or `⌊(p+1)x/p⌋` on the positive integers and fixes `1, …, p − 1`.
     * For `p = 2`, `F_2` is the `3x − 1` map. Its positive-integer cycles `{1}`, `{5, 7, 10}` and `{17, …, 34}` are the negative-integer cycles of `3x + 1` (Kontorovich–Lagarias §1.1; Gonçalves–Greenfeld–Madrid Remark 1.7; THM-4590's negative basins).
     * So no argument built only from these properties can prove that every positive integer reaches the trivial cycle, for any `p`, including `p = 2`. Such an argument must distinguish `3x + 1` from `3x − 1`. This is the contracting, same-family form of Akin's observation (2004, §3) that the 2-adic conjugacy method "is too general": it works for every odd multiplier, and `5x + 1` has extra cycles.
     * In the base-`p` family the unconjugated maps show it too. `C_p` itself has extra positive-integer cycles for `p = 3, 11, 17, 23, 29, 31`, such as minimum 7 and length 9 for `p = 3`, and minimum 642 and length 57 for `p = 11` (Carnielli's Table 1; re-verified here for starts `≤ 10^6`, results note §3).
6. **Numerics on `Z_7`, `Z_11`, `Z_13` (NUMERICAL).** See the results note: merge-time tails, cycles of `C_p` on positive integers, and the repunit lines.

## Proofs

**1.**
* (a), (b) come from writing `e = A P^(−max(0,−k))` and expanding each branch of the table.
  * **`k ≥ 1`.**
    * `(≠, ≠)`: `f′ = (P/p) f + [(p − j)/P^k − (p − i)]/p`.
    * `(0, ≠)`: away, `f′ = f/p + (p − j)/(p P^(k+1))`.
    * `(≠, 0)`: toward, `f′ = (P/p) f − (p − i)/p`.
    * `(0, 0)`: `f′ = f/p`.
  * **`k ≤ −1`** (`f = e`; the reference digit is `u`'s digit `j`).
    * `(≠, ≠)`: `e′ = (P/p) e + [(p − j) − (p − i)/P^h]/p`.
    * `(0, ≠)`: toward, `e′ = (P/p) e + (p − j)/p`.
    * `(≠, 0)`: away, `e′ = e/p − (p − i)/(p P^(h+1))`.
    * `(0, 0)`: `e′ = e/p`.
  * **`k = 0`.**
    * Departures give `f′ = e/p ± (p − ·)/(pP)`.
    * The non-flip `(≠, ≠)` gives `e′ = (P/p) e + (i − j)/p`.
  * Every additive term is at most `(p − 1)/p` in absolute value.
* (c) Given `F_n`, `β_n` is uniform and `A_n` is known, so `j_n = (β_n + A_n) mod p` is uniform too.
  * If `p ∤ A_n`, exactly one value of `β_n` makes `i = 0` (a flip, `j ≠ 0`) and exactly one makes `j = 0` (a flip, `i ≠ 0`).
  * The other `p − 2` values give both digits nonzero.
* (d) If `p | A` and `k = 0`, then `i = j`. The branch is `e/p` or `(Pe + i − i)/p = Pe/p`, and `v_p(P) = 0`.
* (e) At `k ≥ 1` with `p | A` we have `i = j`.
  * The nonzero digits give `A′ = [PA − (p − i)(P^k − 1)]/p`, and the zero digit gives `A/p`.
  * Since `p − i` is a unit, `v_p((p − i)(P^h − 1)) = μ_h`.
  * **`w < μ_h`.** The first term dominates, so `v_p(A′) = w − 1` for every digit.
  * **`w > μ_h`.** The `p − 1` nonzero digits give `μ_h − 1`, while the zero digit gives `w − 1 ≥ μ_h`.
  * **`w = μ_h`.** Write `A = p^w a` and `P^h − 1 = p^w b` with `a`, `b` units.
    * Valuation `≥ μ_h` survives iff `p − i ≡ a b^(−1) mod p`, which happens for exactly one nonzero digit.
    * The zero digit gives `w − 1 < μ_h`. So again exactly `p − 1` digits leave.
  * `k ≤ −1` is the same computation with `[PA + (p − i)(P^h − 1)]/p`.
  * `μ_h = 1 + v_p(h)` is the lifting-the-exponent lemma for `(1 + p)^h − 1`, `p` odd. ∎

**2.**
* **Level `h ≥ 1`, `p ∤ A`.** By 1(a)–(c) and subadditivity of `x ↦ x^θ`,
  `E[V′] ≤ [(1/p) p^(−θ) s^(h+1) + (1/p)(P/p)^θ s^(h−1) + ((p−2)/p)(P/p)^θ s^h] |f|^θ + c^θ s^(h−1) = g(s) V + c^θ s^(h−1)`.
  * The additive weight is at most `s^(h−1)` because the next level is at least `h − 1`.
* **Level `h ≥ 1`, `p | A`.** No flips occur, so `E[V′] ≤ κ(θ) V + c^θ s^h`.
* **Level 0, `p ∤ A`.** Departures (probability `2/p`) pay `1/p` and gain the weight `s`; non-flips (probability `(p − 2)/p`) cost `P/p`. Also `(2/p) p^(−θ) s + ((p − 2)/p)(P/p)^θ ≤ κ(θ) − (1/p)((P/p)^θ − p^(−θ)) < κ(θ)`.
* **Level 0, `p | A`.** This is exact by 1(d): `E[V′] = κ(θ) V`.
* **Return drift.**
  * Fix `j` and let `σ = τ_(j+1)`. The first step gives `E[V_(τ_j + 1)] ≤ λ_0 V_(τ_j) + c^θ`.
  * After it, `W_n = V_n − Σ_(m < n) a_m`, with `a_m` the additive term of step `m`, is a supermartingale until `σ`. It need not be nonnegative; `V` is.
  * `σ < ∞` almost surely: in the stay at level 0, every step with `p ∤ A` departs with probability `2/p`, and runs are finite by 1(d); the excursion is a recurrent SRW (statement 3(i) or directly).
  * Optional stopping at `σ ∧ N` gives `E[V_(σ∧N)] ≤ λ_0 V_(τ_j) + c^θ + E[Σ_(m < σ) a_m]`. Fatou, applied to `V ≥ 0` as `N → ∞`, gives `E[V_σ] ≤ λ_0 V_(τ_j) + c^θ + E[Σ additive terms]`.
  * **Additive terms of the stay at level 0.** At most `c^θ` per step with `p ∤ A`, and the number of such steps before departure is `Geom(2/p)`; runs at level 0 carry no additive term.
  * **Additive terms of the excursion.**
    * Its flip skeleton is a fair simple random walk from 1, killed at 0: by 1(b, c), given a flip, its direction is decided by a fresh digit with probability 1/2 each (optional skipping). So the expected number of visits to level `h ≥ 1` is `G(1, h) = 2`.
    * A visit lasts at most `L_h = (p/2 + 1)(μ_h + p/(p − 1))` steps in expectation: there are `Geom(2/p)` steps with `p ∤ A`, each followed by at most one valuation run of expected length at most `μ_h − 1 + p/(p − 1)` (1(e)).
    * So the expected additive total is at most `C_exc = Σ_(h≥1) 2 L_h c^θ s^(h−1) < ∞`, because `s < 1` and `μ_h ≤ 1 + log_p h`.
  * At `σ` the level is 0, so `V_σ = |e_σ|^θ`. ∎

**3.**
* **(i) Recurrence.**
  * `k` is a martingale with increments in `{−1, 0, 1}` (1(c)). Almost surely it either converges or has `limsup = +∞` and `liminf = −∞` (Durrett, Thm 4.3.1).
  * If it converges, there are finitely many flips. By Lévy's extension of the Borel–Cantelli lemma (Durrett, Thm 4.3.4) applied to `{flip at n}`, whose conditional probability is `(2/p) 1{p ∤ A_n}`, it follows that `p | A_n` for all `n ≥ N` for some `N`.
  * Then `j_n = β_n` for all `n ≥ N`. So `u_N` and `v_N` have the same digit sequence, and `u_N = v_N` by injectivity.
  * If `k_N = 0` this is absorption. Otherwise `v_N = −e_N/(P^(k_N) − 1)`, an `F_N`-measurable value, which `v_N` (Haar given `F_N`) takes with probability 0.
  * So off absorption, `k` returns to 0 infinitely often.
* **(ii) Tightness at returns.**
  * If `k_0 ≠ 0`, first condition on the first arrival state at `τ_0` (finite by (i)).
  * By 2 and induction, `E[|e_(τ_j)|^θ | F_(τ_0)] ≤ λ_0^j |e_(τ_0)|^θ + C_p(θ)/(1 − λ_0) ≤ |e_(τ_0)|^θ + C_p(θ)/(1 − λ_0) < ∞` for every `j`. Conditional Fatou gives `liminf_j |e_(τ_j)| < ∞` almost surely.
  * At `k = 0`, `e = A` is an integer. So almost surely some bounded set of states `(0, e)`, `|e| ≤ M`, is visited infinitely often.
* **(iii) Accessibility.** By 4, every state `(0, e)` with `|e| ≤ M` is absorbed within a fixed number `N_M` of steps with probability at least `δ_M = p^(−N_M) > 0`.
* **(iv) Conclusion.**
  * By Lévy's 0–1 law, `P(absorption | F_n) → 1_absorption` almost surely.
  * On `{|e_(τ_j)| ≤ M infinitely often}` the left side is at least `δ_M` infinitely often, so absorption occurs almost surely there.
  * Taking the union over `M` proves almost-sure absorption. ∎

**4.** Use the table.
* **`e ≥ 2`, `p ∤ e`.**
  * The digit `i = p − ē` gives `j = 0`, a departure to `k = −1` with `A′ = (Pe − ē)/p ≥ (2P − p + 1)/p > 0`.
  * At `k = −1` with `p | A`, digit 0 gives `(0, 0)`-type steps `A ↦ A/p`. Let `A* ∈ [1, A′]` be the first value with `p ∤ A*`.
  * Digit 0 then has `j = A* mod p ≠ 0`, a toward step to `k = 0` with `A″ = (A* + p − j)/p = ⌈A*/p⌉`.
  * So `1 ≤ A″ ≤ (A′ + p − 1)/p ≤ (Pe − 1 + p² − p)/p²`. This is `< e` iff `e(p² − p − 1) > p² − p − 1`, i.e. `e ≥ 2`.
* **`e ≤ −2`.**
  * Digit 0 has `j = ē ≠ 0`, a departure to `k = 1` with `A′ = (Pe + p − ē)/p < 0`.
  * Digit 0 then divides by `p` while `p | A`, ending at `A* ∈ [A′, −1]`.
  * The digit `i = (−A*) mod p ≠ 0` gives `j = 0`, a toward step with `A″ = (A* − (p − i))/p = ⌈A*/p⌉ − 1 ≤ −1`.
  * Also `A″ ≥ A′/p − 1 ≥ (Pe + 1)/p² − 1`. This is `> e` iff `|e|(p² − p − 1) > p² − 1`, which holds for `|e| ≥ 2` when `p ≥ 3`.
* **`e = ±1`** is direct from the table.
* Strong induction on `|e|` finishes the proof. ∎

**5.**
* (a) is 3 at `(k_0, A_0) = (0, e)`. The density statement uses Fact D: the event depends only on `y mod p^K`, and `C` restricted to `Z` is the integer map.
* (b) For `p` odd, `n ↦ (1 + p)^n` is a measure isomorphism `Z_p → 1 + pZ_p` (the `p`-adic exponential). So `x_E = p P^(E−1) − 1` is Haar on the coset `p − 1 + p²Z_p`, which has positive measure, and 3 applies there (conditioning on a positive-measure set preserves "almost surely"). Absorption by time `n` depends only on the reference digits `β_0, …, β_(n−1)` of `x_E`, hence on `x_E mod p^n`, hence on `P^(E−1) mod p^(n−1)`, hence on `E mod p^(n−2)`, since `1 + p` has order `p^(n−2)` modulo `p^(n−1)`.
* (c) If `k_0 ≠ 0`, the chain starts away from level 0.
  * If `k_0 = 0` and `e_0 ≠ 0`: after the exact run of 1(d) the offset is nonzero and `p ∤ A`, so the next step departs with probability `2/p`. Absorption needs `e = 0` at level 0, so it cannot happen first.
  * From level `±1`, absorption needs the flip skeleton to return to 0. A simple random walk from 1 avoids 0 for `T` flips with probability at least `c′ T^(−1/2)`, and `T` steps contain at most `T` flips.
  * So `P(no merge by T) ≥ (2/p) c′ T^(−1/2)`.
* (d) At `p = 2`: `(p − 2)/p = 0`, `g(s) = ½((1/2)^θ s + (3/2)^θ s^(−1))`, so `g(s) = 1` is THM-4581's equation for `s`, and `c = 1/2`. The valuation-run lemma 1(e) holds with `μ_h = v_2(3^h − 1)` (lifting the exponent at 2); `λ_0` comes out weaker than THM-4581's, which does not matter for 3.
* (e) If `C^t(n) = C^t(n + e)` for some `t`, then `n` and `n + e` lie in the same grand orbit, so in `B` or outside it together.
  * Whether this happens by time `K` depends only on `n mod p^K` (Fact D; `C` on `Z` is the restriction).
  * So the density of `{n : n, n + e not merged by time K}` is exactly `q_e(K) = P(no merge by K)` from the start `(0, e)`.
  * Hence the upper density of the boundary is at most `q_e(K)` for every `K`, and `q_e(K) → 0` by 3.
* **Remark.** `x ↦ −x` preserves Haar measure, commutes with the pair relation up to `e ↦ −e`, and conjugates `C_p` to `F_p`; `F_p(i) = ⌊i + i/p⌋ = i` for `1 ≤ i ≤ p − 1`. The cycles of `3x − 1` and of `C_p` are verified by direct iteration (`zp_basins_survey.py` (S1), (S4)). ∎

## What this does and does not say

* It proves the almost-sure part of HYP-9244's requirement 1 (`debt_rank_least_requirements_20261008.md` §6 item 1) for the maps `C_p` on prime bases. The `T^(−1/2)` upper rate and other rank-one maps remain open.
* **The same proof for general two-valued maps (sketch only).** It applies to any map with multipliers `a` on a residue set `A ∋ 0` (with `r_0 = 0`) and `b` elsewhere, `a ≡ b ≡ 1 mod p`, contracting, *given accessibility*.
  * Replace `1/p, P/p` by `a/p, b/p`.
  * The fair skeleton, the factor law and the valuation runs (with `μ_h = v_p(b^h − a^h)`) go through.
  * Accessibility is the only map-specific input, and the obstruction lemma of HYP-9244 shows it can fail.
* It is a measure statement on `Z_p`. Its integer forms are density-one statements. It says nothing about particular integers, nor about whether integer orbits of `C_p` reach a cycle (Möller's conjecture for `m = p + 1`).
* For `p = 11` a second cycle exists on positive integers (minimum 642, length 57; Carnielli's Table 1, results note), so "every integer reaches 1" is false there. Yet translates still coalesce almost surely in the Haar sense.
* **Prior art.**
  * The maps are Carnielli's `T_d` (Applied Mathematics E-Notes 15 (2015) 207–215, with a table of cycles for `d ≤ 150`, `|x| ≤ 50000` computed by K. R. Matthews). Bruschi (arXiv:0810.5169) studies the undivided form `p C_p`. Tremblay (arXiv:1907.13302) locates their cycles. They lie in Hasse's class and Möller's family (Möller 1978: eventual periodicity conjectured for `m < d^(d/(d−1))`), and in the Matthews–Watts class, whose maps preserve Haar measure and are strongly mixing on `Z_p` (Matthews–Watts 1984, Corollary 1, as quoted in Gonçalves–Greenfeld–Madrid Remark 1.6).
  * arXiv:2111.06170 (Gonçalves–Greenfeld–Madrid) extends Tao's almost-bounded-orbits theorem to a class containing them.
  * No coalescence (merging) statement for them was found in two searches.
  * The proof follows THM-4581's (mac-mini) architecture, along the route of HYP-9244's requirement-1 sketch and audit G's drift check. mac-mini's `rank1_two_valued.out` row `Z_3` (1, 4, 4) is `C_3`; its tail (`√T q` = 7.85, 9.09, 8.11 at `T` = 1024, 4096, 16384) agrees with statement 6's.
* **Audit (2026-10-08, the session's independent audit).** Statements 1–5: proofs sound and complete as elementary probability. Corrected after it: the reason for `κ(1) = 1`, the supermartingale wording, 5(b) stated for `x_E` with modulus `p^(n−2)`, the `p = 2` specialization, 5(f) re-typed as the known `3x − 1` barrier, citations (Carnielli, Bruschi), credits (HYP-9244's sketch, audit G, mac-mini's `C_3` row), and the scope of the claim (MISTAKE-649).
