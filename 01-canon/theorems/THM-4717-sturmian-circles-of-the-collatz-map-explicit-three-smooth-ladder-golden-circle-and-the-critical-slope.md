---
id: THM-4717
title: "Sturmian circles of the Collatz map: (1) the consecutive differences of the sorted balanced free cycle of 3x + (2^q - 3^p) are exactly d_i = 2^i 3^(p-1-floor((i+1)p/q)), i = 0..q-2 (every power of 2 once; the 3-exponents are a Beatty sequence), the step from intercept t to t+1 moving the one at k_1 = (q-1-t) p^-1 mod q; (2) for every irrational slope a in (0,1) with a != log_3 2, the balanced RATIONAL 3x+1 cycles of slopes p/q -> a converge elementwise to a circle x_a(theta), theta in R/(1/a)Z: x_a(theta) = (1/2) sum_{r>=0} 3^r 2^(-floor(r/a + theta)) > 0 if a < log_3 2 and x_a(theta) = - sum_{m>=1} 2^(ceil(m/a - theta) - 1) 3^(-m) < 0 if a > log_3 2, on which the Sturmian-coded Collatz step (letter 1 iff theta >= 1/a - 1) is the rotation theta -> theta + 1, i.e. a minimal rotation of rotation number a; the circle is a strictly monotone pure-jump image, a Cantor set of measure zero; at the critical slope a = log_3 2, the slope of every integer cycle (THM-4633), neither series converges; (3) the golden circle a = 1/phi has c_min = sum 3^r 2^(-ceil((r+1)phi)) = 10.3701271417470491618... and c_max = (1/2) sum 3^r 2^(-floor(r phi)) = 31.6103814252411474854..., the Fibonacci-shape cycles (F_(n+1), F_n) converging alternately (x_min - c_min = 1.4e-7 at (987, 610)), and it exists because log_2 3 = 1.58496 < phi = 1.61803 (ratio of terms 3/2^phi = 0.97734); in THM-4528's coding (1 -> 10) the golden Terras circle is the Fibonacci word of slope 1/phi^2"
status: >
  PROVED: statements 1, 2, 3 (statement 1 from THM-4716's intercept-shift lemma with the explicit index computation; statement 2
  by exact summation of the carry formula and a telescoping identity; the Cantor claim from the pure-jump structure). FINITE-EXACT:
  statement 1 verified for all 1101 coprime (q,p) with q <= 60 including the intercept index; statement 2 checked at the slopes
  sqrt2 - 1, 1/phi, pi - 2.5, 1/sqrt2 along their convergents to q ~ 240 and the rotation conjugacy at 40 random points per slope
  (sturmian_circles.out, golden_limit.out); the golden conjugacy at 300 random points to 1e-40, the constants to 45 digits
  (c_min = 10.3701271417470491618232345771410338827702286, c_max = 31.6103814252411474854697037314231016483106859) and the
  pure-jump identity (total variation 21.240254283494098324 = sum of the term drops) in golden_circle.out. New sequences (not in OEIS 2026-10-19): the
  Fibonacci clocks 2^(F_(n+1)) - 3^(F_n) = 1, -1, 5, 13, 1631, 502829, 6719515981, ... and the golden Collatz minima 1, 5, 23, 319,
  14501, 5805215, 68498816767, .... Not claimed: anything about integer 3x+1 cycles beyond THM-4716 (the critical slope is exactly
  where the construction fails), or priority over Halbeisen-Hungerbuehler's rotation cycles (the real limit circles and the
  constants were not found in the literature by this session). Pre-audit.
session: mac-mini-2026-10-19-duality
source: 05-knowledge/results/duality_20261019_multiplication_addition_golden_circle_gilbreath.md
scripts:
  - 04-computation/experiments/duality_20261019/ladder_formula.py (+ .out, .json), golden_limit.py (+ .out), golden_circle.py (+ .out), sturmian_circles.py (+ .out), critical_slope.py (+ .out)
related:
  - THM-4716 (balanced cycles are 3-smooth ladders; intercept-shift lemma), THM-4715, THM-4714 (Fibonacci bulbs; the golden pair of Farey addresses), THM-4484 (free cycles), THM-4633 (integer cycles have slope at the convergents of log_3 2)
  - THM-4528 (parity words in base phi; golden beta-map; Fibonacci word), THM-4590 (golden Collatz resonance), THM-4656 (golden Wieferich)
  - Halbeisen-Hungerbuehler, Acta Arith. 78 (1997) (rotation cycles), Lagarias 1985 (parity vectors), Denjoy (exceptional minimal sets of circle maps)
---

# THM-4717 — Sturmian circles of the Collatz map

## Setting

* `T(x) = x/2` (even), `(3x+1)/2` (odd); words are Terras parity words; for coprime `(q, p)` the balanced necklace of slope `p/q` has the `q` rotations `s^(t)`, `s^(t)_k = ⌊((k+1)p + t)/q⌋ − ⌊(kp + t)/q⌋` (intercept `t/q`, `0 ≤ t < q`), with carries `B_t` (THM-4716). The clock is `D = 2^q − 3^p`. The **rational cycle** of `3x+1` with word `s^(t)` has elements `x_t = B_t/D`; the free cycle of `3x + D` has elements `B_t` (THM-4484).
* A **slope** is an irrational `a ∈ (0,1)`; its convergents `p_n/q_n` are the Fibonacci ratios `F_n/F_(n+1)` for `a = 1/φ`.

## Statements

1. **Explicit ladder (PROVED).** The consecutive differences of the sorted carries `B_(q−1) < … < B_0` are exactly the `q − 1` numbers

       d_i = 2^i · 3^(p − 1 − ⌊(i+1)p/q⌋),      i = 0, 1, …, q − 2,

   each power `2^i` occurring once; the step from intercept `t` to `t + 1` (`0 ≤ t ≤ q − 2`) realises `d_i` with `i = k_1 − 1`, `k_1 = (q − 1 − t)·p^(−1) mod q`. The 3-exponents `p − 1 − ⌊(i+1)p/q⌋` form the (reversed) lower Beatty sequence of slope `p/q`.

2. **Sturmian circles (PROVED).** Let `a ∈ (0,1)` be irrational, `a ≠ log_3 2`, and let `p_n/q_n → a` be its convergents. Write `θ = t/p ∈ [0, q/p)` for the intercept. Then, as `n → ∞` with `t/p_n → θ`, the rational cycle elements converge:
   * if `a < log_3 2` (so `3 < 2^(1/a)`): `x_t → x_a(θ) = (1/2) Σ_(r≥0) 3^r 2^(−⌊r/a + θ⌋)`, positive;
   * if `a > log_3 2`: `x_t → x_a(θ) = − Σ_(m≥1) 2^(⌈m/a − θ⌉ − 1) 3^(−m)`, negative.

   In both cases, with the rotation `R(θ) = θ + 1 mod 1/a` of the circle `R/(1/a)Z` and the letter `s(θ) = 1` iff `θ ≥ 1/a − 1`,

       x_a(R θ) = (3 x_a(θ) + 1)/2  if s(θ) = 1,      x_a(R θ) = x_a(θ)/2  if s(θ) = 0,

   so the Sturmian-coded Collatz step on the circle `X_a = x_a(R/(1/a)Z)` is conjugate to the irrational rotation of rotation number `a` (minimal, uniquely ergodic). `x_a` is strictly monotone (decreasing in the subcritical case) and a pure-jump function with jumps at the dense set `{j − r/a}`, so `X_a` is a Cantor set of Lebesgue measure zero in `[x_a(1/a^−), x_a(0)]` (resp. the mirror interval), the Collatz counterpart of a Denjoy exceptional minimal set. At `a = log_3 2` the geometric ratios `3/2^(1/a)` and `2^(1/a)/3` are both 1 and neither series converges: the critical slope is exactly the slope of every integer `3x+1` cycle (THM-4633).

3. **The golden circle (PROVED; constants FINITE-EXACT to 40 digits).** For `a = 1/φ` (`1/a = φ`): `c_min = x(φ^−) = Σ_(r≥0) 3^r 2^(−⌈(r+1)φ⌉) = 10.37012714174704916182323457714103388276…` and `c_max = x(0) = (1/2) Σ_(r≥0) 3^r 2^(−⌊rφ⌋) = 31.6103814252411474854697037314231016483…`, with `c_max/c_min = 3.0482154165677630541…`. The Fibonacci-shape cycles `(q, p) = (F_(n+1), F_n)` converge alternately: `x_min − c_min = 14.2, −1.48, 1.18, −0.18, 0.19, −0.0089, 0.026, −3.1e−5, 8.0e−4, −2.5e−11, 1.4e−7` for `n = 5, …, 15`. The series converge because `3 < 2^φ`, i.e. `log_2 3 = 1.58496 < φ = 1.61803`, the term ratio being `3/2^φ = 0.977337…`. Under the substitution `1 → 10, 0 → 0` (Terras word to THM-4528's `n/2, 3n+1` word) the golden Terras word of slope `1/φ` becomes the Fibonacci word of slope `1/φ²`, THM-4528's golden-mean word; `Θ` of that word is its base-`φ` value.

## Proofs

**1.** THM-4716 (1)–(2): the step from intercept `t` to `t+1` moves the one at `k_1` to `k_0 = k_1 − 1` and changes the carry by `2^(k_0) 3^j`, where `j` counts the ones to the right of `k_1`. The changing index is the unique `k_1 ∈ {1, …, q−1}` with `k_1 p + t ≡ q − 1 (mod q)`, i.e. `k_1 = (q − 1 − t) p^(−1) mod q`; as `t` runs over `0, …, q−2`, `k_1` runs over `1, …, q−1` once each, so `i = k_1 − 1` runs over `0, …, q−2` once each. Writing `k_1 p = qu + r` with `1 ≤ r ≤ q−1` gives `t = q − 1 − r` and the number of ones at positions `≤ k_1` equals `⌊((k_1+1)p + t)/q⌋ = ⌊(qu + r + p + q − 1 − r)/q⌋ = u + 1`, so `j = p − 1 − u = p − 1 − ⌊k_1 p/q⌋ = p − 1 − ⌊(i+1)p/q⌋`. ∎

**2.** The `m`-th one of `s^(t)` sits at position `⌈(mq − t)/p⌉ − 1`, so `B_t = Σ_(m=1)^p 3^(p−m) 2^(⌈(mq − t)/p⌉ − 1)`. With `r = p − m` and `θ = t/p`, `⌈((p − r)q − t)/p⌉ = q − ⌊(rq + t)/p⌋ = q − ⌊r q/p + θ⌋` (non-integrality for `1 ≤ r ≤ p−1`; the term `r = p` is `2^(q−1)` only when `t = 0` and is otherwise absorbed), hence

    x_t = B_t/D = (1/2) Σ_(r=0)^(p−1) 3^r 2^(−⌊r q/p + θ⌋) / (1 − 3^p 2^(−q)).

For `a < log_3 2`, `3^p 2^(−q) = 2^(−q(1 − a·log_2 3)(1+o(1))) → 0`, and the terms converge for each fixed `r` (`q/p → a^(−1)`, `r/a + θ ∉ Z` for irrational `a` and generic `θ`) with the dominating bound `3^r 2^(−r/a + 1)` (summable as `3 < 2^(1/a)`), so `x_t → x_a(θ)`. For `a > log_3 2` divide numerator and denominator by `3^p` and index from the other end, `m = p − r`: `B_t/3^p = Σ_(m=1)^p 3^(−m) 2^(⌈(mq − t)/p⌉ − 1)` and `D/3^p → −1`, giving the second formula. Conjugacy: shifting the word by one letter sends intercept `t` to `t + p`, i.e. `θ` to `θ + 1`, reduced modulo `q/p → 1/a`; the letter is `s_0 = ⌊(p + t)/q⌋ = ⌊(1 + θ)a⌋ = 1` iff `θ ≥ 1/a − 1`. If `s = 0` then `θ + 1 < 1/a` and `⌊r/a + θ + 1⌋ = ⌊r/a + θ⌋ + 1`, so `x_a(θ + 1) = x_a(θ)/2`. If `s = 1` then `θ' = θ + 1 − 1/a ∈ [0, 1)` and `x_a(θ') = (1/2) Σ_r 3^r 2^(−⌊(r−1)/a + θ⌋ − 1)`; the `r = 0` term is `(1/4)·2^(−⌊θ − 1/a⌋) = 1/2` (as `θ − 1/a ∈ [−1, 0)`), and the rest is `(3/2)·x_a(θ)`: total `(3x_a(θ) + 1)/2`. (The supercritical case is the same computation on the second series.) The function `θ ↦ ⌊r/a + θ⌋` is a step function with one jump per unit; the sum of the jump sizes over a period equals the total variation, so `x_a` is a strictly decreasing pure-jump function whose range misses the union of its jump gaps (open intervals of total length equal to the variation): a measure-zero Cantor set. The rotation by `1` on `R/(1/a)Z` has rotation number `a`, irrational, hence minimal. At `a = log_3 2` the term ratios are `3/2^(1/a) = 1` and `2^(1/a)/3 = 1`. ∎

**3.** Specialise 2 to `a = 1/φ`; `x_a` is decreasing, so `c_max = x(0)` and `c_min = x(φ^−) = (1/2) Σ 3^r 2^(−(⌈(r+1)φ⌉ − 1))`. The numerical values use 6000 terms (ratio `0.977`, error `< 10^(−55)`). The substitution statement is the definition of THM-4528's coding: an odd Terras step is `3n+1` followed by a halving. ∎

## Remarks

* **Duality.** Statement 1 says that the additive structure (consecutive differences) of a multiplicative-looking cycle is a window of the multiplicative semigroup `⟨2, 3⟩`, indexed additively by a Beatty sequence; statement 2 says the limit of these additive ladders is a circle on which multiplication by `3/2` and `1/2`, corrected by the carry `+1/2`, acts as a pure rotation. The constant `c_min` is a Lambert-type series over the lower Wythoff sequence, `(1/2) Σ 3^r 2^(−⌊(r+1)φ⌋)`.
* **Why the golden ratio.** The golden slope is subcritical by `log_3 2 − 1/φ = 0.0129`; the Fibonacci clocks `2^(F_(n+1)) − 3^(F_n)` grow like `2^(F_(n+1))(1 − 2^(−0.02 F_(n+1)))`, positive from `n = 3` on, which is why the Fibonacci bulbs of THM-4714 beyond the `2/3`-bulb carry positive free cycles (`3x + 5`, `3x + 13`, `3x + 1631`, …).
