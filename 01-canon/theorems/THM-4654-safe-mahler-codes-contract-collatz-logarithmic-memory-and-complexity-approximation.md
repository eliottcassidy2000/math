---
id: THM-4654
title: "Safe Mahler codes contract the Collatz maps: the slack automaton, the block law (no 11; a gap 2 forces a gap >= 6), at most floor(l/3)+1 odd steps, slopes <= 3(3^(1/3)/2)^l, logarithmic memory of rational shadows, and a complexity-approximation inequality for the irrational Collatz sources selected by positive Mahler integers"
status: >
  PROVED: statements 1-7 (elementary; statement 2's general form is Parry's criterion, CITED). FINITE-EXACT: statement 2 to
  depth 60 and the exact counts; statement 7 on 4,980,199 pairs (A, u/v) with 0 violations; sharpness of statement 6 on the
  family A_t = (2^(18t+8) - 9)/19. The Mahler Z-number problem and Collatz remain OPEN.
source: mac-mini-2026-10-14-wilson
depends_on: [THM-2228 (strict carry tails and integral stabilization), mahler_carry_interface_20261014i.md (the isometry Psi with Psi M = T Psi, eq. (8)), THM-4072 (the safe follower graph)]
verification:
  - 04-computation/experiments/mahler_collatz_sources_20261014/safe_language.py (+ .out, .json)
  - 04-computation/experiments/mahler_collatz_sources_20261014/collatz_sources_check.py (+ .out)
note: 05-knowledge/results/mahler_collatz_sources_logarithmic_memory_20261014.md
related: [THM-2228, THM-4072, THM-4469, THM-4627, mahler_carry_interface_20261014i.md, recursion_carry_bridges_20261014i.md, MAHLER-THREE-HALVES-FRONTIER-2026-08-23.md]
---

# THM-4654: safe Mahler codes contract the Collatz maps

## Setting

- `M(x) = (3x + (x mod 2))/2` is the Mahler ceiling map, and `T(x) = x/2` or `(3x+1)/2` is the Collatz map, both acting on `Z₂` and on rationals of odd denominator.
- A finite carry word `w` is **safe** if every suffix sum `Σ_j w_(n+j) (2/3)^(j+1)` is `< 1`. This is the finite language of the strict safe shift of THM-2228 / THM-4072.
- An infinite code is strictly safe if every tail `Y_n < 1`.
- A hypothetical Mahler Z-number has an integer part `A >= 1` whose Mahler code `e` is strictly safe at all depths (THM-2228).
- `Ψ` is the isometry of `Z₂` with `ΨM = TΨ` that preserves codes, and `z_A = Ψ(A)` is the **Collatz source selected by `A`**.
  - In closed form, if the ones of `e` sit at positions `n₁ < n₂ < …`, then `A = −(1/3) Σ_m (2/3)^(n_m)` and `z_A = −Σ_m 2^(n_m) 3^(−m)`, both 2-adically.
  - The real value of the same series is `(3/2)Y₀`, so a Z-number is `ξ = A + Y₀/2`.
- Put `λ = 3^(1/3)/2 = 0.72112…` and `c = 1 − (log₂3)/3 = 0.47168…`, so that `λ = 2^(−c)`.

## Statement

1. **(Slack automaton.)** A word is safe iff it is accepted by
   `m ↦ min((3/2)m − e, 1)`, where each step is valid iff `(3/2)m − e > 0`, starting from `m₀ = 1`.

2. **(Block law.)**
   - No safe word contains `11`.
   - A gap 2 between consecutive ones (`101`) must be followed by a gap of at least 6.
   - A pair of gaps `(2, 6)` must be followed by a gap of at least 3.
   - In general, the gap sequence of every suffix that starts at a 1 must exceed, at the first place where they differ, the gap sequence `(2, 6, 3, 3, 2, 10, 7, 5, 7, 3, 3, 2, 12, …)` of the greedy base-3/2 expansion `1 = Σ d_j (2/3)^j`, `d = 1010000010010010100000000010…`.

3. **(Density.)** A safe word of length `ℓ >= 1` has at most `⌊ℓ/3⌋ + 1` ones, and this is sharp for every `ℓ`. The number of safe words of length `ℓ` is `~1.551 (3/2)^ℓ`.

4. **(Uniform contraction of the Collatz affine maps.)** Along a safe word `w` of length `ℓ` with `k` ones, `T^ℓ(y) = (3^k y + C_w)/2^ℓ`. Then
   - `3^k/2^ℓ <= 3^(⌊ℓ/3⌋+1)/2^ℓ <= 3λ^ℓ`;
   - `0 <= T_w(0) = C_w/2^ℓ < 11/10`, and `11/10` is the supremum.

5. **(Entry window and logarithmic memory.)** Let `x` be rational with odd denominator `v` whose `T`-orbit follows a safe word for `n` steps.
   - **Entry window.** `T^n(x)` lies in `[−3λ^n|x|, 3λ^n|x| + 11/10)`. Hence `T^n(x) ∈ (1/v)Z ∩ [−1/2, 8/5)` as soon as `n >= n₁(x)`, where `n₁(x) = ⌈log₂(6|x|)/c⌉` (and `0` if `|x| <= 1/6`).
   - **Logarithmic memory.** Let `x, x' ∈ (1/v)Z` with `|x − x'| <= R` follow a common safe word of length `ℓ` with `3λ^ℓ R < 1/v`. Then `T^ℓ(x) = T^ℓ(x')`.
   - Inside the window this holds as soon as `ℓ >= ℓ(v) := ⌊log₂(6.3v)/c⌋ + 1 ≈ 2.12 log₂ v`. So a rational shadow's state is a **sliding-block function of the last `O(log v)` code bits**.

6. **(Periodic runs; the interface's eq. (8), recorded in this form.)** Suppose the Mahler code of `A >= 1` agrees on `[i, i + h)` with a purely periodic word `w^∞` with `|w| = p`. Then

   `h < p log₂3 + i log₂(3/2) + log₂(A + 1)`.

   This is sharp to within `log₂(27/19) = 0.507` bits: for `A_t = (2^(18t+8) − 9)/19` against `1/5`, the agreement is exactly `18t + 8`.

7. **(Complexity–approximation inequality.)** Let `A >= 1` and `x = u/v`, `v` odd. Suppose the Mahler code `e` of `A` and the Collatz code of `x` agree on `[0, h)`, and `e[0, h)` is safe. Let `N_e(ℓ)` be the number of distinct length-`ℓ` factors of `e[0, h)`. Then

   `h < log₂3 · (n₁(x) + min(⌊2.1v⌋ + 1, ℓ(v) + N_e(ℓ(v)))) + log₂(A + 1)`.

   For a hypothetical Z-number's integer part this holds for every `h`, so for every rational `u/v`:

   `v₂(z_A − u/v) < log₂3 · (n₁(u/v) + min(⌊2.1v⌋ + 1, ℓ(v) + N_e(ℓ(v)))) + log₂(A + 1)`.

### Consequences

- **(a) Unconditional.** `v₂(z_A − u/v) < 3.33v + 3.36 log₂(6|u|/v + 1) + log₂(A + 1) + 3.2`. This improves the interface's `14·max(|u|, v) + 2 + bitlength(A + 1)`: the numerator now enters only logarithmically.
- **(b) Polynomial complexity.** If the code of `A` has factor complexity `N_e(ℓ) <= Kℓ^s`, then `|z_A − u/v|₂ >= c_A exp(−C(log H)^s)` with `H = max(|u|, v)`.
- **(c) Linear complexity.** This covers Sturmian, automatic and morphic codes. If `N_e(ℓ) <= Kℓ`, then `|z_A − u/v|₂ >= c_A H^(−μ)` with `μ <= 3.36(K + 2)`. So **`z_A` has a finite 2-adic irrationality exponent**.
- **(d) Dichotomy.** For a hypothetical Mahler Z-number, either the carry code has super-polynomial factor complexity, or the irrational Collatz source it selects is badly approximable by rationals to polynomial order.

## Proof

**1.**
- For a word of length `ℓ` and a suffix starting at `n` (the empty suffix included), let the slack be `r_n = (1 − S_n)(3/2)^(ℓ−n)`. Safety is `r_n > 0` for all `n`.
- Appending `e` sends each slack to `r ↦ (3/2)r − e` and adds a fresh empty suffix with slack 1.
- The minimum slack `m` is therefore updated by `m ↦ min((3/2)m − e, 1)`, because the map `r ↦ (3/2)r − e` is increasing. The word stays safe iff the new minimum is positive.
- `safe_language.py` checks this equivalence on all 131,071 words of length at most 16.

**2.**
- A 1 is allowed iff `m > 2/3`. After a 1, `m <= 1/2`, so `11` never occurs.
- A gap 2 needs `(3/2)((3/2)m − 1) > 2/3`, i.e. `m > 26/27` before the first 1. After the second 1, `m <= (9/4)·(1/2) − 1 = 1/8`.
- The next 1 needs `(3/2)^(g−1)/8 > 2/3`, i.e. `(3/2)^(g−1) > 16/3`, so `g − 1 >= 5`, i.e. `g >= 6`.
- The `(2, 6)` rule is the same computation one step further (exact output of `safe_language.py`).
- The general form is Parry's lexicographic criterion for β-expansions with `β = 3/2`. Because `d(1)` is infinite, the greedy and quasi-greedy expansions coincide.

**3.**
- Take a safe word with `k` ones and internal gaps `g₁, …, g_(k−1)`. Each `g_i >= 2`.
- A gap `g_i = 2` with `i < k − 1` is followed by `g_(i+1) >= 6`, so the pair contributes at least `8 >= 3 + 3 + 2`. Hence `Σ g_i >= 3(k − 1) − 1`.
- The word has length `ℓ >= 1 + Σ g_i >= 3k − 3`, so `k <= ℓ/3 + 1`.
- Sharpness: the words `(100)^r 101` and their truncations attain the bound. The DP of `safe_language.py` confirms `k_max(ℓ) = ⌊ℓ/3⌋ + 1` exactly for `ℓ <= 60`.

**4.**
- The slope bound is statement 3.
- For the offset, `T_w(0) = Σ_r 3^(r−1) 2^(−t_r)`, where `t₁ < t₂ < …` are the distances of the ones from the end of `w`.
- By statement 2, `t₁ >= 1` and `t_r >= 3r − 3` for `r >= 2`: every gap is at least 2, and a gap 2 at depth at least 2 is paired with the gap of length at least 6 that follows it in reading order.
- So `T_w(0) <= 1/2 + Σ_(s>=1) (3/8)^s = 1/2 + 3/5 = 11/10`, and the inequality is strict for finite words. The words `(100)^r 101` approach the supremum.

**5.**
- `T^n(x) = slope·x + T_w(0)` with the bounds of statement 4.
- `T` preserves `(1/v)Z` for odd `v`.
- Two lattice points following the same word move by the same affine map, so their difference is multiplied by `slope <= 3λ^ℓ`. A nonzero difference in `(1/v)Z` has absolute value at least `1/v`.

**6.**
- Put `a = M^i(A)`. Then `a + 1 <= (3/2)^i (A + 1)`, because `M(y) + 1 <= (3/2)(y + 1)`.
- Let `z_w = −C_M(w)/D` with `D = 3^p − 2^p` and `0 <= C_M(w) <= D`.
- Codes determine points and `M` is a 2-adic isometry on each parity class. So `h <= v₂(a − z_w) = v₂(Da + C_M(w))`, and `Da + C_M(w)` is a positive integer at most `D(a + 1) < 3^p (3/2)^i (A + 1)`.
- For `A_t` we have `19A_t + 9 = 2^(18t+8)`, so `h = log₂(19A_t + 9)`, against the bound `log₂(27(A_t + 1))`.

**7.**
- Suppose `h` exceeds the bound's integer part. Inside the agreement window, the `T`-orbit `x_n` of `x` follows `e`.
- By statement 5 (entry window), `x_n ∈ (1/v)Z ∩ [−1/2, 8/5)` for `n >= n₁`. That set has at most `⌊2.1v⌋ + 1` points.
- **State count.** Among the `⌊2.1v⌋ + 2` times in `[n₁, n₁ + ⌊2.1v⌋ + 1]`, two states coincide: `x_m = x_(m')` with `m < m' <= n₁ + ⌊2.1v⌋ + 1`.
- **Memory.** Among the `N_e(ℓ) + 1` times `n ∈ [n₁ + ℓ, n₁ + ℓ + N_e(ℓ)]`, with `ℓ = ℓ(v)`, two have the same preceding `ℓ`-word. By logarithmic memory (statement 5, window length 2.1) their states coincide.
- Either way, the rational orbit is periodic from some `m` with period `p = m' − m`. So `e[m, h)` is `p`-periodic, and `e[m, h)` equals `w^∞` with `w = e[m, m')`.
- Statement 6 with `i = m` gives `h − m < p log₂3 + m log₂(3/2) + log₂(A + 1)`, i.e. `h < m' log₂3 + log₂(A + 1)`. This contradicts the bound.
- The constants in (a) are `log₂3 · 2.1 = 3.33` and `log₂3/c = 3.36`.
- For (c): `ℓ(v) + N_e(ℓ(v)) <= (K + 1)ℓ(v) <= 2.12(K + 1) log₂(6.3v) + K + 1`. Exponentiate. ∎

## Checks (FINITE-EXACT)

| Check | Scope |
|---|---|
| Automaton equal to the suffix-sum definition | all 131,071 words of length `<= 16` |
| `N(ℓ)` (1, 2, 3, 5, 8, 12, 18, 27, 40, 60, 90, 134, …), with `N(ℓ)/(3/2)^ℓ → 1.551` | `ℓ <= 60` |
| `k_max(ℓ) = ⌊ℓ/3⌋ + 1` | `ℓ <= 60` |
| Greedy expansion of 1 and its gap sequence | as listed in statement 2 |
| Gap rules | as listed in statement 2 |
| Offset `< 11/10` and slope `<= 3λ^ℓ` | all 34,805 safe words of length 1..22 (max offset 1.099374) |
| Entry window | 15,615 checks on rationals of height `<= 60` |
| Statement 7 | all 981 positive integers `A <= 200,000` with safe depth `>= 20` (deepest 38: `A = 8101, 105109, 163493, 192293`), against every `u/v` with `v <= 45` odd and `|u| <= 12v`: 4,980,199 pairs, 0 violations, tightest ratio `h`/bound `= 0.629` |
| Statement 6 equality family `A_t` | `t = 0..5`: agreement 8, 26, 44, 62, 80, 98 against bounds `+0.51` |

## Reading

- **Real–2-adic duality.** Mahler's condition is real (all tails below 1). The Collatz source is 2-adic. This theorem is the bridge between them.
  - Real safety forces the Collatz affine maps along the code to contract in the **real** metric (slope `<= 3λ^ℓ`).
  - Real contraction turns every rational 2-adic shadow into a finite-memory object: logarithmic memory.
  - The positive-integer premise then cuts the shadow off by the periodic-run inequality.
- **What it does not do.** It neither proves nor refutes the existence of Z-numbers. Statement 7 needs only finite safety, so it is a terminating, checkable inequality.
- **The open dichotomy (d).** This is the new restriction. To exclude Z-numbers whose codes have low complexity, one now needs only a 2-adic transcendence-type input for `z_A`. Conversely, a Z-number with a Sturmian-like code would select a Collatz source with finite 2-adic irrationality exponent.
