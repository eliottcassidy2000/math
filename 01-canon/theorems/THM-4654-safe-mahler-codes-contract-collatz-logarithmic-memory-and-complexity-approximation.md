---
id: THM-4654
title: "Safe Mahler codes and their Collatz sources: the slack automaton, the block law (no 11; a gap 2 forces a gap >= 6), at most floor(l/3)+1 odd steps, Collatz slopes <= 3(3^(1/3)/2)^l (real entry into a bounded window), 2-adic memory of rational shadows, and a complexity-approximation inequality for the irrational Collatz sources selected by positive Mahler integers"
status: >
  PROVED: statements 1-7 (elementary; statement 2's general form is Parry's criterion, CITED). FINITE-EXACT: statements 1-4 to
  length 400 (audit D), statement 7 on 4,980,199 pairs (session) and on 2.1 million worst-case (A, u/v, h) cases (audit D,
  including |u| up to 10^15 v and adversarial codes), with 0 violations; sharpness of statement 6. KNOWN in framework: Mahler
  Z-numbers via the beta = 3/2 shift (Flatto 1992); the maximal 1-frequency 1/3 in principle (Boyland-de Carvalho-Hall 2016);
  Mahler codes of positive integers have factor complexity liminf N(l)/l >= log 2/log(3/2) (Dubickas 2009). Audited by audit D
  (MISTAKE-612): the pre-audit "logarithmic memory by real contraction" lemma was vacuous (the memory is 2-adic), and the
  pre-audit dichotomy was mis-stated. The Mahler Z-number problem and Collatz remain OPEN.
source: mac-mini-2026-10-14-wilson
depends_on: [THM-2228 (strict carry tails, integral stabilization; Mahler cylinders are residue classes), THM-4072 (the safe follower graph), THM-3848 (renewal law), mahler_carry_interface_20261014i.md (the isometry Psi with Psi M = T Psi, eq. (8)), Terras/Everett (Collatz parity vectors determine residues)]
verification:
  - 04-computation/experiments/mahler_collatz_sources_20261014/safe_language.py (+ .out, .json)
  - 04-computation/experiments/mahler_collatz_sources_20261014/collatz_sources_check.py (+ .out)
  - 04-computation/experiments/audit_D_20261014/ (independent checks of statements 1-7)
note: 05-knowledge/results/mahler_collatz_sources_logarithmic_memory_20261014.md
related: [THM-2228, THM-3848, THM-4072, THM-4469, THM-4627, mahler_carry_interface_20261014i.md, recursion_carry_bridges_20261014i.md, MAHLER-THREE-HALVES-FRONTIER-2026-08-23.md]
---

# THM-4654: safe Mahler codes and the Collatz sources they select

## Setting

- `M(x) = (3x + (x mod 2))/2` is the Mahler ceiling map. `T(x) = x/2` or `(3x+1)/2` is the Collatz map. Both act on `Z₂` and on rationals of odd denominator.
- A finite carry word `w` is **safe** if every suffix sum `Σ_j w_(n+j) (2/3)^(j+1)` is `< 1`. This is the finite language of the strict safe shift (THM-2228, THM-4072).
- A hypothetical Mahler Z-number has an integer part `A >= 1` whose Mahler code `e` is strictly safe at all depths (THM-2228).
- `Ψ` is the code-preserving isometry of `Z₂` with `ΨM = TΨ`. `z_A = Ψ(A)` is the Collatz source selected by `A`.
- With the ones of `e` at positions `p₁ < p₂ < …`:
  - `A = −(1/3) Σ_m (2/3)^(p_m)` and `z_A = −Σ_m 2^(p_m) 3^(−m)`, both 2-adically;
  - the real value of the same series gives the Z-number `ξ = A + (1/3) Σ_m (2/3)^(p_m)`.
- `λ = 3^(1/3)/2 = 0.72112…`, `c = 1 − (log₂3)/3 = 0.471677…`, so `λ = 2^(−c)`, `1/c = 2.12009…` and `log₂3/c = 3.36026…`.

## Statement

1. **(Slack automaton.)** A word is safe iff it is accepted by `m ↦ min((3/2)m − e, 1)`, valid iff `(3/2)m − e > 0`, from `m₀ = 1`.
   - This is THM-4072's safe follower graph in slack coordinates.
   - The number of safe words of length `ℓ` is `N(ℓ) ~ (3/Σ_j j d_j (2/3)^j)·(3/2)^ℓ = 1.55104519…·(3/2)^ℓ` (THM-3848's renewal law; Rényi–Parry entropy `log(3/2)`).

2. **(Block law.)**
   - No safe word contains `11`.
   - A gap 2 between consecutive ones must be followed by a gap of at least 6.
   - The gaps `(2, 6)` must be followed by a gap of at least 3.
   - In general (Parry, CITED), the gap sequence of every suffix starting at a 1 must exceed, at the first place where they differ, the gap sequence `(2, 6, 3, 3, 2, 10, 7, 5, 7, 3, 3, 2, 12, …)` of the greedy expansion `1 = Σ d_j (2/3)^j`.

3. **(Density.)** A safe word of length `ℓ >= 1` has at most `⌊ℓ/3⌋ + 1` ones, and this is sharp for every `ℓ`. (Maximal frequency 1/3 is known in principle; we did not find the exact finite bound stated.)

4. **(Contraction of the Collatz affine maps.)** Along a safe word `w` of length `ℓ` with `k` ones, `T^ℓ(y) = (3^k y + C_w)/2^ℓ`, and
   - `3^k/2^ℓ <= 3^(⌊ℓ/3⌋+1)/2^ℓ <= 3λ^ℓ`;
   - `0 <= T_w(0) < 11/10`, where 11/10 is the supremum.

5. **(Bounded entry; 2-adic memory.)** Let `x` be rational with odd denominator `v`, and suppose its `T`-orbit follows a safe word for `n` steps.
   - **Entry (real contraction).** `T^n(x) ∈ [−3λ^n|x|, 3λ^n|x| + 11/10)`. So `T^n(x) ∈ W_v := (1/v)Z ∩ [−1/2, 8/5)` for all `n >= n_in(x)`, where `n_in(x) = 0` if `|x| <= 1/6` and `n_in(x) = ⌈log₂(6|x|)/c⌉` otherwise. `W_v` has at most `⌊2.1v⌋ + 1` points.
   - **Memory (2-adic; Terras–Everett).** Two points of `(1/v)Z` whose `T`-codes agree for `ℓ` steps are congruent mod `2^ℓ` in `Z_(2)`, so they are equal or at real distance at least `2^ℓ/v`. Hence a state in `W_v` is determined by its next `ℓ(v) := ⌊log₂(2.1v)⌋ + 1` code bits.
   - Real contraction supplies only the entry into the bounded window. The memory is 2-adic.

6. **(Periodic runs; the interface's eq. (8).)** Suppose the Mahler code of `A >= 1` agrees on `[i, i + h)` with a purely periodic `w^∞`, `|w| = p`. Then

   `h < p log₂3 + i log₂(3/2) + log₂(A + 1)`.

   The additive constant is optimal. For `w = 0^(p−1)1` and `A = (2^h − 2^(p−1))/(3^p − 2^p)` (for those `h` that make `A` a positive integer), the gap tends to `log₂(3^p/(3^p − 2^p)) → 0` (0.011 at `p = 12`). For `A_t = (2^(18t+8) − 9)/19` against `1/5` (period `100`), the agreement is exactly `18t + 8`. The bound exceeds it by 0.562 at `t = 0`, decreasing to `log₂(27/19) = 0.507`.

7. **(Complexity–approximation inequality.)** Let `A >= 1` and `x = u/v`, `v` odd. Suppose the Mahler code `e` of `A` and the Collatz code of `x` agree on `[0, h)`, and `e[0, h)` is safe. Let `n_in = n_in(x)` and `ℓ = ℓ(v)`, and let `N` be the number of distinct factors `e[n, n+ℓ)` with `n_in <= n` and `n + ℓ <= h`. Then

   `h < min( log₂3·(n_in + ⌊2.1v⌋ + 1) + log₂(A+1),  max( n_in + N + ℓ,  log₂3·(n_in + N) + log₂(A+1) ) )`.

   For a hypothetical Z-number's integer part this holds for every `h`. So `v₂(z_A − u/v)` obeys the same bound with `N` replaced by the factor complexity `N_e(ℓ(v))` of the code.

### Consequences

- **(a) Unconditional.** `v₂(z_A − u/v) < 3.33v + 3.361·log₂(6|u|/v + 1) + log₂(A + 1) + 3.2`. Against the interface's `14·max(|u|, v) + 2 + bitlength(A+1)`, the numerator now enters only logarithmically.
- **(b) Polynomial complexity.** If `N_e(ℓ) <= Kℓ^s`, then `|z_A − u/v|₂ >= c_A exp(−C(log H)^s)` with `H = max(|u|, v)`. This is a quasi-polynomial bound, not a polynomial one, when `s > 1`.
- **(c) Linear complexity.** This covers automatic codes and primitive morphic codes (Pansiot), whose complexity is linear. Sturmian codes are impossible, since Dubickas (2009) gives `liminf N_e(ℓ)/ℓ >= log 2/log(3/2) ≈ 1.71` for the Mahler code of every positive integer. If `N_e(ℓ) <= Kℓ`, then `z_A` has finite 2-adic irrationality exponent: `|z_A − u/v|₂ >= c_A H^(−μ)` with `μ <= log₂3·(K + 1) + 3.361`.
- **(d) What this does and does not restrict.** For a hypothetical Z-number, either `limsup N_e(ℓ)/ℓ = ∞`, or the selected Collatz source `z_A` has a finite 2-adic irrationality exponent.
  - This is not by itself a restriction: a finite irrationality exponent is the generic, expected behaviour.
  - To exclude low-complexity Z-number codes one would need an infinite-exponent (Liouville-type) input for `z_A`, or an argument on the Mahler side, where the 2-adic value `−(1/3) Σ (2/3)^(p_m)` must be a positive integer.

## Proof

**1.**
- For a word of length `ℓ` and a suffix starting at `n` (the empty suffix included), the slack is `r_n = (1 − S_n)(3/2)^(ℓ−n)`.
- Appending `e` sends each slack to `(3/2)r − e` and adds a fresh slack 1.
- The map is increasing, so the minimum slack updates as stated.
- The equivalence was checked on all words of length at most 18 by four independent acceptors (audit D). The count is THM-3848's renewal law.

**2.**
- A 1 needs `m > 2/3`. After a 1, `m <= 1/2`, so `11` is impossible.
- A gap 2 needs `m > 26/27` before the first 1, and leaves `m <= 1/8` after the second.
- The next 1 then needs `(3/2)^(g−1) > 16/3`, so `g >= 6`.
- The `(2, 6)` rule is the next step of the same computation. The general form is Parry's criterion for `β = 3/2`; `d(1)` is infinite, so the greedy and quasi-greedy expansions coincide.

**3.**
- Internal gaps are all at least 2.
- A gap 2 that is not last is followed by a gap of at least 6, so each such pair contributes at least `8 >= 6 + 2`.
- Hence `ℓ >= 1 + Σ g_i >= 3k − 3`.
- The words `(100)^r 1`, `(100)^r 10` and `(100)^r 101` attain the bound.

**4.**
- The slope bound is statement 3.
- `T_w(0) = Σ_r 3^(r−1) 2^(−t_r)`, where `t_r` is the distance of the `r`-th one from the end. Statement 2 gives `t₁ >= 1` and `t_r >= 3r − 3` for `r >= 2`, with equality attained for `r = 2..8` (audit D).
- So `T_w(0) <= 1/2 + Σ_(s>=1)(3/8)^s = 11/10`, strictly for finite words. The words `(100)^r 101` approach this.

**5.**
- Entry: `T^n(x) = slope·x + T_w(0)`, with statement 4.
- Memory: the length-`ℓ` Collatz parity vector of `y ∈ Z_(2)` determines `y mod 2^ℓ` (Terras; Everett). The difference of two points of `(1/v)Z` is `(u − u')/v` with `v` odd.

**6.**
- Mahler cylinders are residue classes (THM-2228 §1), so agreement for `h` steps gives `v₂(a − z_w) >= h`, where `a = M^i(A)` and `z_w = −C_M(w)/D`, `D = 3^p − 2^p`, `0 <= C_M(w) <= D`.
- `Da + C_M(w)` is a positive integer at most `D(a + 1) < 3^p (3/2)^i (A + 1)`.
- The sharpness families are explicit.

**7.**
- **State count.** Among the `⌊2.1v⌋ + 2` times in `[n_in, n_in + ⌊2.1v⌋ + 1]`, two states in `W_v` coincide.
- **Factor count.** If `h >= n_in + N + ℓ`, then among the `N + 1` times `n ∈ [n_in, n_in + N]` two have the same factor `e[n, n+ℓ)`. Inside the agreement window these factors are the rational's own codes, so the two states coincide by statement 5's memory. Otherwise `h < n_in + N + ℓ`.
- In the first two cases the rational orbit is periodic from some `m` with period `m' − m`, so `e[m, h)` is periodic. Statement 6 gives `h < m' log₂3 + log₂(A + 1)`.
- For the constants: `log₂3·2.1 = 3.3284`, `log₂3/c = 3.36026`, and `2·log₂3 = 3.17`.
- For (c): `N_e(ℓ(v)) <= K(log₂(2.1v) + 1)`. Exponentiate. ∎

## Checks (FINITE-EXACT)

**Session.**
- Automaton equals definition on all words of length `<= 16`.
- `N(ℓ)` and `k_max(ℓ)` to `ℓ = 60`.
- Offset and slope on all 34,805 safe words of length `<= 22`.
- Entry window: 15,615 checks.
- Statement 7 (the corrected form) on all 981 positive `A <= 200,000` with safe depth `>= 20`, against `u/v` with odd `v <= 45` and `|u| <= 12v`: 4,980,199 pairs, 0 violations (`collatz_sources_check.out`).

**Audit D (independent).**
- Four acceptors agree on all 524,287 words of length `<= 18`.
- Statements 2–4 hold on all 396,551 safe words of length `<= 28`, with `N(ℓ)` and `k_max` to `ℓ = 400`.
- `sup T_w(0) = 11/10 − 4.6·10^(−29)` at length 200.
- Entry window: 214,575 rationals and 1.46 million orbit points.
- Statement 6: `10⁶` periodic runs, and `A_t` for `t <= 11`.
- Statement 7 for every positive `A` at once (the smallest `A` in each cylinder is the worst case):
  - 592,561 cases with `v <= 151`;
  - 1.45 million adversarial cases (`|u|` up to `10^15 v`, long zero runs, 1,953 periodic safe cycles).
  - 0 violations; minimum slack 4.75 bits.

## Reading

- **Three ingredients, three metrics.**
  - Real safety contracts the Collatz maps along the code in the **real** metric. That gives **boundedness**: every rational shadow enters a fixed window.
  - The **2-adic** parity bijection then makes the shadow's state a function of `O(log v)` future bits. That gives **finite memory**.
  - The **positive-integer** premise (statement 6) cuts off every periodic shadow.
  - The pre-audit text credited real contraction with the memory. That was wrong (MISTAKE-612).
- **What it does not do.** It neither proves nor refutes the existence of Z-numbers. Statement 7 needs only finite safety, so it is a terminating, checkable inequality.
