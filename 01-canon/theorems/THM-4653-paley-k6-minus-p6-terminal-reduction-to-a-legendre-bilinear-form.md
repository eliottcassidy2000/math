---
id: THM-4653
title: "The Paley six-vertex terminal K6 - P6 reduces exactly to one Legendre-family bilinear form M_q, equivalently the 2F1 x 3F2 beta moment -q^3 sum phi(x)phi(1-x) 2F1(x) 3F2(x)"
status: PROVED (elementary character-sum identities, valid over every finite field of odd characteristic); FINITE-EXACT checks listed below. Audited (audit C): CORRECT; every step re-derived by hand, V_q brute-forced with injectivity enforced on all 15 pairs at 24 fields (including fully unnormalised runs to q = 31), the proof internals checked at the prime powers 9..125, and eq. (4)'s sign confirmed at q = 11 by an independent subset DP
source: mac-mini-2026-10-14-wilson
depends_on: [orthogonality of the quadratic character; sum_x chi(x-u)chi(x-v) = -1 for u != v; Greene's finite-field 2F1 and 3F2 (statement 4 only)]
verification:
  - 04-computation/experiments/paley_k6p6_20261014/k6p6_values.py (+ .out): O(q^3) formula = O(q^4) brute force for q <= 43; recorded V_7, V_11, V_19, V_31 reproduced
  - 04-computation/experiments/paley_k6p6_20261014/k6p6_prime_powers.py (+ .out): statements 1-2 at 25 fields (q = 9, 25, 27, 49, 81, 121, 125, 169, 243, 289, 343, 361 and the primes 5..47), brute force at the 17 fields with q <= 49
  - 04-computation/experiments/paley_k6p6_20261014/k6p6_fft.py (+ .out): O(q log q) evaluator of M_q, equal to the O(q^3) values at all 237 primes below 1500
  - 04-computation/experiments/paley_k6p6_20261014/k6p6_proof_steps.py (+ .out): all 20 intermediate identities of the proof (and Greene's 2F1, 3F2 as exact rationals) at the 24 primes 5..101
note: 05-knowledge/results/paley_k6p6_adjoint_weight3_level32_20261014.md
related: [HYP-9276 (the modular evaluation of M_q), THM-4652, HYP-9274, HYP-9275, paley_next_cores_20261014i.md (the terminal and eq. (4)), wheel_apery_extension_20261014h.md]
---

# THM-4653: the K6 − P6 terminal is one Legendre bilinear form

## Setting

- `q` is any odd prime power, `χ` is the quadratic character of `F_q` (`χ(0) = 0`), and `ε = χ(−1)`.
- `V_q = I_(K6−P6)(q)` is the terminal of [paley_next_cores_20261014i.md](../../05-knowledge/results/paley_next_cores_20261014i.md), eq. (4):
  - take vertices `0, …, 5`;
  - delete the five path edges `01, 12, 23, 34, 45` from `K_6`;
  - orient each of the ten remaining edges from the lower label to the higher label;
  - `V_q = Σ_x Π_(edges ab, a<b) χ(x_b − x_a)`, summed over injective `x : {0..5} → F_q`.
- The **Legendre trace** is `t(a) = Σ_y χ(y)χ(y−1)χ(y−a)`.
  - For `a ∉ {0, 1}` it equals `−a_q(E_a)`, where `E_a : Y² = X(X−1)(X−a)`.
  - At the boundary, `t(0) = −ε` and `t(1) = −1`.
- The **bilinear form** is `M_q = −Σ_(a, x ∈ F_q) χ(a) χ(x) χ(x−1) χ(x−a) t(a) t(x)`.

## Statement

1. **(Affine normalisation.)** `V_q = q(q−1) W_q`, where
   `W_q = Σ χ(x₂) · χ(x₃)χ(x₃−1) · χ(x₄)χ(x₄−1)χ(x₄−x₂) · χ(x₅)χ(x₅−1)χ(x₅−x₂)χ(x₅−x₃)`.
   The sum runs over `(x₂, …, x₅) ∈ F_q⁴` subject to the four non-edge conditions `x₂ ≠ 1`, `x₃ ≠ x₂`, `x₄ ≠ x₃` and `x₅ ≠ x₄`.

2. **(Reduction.)** `W_q = −ε M_q + 7 + 13ε − 2εq`. Explicitly:
   - `W_q = M_q + 2q − 6` when `q ≡ 3 (mod 4)` (the Paley case);
   - `W_q = −M_q − 2q + 20` when `q ≡ 1 (mod 4)`.

   Moreover, `−ε M_q` is the sum in statement 1 with all four non-edge conditions dropped.

3. **(Convolution form.)** `M_q = −Σ_x χ(x(x−1)) t(x) h(x)`, where `h = (χ t) ∗ χ` is the additive convolution:
   `h(x) = Σ_a χ(a)χ(x−a) t(a) = Σ_y χ(y(1−y)) t(xy)` for `x ≠ 0`.
   Two FFTs therefore give `M_q`, and hence `V_q`, in `O(q log q)` operations when `q` is prime.

4. **(Hypergeometric form.)** In Greene's normalisation, with `φ = χ` and every function taken with quadratic parameters, the following hold for all `x ≠ 0`. (Here the sign `φ(−1)` is written out, and `ε` inside the hypergeometric symbols is Greene's trivial character.)
   - `t(x) = φ(−1) q · ₂F₁(φ, φ; ε | x)`;
   - `h(x) = q² · ₃F₂(φ, φ, φ; ε, ε | x)`.

   Hence

   `M_q = −q³ Σ_(x ∈ F_q) φ(x) φ(1−x) ₂F₁(x) ₃F₂(x)`.

   Greene's own recursion reads `₄F₃(1) = (φ(−1)/q) Σ_x φ(x)φ(1−x) ₃F₂(x)`. So `M_q` is that `₄F₃(1)` sum with one extra `₂F₁(x)` inserted. That `₄F₃(1)` is the Ahlgren–Ono value behind the wheel `W_4` in HYP-9274 (B).

5. **(Consequence for the run coefficient.)** Combine statement 2 with eq. (4) of the source note. For `q ≡ 3 (mod 4)` and `q >= 11`:

   `[t¹⁰] R_q(t) = (P₁₀(q) + (7q² − 173q + 1166) U_q + 64 q(q−1)(M_q + 2q − 6)) / (q)₁₀`.

   Every term is known in closed form except the wheel value `U_q` and the single bilinear form `M_q`. For primes, `U_q = −q(q−1)(a_q + 4)` with the Apéry form (PROVED + CITED in HYP-9274 (B)). The evaluation of `M_q` is HYP-9276.

## Proof

**1.** There are ten edges, so the summand is invariant under `x ↦ λx + c` with `λ ∈ F_q^×`: it changes by `χ(λ)¹⁰ = 1`.
- The affine group acts freely on injective tuples.
- Because `{0, 1}` is a non-edge, `x₀ ≠ x₁` holds for every injective tuple, so each orbit has exactly one point with `x₀ = 0`, `x₁ = 1`.
- After this normalisation, the edges at vertex `0` give `χ(x_j)` for `j = 2..5`, and the edges at vertex `1` give `χ(x_j − 1)` for `j = 3, 4, 5`. The remaining edges are `24`, `25` and `35`.
- Every pair joined by an edge is automatically distinct, since `χ(0) = 0`. The only distinctness conditions left are therefore the non-edges `12`, `23`, `34` and `45`.

**2.** Write
- `f₄(x; a) = χ(x)χ(x−1)χ(x−a)`, so that `F₄(a) = Σ_x f₄(x; a) = t(a)`;
- `F₅(a, b) = Σ_x χ(x)χ(x−1)χ(x−a)χ(x−b)`;
- `w(a, b) = χ(a)χ(b)χ(b−1)`.

For fixed `x₂ = a ∉ {0, 1}` and `x₃ = b ∉ {0, 1, a}`, sum over `x₅ ≠ x₄` first and then over `x₄ ≠ b`. On `x ∉ {0, 1, a}`, the product `f₄(x)·χ(x)χ(x−1)χ(x−a)χ(x−b)` equals `χ(x−b)`, which gives

`W_q = Σ_D w(a, b) [F₅(a, b)(F₄(a) − f₄(b; a)) + χ(−b) + χ(1−b) + χ(a−b)]`,

with `D = {a, b ∉ {0, 1}, a ≠ b}`. This is the identity coded in `k6p6_values.py`. Split it as `TermA + TermB + TermC` in the order written.

*Main term.* Over all `(a, b)` we have `Σ_b χ(b)χ(b−1)χ(x−b) = ε t(x)`. Therefore

`Σ_(a,b) w F₅ F₄(a) = ε Σ_(a,x) χ(a) f₄(x; a) t(a) t(x) = −ε M_q`.

This sum is the same sum as in statement 1 with the four distinctness conditions dropped.

*TermA − main.* Subtract the pairs outside `D` that `w` does not kill: `a = 1`, and `a = b`. Two evaluations are needed:
- `F₅(1, b) = −1 − χ(1−b)`, so the `a = 1` pairs contribute `−1 − ε`;
- `F₅(a, a) = −1 − χ(a)χ(a−1)` together with `Σ_(a∉{0,1}) χ(a−1) t(a) = Σ_(a∉{0,1}) χ(a) t(a) = 1 + ε`, so the `a = b` pairs contribute `−2 − 2ε`.

Hence `TermA = −ε M_q + 3 + 3ε`.

*TermC.* For `b ∉ {0, 1}`:
- `w χ(−b) = ε χ(a)χ(b−1)`;
- `w χ(1−b) = ε χ(a)χ(b)`.

Summing these and the third term over `D` with `Σ_x χ(x−u)χ(x−v) = −1` gives `(1 + ε) + ε(3 − q) + (1 + ε)`. So `TermC = 2 + 5ε − εq`.

*TermB.* Use `w(a, b) f₄(b; a) = χ(a)χ(b−a)` and exchange the sums. The inner sum over `b ∉ {0, 1, a}` equals `−ε − εχ(a)χ(x) − χ(1−a)χ(x−1)`. The three pieces then reduce, via `Σ_(a∉{0,1}) χ(a)χ(x−a) = −ε − χ(x−1)` and `Σ_a χ(a)χ(1−a)χ(x−a) = t(x)`, to

`TermB = −[(1 + ε)(−1) − ε + ε(q − 2) − (1 + ε)] = 2 + 5ε − εq`.

The last step uses `Σ_(x∉{0,1}) χ(x) t(x) = 1 + ε`.

Adding: `W_q = −εM_q + (3 + 3ε) + 2(2 + 5ε − εq) = −εM_q + 7 + 13ε − 2εq`.

**3.** Sum the definition of `M_q` over `a` first. The second form of `h` is the substitution `a = xy`, which gives `χ(xy)χ(x − xy) = χ(y)χ(1−y)` for `x ≠ 0`.

**4.** Greene defines `₂F₁(φ, φ; ε | x) = ε(x)(φ(−1)/q) Σ_y φ(y)φ(1−y)φ(1−xy)`.
- The substitution `y ↦ 1/y` in `t(x)` gives `t(x) = Σ_y φ(y)φ(1−y)φ(1−xy)`. So `t = φ(−1) q ₂F₁` for `x ≠ 0`, including `x = 1`.
- Greene's recursion gives `₃F₂(x) = (φ(−1)/q) Σ_y ₂F₁(xy) φ(y)φ(1−y) = q⁻² h(x)`. The `y = 0` term is killed by `φ(0) = 0`.
- Substituting both into statement 3 and using `φ(x−1) = φ(−1)φ(1−x)` proves statement 4.

**5.** Substitute statement 2 into eq. (4) of the source note. ∎

## Checks (FINITE-EXACT)

| Check | Scope |
|---|---|
| O(q³) formula of statement 2 against the O(q⁴) brute force of statement 1 | all primes `q <= 43`, `q ≡ 3 (mod 4)`; then 17 fields with `q <= 49`, both classes, including `q = 9, 25, 27, 49` |
| Reduction `W_q = −εM_q + 7 + 13ε − 2εq` | 25 fields up to `q = 361`: every listed prime power and the primes 5..47 |
| Recorded values `V_7 = 336`, `V_11 = −12320`, `V_19 = −76608`, `V_31 = 52080` | reproduced |
| FFT values of `M_q` against the exact O(q³) values | all 237 primes below 1500 |
| Every intermediate identity of the proof, plus statements 3–4 with Greene's functions as exact rationals | 20 identities at each prime 5..101 |
