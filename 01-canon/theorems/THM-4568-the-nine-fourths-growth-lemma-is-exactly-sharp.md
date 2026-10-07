---
id: THM-4568
title: "The growth lemma behind omega <= 9/4 (openai/math #107) is exactly sharp: the symmetric, separately concave profiles with P(1,b) = b and the shifted tripling P(a, 3h+a-1) >= 3 P(a,h) have a pointwise least element P_min(a,b) = sigma_m (2M+m-1)/2, sigma_m = Gamma(m+1/3)/(Gamma(4/3) Gamma(m)), which meets the rank bound exactly at t = 3/4; the critical exponent 4/3 = 2/(1 - s*) with s* = -1/2 the fixed point of s -> 3s+1; finite-size corollary omega <= omega_a = 3 ln(2a-1)/ln P_min(a,a) = 9/4 + 0.6843/ln a + O(1/ln^2 a)"
status: "PROVED (elementary; the session's mm94 reader; admissibility, the diagonal values and the rank bound re-checked independently here in exact rationals for a < 40; minimality checked by the reader's LP to 1e-12 for n <= 25). The finite-size bounds are PROVED modulo openai/math #107 (accepted per owner directive 2026-10-07). The Collatz links are NUMEROLOGY except the fixed-point linearization (ANALOGY with an exact identity)."
session: mac-mini-2026-10-07-oaimath3 (owner: "think of how 9/4 relates to our collatz and dynamical systems work")
source: 05-knowledge/results/oai3_two_orbits_twos_and_threes_20261007.md
scripts:
  - 04-computation/experiments/oai3_20261007_readers/mm94/ (mm94_least_profile_exact.py (ALL CHECKS PASSED), mm94_lp.py, mm94_diagonal_gadget_lp.py, mm94_quantum_functionals.py, + .out)
related:
  - openai/math #107, An Upper Bound of 9/4 for the Matrix Multiplication Exponent (2026-10-02)
  - THM-4556 (T(n) + 1 = (3/2)(n + 1): the Mersenne chain), THM-4555 (the root collision 1 + 2 = 3; -1/2 is the fixed point of 3x + 1), THM-4504 (Moran function g), HYP-9210 (tau*(2) = 4/9)
  - Ambainis-Filmus-Le Gall 2015; Alman-Vassilevska Williams 2018; Christandl-Vrana-Zuiddam 2019 (limits of methods)
---

# THM-4568 — the 9/4 growth lemma is exactly sharp

## Where 9/4 comes from (#107)

* A character `λ` has `λ(T_n) = n^(3t)`, where `t` is the mean of its three dot-product exponents.
* The symmetrized profile `P(a,b)` of the polynomial-multiplication tensors `C(a,b)` satisfies four properties:
  * separate concavity (Lemma 4.1);
  * the shifted tripling `P(a, 3h + a − 1) ≥ 3P(a,h)` (Lemma 4.2);
  * symmetry;
  * `P(1,b) = b`.
* These force `P(a,a) ≥ a^(4/3)`. The rank bound `P(a,a) ≤ (2a−1)^(1/t)` then gives `t ≤ 3/4`, hence `ω ≤ 3t ≤ 9/4`.
* So `9/4 = 3 · (3/4) = (3/2)·(3/2)`.

## Statements (PROVED)

1. **Least element.**
   * The admissible set has the pointwise least element `P_min(a,b) = σ_m (2M + m − 1)/2`, with `m = min(a,b)` and `M = max(a,b)`.
   * Here `σ_m = Π_(j<m)(1 + 1/(3j)) = Γ(m + 1/3)/(Γ(4/3)Γ(m))`, which is OEIS A004991(m−1)/9^(m−1).
   * Diagonal values: `1, 10/3, 56/9, 770/81, 3185/243, …`, with generating function `(1+x)(1−x)^(−7/3)`.
   * `P_min(a,a) ~ (3/(2Γ(4/3))) a^(4/3)`.
2. **Exact ceiling.** `P_min(a,b) ≤ (a+b−1)^(4/3)` everywhere, with equality only at `(1,1)`. So the lemma's hypotheses together with the rank bound are consistent exactly for `t ≤ 3/4`: **9/4 is the ceiling of the argument**. Every inequality in the paper's chain up to `H_a ≥ Π_(m<a)(1 + 1/(3m))` is an equality on `P_min`. The two final simplifications, `(1 + 1/(3m))³ ≥ 1 + 1/m` and `(3a − 1)/2 ≥ a`, are strict (e.g. `P_min(2,2) = 10/3 > 2^(4/3)`). They cost only the constant `3/(2Γ(4/3)) ≈ 1.68`, not the exponent. (Corrected after audit C, 2026-10-07; "everywhere" holds by Wendel's inequality `Γ(m+1/3)/Γ(m) ≤ m^(1/3)`.)
3. **Gadget barrier.**
   * Shared-leg gadgets `P(a, mh + c_a) ≥ mP(a,h)` that are valid for all characters need `c_a ≥ (m−1)(a−1)/2`, so the paper's shift `a − 1` is minimal. All such gadgets are pinned at 9/4.
   * Diagonal-only gadgets give at best `κ = 2(m−1)/(u−1)` subject to `u(u+1) ≥ 2m²`, and so can never get below about 2.058.
4. **Fixed point.**
   * Rescaled to aspect ratios, tripling is `s ↦ 3s + 1 − 1/a`; at `a = 2` it is literally `h ↦ 3h + 1`.
   * `b ↦ 3b + a − 1` is conjugate to `ũ ↦ 3ũ` through its repelling fixed point `−(a−1)/2`, and `κ* = 4/3 = 2/(1 − s*)` with `s* = −1/2`.
   * (ANALOGY) The unhalved odd step satisfies `3h + 1 + 1/2 = 3(h + 1/2)`, with fixed point `−1/2` (THM-4555). In the halved normalization this is `T(n) + 1 = (3/2)(n + 1)`, with fixed point `−1`, which drives THM-4556's Mersenne chain. (Corrected after audit C: the earlier text paired `s* = −1/2` with the identity whose fixed point is `−1`.)
   * The difference: #107's gadgets commute and share one fixed point, which is why a closed-form least element exists. Collatz words have word-dependent fixed points `c_w/(2^L − 3^w)` (the rational cycles), so no simultaneous linearization exists.
5. **Finite size (PROVED modulo #107).**
   * `ω ≤ ω_a := 3 ln(2a−1)/ln P_min(a,a)`, using only `C(a′, b′)` with `a′ ≤ a` and `b′ ≤ 4a − 1`.
   * `ω_2 = 2.7375 < log₂ 7`.
   * `ω_188 < 2.371339` (Alman–Duan–Vassilevska Williams–Xu–Xu–Zhou, SODA 2025), and `ω_190 < 2.371177` (Dupont et al. 2026, the latest pre-#107 bound that #107 cites). (Corrected after audit C: 2.371339 is not the bound immediately preceding #107.)
   * `ω_a < 2.30` from `a = 595,975` on.
   * `ω_a = 9/4 + 0.6843/ln a + O(1/ln² a)`.

## Collatz-side numerology (registered so no later session re-derives it)

`(3/2)²` is also:
* the two-odd-step multiplier;
* the tube slope `(9/4)ξ` (S19);
* a point on the upper root of THM-4504's Moran function, where `1/4 + (1/3)(9/4) = 1`.

Related numbers:
* `4/3` is also the inverse tree's mean offspring `g(0)`;
* `3/4` is also the per-odd-step drift factor;
* `4/9` is also the latest first lonely time `τ*(2)` (HYP-9210).

All agree only because `3 = 1 + 2` in both settings. NUMEROLOGY.

**Audit (2026-10-07, independent audit C).** Statements 1 and 3, the identities (OEIS A004991), the constants (`K = 0.684348`, `ω_188 = 2.371335`, `ω_190 = 2.371111`) and the numerology register are CONFIRMED. The auditor added a proof of minimality off the diagonal, by induction in `h` from the diagonal LP bound, and a tangent-line proof that minimal-shift gadgets are pinned at 9/4 when `mh + c ≥ a − 1`. Statements 2, 4 and 5 are corrected above (MISTAKE-583).
