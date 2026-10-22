---
id: HYP-9366
title: "The cosine twin and the alternating end of the forest series. PROVED parts: (1) the forward-step polynomial in forest form, F(T,x) = sum over even L of N_L(S) (n-L)! ((x-1)/2)^L ((x+1)/2)^(n-1-L); odd N_L vanish (path reversal), so F is palindromic and F(T,-1) = N_(n-1) (odd n); (2) the alternating count as a tangent-number OCF, F(T,-1) = (-1)^((n-1)/2) sum_S prod_(c in S) [(-1)^((|c|-1)/2) 2^|c|] E_(n-phi(S)) over collections S of disjoint odd directed cycles; (3) for odd n the rotational R_n and the transitive tournament are cospectral twins, both with skew spectrum {i cot(pi(2m+1)/(2n))}, so they have the same blow-up gain rho, and their H differ only through the tree level. CONJECTURE: |F(R_n,-1)|/n! ~ C (2/pi)^n / n, with 2/pi the reciprocal of the first zero of cos (the R_n limit spectrum is the set of reciprocal zeros of cos), and the n mod 4 wobble of H(R_n) is proportional to it."
status: >
  (1)-(3) PROVED (elementary; (2) via the Redei-Berge power-sum expansion, the OCF form used by the repo since THM-062/063
  and Grinberg-Stanley), verified on random tournaments n <= 7 (fneg1_ocf.out). The rate statements are FINITE-EXACT fits
  plus conjecture (mac-mini-2026-10-22-sierp, pre-audit). The cospectrality gives a clean example that the one-loop gain
  cannot see regularity: rho(R_n) = rho(transitive), while kappa = 2.86 against 2^(n-1)/n!. A first "resurgence" guess
  (non-perturbative rate = reciprocal of the first zero of prod(1 - z^2 x_k^2) for every regular family) FAILS for Paley.
  PRIOR ART: (1) and (2) are the repo's W-polynomial theory in other coordinates. THM-059 (kind-pasteur, 2026-03) defines
  W(r) = sum_P prod (r + s_i), s_i = A - 1/2, so [r^(n-1-L)] W(r) = 2^(-L) (n-L)! N_L: the forest numbers ARE the
  W-coefficients. THM-059 also proves W(r) = sum over independent S of 2^|S| F_(n-1-phi(S))(r) (a deformed independence
  polynomial); THM-061 gives W(-1/2) = (-1)^(n-1) H; THM-062 gives the deformed-Eulerian closed form. (2) is THM-059/062
  at the alternating point; the tangent-number reading, the cospectral twins (3) and the rate statements are new here.
session: mac-mini-2026-10-22-sierp
source: 01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md
scripts:
  - 04-computation/experiments/sierpinski_ocf_20261022/fneg1_ocf.py (+ .out): (1), (2)
  - 04-computation/experiments/sierpinski_ocf_20261022/fneg1.c (+ fneg1.out, fneg1_rate.out), fneg1_gen.c (+ fneg1_paley.out)
  - 04-computation/experiments/sierpinski_ocf_20261022/forest_top.c (+ .out): top of the forest series for R_n
related:
  - 05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md
  - 01-canon/theorems/THM-4756-forward-step-polynomial-mod-3-is-h-times-eulerian.md
  - 01-canon/theorems/THM-059-universal-master-polynomial.md
  - 01-canon/theorems/THM-061-anti-evaluation.md
  - 01-canon/theorems/THM-062-forward-edge-distribution.md
---

# HYP-9366 — the cosine twin and the alternating end of the forest series

## 1. The forward-step polynomial in forest form (PROVED)

For a permutation `w` with step signs `s_i = S_(w_i w_(i+1))`, a step is forward exactly when `s_i = +1`. Hence `x^[forward] = (x+1)/2 + s_i (x−1)/2`. Multiply over the steps, expand, and use THM-4757 (1):

    F(T, x) = Σ_L N_L(S) (n − L)! ((x − 1)/2)^L ((x + 1)/2)^(n−1−L).

Reversing every path of a linear forest negates each arc, so `N_L = (−1)^L N_L`, and `N_L = 0` for odd `L`. Consequences:
* `F(T, x) = x^(n−1) F(T, 1/x)` (palindromic), and `F(T) = F(T^op)`;
* `F(T, −1) = N_(n−1)`, the signed count of Hamiltonian paths of `K_n`, for odd `n`, and 0 for even `n`;
* the Taylor coefficients of THM-4756 are `c_m = 2^(m−n) Σ_L N_L (n−L)! C(n−1−L, m−1)`.

## 2. The alternating count as a tangent-number OCF (PROVED)

The Rédei–Berge function has the power-sum form `U_T = Σ_S 2^|S| p_type(S)`, over collections `S` of disjoint odd directed cycles (the OCF bridge). Combine it with the specialisation `Σ_w t^|D(w)| = (1−t)^(n+1) Σ_(m≥1) U(1^m) t^(m−1)`. This gives

    F(T, t) = Σ_S 2^|S| (1 − t)^φ(S) Ã_(n−φ(S))(t),     φ(S) = Σ_(c∈S) (|c| − 1),

with `Ã_k(t) = Σ_(w∈S_k) t^des(w)`. At `t = 0` this is the OCF `H = I(Ω, 2)`. At `t = −1`, using `Ã_k(−1) = (−1)^((k−1)/2) E_k` for odd `k` (tangent numbers 1, 2, 16, 272, 7936, …):

    F(T, −1) = (−1)^((n−1)/2) Σ_S Π_(c∈S) [(−1)^((|c|−1)/2) 2^|c|] E_(n−φ(S)).

So the top term of the forest series is an OCF with tangent-number weights. This is checked exactly on random tournaments with `n = 3, 5, 7` (`fneg1_ocf.out`).

## 3. Cospectral twins (PROVED)

For odd `n`, let `T_n` be the transitive tournament, with skew matrix `sign(j − i)`. This matrix is a skew-circulant, with first row `(0, 1, …, 1)` and `c_(−k) = −c_(n−k)`. So it is diagonalised by `(ζ^j)` with `ζ^n = −1`, and its eigenvalues are

    Σ_(k=1)^(n−1) ζ^k = (1 + ζ)/(1 − ζ) = i cot(π(2m+1)/(2n)),   m = 0..n−1.

The rotational circulant `R_n` has eigenvalues `i cot(πk/2n)` for odd `k` and `−i tan(πk/2n) = i cot(π(n+k)/(2n))` for even `k`: the same set. Both are normal, so the two skew matrices are orthogonally similar (checked numerically to `n = 21`). Hence:
* `ρ(R_n) = ρ(T_n)`. Wormald's gain cannot tell the best-known Hamiltonian tournament from the worst one.
* `κ(R_n) → 2.856` while `κ(T_n) = 2^(n−1)/n!`. The difference is entirely the tree level, `S·1 ≠ 0` for `T_n`, which the regular loop expansion of THM-4757 excludes.
* In the limit the normalised spectrum is `{±2/(πj) : j odd}`, the reciprocal zeros of `cos`. Hence `Π(1 − z²x²) = cos² z` and `Π(1 + z²x²) = cosh² z` over the spectrum, so `ρ(R; z) = cosh z/cos z`, `tr G = z tan z` and `tr H = −z tanh z`. These are the constants of HYP-9365 (b).

## 4. The alternating end for `R_n` (FINITE-EXACT and conjecture)

Exact `F(R_n, −1) = 6, 0, 112, −2304, 84480, −4499456, 330332160, −31989497856, 3951688548352, −606498931605504` (`n = 3..21`, `fneg1.out`). Neither this sequence nor `F/n` nor `F/(n 2^((n−1)/2))` is in OEIS (2026-10-22).
* **Rate.** With the rate fixed at `2/π`, `log|F/n!| = n log(2/π) + 1.29 − 1.03 log n` fits `n = 11..21` with residuals ≤ `6·10^(−4)` (≤ `10^(−4)` from `n = 13`). A free-rate fit gives 0.6359, against `2/π = 0.6366`.
* **Comparison.** The transitive twin has `|F(T_n, −1)|/n! = E_n/n! ~ (4/π)(2/π)^n`: the same rate, without the `1/n`.
* **Conjecture (a).** `|F(R_n, −1)|/n! ~ C (2/π)^n / n`, with `C ≈ 3.6` and sign `(−1)^((n+1)/2)`.
* **Conjecture (b).** The `n mod 4` wobble of `κ(R_n)` (HYP-9365 (e)) equals `c·F(R_n, −1)/n!·(1 + o(1))`, with measured `c ≈ 1.6–1.7` (`n = 13–21`).
* The top terms `N_(n−1−j)/(n)_(n−1−j)` (`forest_top.out`) show that the alternating part is spread over the top of the series. The `j = 2` term also alternates (about −20% of `j = 0`), and the `j ≥ 4` terms are positive with an alternating modulation. So (b) is a proportionality, not an identity term by term.

## 5. What fails

The natural "resurgence" generalisation does not hold. The guess is that, for every regular family, the alternating end decays at `1/z*`, where `z*` is the first zero of `Π(1 − z² x_k²)`.
* For Paley `z* = √q`, which would make `F(P_q,−1)/q!` about `q^(−q/2)`. The exact values are 0.2, 1.2e-2, 9.0e-6 at `q = 7, 11, 19`, against `q^(−q/2) = 1.1e-3, 1.9e-6, 7.1e-13` (`fneg1_paley.out`).
* So the `2/π` for `R_n` is a property of the cosine twin pair, not a general law.
* The doubly regular `D_19` has `F = 1088644055040` against `F(P_19, −1) = 1094419611648`, a difference of 0.53%. The alternating end is far less universal than `H` (3.2e-5), as expected for the far end of the series.
