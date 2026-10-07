---
id: THM-4563
title: "Tube theorem for real extensions of the Collatz map: Lygeros-Rozier's tube lemma for Chamberland's C(x) = x + 1/4 - (2x+1)/4 cos(pi x) (2014; re-proved here with interval arithmetic for all integers at once, a = 0.8, 0.9) and its analogue for the Dumont-Reiter 3-power extension D(x) = (3^s x + s)/2 (a = 0.6, 0.7), which proves Dumont-Reiter's Odd Critical Point Conjecture (ii) and (iii) (real sense) for every odd n and shows (i) is equivalent to n reaching 1; on the Chamberland side, rigorous c_1, c_3, c_5 -> A2 and an even-side flip lemma"
status: "PROVED (computer-assisted: mpmath interval arithmetic over the continuum h = 1/(2m+1) in [0, 1/3]). For C: KNOWN (Lygeros-Rozier 2014, Lemma 2.4, Theorem 3.3, Corollary 3.4; small odd n there checked in floating point), re-proved with full interval certification. For D: no prior proof found (Dumont-Reiter 2003 conjectured it; Lygeros-Rozier footnote 4 mention it), proved here by Lygeros-Rozier's method. Singer corollary uses Singer 1978 (CITED). INDEPENDENTLY AUDITED (2026-10-07; corrections in MISTAKE-580)."
session: opus-2026-10-06-S19
source: 05-knowledge/results/collatz_cycles_tubes_debt_walk_openai_20261006.md (section 2)
scripts:
  - 04-computation/experiments/chamberland_tubes_dumont_reiter_20261006.py (+ .out, ALL CHECKS PASSED, about 30 s)
  - 04-computation/experiments/chamberland_dumont_reiter_critical_census_20261006.py (+ .out, NUMERICAL census)
related:
  - N. Lygeros, O. Rozier, Dynamique du probleme 3x+1 sur la droite reelle, Ratio Mathematica 26 (2014) 77-94, arXiv:1402.1979 (Lemmas 2.3, 2.4, Theorem 3.3, Corollary 3.4, (5.2)-(5.3), census n <= 2000)
  - M. Chamberland, A continuous extension of the 3x+1 problem to the real line, Dynam. Contin. Discrete Impuls. Systems 2 (1996) 495-509
  - J. P. Dumont, C. A. Reiter, Real dynamics of a 3-power extension of the 3x+1 function, Dyn. Contin. Discrete Impuls. Syst. Ser. A 10 (2003) 875-893 (Conjecture 1)
  - D. Singer, Stable orbits and bifurcation of maps of the interval, SIAM J. Appl. Math. 35 (1978) 260-267
  - 05-knowledge/results/collatz_procgen_20260924_fixed_points_kawasaki_audit.md (Theorem K(b), K(c))
---

# THM-4563 — tubes for real extensions of the Collatz map

**Setting.**
* `T(n) = n/2` for even `n` and `(3n+1)/2` for odd `n`.
* The two extensions agree with `T` on `Z`:
  * `C(x) = x + 1/4 − (2x+1)/4 · cos(πx)` (Chamberland);
  * `D(x) = (3^s x + s)/2` with `s = sin²(πx/2)` (Dumont–Reiter).
* For an integer `m ≥ 1` and `a > 0`, the right tube is `τ_a(m) = [m, m + a/(2m+1)]`.
* Scaled coordinates: `x = m + ξh` with `h = 1/(2m+1)`, and image coordinate `ξ′ = (2T(m)+1)(F(x) − T(m))`.

**(i) Scaled maps (PROVED, exact identities).** `ξ′ = O(ξ, h)` for odd `m` and `ξ′ = E(ξ, h)` for even `m`, where `S = sinc(πξh/2)`, `σ = sin²(πξh/2)`, `φ(y) = (1 − e^(−y))/y` and `ψ(y) = (e^y − 1)/y`:

* `C`, odd `m`: `O = (3+h)/2 · [ξ(1 + cos(πξh)/2) − (π²ξ²/8)S²]`.
* `C`, even `m`: `E = (1+h)/2 · [ξ(1 − cos(πξh)/2) + (π²ξ²/8)S²]`.
* `D`, odd `m`: `O = (3+h)/4 · [−(3(1−h)/2) ln3 (π²ξ²/4) S² φ(σ ln3) + 3^(1−σ) ξ − h(π²ξ²/4)S²]`.
* `D`, even `m`: `E = (1+h)/4 · [((1−h)/2) ln3 (π²ξ²/4) S² ψ(σ ln3) + 3^σ ξ + h(π²ξ²/4)S²]`.

**(ii) Tube invariance (PROVED, computer-assisted).**
* For `C` with `a ∈ {0.8, 0.9}`, and `D` with `a ∈ {0.6, 0.7}`:
  * `0 < O(ξ, h) ≤ a` on `(0, a] × [0, 1/3]`;
  * `0 ≤ E(ξ, h) ≤ a` on `[0, a] × [0, 1/5]`.
* Hence `F(τ_a(m)) ⊂ τ_a(T(m))` for every integer `m ≥ 1`, and `F^k(τ_a(n)) ⊂ τ_a(T^k(n))` for all `k ≥ 0`.
* For `C` this is Lygerōs–Rozier's Lemma 2.4 in other units: their `[n, n + a_LR/(π²n)]` with `a = 2a_LR/π²` asymptotically, their admissible `27/8 < a_LR < 6` matching `6.75/π² < a < 12/π²`.

**(iii) Odd critical points (PROVED, computer-assisted).**
* For odd `n`, `F′(n) = 3/2` and `F′(n + a/(2n+1)) < 0`, and `F″ < 0` on `τ_a(n)`.
* So `F` has a unique critical point `c_n` (a local maximum) in `τ_a(n)`.
  * For `C`, cf. Lygerōs–Rozier, Lemma 2.3.
  * For `D` it is Dumont–Reiter's `c_n` (their Theorem 6, CITED: `μ_n ≤ c_n ≤ μ_(n+1)`; the tube lies in `(μ_n, μ_(n+1))`).
* Therefore `0 ≤ F^k(c_n) − T^k(n) ≤ a/(2T^k(n) + 1)` for all `k ≥ 0`.

**(iv) Dumont–Reiter's Odd Critical Point Conjecture.** Total stopping time is measured by entry into `(μ1, μ2) = (0.3158162, 1.5155526)`.
* For every odd `n ≥ 1`, (ii) holds: `c_n` and `n` have the same total stopping time, finite or infinite.
* (iii) holds in the real sense: `τ(n)` is a connected set containing `n` and `c_n` on which the total stopping time is constant. Dumont–Reiter draw (iii) in the complex plane, where we claim nothing.
* (i), "`c_n` is attracted to `(1,2)`", holds **iff the Collatz orbit of `n` reaches 1**. This uses `W_D(ξ) < ξ` on `(0, 0.6]` for the scaled return map on `τ(1)`. Otherwise `c_n` stays in tubes of integers `≥ 3`.
* So (i) for all odd `n` is equivalent to the 3x+1 conjecture.

**(v) Chamberland's map.** Known parts are Lygerōs–Rozier 2014, Theorem 3.3 and Corollary 3.4, re-proved here with interval arithmetic.
* For every odd `n ≥ 7` whose orbit reaches 1, `c_n → A1 = {1,2}`. The orbit passes through 13, 21, 40 or 64, and those tubes enter `τ(1)` at `ξ ≤ 0.0572 < 0.0705`, where `[0, 0.0705]` lies in `A1`'s basin.
* `c_1`, `c_3` and `c_5` go to `A2`: interval enclosures enter a trap `J = [0.55, 0.60]` with `|W′| ≤ 0.377`. This was numerical in Chamberland and Lygerōs–Rozier.
* `W_C` has a repelling fixed point in `(0.0705, 0.3)` (IVT; numerically `0.0710584`, Lygerōs–Rozier's `x1 = 1.023686`) and the `A2` fixed point `0.5775957` (NUMERICAL).
* Hence the 3x+1 conjecture holds iff every odd critical point `c_n`, `n ≥ 7`, of `C` is attracted to `{1,2}` (Lygerōs–Rozier, Corollary 3.4).

**(vi) Flip lemma (PROVED, computer-assisted; new).** For every even `k ≥ 2` and `ξ ∈ [−1.3, −0.45]`, `0 ≤ E(ξ, 1/(2k+1)) ≤ 0.8`. So a point at scaled position `ξ` left of `k` is mapped into `τ(T(k))`, and from then on shadows the orbit of `T(k)`.

**(vii) Singer corollary (PROVED, using Singer 1978).**
* `C` has negative Schwarzian on `[0, ∞)`: `2C′C‴ − 3C″² = π²Q(π(x+1/2))`, with `Q < 0` for `x ≥ 0`. This re-proves Chamberland's claim, used by Lygerōs–Rozier as (2.3).
* Every immediate-basin component of an attracting cycle in `(0, ∞)` is bounded. This uses Theorem K(c) of the Kawasaki audit note: intermediate-value chains through `[2j+2, 2j+3]` give divergent orbits starting beyond any bound. So Singer's compact-interval argument applies.
* Hence an attracting or neutral cycle of `C` in `(0, ∞)` other than `A1`, `A2` whose immediate basin contains an odd critical point lies in the closed tubes of a nontrivial positive integer cycle of `T`.
* Even critical points are not controlled. `c_54` leaves the near-integer regime, and `c_382`, `c_496`, `c_502` go to `A2` (NUMERICAL, as in Lygerōs–Rozier).

## Proofs

* **(i).** Use `cos(π(m+ε)) = (−1)^m cos(πε)`, `sin²(π(m+ε)/2) = cos²(πε/2)` (odd `m`) or `sin²(πε/2)` (even `m`), `2T(m) + 1 = 3m + 2` or `m + 1`, `1 − cos z = 2 sin²(z/2)`, `(3^(−σ) − 1)/h² = −ln3 (σ/h²) φ(σ ln3)`, `(3^σ − 1)/h² = ln3 (σ/h²) ψ(σ ln3)` and `σ/h² = (πξ/2)² S²`. Checked against direct evaluation in 48 cases.
* **(ii)–(iii) and (vi).**
  * Interval arithmetic on boxes covering the rectangles, with outward rounding. `sinc`, `φ` and `ψ` are enclosed by monotonicity on the ranges used.
  * The tube constants are exact decimals: boxes cover `ξ ∈ [0, A.b]`, and images are compared with `A.a`, where `A` is an interval enclosing `a`.
  * Since `1/(2m+1) ∈ [0, 1/3]`, all `m` are covered at once.
  * For `C`, `C″(m+ε) = −π sin(πε) − (2m+1+2ε)(π²/4)cos(πε) < 0` on `[0, 1/2)` for odd `m`. For `D`, `h·D″ < 0` is interval-checked for both `a`.
* **(iv).**
  * By (ii), `D^k(c_n) ∈ τ(T^k(n))`.
  * `τ(1) ⊂ (μ1, μ2)`, and `τ(m) ⊂ [2, ∞)` for `m ≥ 2`.
  * `W_D(ξ)/ξ < 1` on `(0, 0.6]` is interval-checked. With `W_D ≥ 0`, every orbit in `τ(1)` decreases to the fixed point 0, that is, to `(1,2)`.
  * Conversely, an orbit that never reaches 1 stays in tubes of integers `≥ 3`.
* **(v).**
  * The backward tree: `8 ← {5, 16}`, `5 ← {10, 3}`, `10 ← 20 ← {13, 40}`, `16 ← 32 ← {21, 64}`.
  * Odd `n ≥ 5` never meet a multiple of 3 after the start (`T(x) = 3·2^k` forces `x = 3·2^(k+1)`).
  * Interval pushes of the four tubes along the fixed tails.
  * `W_C(ξ)/ξ < 1` on `(0, 0.0705]` with `W_C ≥ 0`, and the trap `J` with `|W′| ≤ 0.377`.
* **(vii).** In Singer's argument, components of immediate basins are bounded open intervals with endpoints in `(0, ∞)`. If no critical point lies in the immediate basin, the minimum principle for the negative-Schwarzian iterate gives a contradiction. A captured odd `c_n` shadows `n` (by (iii)); then `n` reaches 1 (contradicting (v)), diverges (contradicting convergence), or enters a nontrivial cycle. ∎

## Remarks

* **What is new.**
  * For `D`, everything: Dumont–Reiter conjectured (i)–(iii) with evidence to `n < 180 000`. Lygerōs–Rozier's footnote 4 mentions the conjecture without proving it, and no other proof was found.
  * For `C`: the full interval certification for all `n` (Lygerōs–Rozier used analytic bounds plus floating-point checks at `n = 1, 3, 5, 7, 9`), the rigorous `A2` capture of `c_1, c_3, c_5`, the flip lemma, and the Singer refinement (vii).
* **Scope.** Chamberland's survey calls his two-cycle conjecture "equivalent to the 3x+1 problem". By (vii) it implies the cycle half. The converse would need control of the even critical points, and it says nothing about divergence (Lygerōs–Rozier relate divergence to wandering intervals).
* These statements transport the Collatz question into critical-orbit language exactly. They do not prove that any orbit reaches 1. Collatz is OPEN.
