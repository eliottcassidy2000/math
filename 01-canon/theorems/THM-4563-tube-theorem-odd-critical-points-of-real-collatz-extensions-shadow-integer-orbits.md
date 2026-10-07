---
id: THM-4563
title: "Tube theorem: for Chamberland's C(x) = x + 1/4 - (2x+1)/4 cos(pi x) and the Dumont-Reiter 3-power extension D(x) = (3^s x + s)/2 (s = sin^2(pi x/2)), every right tube tau(m) = [m, m + a/(2m+1)] (m >= 1; a = 0.8 for C, 0.6 for D) is mapped into tau(T(m)), and the odd critical point c_n lies in tau(n); so odd critical orbits shadow Collatz orbits forever. Dumont-Reiter's Odd Critical Point Conjecture (2003): parts (ii)-(iii) PROVED for every odd n, part (i) equivalent to n reaching 1; for C, every odd n >= 7 that reaches 1 sends c_n to {1,2} (c_1, c_3, c_5 go to A2, the satellite of {1,2} in its tubes), so the 3x+1 conjecture holds iff every odd critical point c_n (n >= 7) of C is attracted to {1,2}"
status: "PROVED (computer-assisted: mpmath interval arithmetic over the continuum h = 1/(2m+1) in [0, 1/3], covering all integers at once); Singer corollary PROVED using Singer 1978 (CITED); even critical points NOT covered (NUMERICAL census only). Independent audit: see section 7 of the results note."
session: opus-2026-10-06-S19
source: 05-knowledge/results/collatz_cycles_tubes_debt_walk_openai_20261006.md (section 2)
scripts:
  - 04-computation/experiments/chamberland_tubes_dumont_reiter_20261006.py (+ .out, ALL CHECKS PASSED, about 25 s)
  - 04-computation/experiments/chamberland_dumont_reiter_critical_census_20261006.py (+ .out, NUMERICAL census)
related:
  - M. Chamberland, A continuous extension of the 3x+1 problem to the real line, Dynam. Contin. Discrete Impuls. Systems 2 (1996) 495-509 (via Lagarias's bibliography, entry 35, and Chamberland's 2003 survey)
  - J. P. Dumont, C. A. Reiter, Real dynamics of a 3-power extension of the 3x+1 function, Dyn. Contin. Discrete Impuls. Syst. Ser. A 10 (2003) 875-893 (Conjecture 1, Odd Critical Point Conjecture)
  - D. Singer, Stable orbits and bifurcation of maps of the interval, SIAM J. Appl. Math. 35 (1978) 260-267
  - 05-knowledge/results/collatz_procgen_20260924_fixed_points_kawasaki_audit.md (Theorem K(b): integer cycles of C are attracting iff positive)
---

# THM-4563 — the tube theorem for real extensions of the Collatz map

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

**(iii) Odd critical points (PROVED, computer-assisted).**
* For odd `n`, `F′(n) = 3/2` and `F′(n + a/(2n+1)) < 0`, and `F″ < 0` on `τ_a(n)`.
* So `F` has a unique critical point `c_n` (a local maximum) in `τ_a(n)`. For `D` this is Dumont–Reiter's `c_n`.
* Therefore `0 ≤ F^k(c_n) − T^k(n) ≤ a/(2T^k(n) + 1)` for all `k ≥ 0`: **the odd critical orbit shadows the Collatz orbit of `n` forever.**

**(iv) Dumont–Reiter's Odd Critical Point Conjecture.** Total stopping time is measured by entry into `(μ1, μ2) = (0.3158162, 1.5155526)`.
* For every odd `n ≥ 1`, (ii) holds: `c_n` and `n` have the same total stopping time, finite or infinite.
* (iii) holds in the real sense: `τ(n)` is a connected set containing `n` and `c_n` on which the total stopping time is constant.
* (i), "`c_n` is attracted to `(1,2)`", holds **iff the Collatz orbit of `n` reaches 1**. This uses `W_D(ξ) < ξ` on `(0, 0.6]` for the scaled return map on `τ(1)`.
* So (i) for all odd `n` is equivalent to the 3x+1 conjecture.

**(v) Chamberland's map.**
* On `τ(1)` the scaled return map `W_C = E(·, 1/5) ∘ O(·, 1/3)` has fixed points:
  * `0`, which is `A1 = {1,2}`;
  * `0.0710584`, repelling;
  * `0.5775957`, which is `A2 = {1.1925319, 2.1386563}`.
* `A2` is the satellite of `A1` in `τ(1) ∪ τ(2)`.
* For every odd `n ≥ 7` whose orbit reaches 1, `c_n → A1`. Such orbits pass through 13, 21, 40 or 64, and those tubes enter `τ(1)` at `ξ ≤ 0.0572 < 0.0705`, inside `A1`'s basin.
* `c_1`, `c_3` and `c_5` go to `A2`.
* Hence **the 3x+1 conjecture holds iff every odd critical point `c_n`, `n ≥ 7`, of `C` is attracted to `{1,2}`.**

**(vi) Singer corollary (PROVED, using Singer 1978).**
* `C` has negative Schwarzian on `[0, ∞)`: `2C′C‴ − 3C″² = π²Q(π(x+1/2))`, with `Q < 0` for `x ≥ 0`.
* An attracting or neutral cycle of `C` in `(0, ∞)` other than `A1`, `A2` whose immediate basin contains an odd critical point lies in the closed tubes of a nontrivial positive integer cycle of `T`.
* Even critical points are not controlled. `c_54` leaves the near-integer regime (as in Dumont–Reiter's Table 6 for `D`). All even ones up to 300 end at `A1` (NUMERICAL).

## Proofs

* **(i).** Use `cos(π(m+ε)) = (−1)^m cos(πε)`, `sin²(π(m+ε)/2) = cos²(πε/2)` (odd `m`) or `sin²(πε/2)` (even `m`), `2T(m) + 1 = 3m + 2` or `m + 1`, and `1 − cos z = 2 sin²(z/2)`. Then `(3^(−σ) − 1)/h² = −ln3 (σ/h²) φ(σ ln3)` and `(3^σ − 1)/h² = ln3 (σ/h²) ψ(σ ln3)`, with `σ/h² = (πξ/2)² S²`. The script checks the identities against direct evaluation in 48 cases.
* **(ii)–(iii).** Interval arithmetic on boxes covering the rectangles, with outward rounding. `sinc`, `φ` and `ψ` are enclosed by monotonicity. Since `1/(2m+1) ∈ [0, 1/3]`, all `m` are covered at once. For `C`, `C″(m+ε) = −π sin(πε) − (2m+1+2ε)(π²/4)cos(πε) < 0` on `[0, 1/2)` for odd `m`. For `D`, `h·D″ < 0` is interval-checked.
* **(iv).** By (ii), `D^k(c_n) ∈ τ(T^k(n))`. `τ(1) ⊂ (μ1, μ2)`, and `τ(m) ⊂ [2, ∞)` for `m ≥ 2`. `W_D(ξ)/ξ < 1` on `(0, 0.6]` is interval-checked, and an orbit that never reaches 1 stays in `[2, ∞)`.
* **(v).**
  * The backward tree: `8 ← {5, 16}`, `5 ← {10, 3}`, `10 ← 20 ← {13, 40}`, `16 ← 32 ← {21, 64}`.
  * An odd `n ≥ 5` never meets a multiple of 3 after its start: `T(x) = 3·2^k` forces `x = 3·2^(k+1)`.
  * The whole tubes at the four nodes are pushed by interval arithmetic along the fixed tails.
  * `W_C(ξ) < ξ` on `(0, 0.0705]` with `W_C` increasing there. `J = [0.55, 0.60]` is a trap with `|W′| ≤ 0.377`. The enclosures of `c_1`, `c_3`, `c_5` enter `J`.
* **(vi).** Singer's theorem on `[0, ∞)`:
  * `C` is `C³`, `SC < 0`, and `C([0, ∞)) ⊂ [0, ∞)`.
  * Immediate basins are bounded, because monotone divergent orbits start at arbitrarily large points (Theorem K(c) of the Kawasaki audit note).
  * So the immediate basin contains a critical point.
  * If it is an odd `c_n`, then `n` reaches 1 (contradicting (v)), diverges (contradicting convergence), or enters a nontrivial cycle. ∎

## Remarks

* **Prior work.** Dumont and Reiter (2003) conjectured (i)–(iii) for `D` with evidence to `n < 180 000`. We found no proof in the sources checked (the papers, Lagarias's bibliographies, Chamberland's survey, web searches on 2026-10-06).
  * Chamberland's survey calls his two-cycle conjecture "equivalent to the 3x+1 problem". By (vi) it implies the cycle half; the converse would need control of the even critical points.
* **Scope.** These are statements about real extensions. They transport the Collatz question into critical-orbit language exactly. They do not prove that any orbit reaches 1. Collatz is OPEN.
