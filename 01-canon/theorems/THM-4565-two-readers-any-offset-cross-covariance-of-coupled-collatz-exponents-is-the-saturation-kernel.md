---
id: THM-4565
title: "Two readers at any offset: for the coupled orbits y and x = 3*2^v*y + 1 (Haar model), at every pair of steps (s, k) the reader whose window starts later on the common 2-adic tape is Haar given the joint past; the other reader's exponent minus the offset is min(fresh exponent, M) away from ties, with a past-measurable saturation depth M = v_2((3x_k+1) - 2^lam 3^(k+1-s)(3y_s+1)) - max(lam, 0); the windows overlap iff M >= 1, M is infinite only at k - s = -1, and Cov(A_s, B_k) = E[(2 - 6*2^-M) 1{M >= 1}] exactly; the first re-read bit is the fresh dictator bit XOR 1{M = 1}, the conditional information is 2(1 - 2^-M) bits, and within an offset class the fresh and re-read exponents are independent iff M is Geom(1/2)"
status: "PROVED (elementary, 2-adic): statements 1-5 and the merge remark (almost surely). FINITE-EXACT: the kernel algebra in exact rational arithmetic (m <= 30; Geom mixture = product law for a, r <= 25); the pointwise statements (overlap iff M >= 1, the reading identity, the XOR identity, M = inf only at d = -1) at all 8,121,017 overlapping pairs and 12,966,272 adjacent non-overlapping pairs of 12,000 exact integer pairs (3000-bit sources, 600 steps, |k - s| <= 40); D = 0 at all 7444 merges. VERIFIED: the kernel by truncated enumeration over z mod 2^22 at 20 maps (error < 3e-3). NUMERICAL: the depth law by lag, the covariance tables, the near-synchronous coupling. Independently audited (session subagent; findings in MISTAKE-582); the kernel and covariance identity were also independently confirmed by codex-tiling (collatz_overlap_kernel_integration_20261007.md). Session opus-2026-10-07-S20. Extends THM-4564 (mac-mini, concurrent) from exact alignments (lam = 0, d >= 0) to all offsets and lags."
session: opus-2026-10-07-S20 (owner prompt "figure out how the two orbits' exponents correlate over long gaps")
source: 05-knowledge/results/two_readers_exponent_correlation_20261007.md
scripts:
  - 04-computation/experiments/two_readers_correlation_20261007.py (+ .out; ALL CHECKS PASSED, about 1.5 min on 12 cores)
  - 04-computation/experiments/two_readers_merge_side_20261007.py (+ .out)
related:
  - THM-4564 (mac-mini: the lam = 0, d >= 0 case, i.e. exact tape alignments; lockstep minimum rule with an equality/cancellation exception; one-sided causality; comparison perpetuity). Index translation - their y_t is our y_(t+1); their (s, t) is our (s+1, t); their k is our j.
  - THM-4581 (mac-mini: Haar coalescence with P(no merge by T) of order T^(-1/2) up to logs; with it, statement 4 gives the long-time limit Cov(A_s, B_(s+d)) -> 2 1{d = -1}, see the note's Corollary 6)
  - HYP-9217 (debt recurrence law); HYP-9218 (as first stated, conditioning on the full past, refuted by codex-tiling; a replacement must name a coarser sigma-field and cover offset overlaps)
  - 05-knowledge/results/collatz_overlap_kernel_integration_20261007.md (codex-tiling: independent derivation of statements 3-4, hostile two-depth family)
  - 05-knowledge/results/collatz_cycles_tubes_debt_walk_openai_20261006.md (S19: Theorem D, the debt walk, the two streams)
  - Terras (1976), Everett (1977), Lagarias (1985), Bernstein-Lagarias (1996), Tao (2022): the one-orbit parity-vector bijection that statement 1 extends to two orbits
---

# THM-4565 — two readers of one tape, at any offset

## Setting

* `e(z) = v_2(3z + 1)` and `U(z) = (3z + 1)/2^e(z)` on odd 2-adic integers.
* `y` is Haar on the odd 2-adic integers, and `v ≥ 2` with `P(v = k) = 2^-(k-1)`, independent of `y`. Set `x = 3·2^v·y + 1`.
  * Then `x` is Haar on `1 + 4Z_2`.
  * For general `x = 2^v·u·y + c` (`u` a 2-adic unit), replace `3^j` below by `u·3^(k−s)`.
* `y_s = U^s(y)` and `x_k = U^k(x)`, with exponents `A_s = e(y_s)` and `B_k = e(x_k)`.
* Partial sums: `S_y(s) = A_0 + … + A_(s−1)` and `S_x(k) = B_0 + … + B_(k−1)`.
* **The tape.** In the digit frame of `x`:
  * `y_s` reads the window `(S_y(s) + v, S_y(s+1) + v]`;
  * `x_k` reads the window `(S_x(k), S_x(k+1)]`.
* **Quantities for a pair `(s, k)`.**
  * Offset `λ = v + S_y(s) − S_x(k)`, `j = k + 1 − s`, and lag `d = k − s`.
  * Comparison constant `E = (3x_k+1) − 2^λ 3^j (3y_s+1)`, computed 2-adically.
  * `μ = v_2(E)`, and saturation depth `M = μ − max(λ, 0)`.
* `G = σ(v, A_(<s), B_(<k))` is the joint past.
* **Fresh and re-read exponents.**
  * If `λ ≥ 0`: `F = A_s`, `R = B_k − λ`, `w = 3^j(3y_s+1)/2`, `c = −E/2^(λ+1)`.
  * If `λ < 0`: `F = B_k`, `R = A_s + λ`, `w = (3x_k+1)/2`, `c = E/2`.
  * In both cases `F = 1 + v_2(w)`, `R = 1 + v_2(w − c)` and `v_2(c) = M − 1`.

## Statements (PROVED)

1. **The later reader is fresh.** Given `G`, the later-starting reader (`y_s` if `λ ≥ 0`, `x_k` if `λ < 0`) is Haar on the odd 2-adic integers. Also `w` is Haar on `Z_2`, and `E` is `G`-measurable:
   `E = (3^(k+1) + c'_(k+1) − 2^v 3^j c_(s+1))/2^(S_x(k))`.
2. **Reading identity.**
   * `R = min(F, M)` if `F ≠ M`; if `F = M`, then `R − M ~ Geom(1/2)` given `G`.
   * The windows overlap iff `M ≥ 1`. If `M ≤ 0`, the earlier reader's exponent is the `G`-measurable number `μ + |λ|·1{λ<0}`.
   * `M = ∞` (that is, `E = 0`) only at `j = 0`, i.e. `d = −1`. Proof: mod 3, `c'_(k+1) ≡ 2^(S_x(k))` and `c_(s+1) ≡ 2^(S_y(s))`, and both are nonzero.
3. **Memoryless competition.** On `{M ≥ 1}`, given `G`:
   * (a) `Cov(F, R | G) = κ(M) = 2 − 6·2^−M`, with `κ(1) = −1`, `κ(∞) = 2`.
   * (b) `1{R = 1} = 1{F = 1} XOR 1{M = 1}`.
   * (c) `I(F; R | G) = 2(1 − 2^−M)` bits and `P(F = R | G) = 1 − 2^(1−M)`. The three summaries (a) and (c) are affine in `E[2^−M]`.
   * (d) Given the overlap event and `λ`, `(F, R)` is `Geom(1/2) ⊗ Geom(1/2)` iff `M` is `Geom(1/2)` given the overlap event and `λ`.
4. **Covariance at every pair.** For all `s, k ≥ 0`:
   `Cov(A_s, B_k) = E[κ(M_(s,k)); M_(s,k) ≥ 1] = 6·E[(1/3 − 2^−M); overlap]`, and `−P(overlap) ≤ Cov ≤ 2·P(overlap)`.
   * If `P(overlap) > 0`, the covariance vanishes iff `E[2^−M | overlap] = 1/3`. That condition is necessary but not sufficient for independence; see 3(d).
   * Since `Var(A_s) = Var(B_k) = 2` and both streams are independent sequences, `Var(L_(t+K) − L_t) = 4K − 2 Σ Cov(A_s, B_k)` exactly.
5. **Static form (elementary, likely folklore).** For `n` Haar on `Z_2` and fixed `h`: `Cov(v_2(n), v_2(n+h)) = 2 − 3·2^−v_2(h)`.

**Remark (S19 Theorem D, both signs; almost surely).**
* In the Haar model, almost surely `D_(τ−1) = 0` at the merge time `τ` (the first `t` with `x_t = y_(t+1)`). One step before the merge the relation is then a sibling relation:
  * `x_t = 4^k y_(t+1) + (4^k − 1)/3` (`L_t = 2k > 0`), or
  * `y_(t+1) = 4^|k| x_t + (4^|k| − 1)/3` (`L_t = −2|k| < 0`).
* Integer value coincidences without `D = 0` exist. Example: `y = 1, v = 2`.

## Numbers (NUMERICAL unless stated)

12,000 exact pairs; script `two_readers_correlation_20261007.py`.

* **Depth law on overlaps, single lags.**
  * Pearson `χ² < 28` (6 dof) at every lag with `d ≥ 6` or `d ≤ −9`, out to `|d| = 40`.
  * Pooled over `13.5 ≤ |d+½| ≤ 40.5`: `P(M=1..4) = .4993/.2502/.1252/.0625` (2.03M overlaps, `χ² = 8.2`).
  * The near-synchronous lags `d ∈ [−8, 5]` deviate, with effect `√(χ²/n)` about halving per lag (`.17, .085, .043, .030, .020` at `d = −4…−8`).
  * This is an ensemble average at `s ≥ 10`. On explicit cylinders (`v = 4k`) codex-tiling exhibits non-Geom depths at arbitrarily large gaps.
* **Offset overlaps** are two thirds of pre-merge overlaps. They carry the same conditional coupling (`E|κ| ≈ 1`), and at long lags the same depth law as exact alignments.
* **Covariances.** Below detection at bit-gap `|L| ≥ 9`. Near synchrony they peak at the re-reading lag with parity-dependent sign, up to `+0.18`.
* **Merges.** All 7444 have `D = 0` (FINITE-EXACT), at even nonzero `L`. The `−2` side dominates early (1128 : 24 for `τ ≤ 5`); later ratios are 1.14–1.27.

## Not claimed

* A geometric law for `M` conditional on the full joint past. `M` is measurable in that past, and HYP-9218 as first stated is refuted. A replacement must specify a coarser conditioning or an averaged offset class.
* Independence of the raw pair `(A_s, B_k)` from zero covariance; independence is asserted only for `(F, R)` within an offset class. The whole processes are totally dependent: `B` is a function of `(v, A)`.
* Anything beyond the Haar model.
* The long-time limit `Cov(A_s, B_(s+d)) → 2·1{d=−1}` holds only conditionally on THM-4581 (note, Corollary 6).
