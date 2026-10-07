---
id: THM-4565
title: "Two readers at any offset: for the coupled orbits y and x = 3*2^v*y + 1 (Haar model), at every pair of steps (s, k) the reader whose window starts later on the common 2-adic tape is Haar given the joint past, the other reads min(fresh exponent, M) with a past-measurable saturation depth M = v_2((3x_k+1) - 2^lam 3^(k+1-s)(3y_s+1)) - max(lam, 0), the windows overlap iff M >= 1, and Cov(A_s, B_k) = E[(2 - 6*2^-M) 1{M >= 1}] exactly; the first re-read bit is the fresh dictator bit XOR 1{M = 1}, the conditional information is 2(1 - 2^-M) bits, and the fresh and re-read exponents are independent iff M is Geom(1/2)"
status: "PROVED (elementary, 2-adic). FINITE-EXACT: the kernel by exact enumeration over z mod 2^22 (20 affine maps); the pointwise statements (overlap iff M >= 1, the reading identity, the XOR identity) at all 8,121,017 overlapping pairs and 12,966,272 adjacent non-overlapping pairs of 12,000 exact integer pairs (3000-bit sources, 600 steps, |k - s| <= 40). NUMERICAL: the law of M by time lag (Geom(1/2) to 4 decimals at |d + 1/2| >= 6), the covariance tables, the near-synchronous coupling. Session opus-2026-10-07-S20. Extends THM-4564 (mac-mini, concurrent, same prompt) from exact alignments (lam = 0) to all offsets."
session: opus-2026-10-07-S20 (owner prompt "figure out how the two orbits' exponents correlate over long gaps")
source: 05-knowledge/results/two_readers_exponent_correlation_20261007.md
scripts:
  - 04-computation/experiments/two_readers_correlation_20261007.py (+ .out; ALL CHECKS PASSED, about 85 s on 12 cores)
related:
  - THM-4564 (mac-mini: the lam = 0 case, i.e. exact tape alignments; lockstep min-depth rule with an equality/cancellation exception; one-sided causality; Kesten perpetuity) - statements (2)-(4) below reduce to its statements 2-3 at lam = 0
  - HYP-9217 (debt recurrence law), HYP-9218 (asymptotic Haar alignment; should be read for all overlaps, see the note)
  - 05-knowledge/results/collatz_cycles_tubes_debt_walk_openai_20261006.md (S19: Theorem D, the debt walk, the two streams)
---

# THM-4565 — two readers of one tape, at any offset

## Setting

* `e(z) = v_2(3z + 1)` and `U(z) = (3z + 1)/2^e(z)` on odd 2-adic integers.
* `y` is Haar on the odd 2-adic integers, and `v ≥ 2` with `P(v = k) = 2^-(k-1)`, independent of `y`. Set `x = 3·2^v·y + 1`.
  * Only the form `x = 2^v·u·y + c` with `u` a 2-adic unit is used.
* `y_s = U^s(y)` and `x_k = U^k(x)`, with exponents `A_s = e(y_s)` and `B_k = e(x_k)`.
* Partial sums: `S_y(s) = A_0 + … + A_(s−1)` and `S_x(k) = B_0 + … + B_(k−1)`.
* **The tape.** Both orbits are functions of the 2-adic digits of `x`. In that frame:
  * `y_s` reads the window `W_y(s) = (S_y(s) + v, S_y(s+1) + v]`;
  * `x_k` reads the window `W_x(k) = (S_x(k), S_x(k+1)]`.
* **Quantities for a pair `(s, k)`.**
  * Offset: `λ = v + S_y(s) − S_x(k)`.
  * `j = k + 1 − s`.
  * Comparison constant: `E = (3x_k + 1) − 2^λ 3^j (3y_s + 1)`, computed 2-adically (`3^j` is a unit).
  * `μ = v_2(E) ∈ Z ∪ {∞}`.
  * **Saturation depth:** `M = μ − max(λ, 0)`.
* `G_(s,k) = σ(v, A_0, …, A_(s−1), B_0, …, B_(k−1))` is the joint past of the pair.
* Define the fresh and re-read exponents:
  * `F` is the exponent of the reader whose window starts later: `A_s` if `λ ≥ 0`, `B_k` if `λ < 0`.
  * `R` is the other exponent minus `|λ|`: `R = B_k − λ` if `λ ≥ 0`, `R = A_s + λ` if `λ < 0`.

## Statements (PROVED)

1. **The later reader is fresh.**
   * Given `G_(s,k)`, `y_s` is Haar on the odd 2-adic integers if `λ ≥ 0`, and `x_k` is if `λ < 0`.
   * `E`, and hence `μ` and `M`, is `G_(s,k)`-measurable.
2. **Reading identity.** `R = min(F, M)` if `F ≠ M`, and `R ≥ M + 1` if `F = M`. In the tie case, `R − M ~ Geom(1/2)` given `G_(s,k)`, independently of `F`.
   * The windows `W_y(s)` and `W_x(k)` overlap iff `M ≥ 1`.
   * If `M ≤ 0`, the earlier reader's exponent equals `μ + |λ|·1{λ<0}`, which is `G`-measurable.
3. **Memoryless competition.** On `{M ≥ 1}`, given `G_(s,k)`:
   * (a) `Cov(F, R | G) = κ(M) = 2 − 6·2^−M`, with `κ(1) = −1`, `κ(2) = 1/2`, `κ(3) = 5/4` and `κ(∞) = 2`.
   * (b) `1{R = 1} = 1{F = 1} XOR 1{M = 1}`. The first re-read bit is the fresh dictator bit through a binary symmetric channel with crossover `P(M = 1 | ·)`.
   * (c) `I(F; R | G) = 2(1 − 2^−M)` bits. `R` determines exactly `min(F, M + 1)`.
   * (d) Given the overlap event and the offset, `(F, R)` is independent `Geom(1/2) ⊗ Geom(1/2)` **iff** `M` is `Geom(1/2)` given the overlap event and the offset.
   * Proof idea: `R` is the first position where the fresh digit string disagrees with the past-measurable string `−E/2^(λ⁺+1)`, and `F` is its first disagreement with `0`.
4. **Covariance at every pair.** For all `s, k ≥ 0`:
   `Cov(A_s, B_k) = E[κ(M_(s,k)); M_(s,k) ≥ 1] = 6·E[(1/3 − 2^−M); overlap]`.
   * Non-overlapping windows never contribute.
   * `−P(overlap) ≤ Cov(A_s, B_k) ≤ 2·P(overlap)`.
   * If `P(overlap)>0`, the covariance vanishes iff `E[2^−M | overlap] = 1/3`, the `Geom(1/2)` value. If overlap has probability zero, covariance is zero without defining that conditional law.
5. **Static form.** For `n` Haar on `Z_2` (or uniform on `[1, N]`, `N → ∞`) and fixed `h`: `Cov(v_2(n), v_2(n + h)) = 2 − 3·2^−v_2(h)`. Averaging over a Haar shift `h` gives 0. Statements 2–4 are this two-point function of `v_2`, read along two affine forms of one Haar variable whose shift `E` is written by the past.

**Remark (S19 Theorem D, both signs).** One step before the merge (the first `t` with `x_t = y_(t+1)`), the relation is a sibling relation:
* `x_t = 4^k y_(t+1) + (4^k − 1)/3` if `L_t = 2k > 0`;
* `y_(t+1) = 4^k x_t + (4^k − 1)/3` if `L_t = −2k < 0`.

So `L_t` is even and nonzero there.

## Numbers

12,000 exact pairs, 3000-bit sources, 600 steps; script `two_readers_correlation_20261007.py`.

* **FINITE-EXACT.**
  * The kernel at 20 maps, including `m = ∞` for the sibling maps `z ↦ 4^k z + (4^k − 1)/3`.
  * Statement 2's iff, the reading identity, and 3(b), at every overlapping pair (8,121,017) and its non-overlapping neighbours (12,966,272).
  * The leader's prediction of the lagger in every step.
* **NUMERICAL: the law of `M` on overlaps by time lag `d = k − s`.** Geom: `0.5, 0.25, 0.125, 0.0625`, with `E[2^−M] = 1/3`.

  | lag class | n | `P(M=1)` | `P(M=2)` | `P(M=3)` | `P(M=4)` | `E[2^−M]` | `E[κ]` |
  |---|---|---|---|---|---|---|---|
  | `6 ≤ \|d+1/2\| ≤ 12` | 1,472,757 | 0.5012 | 0.2488 | 0.1241 | 0.0633 | 0.3336 | −0.0016 |
  | `13 ≤ \|d+1/2\| ≤ 40` | 2,203,461 | 0.4993 | 0.2502 | 0.1252 | 0.0625 | 0.3331 | +0.0016 |

  * In the long-gap classes the joint law of `(F, R)` is within 0.0006 of the product. `M` is never infinite at `d ≠ −1`.
  * Near synchrony, `E[κ]` is −0.044 (`d = −1`, pre-merge), +0.022 (`d = 0`), +0.021 (`d ∈ {−2, 1}`) and −0.022 (`2 ≤ |d + 1/2| ≤ 5`).
* **NUMERICAL: offsets.**
  * Two thirds of all pre-merge overlaps have `λ ≠ 0`: shares 0.669 and 0.667. These are outside THM-4564's alignment law.
  * Both kinds carry the same conditional coupling, `E|κ(M)| = 1.00` and `E[I(F;R | past)] = 1.33` bits (`4/3` for Geom `M`), and both average to zero covariance at long gaps.
* **NUMERICAL: near-synchronous coupling.**
  * It is confined to `|L_t| ≤ 8` (max `|Cov|/s.e.` from 2.2 down to 1.2 at `L = 9..12`). It peaks at the re-reading lag `d ≈ L/2` with parity-dependent sign, e.g. `+0.181` (`L = −2, d = −1`), `+0.149` (`L = 2, d = 1`), `−0.087` (`L = 6, d = 4`).
  * Merging steps occur only at even `L ≠ 0`: 4548 at `L = −2` against 2553 at `L = +2`.

## Not claimed

* A geometric law for `M` conditional on the full joint past: `M` is measurable in that past. The original HYP-9218 formulation is refuted; a replacement must specify a coarser conditioning or an averaged offset class. See [the integration audit](../../05-knowledge/results/collatz_overlap_kernel_integration_20261007.md).
* Anything beyond the Haar model.
* That zero covariance of the raw pair `(A_s, B_k)` implies its independence. The offset class is past-measurable but correlated with both exponents; independence is asserted only for `(F, R)` within an offset class.
