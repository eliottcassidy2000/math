---
id: THM-4564
title: "Two orbits read one 2-adic tape: the exponent streams of two Collatz orbits related by x = 2^L y + Delta have an exact law at tape alignments, where x_t = 3^k y_s + kappa with kappa fixed by the past and the exponents obey an ultrametric law of depth delta = v_2(3 kappa + 1 - 3^k); averaging this law over a Haar-distributed depth gives exactly independent Geom(1/2) pairs; lockstep continuations obey E' = 3E/2^a - (3^k - 1), with the min-depth rule away from equality, unbounded cancellation possible at equality, and conditional preservation of a Haar depth input; the normalized debt has an exact forced recurrence, and its unforced comparison perpetuity has the Moran moment function of THM-4554"
status: "PROVED (elementary, 2-adic). FINITE-EXACT: every identity checked on exact integer orbits of S19's lag-1 Mersenne debt pairs: the coupling law at all 1,659,753 exact alignments, the lockstep recursion at all 554,893 continuations, past-measurability of kappa, the rho recursion; 3000 pairs, 20,000-bit sources, 5.58M steps. NUMERICAL: the law of delta at alignments, cross-covariances, block variances, the decay profile, the perpetuity law. Session mac-mini-2026-10-07-oaimath3. Independent audit: see the results note."
session: mac-mini-2026-10-07-oaimath3 (owner prompt "figure out how the two orbits' exponents correlate over long gaps")
source: 05-knowledge/results/oai3_two_orbits_twos_and_threes_20261007.md
scripts:
  - 04-computation/experiments/twoorbit_exponents_20261007.py (+ .out (400 pairs), _large.out (3000 pairs); ALL CHECKS PASSED)
  - 04-computation/experiments/twoorbit_decay_perpetuity_20261007.py (+ .out)
  - 04-computation/experiments/twoorbit_conditional_haar_20261007.py (+ .out)
related:
  - HYP-9217 (the debt recurrence law; this theorem isolates the only coupling channel); HYP-9218 (asymptotic Haar alignment)
  - 05-knowledge/results/collatz_cycles_tubes_debt_walk_openai_20261006.md (S19: Theorem D, the debt walk, the two exponent streams)
  - THM-4554 (backward sieve; Moran function rho(theta) = 3^(theta-1)/(2^theta - 1)); THM-4555, THM-4556 (Mersenne debt states)
  - H. Kesten (1973), C. M. Goldie (1991): tails of perpetuities (lattice case)
---

# THM-4564 — two orbits read one tape

**Correction, 2026-10-07:** [the integration audit](../../05-knowledge/results/collatz_overlap_kernel_integration_20261007.md) repairs the lockstep equality branch, numerical-to-exact summary, full-past conditioning, and actual-debt/perpetuity distinction. The exact local identities survive.

## Setting

* `U(x) = (3x + 1)/2^(v_2(3x+1))` on odd 2-adic integers.
* Two orbits `x_t = U^t(x)`, `y_t = U^t(y)`, with exponents `b_t = v_2(3x_t + 1)` and `a_t = v_2(3y_t + 1)`, and partial sums `B_t`, `A_t`.
* They are related by `x_t = 2^(L_t) y_t + Δ_t`, with `L_t = L_0 + A_t − B_t` and `Δ_t ∈ Z[1/2]` (S19, Theorem D). Both orbits are functions of the same 2-adic digits of `y`, the **tape**.
* The `y`-orbit has consumed `A_t` digits, and the `x`-orbit's head sits `L_t` digits behind.
* **Exact alignment** `(s, t)`, `s < t`: `B_t − B_s = L_s`. At time `t` the `x`-orbit starts reading exactly where the `y`-orbit started at time `s`. Write `k = t − s`.

## Statements (PROVED)

1. **Alignment identity.** At an exact alignment, `x_t = 3^k y_s + κ`, with `κ = (3^k Δ_s + c)/2^(L_s)` a 2-adic integer (an integer for integer orbits). Here `c` is the constant of the `x`-word on `[s, t)`.
   * `κ` is a function of the tape before position `A_s` only (**past-measurable**).
   * In the Haar model, `y_s` is Haar on the odd 2-adic integers and independent of `κ`.
2. **Ultrametric coupling law.** Put `E = 3κ + 1 − 3^k` and `δ = v_2(E)`. Since `3x_t + 1 = 3^k(3y_s + 1) + E`:
   * `b_t = a_s` if `a_s < δ`;
   * `b_t = δ` if `a_s > δ`;
   * `b_t ≥ δ + 1` if `a_s = δ`.
   Given the past, `a_s ~ Geom(1/2)` and, in the tie case, `b_t − δ ~ Geom(1/2)`. Hence `Cov(a_s, b_t | past) = 2 − 3·2^(1−δ)`, and `P(b_t = a_s | past) = 1 − 2^(1−δ)`.
3. **Haar depth ⟺ independence.** Averaging the law in 2 over `P(δ = d) = 2^(−d)` gives exactly `P(a_s = i, b_t = j) = 2^(−i−j)`. Conversely the aligned covariance is `2 − 3 E[2^(1−δ)]`, which is zero iff `E[2^(−δ)] = 1/3`.
4. **Lockstep and the 2-adic clock of 3.** If `a_s < δ`, then `(s+1, t+1)` is again aligned with the same `k`, and
   `E' = 3E/2^(a_s) − (3^k − 1)`.
   So `δ' = min(δ − a_s, ν_k)` when `δ − a_s ≠ ν_k`, and `δ' > ν_k` when they are equal, where
   `ν_k = v_2(3^k − 1) = 1` (`k` odd) and `= 2 + v_2(k)` (`k` even).
   * The bound by `ν_k` holds only when `δ - a_s != ν_k`. At equality, cancellation can give `δ' > ν_k`; there is no unconditional depth cap. The lawful states `y_s=17, x_t=49`, `k=1`, give `E=-8, δ=3, a_s=2` and `E'=-8, δ'=3>ν_1=1`. See the integration audit linked below.
   * If `E/2^(a_s)` is Haar on `2Z_2`, so is `E'`. Lockstep preserves Haar depth.
5. **One-sided causality.** If `b_t ≤ L_t` (the `x`-orbit's next read stays inside the digits the `y`-orbit has already read), then `b_t` is a function of `a_0, …, a_(t−1)`. Hence `Cov(a_(t'), b_t 1{b_t ≤ L_t}) = 0` for every `t' ≥ t`.
6. **The forced normalized-debt recurrence and its unforced comparison perpetuity.** The normalized debt `ρ_t = D_t/2^(L_t)`, with `D_t = 3Δ_t + 1 − 2^(L_t)`, satisfies exactly
   `ρ_(t+1) = (3/2^(a_t)) ρ_t − 1 + 2^(−L_(t+1))`.
   * For `M = 3/2^a` with `a ~ Geom(1/2)`, `E[M^θ] = 3^θ/(2^(θ+1) − 1)`. This is THM-4554's Moran function `ρ(θ+1)`.
   * It equals 1 at `θ = 0` and at `θ = 1`, because `3 = 2² − 1`. So `Π_(i<j) M_i = 3^j/2^(A_j)` is a mean-one martingale.
   * For `u>=1`, by Ville's inequality and optional stopping with overshoot at most `3/2`: `2/(3u) ≤ P(sup_j 3^j/2^(A_j) ≥ u) ≤ 1/u`.
   * The stationary perpetuity `Y = Σ_j 3^j 2^(−A_j)` (`ρ ≈ −Y` for `L ≫ 0`) has tail index exactly 1: `P(Y > u) ≍ 1/u` (Kesten–Goldie, lattice case).

## Numbers (S19's lag-1 Mersenne pairs; 3000 pairs, 20,000-bit sources)

* **FINITE-EXACT.** Statements 2 and 4 hold at every alignment (1,659,753) and continuation (554,893). `κ` is unchanged when the tape above the alignment position is replaced (104 cases). The `ρ` recursion is exact.
* **The depth law** (NUMERICAL).
  * Away from merges (`L_s ≥ 8`), `P(δ ≥ d)/2^(1−d) = 1.000 ± 0.001` for `d ≤ 6`, for both fresh alignments (n = 1.06M) and continuations (n = 0.53M).
  * The joint law of `(a_s, b_t)` against independent `Geom × Geom`: `χ² = 66.2` on 48 dof (n = 1.59M), and `P(b_t = a_s) = 0.33343` (independent: 1/3).
  * Successive fresh depths of one pair are jointly Haar × Haar: `χ² = 32.6` on 35 dof (n = 638,702).
* **Decay with the gap** (NUMERICAL).
  * The excess `ε(l) = E[2^(1−δ) | L_s = l] − 2/3` is `+0.043, −0.103, +0.029, −0.094` for `l = 1, 2, 3, 4`. It alternates with the parity of `l`: even `l` phase-locks (deep `δ`), matching Theorem D's rule that merges need `L` even.
  * It is at most 0.024 in size for `5 ≤ l ≤ 8`, and within about 1σ (0.006) for `9 ≤ l ≤ 40`.
* **No long-gap correlation** (NUMERICAL).
  * `|Cov(a_s, b_(s+k))| ≤ 0.0014` for all `k ≤ 160` (5.58M steps). Conditioned on `L_s` bins it is within 2σ, with no peak at the surfacing lag `k ≈ L/2`; the mean `k/L_s` at alignment is 0.511.
  * Block variances `Var(L_(t+K) − L_t)/K = 4.002, 4.002, 4.017, 4.066, 3.998` for `K = 1, 4, 16, 64, 256`.
* **Perpetuity** (NUMERICAL). At `L ≥ 24`, `−ρ ≥ 1` always. Its quantiles are within 2–14% of `Y`'s, and `u P(−ρ > u) ≈ 4–5` for `10 ≤ u ≤ 10⁴`.

## Answer to "how do the two orbits' exponents correlate over long gaps"

* The two exponent streams interact only where one orbit re-reads digits the other has already read.
  * That happens at gap `k ≈ L/2`.
  * There the interaction is the ultrametric law with depth `δ`.
  * Lockstep follows the min-depth rule away from equality; the equality branch permits cancellation beyond `v_2(3^k − 1)`.
* The reported finite tables agree with a geometric depth law to about 0.1%. An exact geometric mixing law would give independence at the corresponding aligned pair; these tables do not prove that law, all-lag independence, a variance-four limit, or recurrence.
* The recorded numerical deviations are concentrated near small `L`, with an even-gap phase pattern. No exact cutoff at `L=8` has been proved.
* The unforced comparison perpetuity has the stated tail law. Applying it to the actual debt requires control of the retained `2^(-L_next)` forcing; the finite comparison is numerical.

## Not claimed

* A geometric depth law after conditioning on the full joint past: the depth is already measurable there, so that formulation of HYP-9218 is refuted. A useful replacement requires a specified coarser conditioning; recurrence and per-visit success remain separate obligations.
* Anything about actual integers beyond the Haar model.

---

## Companion (2026-10-07): the Terras clock, THM-4569

The same pair, advanced one Terras step at a time, has relation `x_s = 3^(j_s) y_s + c_s`. The two parity streams then differ exactly at the predictable times where `c_s` is odd.
* So the odd-step difference `j` is an exact simple random walk on that clock, and its recurrence is unconditional.
* The coupling studied here (alignment depth, odd-step clock) governs the *density* of that clock: disagreement density 1/2 is equivalent to variance 4 per odd step.

---

## Update (2026-10-07, same session): scope correction from S20 (THM-4565), and almost-sure merging (THM-4581)

* **Scope.** Concurrent work by opus-2026-10-07-S20 (THM-4565) shows that the coupling sits at every **window overlap**, not only at exact alignments.
  * The later-starting reader is fresh given the joint past. The other reads `min(fresh, M)` with a past-written saturation depth `M`. The covariance is `2 − 6·2^(−M)` at every overlapping pair.
  * Statements (2)–(3) here are its `λ = 0` case. Offset overlaps are 2/3 of all coupling events.
  * Read the title's "coupled only at tape alignments" as "coupled only at window overlaps". The lockstep recursion (4) (min-depth rule off the equality branch, per the integration audit), causality (5) and the comparison perpetuity (6) are not in THM-4565.
* **The merge is almost sure** (THM-4581 (3), in the Terras clock). After it, the two exponent streams coincide exactly, shifted by the merge lag. So at long gaps the cross-covariance vanishes away from overlaps, and at the merge lag the correlation tends to 1, with a `T^(−1/2+o(1))` deficit (upper bound at sketch level).
