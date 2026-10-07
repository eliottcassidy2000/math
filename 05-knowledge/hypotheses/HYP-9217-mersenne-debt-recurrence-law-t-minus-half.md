---
id: HYP-9217
title: "Debt recurrence law: in the 2-adic Haar model the lag-1 Mersenne debt merges almost surely with P(no merge by template total T) = c1 T^(-1/2)(1 + o(1)), c1 ~ 14-17, and the any-lag Mersenne non-switch probability is q(T) = T^(-alpha + o(1)) with alpha ~ 0.7 +- 0.1; hence mu_2(S) = 1 and HYP-9213 (o(A) new sigma-levels)"
status: >
  PARTLY RESOLVED 2026-10-07 (THM-4581; THM-4593: every nonempty set of odd lags D <= 61 has alpha = 1/2, lower bound PROVED via one witness covering all subsets, upper sketch level; for general finite lag sets alpha >= 1/2 (sketch) and = 1/2 is a CONJECTURE; the measured 0.69 was a coalescence transient). PROVED: almost-sure lag-1 merging, hence mu_2(S) = 1 and HYP-9213; the decay exponent
  is exactly 1/2, with the upper bound at sketch level: c T^(-1/2) <= q1(T) <= C T^(-1/2) (log T)^2 in Terras time. Template total =
  Terras merge time + v_2(merge value), with a Geom(1/2) overshoot, so the bounds transfer (audit A). If the any-lag exponent
  alpha exists then alpha >= 1/2 (sketch level). OPEN: q1(T) ~ c1 T^(-1/2) with c1 ~ 16.7 (NUMERICAL; HEURISTIC one-big-jump form
  c1 = E[J] sqrt(4/pi), THM-4581 (7)) and the value of alpha for unbounded lag sets (finite witnessed sets: alpha = 1/2, THM-4593; the earlier NUMERICAL 0.66-0.78 was a coalescence transient). Earlier evidence: two committed seeded Haar
  Monte Carlo runs (N = 3000 to total 19 924, N = 1200 to 9600), audited twice (MISTAKE-580).
source: opus-2026-10-06-S19, 05-knowledge/results/collatz_cycles_tubes_debt_walk_openai_20261006.md, section 1
related:
  - 05-knowledge/hypotheses/HYP-9213-mersenne-collatz-trajectories-coalesce-plateaus-of-odd-step-time.md
  - 05-knowledge/hypotheses/HYP-9214-reset-two-debt-resolves-with-probability-tending-to-one.md (mac-mini; the walk model)
  - 01-canon/theorems/THM-4556-the-mersenne-line-is-a-chain-of-debt-states-odd-shift-distance-2-adic-periodicity.md ((iv) 2-adic periodicity)
scripts:
  - 04-computation/experiments/mersenne_haar_debt_walk_20261006.py (+ .out default run, + _large.out for "3000 20000 61 7"; ALL CHECKS PASSED)
---

# HYP-9217 — the debt recurrence law

**Haar model.** For odd `a`, `X = 3^(a−1)` is Haar on `1 + 8Z_2`. The Mersenne switch at lag `D` is the event `U^i(2X − 1) = U^(i+D)(2X·3^(−D) − 1)` for some `i`. By THM-4556 (iv), a switch of template total `T` is a clopen event: it depends only on `a mod 2^(T−2)`. The certified share at `K` is the Haar measure of "switch with total `≤ K`".

**Conjecture.**
1. Let `q₁(T)` be the probability of no lag-1 merge with total `≤ T`. Then `q₁(T) = c₁ T^(−1/2)(1 + o(1))`.
   * Large run: `√T q₁(T) = 15.6, 16.1, 16.2, 16.3` at `T = 3200, 6400, 12 800, 19 000`.
   * Small run: `14.7, 14.5, 14.4` at `3200, 6400, 9600`, with bootstrap interval `[12.4, 16.3]` at 9600.
   * Tail exponents on `T ≥ 3200` are 0.477 (large, `[0.435, 0.520]`) and 0.521 (small, `[0.435, 0.609]`); `√T q₁` at 19 000 is 16.27 (`[14.70, 17.87]`).
2. The any-lag non-switch probability satisfies `q(T) = T^(−α + o(1))` with `α ≈ 0.7 ± 0.1`: least squares on `T ≥ 400` gives 0.679 (large run, bootstrap `[0.625, 0.745]`) and 0.782 (small run, bootstrap `[0.689, 0.900]`).

**Consequences.**
* (2) gives `μ_2(S) = 1`, the switching set having full 2-adic measure. By S18 Proposition 6 this gives HYP-9213: `o(A)` new `σ`-levels.
* The finer count `#{σ(M_a) : a ≤ A} ≈ A^(1−α)` is HEURISTIC. It needs the Haar model to govern actual exponents at `T_post(a)`, which holds only approximately (source note §1.7).

**Mechanism (PROVED pieces + HEURISTIC).**
* **Debt.** With `x = 2^L y + Δ` and `D = 3Δ + 1 − 2^L`: `D = 0` iff the next relation is the identity, and then `x = 4^k y + (4^k − 1)/3` with `k = L/2`. Here `k ≠ 0` unless the current relation is already the identity. Merges with `D ≠ 0` are value coincidences of Haar measure 0. All 834 sampled merges pass through `D = 0`.
* **Integer debts.** `D′_odd = U(D_odd) − 2^(L′ − w′)` with `w′ = v_2(3D_odd + 1)`, and `L′ = L + a − v_2(D)`.
* **The two streams.** `L_t − L_0` is the difference of the bit-consumption counts (exponent sums) of the two coupled orbits `y` and `x = 3·2^v y + 1`.
  * Each stream is marginally i.i.d. `Geom(1/2)` exactly (each orbit is Haar on a coset, and `U` preserves Haar). In generic steps the "debt exponent" `w = v_2(D)` is just `b`, the `x`-orbit's own exponent, and same-step independence of `a` and `w` is also exact.
  * The 672 443-step table (`P(w = k)`, `χ² = 7.4`, autocorrelation `−0.002`) is therefore a consistency check, not evidence.
  * The evidence on the open coupling is the block variances `Var(L_(t+K) − L_t)/K = 4.013, 4.014, 4.011, 4.016, 4.12 ± 0.14` for `K = 1, 4, 16, 64, 256` (blocks starting at `L ≥ 40`), and the lagged cross-covariances `|Cov(a_s, b_(s+k))| ≤ 0.006` for `k ≤ 64` (s.e. 0.0026).
* **Heuristic.** This is a recurrent walk. Returns to `L ∈ {±2, ±4}` bring the debt to small values, where `D = 0` has positive chance, which gives the `T^(−1/2)` law.
* **Unproved.** The long-lag joint law of the two exponent streams, i.e. recurrence of the difference of the bit-consumption counts of `y` and `3·2^v y + 1`. The marginal laws are exact.

**Evidence against actual exponents (NUMERICAL).**
* **Any partner.** For odd `a ∈ [1001, 2001]` the share with a smaller `σ`-partner is 0.970, against a Haar prediction of 0.968–0.974.
* **Level counts.** `σ`-levels for `a ≤ 100, 400, 1000, 2001` are `23, 37, 54, 69`, against predicted `22–24, 40–41, 56, 69–72`.
* **Lag 1.** The share is over-predicted by 0.02–0.04 (0.846–0.864 predicted, 0.824 observed). Actual orbits end, so their effective time is shorter.

**Not covered.** The reset-2 debt of HYP-9214 (general sources) is a natural extension; it was not tested here.

---

## Update (2026-10-07, mac-mini-2026-10-07-oaimath3; THM-4569, THM-4564)

* **The recurrence half is PROVED (THM-4569).** In the Terras clock the pair relation is `x_s = 3^(j_s) y_s + c_s`. The odd-step difference `j` is exactly a simple random walk run on the predictable disagreement clock. So almost surely the pair either merges or `j` visits every integer infinitely often. This needs no input on the long-lag joint law.
* **The exponent cannot exceed 1/2 (PROVED).** `q₁(T) ≥ (1.59 + o(1)) T^(−1/2)`.
* **Box reduction (PROVED).** Almost-sure merging is equivalent to box recurrence, an archimedean condition.
  * Computer-assisted: the lag-1 merge probability is at least 0.3853, so the lower density of odd `a` with a lag-1 switch is at least 0.3853.
* **Rate (NUMERICAL).** An exact Terras-clock sampler (20,000 paths) reproduces this hypothesis's `q₁(T)` to within about 3% and extends it: `√T q₁ = 16.77, 16.74, 16.68` at `T = 5·10⁴, 10⁵, 2·10⁵`, with tail exponent 0.476.
* **What is still needed for (1).** The `T^(−1/2)` rate: the disagreement clock must run at density 1/2 (THM-4564; HYP-9218), with uniform per-visit success.

---

## Update 2 (2026-10-07, mac-mini-2026-10-07-oaimath3): almost-sure merging and the exponent PROVED (THM-4581)

* **Almost-sure merging (THM-4581 (3)).** The Terras-clock pair chain is absorbed almost surely from every admissible start.
  * The proof combines the recurrence of `k` (THM-4569) with tightness of `|e|` at returns. Tightness comes from a Lyapunov weight `|f|^θ s^|k|`, a martingale on flips, plus fresh-coin control of runs.
  * This is the "uniform per-visit success" that the update above asked for. It holds in the averaged form `E|e_return|^θ ≤ 0.634 |e|^θ + C`.
* **Exponent (THM-4581 (4)).** `c T^(−1/2) ≤ q₁(T) ≤ C T^(−1/2) (log T)^2`.
  * The upper bound uses geometric tails for the number of excursions to absorption, plus excursion-length tails `≤ C (log t / t)^(1/2)`.
  * HYP-9218 is not needed.
* **The constant (HEURISTIC + NUMERICAL).** One big jump: `q(T) ~ E[J] √(4/(πT))`.
  * For `y` vs `y + 1`: `E[J] ≈ 9.93` (truncation-corrected) gives 11.2, against 11.06–11.18 measured (audit A).
  * For S19's lag-1 start, the post-prefix state is `(2, 1)` with `E[J] ≈ 12.8`. The constant is `(|k_0| + E[J])·√(4/π) = (2 + 12.8)·1.128 = 16.7`, which is S19's value (audit A). The earlier "`E[J] ≈ 14.8`" ignored the `|k_0|` term.
  * Proving the `(1 + o(1))` form needs a subexponential renewal argument for the Markov-modulated excursion sums.

---

## Update 3 (2026-10-07, mac-mini-2026-10-07-golden; THM-4593): every nonempty set of odd lags D ≤ 61 has exponent 1/2, and 0.69 is a coalescence transient

* **Lower bound (PROVED, THM-4593 (2)).** For the lag sets `D ≤ 7` and `D ≤ 61` there are explicit witness classes on which the partners merge with each other before the source merges with any of them. From there the source faces one cluster, i.e. one adjacent `±1` walk, so `q(T) ≥ c·T^(−1/2)`.
  * A witness covers every subset, so this holds for every nonempty set of odd lags `D ≤ 61`, including S19's `D ≤ 41` and `D ≤ 61` runs.
  * With THM-4581 (4), `α = 1/2` for these sets (upper bound at sketch level).
  * For general finite lag sets `α ≥ 1/2` (sketch level); `α = 1/2` is a CONJECTURE (scope corrected after audit A2).
* **NUMERICAL.** `D ≤ 61` has local slopes 0.70, 0.69, 0.64, 0.59, 0.51, then 0.48 ± 0.03 across `[4·10^2, 10^6]` (2·10^5 paths), with plateau `√T q ≈ 2.6`. Lag 1 alone gives `√T q_1 → 16.7`. So the measured `α ≈ 0.69` (and the earlier 0.66–0.78) was the pre-asymptotic transient while partners coalesce.
* **OPEN.** Unbounded lag sets: `D ≤ 241` still has slope 0.59 ± 0.06 on `[10^5, 10^6]`, and `√T q` (1.8 at `10^6`) is still falling. A coalescing-front heuristic predicts 1/2.
* **Consequence for HYP-9213's finer count.** The heuristic `#σ`-levels `≈ A^(1−α)` would then be `≈ A^(1/2)` asymptotically; the empirical `A^0.37` is likewise pre-asymptotic. (HYP-9213 itself, the `o(A)` statement, is PROVED by THM-4581.)
