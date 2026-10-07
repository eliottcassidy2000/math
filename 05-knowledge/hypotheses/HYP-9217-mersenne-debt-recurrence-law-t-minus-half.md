---
id: HYP-9217
title: "Debt recurrence law: in the 2-adic Haar model the lag-1 Mersenne debt (and the reset-2 debt of HYP-9214) merges almost surely with P(no merge by template total T) ~ c1 T^(-1/2) (c1 ~ 16.3), and the any-lag Mersenne non-switch probability is q(T) = T^(-alpha + o(1)) with alpha ~ 0.68 (bootstrap [0.63, 0.75]); hence mu_2(S) = 1 and HYP-9213, with level-count exponent 1 - alpha"
status: >
  OPEN. NUMERICAL (seeded Haar Monte Carlo, N = 3000, totals to 19 924; committed rerun N = 1200) + PROVED structure
  (Theorem D of the source note: merge iff the debt D = 3 Delta + 1 - 2^L vanishes, forcing a sibling pair;
  the debt's odd part runs the Collatz map kicked by 2^(L'-w'); the walk L moves by the partner's exponent minus the debt's).
source: opus-2026-10-06-S19, 05-knowledge/results/collatz_cycles_tubes_debt_walk_openai_20261006.md, section 1
related:
  - 05-knowledge/hypotheses/HYP-9213-mersenne-collatz-trajectories-coalesce-plateaus-of-odd-step-time.md
  - 05-knowledge/hypotheses/HYP-9214 (reset-2 debt resolves almost surely)
  - 01-canon/theorems/THM-4556-the-mersenne-line-is-a-chain-of-debt-states-odd-shift-distance-2-adic-periodicity.md ((iv) 2-adic periodicity)
scripts:
  - 04-computation/experiments/mersenne_haar_debt_walk_20261006.py (+ .out, ALL CHECKS PASSED)
---

# HYP-9217 — the debt recurrence law

**Haar model.** For odd `a`, `X = 3^(a−1)` is Haar on `1 + 8Z_2`. The Mersenne switch at lag `D` is the event `U^i(2X − 1) = U^(i+D)(2X·3^(−D) − 1)` for some `i`. By THM-4556 (iv), a switch of template total `T` is a clopen event (it depends only on `a mod 2^(T−2)`), and the certified share at `K` is the Haar measure of "switch with total `≤ K`".

**Conjecture.**
1. Let `q₁(T)` be the probability of no lag-1 merge with total `≤ T`. Then `q₁(T) ~ c₁ T^(−1/2)`. Numerically `√T q₁(T) = 15.6, 16.1, 16.2, 16.3` at `T = 3200, 6400, 12 800, 19 000`. The over-`T ≥ 400` exponent 0.40 is pre-asymptotic, matching HYP-9214's `B^(−0.4)`.
2. The any-lag non-switch probability satisfies `q(T) = T^(−α + o(1))` with `α ≈ 0.68` (least squares over `[400, 19 000]`, bootstrap 95% `[0.63, 0.75]`).

**Consequences.** `μ_2(S) = 1` (the switching set has full 2-adic measure). By S18 Proposition 6 this gives HYP-9213, with `#{σ(M_a) : a ≤ A} ≈ A^(1−α + o(1))`, consistent with the observed exponent 0.36.

**Mechanism (PROVED pieces + HEURISTIC).**
* **Debt.** Write `x = 2^L y + Δ` and `D = 3Δ + 1 − 2^L`. The two orbits merge at the next step iff `D = 0`, which forces `L = 2k ≠ 0` and `x = 4^k y + (4^k − 1)/3`.
* **Generic steps.** Here `D′_odd = U(D_odd) − 2^(L′ − w′)` with `w′ = v_2(3D_odd + 1)`, and `L′ = L + a − v_2(D)`.
* **The walk.** `L` moves by the partner's exponent minus the debt's exponent. Measured mean `+0.0004`, variance `4.002 = Var(Geom − Geom)`.
* **Heuristic.** This is a recurrent walk. Returns to `L ∈ {±2, ±4}` bring the debt to small integers, where `D = 0` has positive chance. That gives the `T^(−1/2)` law.
* **Unproved.** That the debt's exponent stream (a deterministic Collatz orbit randomized only through the kicks) is asymptotically fair and independent.

**Evidence against actual exponents (NUMERICAL).**
* **Any partner.** For odd `a ∈ [1001, 2001]` the share with a smaller `σ`-partner is 0.970, against a Haar prediction of 0.968–0.974.
* **Level counts.** `σ`-levels for `a ≤ 100, 400, 1000, 2001` are `23, 37, 54, 69`, against predicted `22–24, 40–41, 56, 69–72`.
* **Lag 1.** The share is over-predicted by about 0.04 (0.863 predicted, 0.824 observed). Actual orbits end, so their effective time is shorter.
