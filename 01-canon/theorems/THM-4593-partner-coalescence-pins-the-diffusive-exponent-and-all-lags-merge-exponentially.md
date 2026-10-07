---
id: THM-4593
title: "Partner coalescence pins the diffusive exponent: for any finite set of translation partners y+1, ..., y+R (R <= 32, and S19's Mersenne lag sets D <= 7, D <= 61) the probability that y merges with none of them by Terras time T is >= c_R T^(-1/2), because on explicit residue classes the partners merge with each other first (y = 21 mod 32 for R = 2); with THM-4581 (4) the exponent is exactly 1/2 (sketch level), so extra 'nonadjacent lags' improve only the constant (sqrt(T) q_R ~ 11-14/R); with ALL translation lags at once, P(no merge with any y + r by time T) = N_T/2^T <= (T+1) 2^(-0.0346 T) (pigeonhole on (a, T^T(y) mod 3^a)); near y = -1 the relations fail together only on sets of measure 2^(-L)"
status: >
  PROVED: the witness lower bounds (exact chain verification of the witness classes for every R <= 32 and for the Mersenne lag
  sets D <= 7, D <= 61, uniform on random lifts); the all-lags identity q_inf(T) = N_T/2^T and its bound; the -1 shadow (A4).
  PROVED at sketch level: alpha_R = 1/2, via the upper bound of THM-4581 (4). NUMERICAL: sqrt(T) q_R plateaus over 2e5 paths to T = 1e6
  (bignum chain matched Python on 333 long runs). Unbounded Mersenne lag sets: OPEN, slope 0.59 +- 0.06 on [1e5, 1e6] and falling.
  Found by the session's core reader; the R = 2 witness and q_inf(T) for T <= 20 were re-checked independently here.
session: mac-mini-2026-10-07-golden
source: 05-knowledge/results/golden_collatz_resonance_20261007.md
scripts:
  - 04-computation/experiments/golden_20261007_readers/core/ (qsurv.c, qwit.c, witness_search*.py, witnesses_full.txt, allags.c, check_qsurv_chain.py, slopes.py, fit_rates.py, A/*.txt)
related:
  - THM-4581 (the pair chain; (4) the rate bounds), HYP-9217 (2) (the any-lag exponent: its 0.69 is this theorem's coalescence transient)
  - CrocSwap/integer-mult-bounds (owner-linked inspiration: nonadjacent swaps O(d) against adjacent O(d^2))
---

# THM-4593 — partner coalescence pins the diffusive exponent

## Setting

* `y` is Haar on `Z_2`. The partners `y + r` for `r ∈ P` are a finite set; `P = {1, …, R}` for translations, or S19's Mersenne lags.
* `q_P(T)` is the probability that `y` merges with none of its partners by Terras time `T`.
* Each relation `(y + r, y)` is a pair chain of THM-4581, and all chains are driven by the same parity bits `β` of `y`.

## Statements

1. **One coin per step (PROVED).** At every flip each relation moves `k ↦ k + 1 − 2β`. So all relations move in the same direction whenever they move.
2. **Coalescence lower bound (PROVED, given a witness).** Suppose on a cylinder `{y ≡ ρ mod 2^(t_0)}` all partners have merged with each other by time `t_0` while `y` has merged with none. Then `q_P(T) ≥ 2^(−t_0)·c·T^(−1/2)` for `T ≥ t_0`. From then on `y` faces a single cluster, i.e. one adjacent `±1` walk, and the SRW-skeleton bound of THM-4581 (4) applies.
   * Witnesses exist for every `R ≤ 32` and for the Mersenne lag sets `D ≤ 7` and `D ≤ 61`. They were found by search and verified exactly.
   * `R = 2`: for `y ≡ 21 (mod 32)`, at time 5 both `y + 1` and `y + 2` equal `27q + 20` while `y` is at `3q + 2`. The chain state is `(2, 2)`, so `q_2(T) ≥ 2^(−5)·2√(2/(πT)) ≈ 0.05 T^(−1/2)`. Re-checked here for 2000 lifts.
   * A witness for a lag set is also a witness for every subset, so `D ≤ 61` covers every nonempty set of odd lags `D ≤ 61`, including S19's `D ≤ 41` run (audit A2).
   * With THM-4581 (4), `α_P = 1/2` for every `P` that has a witness (sketch level). For general finite `P`, `α_P ≥ 1/2` holds at sketch level, and `= 1/2` is a CONJECTURE. The claim "nonadjacent lags give `α > 1/2`" is false for every witnessed `P`.
3. **All translation lags (PROVED).** Let `w` be the length-`T` parity word of `y`, with weight `a`, and put `c(w) = 2^T T^T(y) − 3^a y`.
   * `y` merges with some `y + r` (`r ≥ 1`) by time `T` iff some word `w′` of weight `a` has `c(w′) ≡ c(w) (mod 3^a)` and `c(w′) < c(w)`. The partner is then `r = (c(w) − c(w′))/3^a`. Equal-time merges with unequal odd counts are Haar-null.
   * Hence `q_∞(T) = N_T/2^T`, with `N_T = #{(a, c mod 3^a)} ≤ Σ_a min(C(T,a), 3^a) ≤ (T+1)·2^(H(x*)T)`.
   * Here `x* = 0.60909` solves `H(x) = x log_2 3`, so the decay rate is `1 − H(x*) = 0.0346` bits per step.
   * By symmetry the same count governs merging with some *smaller* translate `y − r`.
   * Exact values: `q_∞(T) = 0.4961, 0.3748, 0.2941, 0.2355` at `T = 8, 12, 16, 20` (re-computed here), and `0.1915` at `T = 24` (core reader).
   * The bound `(T+1)2^(−0.0346T)` is vacuous for `T ≲ 230`; it is an asymptotic statement.
4. **The −1 shadow (PROVED).** During a run of `L` odd steps of `y`, no chain is absorbed and every flip lowers `k`. Afterwards the partners' levels are exactly `O_L(r−1) − L`.
   * There is no positive-measure set of permanent joint failure, since each chain is absorbed almost surely.
   * Near −1 the relations fail together only on sets of measure `2^(−L)`. This affects the constant, not the exponent.

## Numbers (NUMERICAL, `T ∈ [10^5, 10^6]`, 2·10^5 paths)

| `R` | 1 | 2 | 3 | 4 | 8 | 16 | 32 | 64 |
|---|---|---|---|---|---|---|---|---|
| `√T q_R` | 11.20 | 6.64 | 4.60 | 3.31 | 1.47 | 0.82 | 0.28 | 0.13 |
| `R·√T q_R` | 11.2 | 13.3 | 13.8 | 13.3 | 11.7 | 13.1 | 9.1 | 8.1 |
| local slope, `[10^2, 10^3]` | 0.24 | 0.30 | 0.35 | 0.38 | 0.49 | 0.59 | 0.72 | 0.82 |
| local slope, `[10^5, 10^6]` | 0.49 | 0.50 | 0.52 | 0.51 | 0.49 | 0.43 | 0.51 | 0.48 |

* **Mersenne lags (S19).**
  * Lag 1 alone gives `√T q_1 → 16.7`.
  * `D ≤ 61`: local slopes 0.70, 0.69, 0.64, 0.59, 0.51, then 0.48 ± 0.03 on `[4·10^2, 10^6]`. So HYP-9217 (2)'s `α ≈ 0.69` is the coalescence transient.
  * `D ≤ 241`: the slope is still 0.59 ± 0.06 on `[10^5, 10^6]` (OPEN).

## Reading

* The CrocSwap analogy holds in a precise sense. A bounded number of extra partners ("nonadjacent lags") buys only a constant factor of about `R`, because the partners coalesce among themselves (Arratia-type).
* Only exponentially many lags, `r` up to `2^T`, change the decay from polynomial to exponential.
* That is the sieve regime of THM-4594.

**Audit (2026-10-07, independent audit A2).** CONFIRMED throughout:
* all 31 translation witnesses and both Mersenne witnesses;
* the all-lags identity against brute force for `T ≤ 10` and the values to `T = 24`;
* `x* = 0.60909` and the rate 0.0346;
* the −1 shadow;
* Monte Carlo for `R = 1, 2, 4, 8`.

Added: the subset property, and the precise scope of `α_P = 1/2`.
