---
id: THM-4581
title: "Haar coalescence: two Collatz (Terras) orbits related by u = 3^k v + e with e in Z[1/3] merge almost surely, at equal time and with equal odd-step counts; P(no merge by T) lies between c T^(-1/2) and C T^(-1/2) (log T)^2; hence y and y+1 merge for almost every 2-adic y (HYP-9220), the Collatz grand-orbit relation on Z_2 equals the orbit relation of Z[1/6] x| <2,3> up to null sets, the Mersenne switching set has full measure, HYP-9213 and HYP-9214 hold, and HYP-9217's decay exponent is 1/2"
status: >
  PROVED (elementary probability: a bounded-increment martingale, an explicit Lyapunov weight s^|k|, optional stopping,
  Levy's 0-1 law): statements 1-3, 5, 6. PROVED at sketch level (standard Foster-Lyapunov and first-passage estimates):
  the rate (4). Local lemmas checked exhaustively on all 336,040 (state, bit) pairs with |k| <= 10, |3^max(0,-k) e| <= 4000.
  The chain was matched against direct 2-adic orbits on 450,000 steps with 0 mismatches.
  FINITE-EXACT (6c): for every residue n mod 2^K and every even K in [2, 20], chain absorption by step K equals
  an actual integer merge of n and n+1 (n = 2^(K+40) + r): 0 disagreements; unmerged fraction 626933/2^20 = 0.598 at K = 20. NUMERICAL: the return drift is 0.59-0.61
  against the proved 0.634. Direct big-integer orbits for six relations merge within 8000 steps in 80-89% of 120 trials each.
  HEURISTIC + NUMERICAL: the one-big-jump constant (7). Found by the session's pick reader (pair-chain lane); the proof was
  re-derived and re-checked here, and the rate (4) was added here. Independent audit pending at checkpoint.
session: mac-mini-2026-10-07-oaimath3
source: 05-knowledge/results/oai3_two_orbits_twos_and_threes_20261007.md
scripts:
  - 04-computation/experiments/oai3_20261007_coalescence/ (coal_check.py: chain vs direct orbits, run lengths, return drift, survival; lemmas_exhaustive.py; direct_merge.py: chain-free merges; excursions_J.py: E[J] and the one-big-jump constant; residue_exact.py: exact residue census for (6c); + .out)
  - 04-computation/experiments/oai3_20261007_readers/pick/ (the lane's pair_chain_*.py and outputs: integer check, chain vs direct merges on 3000 random 404-bit p, lemma check, exact dyadic values, Monte Carlo, flip autocorrelation, drift test)
related:
  - THM-4569 (the same chain in the notation (j, c); this theorem closes its box recurrence (5) and its index question (7))
  - THM-4564, THM-4565 (S20) (the same pair in the odd-step clock: coupling only at window overlaps, Haar depth = independence)
  - HYP-9220 (now PROVED), HYP-9213 (PROVED), HYP-9214 (PROVED), HYP-9217 (almost-sure part and exponent PROVED; constant open), HYP-9218 (no longer needed for any of these)
  - THM-4556 (ii), (iv) (pairing; 2-adic periodicity of templates); S18 Proposition 6 (mersenne_switch_parity_f21_compression_20261006.md)
  - Terras (1976), Lagarias (1985) (the parity-vector map is a measure-preserving bijection of Z_2); Durrett, Probability, Thm 4.3.1 (bounded-increment martingales)
---

# THM-4581 — Haar coalescence of affinely related Collatz orbits

## Setting

* `T(x) = x/2` for even `x` and `(3x+1)/2` for odd `x`, on `Z_2`. The variable `y` is Haar. Write `v_n = T^n(y)`, `β_n = v_n mod 2` and `F_n = σ(y mod 2^n)`.
  * By Terras, `β_n` is a fair coin independent of `F_n`, and `v_n` is Haar given `F_n`.
* Let `u_0 = 3^(k_0) y + e_0` with `3^max(0,−k_0) e_0 ∈ Z`, `(k_0, e_0) ≠ (0, 0)`, and `u_n = T^n(u_0)`.
* **The pair chain** (THM-4569 (1), there written `(j, c)`). For every `n`, `u_n = 3^(k_n) v_n + e_n`, where `σ_n = e_n mod 2` and `(k_n, e_n)` updates as follows:

  | `σ` | `β` | new `(k, e)` |
  |---|---|---|
  | 0 | 0 | `(k, e/2)` |
  | 0 | 1 | `(k, (3e + 1 − 3^k)/2)` |
  | 1 | 0 | `(k + 1, (3e + 1)/2)` |
  | 1 | 1 | `(k − 1, (e − 3^(k−1))/2)` |

  * The state is `F_n`-measurable.
  * The parity of `u_n` is `β_n ⊕ σ_n`.
  * `(0, 0)` is absorbing. Absorption means `u_n = v_n`, a merge at equal Terras time with equal odd-step counts.
  * At `k = 0` the coordinate `e` is an integer, and `3^max(0,−k) e` is always an integer.
* **Notation.**
  * Level `h = |k|`.
  * Size `|f| = |e|·3^(−max(k,0))`, so `|f| = |e|` when `k ≤ 0`.
  * Coin `c_n = β_n` if `k_n ≥ 0`, and `c_n = β_n ⊕ σ_n` if `k_n < 0`. This is a fair coin independent of `F_n`.
  * A *flip* is a step with `σ = 1`; flips are exactly the steps that change `k`. A *run* is a maximal stretch of steps with `σ = 0`.
  * `m_h = v_2(3^h − 1)`, which is 1 for odd `h` and `2 + v_2(h)` for even `h ≥ 2`.
  * For `0 < θ < 1`: `s(θ) = 2^θ(1 − √(1 − (3/4)^θ))` and `ρ(θ) = 2^(−θ) s(θ) = 1 − √(1 − (3/4)^θ)`. At `θ = 1/2`, `s = 0.8966` and `ρ = 0.6340`.

## Statements

1. **Local lemmas (PROVED; exhaustive check on 336,040 state-bit pairs).**
   * (a) **Size.** `|f'| ≤ A(c)|f| + 1/2`, with `A(0) = 1/2` and `A(1) = 3/2`. The exception is a departure from `k = 0` (a flip at `k = 0`), where `|f'| ≤ |f|/2 + 1/6` for either coin.
   * (b) **Direction.** A flip at level `h ≥ 1` moves to level `h − 1` iff `c = 1`. So moving toward 0 always costs the factor 3/2, and moving away pays 1/2.
   * (c) **Runs at `k = 0`.** An even `e` is multiplied exactly by 1/2 (`c = 0`) or 3/2 (`c = 1`). Its valuation drops by one either way, so the chain stays at `k = 0` for exactly `v_2(e)` steps, with no additive term.
   * (d) **Runs at level `h ≥ 1`.** Let `w = v_2(e)`.
     * If `w < m_h`, the next valuation is `w − 1` for either coin.
     * If `w ≥ m_h`, each step leaves the regime `{w ≥ m_h}` with probability exactly 1/2, landing at `m_h − 1`. (At `w > m_h` it leaves on `c = 1`; at `w = m_h` it leaves on `c = 0`, and `c = 1` resets `w` to some value `≥ m_h`; `e = 0` counts as `w = ∞`.)
     * Hence, conditionally on the past, a run has length `R ≤ m_h − 1 + G` with `G ~ Geom(1/2)`. At odd `h` it is exactly `R ~ Geom(1/2)`.
2. **Return drift (PROVED).** Let `τ_0 < τ_1 < …` be the successive arrival times at `k = 0`: time 0 if `k_0 = 0`, then every step from `k = ±1` into `k = 0`. After absorption, set `τ_(j+1) = τ_j + 1`. Then for `0 < θ < 1`:

       E[ |e_(τ_(j+1))|^θ | F_(τ_j) ]  ≤  ρ(θ) |e_(τ_j)|^θ + C(θ),   C(θ) < ∞.

3. **Theorem H: almost-sure coalescence (PROVED).** The chain is absorbed almost surely, from every admissible start. Equivalently, for Haar-almost every `y ∈ Z_2`, `T^n(3^(k_0) y + e_0) = T^n(y)` for some `n`.
4. **Rate (PROVED at sketch level).** For every non-absorbed admissible start there are `0 < c ≤ C` with

       c T^(−1/2)  ≤  P(no merge by Terras time T)  ≤  C T^(−1/2) (log T)^2 .

   So the decay exponent is exactly 1/2.
5. **Long-gap correlation in Terras time (PROVED).**
   * The parity streams `(β_n ⊕ σ_n)` of `u` and `(β_n)` of `v` are pairwise independent at distinct times, exactly: the later bit is a fresh coin.
   * At equal times their correlation is `1 − 2P(σ_n = 1)`. This tends to 1, because `σ_n = 0` forever after absorption.
   * After the merge, the odd-step exponent streams of the two orbits coincide exactly, shifted by `−k_0` odd steps.
6. **Corollaries (PROVED from 3 and the cited bridges).**
   * (a) **HYP-9220.** `y` and `y + 1` merge for Haar-almost every `y`. This is the start `(0, 1)`.
   * (b) **THM-4569 (7).** The index `[R_A : R_C] = 1`. So, up to null sets, two 2-adic integers lie on one Collatz grand orbit iff they are related by an affine map `y ↦ 2^a 3^b y + c` with `c ∈ Z[1/6]`. Every element of `Γ_C` merges almost everywhere; for example `y` and `3y` merge (start `(1, 0)`).
   * (c) **Integers, density one.** For each admissible relation and each `K`, absorption by time `K` depends only on the residue mod `2^K` of the driving integer. Hence:
     * for all but a fraction `q(K) ≤ C K^(−1/2) (log K)^2` of residues `n mod 2^K`, `T^t(n) = T^t(n+1)` for some `t ≤ K`;
     * the integers `n` whose trajectories meet those of `n + 1` (likewise `3n`) at equal Terras time have natural density 1;
     * for `n > 2^(K+1)` the merge happens above the cycle `{1, 2}`.
   * (d) **Mersenne switches.** S19's lag-`D` switch is the chain started at `(D, 3^D − 1)` on `(p, q_D)` with `p = 2X − 1` and `q_D = 2X·3^(−D) − 1`. Absorption at `(0, 0)` is the switch `U^i(p) = U^(i+D)(q_D)`. For `D = 1`, after its forced prefix, the chain is absorbed almost surely. So the switching set `S` has `μ_2(S) = 1`, and **HYP-9213 holds** by S18 Proposition 6 with THM-4556 (ii), (iv).
   * (e) **HYP-9214 holds.**
     * The pair `(n, (n−1)/2)` is the chain from `(0, 1)` on `(n, n − 1)`: the step `n − 1 ↦ (n−1)/2` is even, so absorption is an equal-odd-step-time merge.
     * The first-reset-2 condition (`n = 2^(r+1) t − 1`, `r ≥ 1`, `v_2(3^r t − 1) ≥ 2`) is a positive-measure union of residue classes, and `B`-bit sources are uniform on residues mod `2^K` for `K < B`.
     * Hence `liminf_B P_B ≥ P(absorbed by K | condition)`, which tends to 1 as `K` grows.
   * (f) **HYP-9217.**
     * (1) Almost-sure merging and `c T^(−1/2) ≤ q₁(T) ≤ C T^(−1/2)(log T)^2` are PROVED in Terras time. The template total is the Terras time of the merge: the number of 2-adic digits read, THM-4556 (iv).
     * The constant `c₁` is still open.
     * (2) If the any-lag exponent `α` exists, then `α ≥ 1/2`, since `q ≤ q₁`.
7. **One big jump (HEURISTIC + NUMERICAL).** `P(no merge by T) ~ E[J]·√(4/(πT))`, where `J` is the number of excursions of `k` before absorption and `√(4/(πT))` is the tail of one excursion: a first-passage tail `√(2/(πn))` at `n ≈ T/2` flips, since flips are half the steps.
   * For `(0, 1)`: `E[J] = 9.83`, giving 11.10, against `√T q(T) = 11.0` measured.
   * For S19's lag-1 constant 16.7 this predicts `E[J] ≈ 14.8`.

## Proofs

**1.** Write `e` as `N/3^max(0,−k)` and check the four branches.
* For `k ≥ 1`:
  * `(σ, β) = (0, 0)`: `f' = f/2`.
  * `(0, 1)`: `f' = 3f/2 − (1 − 3^(−k))/2`.
  * `(1, 0)`: `k + 1`, `f' = f/2 + 1/(2·3^(k+1))`.
  * `(1, 1)`: `k − 1`, `f' = 3f/2 − 1/2`.
* For `k ≤ −1`, with `c = β ⊕ σ`:
  * `(0, ·)`: `e' = e/2` or `3e/2 + (1 − 3^k)/2`.
  * `(1, 0)`: `k + 1`, `e' = 3e/2 + 1/2`.
  * `(1, 1)`: `k − 1`, `e' = e/2 − 3^(k−1)/2`.
* At `k = 0`, a flip gives `f' = e/2 + 1/6` (to `k = 1`) or `e' = e/2 − 1/6` (to `k = −1`).
* For (d): at level `h ≥ 1` with even `e`, `c = 0` halves `e`. `c = 1` gives valuation `v_2(3e − (3^h − 1)) − 1`, which is `w − 1` if `w < m_h`, `m_h − 1` if `w > m_h`, and `≥ m_h` if `w = m_h`. For `k < 0` this is the same chain in the mirror coordinates `(|k|, −3^|k| e)`, with the roles of `u` and `v` swapped. ∎

**2.** Put `V = |f|^θ s^h`.
* **Flips at `h ≥ 1`.** By 1(a, b) and subadditivity of `x ↦ x^θ`, `E[V' | F] ≤ ½((3/2)^θ s^(h−1) + (1/2)^θ s^(h+1))|f|^θ + ε_h = V + ε_h`. The equality holds because `s` solves `½((3/2)^θ s^(−1) + (1/2)^θ s) = 1`. Here `ε_h = 2^(−θ) s^(h−1)`.
* **Run steps at `h ≥ 1`.** `E[V'] ≤ ½((1/2)^θ + (3/2)^θ) V + ε_h ≤ V + ε_h`. The bracket is 0.966 at `θ = 1/2`.
* **Start of the segment.** From `τ_j`, the run at `k = 0` (1c) gives `E|e|^θ ≤ |e_(τ_j)|^θ` with no additive term. The departure gives `V_dep ≤ s(2^(−θ)|e|^θ + 6^(−θ))`.
* **The excursion.** It is a stretch at levels `≥ 1`. Its flip skeleton `|k|` is a simple random walk started at 1 and killed at 0: the coins are fresh at flip times (optional skipping, 1b). So the expected number of visits to each level `h ≥ 1` is `G(1, h) = 2`. Each visit costs one flip plus a run of expected length at most `m_h + 1` (1d). So the expected total additive error is at most `C_exc = Σ_(h≥1) 2(m_h + 2) 2^(−θ) s^(h−1) < ∞`.
* **The return.** It happens in finite time almost surely. Since `V ≥ 0`, optional stopping at `τ ∧ N` and Fatou give `E[V_τ] ≤ E[V_dep] + C_exc`. At the return `V_τ = |e_τ|^θ`. So `ρ = s 2^(−θ)` and `C = s 6^(−θ) + C_exc`. ∎

**3.**
* **(i) Recurrence.** `k` is a martingale with increments in `{−1, 0, 1}`. Almost surely it either converges or has `limsup = +∞` and `liminf = −∞` (Durrett 4.3.1).
  * If it converges, flips stop from some `N` on, so `u_N` and `v_N` have the same parity vector, and `u_N = v_N` by injectivity of the parity map.
  * Then either `k_N = 0` (absorption), or `v_N = −e_N/(3^(k_N) − 1)`. The latter is an `F_N`-measurable value, which `v_N` (Haar given `F_N`) takes with probability 0.
  * So off absorption, `k` returns to 0 infinitely often (THM-4569 (3)).
* **(ii) Tightness at returns.** By 2, `sup_j E|e_(τ_j)|^θ ≤ max(|e_(τ_0)|^θ, C/(1 − ρ)) < ∞`. By Fatou, `liminf_j |e_(τ_j)| < ∞` almost surely. So almost surely there is an `M` with `|e_(τ_j)| ≤ M` for infinitely many `j`. These are finitely many states, since `e` is an integer at `k = 0`.
* **(iii) Every state at `k = 0` reaches `(0, 0)` with positive probability.**
  * For odd `e > 0`, use `c = 0` (to `(1, (3e+1)/2)`), then `c = 0` halvings at level 1, then `c = 1`. This lands at `(0, (U(e) − 1)/2)`, where `U(e) = oddpart(3e+1)` and `0 ≤ (U(e) − 1)/2 ≤ (3e − 1)/4 < e`.
  * `c = 0` halvings at `k = 0` then reach the odd part. So `e` decreases strictly to 0 along a finite bit string.
  * For `e < 0` the event `{T^n(v) = T^n(v + e)}` is the translate by `e` of the case `|e|` with roles exchanged, and translation preserves Haar measure.
* **(iv) Conclusion.** For `M` fixed, `δ_M = min_(|e|≤M) P_(0,e)(absorbed within N_M steps) > 0`.
  * By Lévy's 0–1 law, `P(absorption | F_n) → 1_absorption` almost surely.
  * On `{|e_(τ_j)| ≤ M` infinitely often`}` the left side is `≥ δ_M` infinitely often, so absorption occurs almost surely there.
  * Taking the union over `M` proves almost-sure absorption. ∎

**4 (sketch).**
* **Lower bound.** Absorption needs `k` to reach 0 from `±1` after the last departure. The skeleton SRW avoids 0 for `T` flips with probability `≥ c T^(−1/2)`, and there are at most `T` flips in `T` steps.
* **Upper bound.**
  * **(R1) The number of excursions `J` before absorption has a geometric tail.** Put `Z_j = |e_(τ_j)|^θ` and take `K > 2C/(1 − ρ)`.
    * Off `{Z ≤ K}` the drift gives `E[Z_(j+1) | F] ≤ ρ' Z_j` with `ρ' = (1 + ρ)/2`. So hitting times of `{Z ≤ K}` have tails `≤ (Z/K) ρ'^g`.
    * From `{Z ≤ K}` absorption within `D` excursions has probability `≥ δ > 0`, by (iii).
    * Attempts spaced `D` apart, with tails restarted from `E Z ≤ ρ^D K + C/(1 − ρ)`, give `P(J > r) ≤ C_1 e^(−c_1 r)`.
  * **(R2) Excursion lengths.** Let `X_j = τ_(j+1) − τ_j` and let `F` be the number of flips in the excursion, so `P(F > n) ≤ C n^(−1/2)`. Then

        X_j ≤ v_2(e_(τ_j)) + 1 + F(2 + log_2 F) + Σ_(i≤F) G_i ,

    with `G_i` conditionally dominated by `Geom(1/2)` (1d, and `m_h ≤ 2 + log_2 h ≤ 2 + log_2 F`). Chernoff on the `G_i` gives

        P(X_j > t | F_(τ_j)) ≤ 1{v_2(e_(τ_j)) ≥ t/3 − 1} + C_2 (log t / t)^(1/2).

  * **(R3) Union bound.**

        P(no merge by T) ≤ P(J > r) + Σ_(j<r) [ P(|e_(τ_j)| ≥ 2^(T/(3r) − 1)) + C_2 (r log T / T)^(1/2) ].

    Take `r = ⌈(log T)/(2c_1)⌉`. The middle term is `≤ r C_3 2^(−θ(T/(3r) − 1))` by Markov and 2, so the total is `O(T^(−1/2) (log T)^2)`. ∎

**5.**
* For `m > m'`, `β_m` is fresh given `F_m ∋ (β_(m'), σ_m)`, so `β_m ⊕ σ_m` is fair and independent of `β_(m')`. The case `m < m'` is symmetric.
* At equal times, `Cov(β, β ⊕ σ) = (1 − 2P(σ = 1))/4`, and both variances are `1/4`.
* After absorption `σ ≡ 0`. ∎

**6.**
* (a) is 3 at `(0, 1)`.
* (b) is THM-4569 (7), whose equivalences are already proved there.
* (c) The parity vector of `n` of length `K` is a bijection on residues mod `2^K` (Terras), and `T^t(n) ≥ n/2^t`.
* (d) Absorption means `T^n(p) = T^n(q_D)` with `q_D` making `D` more odd steps. By the odd-part argument this is `U^i(p) = U^(i+D)(q_D)`.
  * `X = 3^(a−1)` for odd `a` is Haar on `1 + 8Z_2`, because `b ↦ 9^b` is a measure isomorphism `Z_2 → 1 + 8Z_2`.
  * So `q_1` is Haar on a coset of `16Z_2`. Four forced parities are followed by fair ones, and 3 applies from the post-prefix state.
  * Proposition 6 then needs only `μ_2(S) = 1` and the periodicity THM-4556 (iv).
* (e) Absorption at time `t` gives `T^t(n) = T^(t−1)((n−1)/2)` with equal odd counts, hence `U^j(n) = U^j((n−1)/2)`. The merge value exceeds 1 when `t ≪ B`.
* (f) follows from (d) and 4. ∎

## What this does and does not say

* **It is a measure statement about 2-adic integers.** Its integer forms are density-one statements. It says nothing about any particular integer, and the Collatz conjecture is untouched.
* **It answers the long-gap question.**
  * In Terras time the two parity streams are exactly pairwise independent at distinct times.
  * They become identical after an almost-sure merge, whose waiting time has a `T^(−1/2+o(1))` tail.
  * In odd-step time the exponent streams are uncorrelated away from window overlaps (THM-4564, THM-4565), and identical, shifted by the merge lag, after the merge.
* **Why it is not adversarially fragile.**
  * The schedule of flips is predictable, and an adversary who could wait at level 1 while `|f|` grows would break the drift.
  * Lemma 1(d) shows that every continuation of a run at level `h ≥ 1` beyond `m_h − 1` steps is bought with a fresh fair coin. So the expected time per level is bounded, and the additive errors sum.
* **Prior art.** No prior statement was found in two short searches. The nearest are Garner (1985) on consecutive heights, Burson's withdrawn numerical study of first coalescence points (arXiv:2005.09456), and Kontorovich–Lagarias (arXiv:0910.1944), who note that `8n+4` and `8n+5` coalesce. This is not a claim of novelty.
