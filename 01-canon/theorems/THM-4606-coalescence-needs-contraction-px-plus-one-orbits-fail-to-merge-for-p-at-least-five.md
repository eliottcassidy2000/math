---
id: THM-4606
title: "Coalescence needs contraction: for the px+1 Terras maps on Z_2 with p >= 5 odd, the orbits of Haar y and y+e (e a nonzero integer) ever meet with probability q_p(e) < 1, and q_p(e) -> 0 as |e| -> infinity; for p = 1 and p = 3 they meet almost surely (p = 3 is THM-4581). So among the px+1 maps, Haar coalescence holds exactly at the contracting multipliers p < 4. Mechanism: away from departures at zero debt, every pair-chain step multiplies the normalized offset f by 1/2 or p/2 on a fresh fair coin, a walk with drift (1/2)ln(p/4) > 0; departures are visits of a simple random walk to 0, so they are o(n); a growth lemma plus explicit short paths carry every integer start into the escape region"
status: >
  PROVED (elementary probability: the pair-chain table of THM-4581 with 3 replaced by p, optional skipping, Hoeffding's
  inequality, the simple-random-walk return tail P(tau > 2m) = C(2m,m)/4^m, induction). FINITE-EXACT (escape_check.py, ALL PASS):
  the growth lemma for every odd p in [5, 31] and 1 <= |e| <= 2000; explicit escape paths from every start (0, e) with
  1 <= |e| <= 4 for p = 5 and p = 7 (the only multipliers whose growth threshold 4/(p-4) exceeds 1); the step-multiplier table on
  20,000 random states. NUMERICAL (qp_estimate.py, 200,000 Haar pairs each): q_5 = 0.00834 +- 0.00020, q_7 = 0.1084 +- 0.0007,
  q_9 = 0.00226 +- 0.00011, q_11 = 2.5e-5, q_13 = 5e-6. Exact first-merge times from (0, 1) (shortest_merge.py): 3, 11, 4, 11, 18,
  23, 5, 13 for p = 3, 5, ..., 17; for p = 2^a - 1 the reset path of length a + 1 gives q_p >= 2^-(a+1), whence the spikes at p = 3, 7, 15.
  Not yet independently audited.
session: mac-mini-2026-10-08-reframes
source: 05-knowledge/results/coalescence_phase_diagram_20261008.md
scripts:
  - 04-computation/experiments/reframes_20261007/escape_check.py (+ escape_check.out)
  - 04-computation/experiments/reframes_20261007/qp_estimate.py (+ qp_estimate.out), shortest_merge.py (+ shortest_merge.out)
  - 04-computation/experiments/reframes_20261007/coalescence_phase.py, coalescence_visits.py (the phase diagram of HYP-9244)
related:
  - THM-4581 (Haar coalescence for p = 3; the pair chain, its table and its corollaries; this theorem is its converse side)
  - HYP-9244 (the Polya trichotomy: contraction plus recurrence of the debt walk; rank of the debt lattice)
  - HYP-9242 (the orphan law: its exponent 1/2 is the first-return exponent of a rank-one debt walk)
  - Lagarias (1985), Terras (1976) (the parity-vector map of any px+1 map, p odd, is a measure-preserving bijection of Z_2)
---

# THM-4606 — coalescence needs contraction

## Setting

* `p ≥ 1` odd. `T_p(x) = x/2` for even `x` and `(px + 1)/2` for odd `x`, on `Z_2`. Haar measure `μ`.
  * For any odd `p` the parity vector of `y` is a fair i.i.d. coin sequence under `μ` (Terras 1976; Lagarias 1985). Each branch maps its residue class onto `Z_2`.
* **The pair chain** is THM-4581's chain with 3 replaced by `p`. Put `v_n = T_p^n(y)`, `u_0 = p^(k_0) y + e_0`, `u_n = T_p^n(u_0)`. Then `u_n = p^(k_n) v_n + e_n`.
  * Here `σ = e mod 2` (2-adic parity) and `β = v_n mod 2` is a fresh fair coin. The update is:

    | σ | β | new `(k, e)` |
    |---|---|---|
    | 0 | 0 | `(k, e/2)` |
    | 0 | 1 | `(k, (pe + 1 − p^k)/2)` |
    | 1 | 0 | `(k + 1, (pe + 1)/2)` |
    | 1 | 1 | `(k − 1, (e − p^(k−1))/2)` |

  * `(0, 0)` is absorbing; absorption is an equal-time merge `u_n = v_n`.
  * Translations `y, y + e` start at `(0, e)`.
* **Normalized offset.** `f = e · p^(−max(k,0))`, so `f = e` for `k ≤ 0`.
* **Departure.** A flip (`σ = 1`) at `k = 0`.
* **Drift.** `Λ_p = ½ ln(p/4)`, which is positive iff `p ≥ 5`.

## Statements

1. **Step table (PROVED).**
   * In every state, `f' = λf + η` with `λ ∈ {1/2, p/2}` and `|η| ≤ 1/2 + 1/(2p) < 1`.
   * Except at departures, the two coin values give the two multipliers, one each.
   * At a departure both coin values give `λ = 1/2`.
   * Flips move `k` by ±1 on a fresh fair coin. So the values of `k` at successive flips form a simple random walk (SRW).
2. **Escape lemma (PROVED, `p ≥ 5`).**
   * For every `η ∈ (0, 1)` there is `A_p(η)` such that, from every state with `|f_0| ≥ A_p(η)`, with probability at least `1 − η`:
     * the chain is never absorbed; and
     * `|f_n| ≥ 4 e^(Λ_p n/3)` for all `n`.
   * This holds for every `k_0`, and for every `e_0 ∈ Z_2 ∩ Q`.
3. **Integer translations (PROVED, `p ≥ 5`).**
   * For every `e ∈ Z ∖ {0}`, `q_p(e) := μ{y : T_p^n(y) = T_p^n(y + e) for some n} < 1`, and `q_p(e) → 0` as `|e| → ∞`.
   * Meetings at unequal times, `T_p^m(y) = T_p^n(y + e)` with `m ≠ n`, form a null set. So `q_p(e)` is also the probability that the two orbits ever meet.
4. **The dichotomy in p (PROVED).**
   * For `p = 1` and `p = 3`, `y` and `y + e` meet almost surely:
     * `p = 1` with a geometric tail;
     * `p = 3` by THM-4581, with tail `≍ T^(−1/2)`.
   * Hence, among the maps `px + 1` (p odd), Haar coalescence of translates holds iff `p < 4`, i.e. iff the map contracts on average.
5. **Integers (PROVED; the transfer of THM-4581 6(c)).**
   * For every `K ≥ 1` and `e ≠ 0`, the natural density of `n ∈ N` with `T_p^t(n) = T_p^t(n + e)` for some `t ≤ K` equals `P(absorbed by K) ≤ q_p(e)`.
   * For `p = 5` this density is below 0.0087 for every `K`, numerically.
6. **Sheet-blindness (PROVED).**
   * `x ↦ rx` (r odd) conjugates `px + 1` to `px + r` on `Z_2` and preserves Haar measure.
   * So statements 1–5 hold verbatim for `px + r`, with the translation `e` replaced by `re`.
7. **Numerics (NUMERICAL / FINITE-EXACT).**

   | p | q_p = P(y, y+1 meet) | first merge time from (0,1) | exact lower bound |
   |---|---|---|---|
   | 3 | 1 (THM-4581) | 3 | — |
   | 5 | 0.00834 ± 0.00020 | 11 | `2^-11` |
   | 7 | 0.1084 ± 0.0007 | 4 | `2^-4` |
   | 9 | 0.00226 ± 0.00011 | 11 | `2·2^-11` |
   | 11 | 2.5·10^-5 (5 of 200,000) | 18 | `3·2^-18` |
   | 13 | 5·10^-6 (1 of 200,000) | 23 | `3·2^-23` |
   | 15 | — | 5 | `2^-5` |
   | 17 | — | 13 | `2^-13` |

   * For `p = 2^a − 1` the reset path has length `a + 1`: `(0,1) → (1, 2^(a−1)) → … → (1, 1) → (0, 0)`, with `a − 1` halvings. It gives `q_p ≥ 2^(−a−1) = 1/(2(p+1))`.
   * This is why `q_p` is not monotone in `p`.
   * Merges for `p = 5` occurred as late as Terras time 188 in the sample; none later.

## Proofs

**1.**
* The algebra of THM-4581's proof of statement 1 goes through with 3 replaced by `p`. For `k ≥ 1`:
  * `(0,0)`: `f/2`.
  * `(0,1)`: `(p/2)f − (1 − p^(−k))/2`.
  * `(1,0)`: `k + 1`, `f/2 + 1/(2p^(k+1))`.
  * `(1,1)`: `k − 1`, `(p/2)f − 1/2`.
* For `k ≤ −1`:
  * `(0,β)`: `e/2` or `(p/2)e + (1 − p^k)/2`.
  * `(1,0)`: `k + 1`, `(p/2)e + 1/2`.
  * `(1,1)`: `k − 1`, `e/2 − p^(k−1)/2`.
* At `k = 0`:
  * `(0,0)`: `e/2`.
  * `(0,1)`: `pe/2`.
  * `(1,0)`: `k = 1`, `f' = e/2 + 1/(2p)`.
  * `(1,1)`: `k = −1`, `f' = e/2 − 1/(2p)`.
* A flip at `k ≥ 1` goes up on `β = 0` and down on `β = 1`. At `k ≤ −1` it goes toward 0 on `β = 0`. At `k = 0` it goes up on `β = 0`. Each is fair.
* At the stopping times of flips the coins are still fresh and fair (optional skipping). So the skeleton of `k` is a SRW. ∎

**2.**
* **The coin walk.** Define `ξ_n = ln λ_n` at non-departure steps. At departures define `ξ_n = ln(1/2)` if `β_n = 0` and `ln(p/2)` if `β_n = 1`.
  * Given the past, `ξ_n` is fair on `{ln(1/2), ln(p/2)}`: off departures the state assigns the two multipliers bijectively to the two coin values, and `β_n` is fair and fresh.
  * So `(ξ_n)` is i.i.d. with mean `Λ_p`.
  * Always `ln λ_n ≥ ξ_n − ln p · 1{departure at n}`.
* **Departures.** The number `D_n` of departures before time `n` is at most the number `V_n` of visits to 0 among the first `n` positions of the skeleton SRW.
  * Excursions of a SRW from 0 are i.i.d., with `P(τ > 2m) = C(2m, m)/4^m ≥ 1/(2√m)`. So `P(V_n ≥ m) ≤ exp(−(m − 1)/(2√n))`.
  * Fix `ε = Λ_p/(3 ln p)`. For `n ≤ n_0`, `V_n ≤ n_0` trivially. Hence
    `P(∃n : D_n > εn + n_0) ≤ Σ_(n > n_0) exp(−εn/(2√n)) =: δ_1(n_0)`, which tends to 0 as `n_0 → ∞`.
* **Drift.** By Hoeffding (increments in an interval of length `ln p`):
  `P(∃n : Σ_(s<n) ξ_s < (2Λ_p/3)n − C') ≤ Σ_n exp(−2(Λ_p n/3 + C')²/(n ln² p)) =: δ_2(C')`, which tends to 0 as `C' → ∞`.
* **Induction.** On the complement `G` of both bad events:
  * While `|f_s| ≥ 4` we have `λ_s|f_s| ≥ 2`. So `|f_(s+1)| ≥ λ_s|f_s|(1 − 1/(λ_s|f_s|))`, and `ln|f_(s+1)| ≥ ln|f_s| + ln λ_s − 2/(λ_s|f_s|)`.
  * Summing, `ln|f_n| ≥ ln|f_0| + (Λ_p/3)n − C' − n_0 ln p − Σ_(s<n) 2/(λ_s|f_s|)`.
  * Take `|f_0| ≥ A := 4 e^(C' + n_0 ln p + 1) · max(1, 4/(1 − e^(−Λ_p/3)))`.
  * Induction on `n` gives `Σ_(s<n) 2/(λ_s|f_s|) ≤ 1` and `|f_n| ≥ 4 e^(Λ_p n/3)` for all `n`.
  * Absorption requires `f = 0`. So on `G` the chain is never absorbed.
  * Choose `n_0` and `C'` with `δ_1 + δ_2 ≤ η` and set `A_p(η) = A`. ∎

**3.**
* **Growth lemma.** From `(0, e)`, `e ∈ Z ∖ {0}`, follow this path:
  * at `k = 0` with `e` even, take `β = 1`: `e ↦ pe/2`, exactly;
  * at the first odd `e`, depart with `β = 1` to `k = −1`, where `e_1 = (e − 1/p)/2`;
  * at `k = −1`, take `β = 1` on runs: `e ↦ (pe + 1 − p^(−1))/2`;
  * at the first flip, take `β = 0` to return to `k = 0`: `e ↦ (pe + 1)/2`.
* The path either stays at `k = −1` with `|f|` growing geometrically, or returns at `(0, e*)` with
  `e* = (p/2)^(j+1) e_1 + A_j`, where `0 < A_j < (p/2)^(j+1)/(p − 2)`.
* So `|e*| ≥ (p/4)|e| − 1` in both signs:
  * `e > 0`: `e* ≥ pe/4 − 1/4`;
  * `e < 0`: `|e*| ≥ (p/2)(|e|/2 + 1/(2p) − 1/(p−2))`.
* The path never visits `(0, 0)`.
* Above the fixed point `x* = 4/(p − 4)` of `x ↦ (p/4)x − 1`, the distance to `x*` grows by the factor `p/4` per lap.
  * For `p ≥ 9`, `x* < 1 ≤ |e|`.
  * For `p = 5` (`x* = 4`) and `p = 7` (`x* = 4/3`), escape_check.py (B) exhibits paths from every start with `|e| ≤ 4` to `|e'| > x*`.
* Hence from every `(0, e)` there is a positive-probability path, avoiding absorption, to a state with `|f| ≥ A_p(1/2)`. By 2, `q_p(e) ≤ 1 − 2^(−N(e))/2 < 1`.
* For `|e| ≥ A_p(η)`, 2 gives `q_p(e) ≤ η`.
* **Unequal times.** For fixed `m ≠ n` and fixed parity vectors, `T^m(y) = T^n(y + e)` is an affine equation in `y` with slopes `p^o 2^(−m) ≠ p^(o') 2^(−n)` (p odd). So it has at most one solution. Equal times with unequal odd counts are excluded likewise. Equal times with equal odd counts are the chain's merges; a merge at `k_n ≠ 0` pins `v_n` to an `F_n`-measurable value, which is null (THM-4581 3(i)). ∎

**4.**
* For `p = 1`, `k` plays no role. The relation is `u = v + e`, with `e ↦ e/2` on equal parities and `e ↦ (e ± 1)/2` on unequal ones, fair.
  * So `|e|` enters `{0, ±1}` after `O(log|e|)` steps.
  * From `±1` the chain is absorbed with probability 1/2 per step.
* For `p = 3`, apply THM-4581 3 and 4. ∎

**5.**
* Absorption by time `K` depends only on `n mod 2^K`, and residues are uniform (statement 1's Terras property). For integers it is literally the merge `T^t(n) = T^t(n + e)`. ∎

**6.**
* `T_(p,r)(rx) = r T_(p,1)(x)`, and parities agree since `r` is odd.
* The relation `u = p^k v + e` becomes `ru = p^k (rv) + re`. ∎

## Reading

* **The difference process is the size process.**
  * At zero debt, `e = u − v` is multiplied by 1/2 or `p/2` whenever the two parities agree. That is, the offset between two orbits is driven by the same multipliers as the orbits' real size.
  * Two orbits coalesce exactly when that common dynamics contracts. For `p = 3` the geometric mean `√3/2 < 1` contracts. For `p ≥ 5`, `√p/2 > 1` expands.
  * THM-4581's Lyapunov weight `|f|^θ s^|k|`, with `ρ(θ) = 1 − √(1 − (3/4)^θ)`, is built on the factor 3/4. Here `p/4 > 1` reverses it, and `|f|^(−θ)` is the escaping quantity.
* **Rank one.**
  * The debt group of `px + 1` is `p^Z`, of rank one. Its walk is a SRW in flip time.
  * For `p = 3` the merge tail `T^(−1/2)` is that walk's first-return tail. The orphan law's exponent 1/2 (HYP-9242) is the same number.
  * HYP-9244 conjectures the general picture:
    * coalescence = contraction × recurrence of the debt walk;
    * the tail is set by the rank of the debt lattice: exponential, `T^(−1/2)`, `1/log T`, then failure from rank 3 on.
  * Simulations of maps on `Z_3` and `Z_5` show all five regimes.
* **What it is not.**
  * A statement about 2-adic Haar measure and its density shadows: no individual orbit is decided.
  * For integers under `5x + 1`, merges at unequal times inside cycles are not excluded, and the theorem says nothing about divergence.
