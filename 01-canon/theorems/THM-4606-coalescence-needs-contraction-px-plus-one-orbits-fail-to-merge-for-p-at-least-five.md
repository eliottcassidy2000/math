---
id: THM-4606
title: "Coalescence needs contraction: for the px+1 Terras maps on Z_2 with p >= 5 odd, the orbits of Haar y and y+e (e a nonzero integer) ever meet with probability q_p(e) < 1, and q_p(e) -> 0 as |e| -> infinity; for p = 1 and p = 3 they meet almost surely (p = 3 is THM-4581). So among the px+1 maps with p >= 1 odd, Haar coalescence holds exactly at the contracting multipliers p < 4 (the Matthews-Watts threshold m_0 m_1 < d^d = 4). Mechanism: away from departures at zero debt, every pair-chain step multiplies the normalized offset f by 1/2 or p/2 on a fresh fair coin, a walk with drift (1/2)ln(p/4) > 0; departures are visits of a simple random walk to 0, so they are o(n); a growth lemma plus explicit short paths carry every integer start into the escape region"
status: >
  PROVED (elementary probability: the pair-chain table of THM-4581 with 3 replaced by p, optional skipping, Hoeffding's
  inequality, the simple-random-walk return tail P(tau > 2m) = C(2m,m)/4^m, induction). The escape lemma is qualitative: its
  constant A_p(1/2) is astronomically large (about 10^981058 for p = 5; audit E). FINITE-EXACT (escape_check.py, ALL PASS):
  the growth lemma for every odd p in [5, 31] and 1 <= |e| <= 2000 (audit E: odd p <= 63, |e| <= 5000, 300,000 starts); explicit
  escape paths from every start (0, e) with 1 <= |e| <= 4 for p = 5 and p = 7 (the only multipliers whose growth threshold
  4/(p-4) is at least 1); the step-multiplier table on 20,000 random states (audit E: 200,000 states, and the chain against
  actual big-integer orbits). Rigorous lower bounds by exact enumeration (audit E): q_5 >= 18401/2^22 = 0.004387,
  q_7 >= 449933/2^22 = 0.10727, q_9 >= 36571/2^24 = 0.002180, q_11 >= 2.356e-5, q_13 >= 1.059e-6. NUMERICAL: the author's
  200,000 chain paths per p (random coins, equal in law to Haar pairs) and audit E's 10.3 million actual orbit pairs plus
  exact-mass-and-tail estimates: q_5 = 0.00833 +- 0.00005, q_7 = 0.10857 +- 0.00005, q_9 = 0.00225 +- 0.00001,
  q_11 = 2.43e-5 +- 0.04e-5, q_13 = (1.6 +- 0.4)e-6. First merge times from (0, 1), confirmed by audit E's rigorously pruned
  search: 3, 11, 4, 11, 18, 23, 5, 13, 8, 18, 6 for p = 3, 5, 7, 9, 11, 13, 15, 17, 21, 23, 31; none up to depth 26 for
  p = 19, 25, 27, 29 (none up to 30 for p = 19). Independently audited 2026-10-08 (audit E): corrections applied (statement 5
  restricted to e = +-1 for the numerical bound, q_11 and q_13, the run argument in proof 3, the p = 3 tail typing, the scope
  of sheet-blindness, labels); MISTAKE-589.
session: mac-mini-2026-10-08-reframes
source: 05-knowledge/results/coalescence_phase_diagram_20261008.md
scripts:
  - 04-computation/experiments/reframes_20261007/escape_check.py (+ escape_check.out)
  - 04-computation/experiments/reframes_20261007/qp_estimate.py (+ qp_estimate.out), shortest_merge.py (+ shortest_merge.out)
  - 04-computation/experiments/reframes_20261007/coalescence_phase.py, coalescence_visits.py (the phase diagram of HYP-9244)
  - 04-computation/experiments/reframes_20261007/audit_E/ (REPORT.md; table, escape, growth, density-transfer and direct-orbit checks; exact enumerations)
related:
  - THM-4581 (Haar coalescence for p = 3; the pair chain, its table and its corollaries; this theorem is its converse side)
  - HYP-9244 (the Polya trichotomy: contraction plus recurrence of the debt walk; rank of the debt lattice)
  - HYP-9242 (the orphan law: its exponent 1/2 is the first-return exponent of a rank-one debt walk)
  - Lagarias (1985), Terras (1976) (the parity-vector map of any px+1 map, p odd, is a measure-preserving bijection of Z_2)
  - Kontorovich and Lagarias, arXiv:0910.1944 (the 2-adic 3x+1 and 5x+1 maps are measure-theoretically conjugate)
  - Matthews and Watts (Acta Arith. 43, 1984) (the contraction criterion); Homburg and Kalle (Adv. Math. 2025) (synchronization of random affine IFS)
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
   * In every state, `f' = λf + η` with `λ ∈ {1/2, p/2}` and `|η| ≤ 1/2` (sharp).
   * For `p ≥ 3`, except at departures, the two coin values give the two multipliers, one each. For `p = 1` the multipliers coincide.
   * At a departure both coin values give `λ = 1/2`.
   * Flips move `k` by ±1 on a fresh fair coin. So the values of `k` at successive flips form a simple random walk (SRW).
2. **Escape lemma (PROVED, `p ≥ 5`; qualitative).**
   * For every `η ∈ (0, 1)` there is `A_p(η)` such that, from every state with `|f_0| ≥ A_p(η)`, with probability at least `1 − η`:
     * the chain is never absorbed; and
     * `|f_n| ≥ 4 e^(Λ_p n/3)` for all `n`.
   * This holds for every `k_0`, and for every `e_0 ∈ Z_2 ∩ Q`.
   * The constant the proof produces is astronomically large (about `10^981058` for `p = 5`, `η = 1/2`). The lemma gives existence, not a numerical bound.
3. **Integer translations (PROVED, `p ≥ 5`).**
   * For every `e ∈ Z ∖ {0}`, `q_p(e) := μ{y : T_p^n(y) = T_p^n(y + e) for some n} < 1`, and `q_p(e) → 0` as `|e| → ∞`.
   * Meetings at unequal times, `T_p^m(y) = T_p^n(y + e)` with `m ≠ n`, form a null set. So `q_p(e)` is also the probability that the two orbits ever meet.
4. **The dichotomy in p (PROVED).**
   * For `p = 1` and `p = 3`, `y` and `y + e` meet almost surely:
     * `p = 1` with a geometric tail;
     * `p = 3` by THM-4581, with `c T^(−1/2) ≤ P(no merge by T) ≤ C T^(−1/2)(log T)²` (the upper bound at sketch level).
   * Hence, among the maps `px + 1` with `p ≥ 1` odd, Haar coalescence of translates holds iff `p < 4`, i.e. iff the map contracts on average. This is the Matthews–Watts contraction threshold `m_0 m_1 < d^d = 4`.
5. **Integers (PROVED; the transfer of THM-4581 6(c)).**
   * For every `K ≥ 1` and `e ≠ 0`, the natural density of `n ∈ N` with `T_p^t(n) = T_p^t(n + e)` for some `t ≤ K` equals `P(absorbed by K) ≤ q_p(e)`.
   * For `p = 5` and `e = ±1` this density is below 0.0087 for every `K`, numerically (`q_5 = 0.00833 ± 0.00005`).
   * For other `e` it can be much larger: `q_5(3) ≥ 1/32`, since the whole class `n ≡ 16 (mod 32)` merges with `n + 3` at `t = 5` (e.g. `T^5(16) = T^5(19) = 3`). At `K = 14` the exact densities are 0.0444 for `e = ±3`, 0.0230 for `e = 6` and 0.00995 for `e = 5` (audit E).
6. **Sheet-blindness (PROVED, with this scope).**
   * `x ↦ rx` (r odd) conjugates `px + 1` to `px + r` on `Z_2` and preserves Haar measure.
   * So statements 1–5 hold for `px + r` after conjugating back: translations by `re` with `e` a nonzero integer behave exactly as translations by `e` for `px + 1`. In the `px + r` coordinates the additive error is `|η| ≤ |r|/2`.
   * Translations not divisible by `r` correspond to non-integer `e`. The growth lemma's integrality does not cover them, so they are not covered here.
     * Audit E's breadth-first search found escape paths for `e = ±1/3, ±2/3` (`p = 5, 7`) and `e = ±1/7, ±2/7` (`p = 5`). The extension is very likely true but unproved.
7. **Numerics (NUMERICAL / FINITE-EXACT).**

   | p | q_p = P(y, y+1 meet) | first merge time from (0,1) | rigorous lower bound (exact enumeration, audit E) |
   |---|---|---|---|
   | 3 | 1 (THM-4581) | 3 | — |
   | 5 | 0.00833 ± 0.00005 | 11 | `18401/2^22 = 0.004387` |
   | 7 | 0.10857 ± 0.00005 | 4 | `449933/2^22 = 0.10727` |
   | 9 | 0.00225 ± 0.00001 | 11 | `36571/2^24 = 0.002180` |
   | 11 | (2.43 ± 0.04)·10^-5 | 18 | `101183/2^32 = 2.356·10^-5` |
   | 13 | (1.6 ± 0.4)·10^-6 | 23 | `9099/2^33 = 1.059·10^-6` |
   | 15 | — | 5 | `2^-5` (reset path) |
   | 17 | — | 13 | — |

   * For `p = 2^a − 1` the reset path has length `a + 1`: `(0,1) → (1, 2^(a−1)) → … → (1, 1) → (0, 0)`, with `a − 1` halvings. It gives `q_p ≥ 2^(−a−1) = 1/(2(p+1))`. Audit E checked it for `a = 2..7`.
   * This is why `q_p` is not monotone in `p`.
   * In audit E's 1.9 million `p = 5` samples, merges occurred as late as Terras time 352, and none after 400.

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
* All additive terms are at most 1/2 in absolute value. The value 1/2 is attained.
* A flip at `k ≥ 1` goes up on `β = 0` and down on `β = 1`. At `k ≤ −1` it goes toward 0 on `β = 0`. At `k = 0` it goes up on `β = 0`. Each is fair; in every case `k' = k + 1 − 2β`.
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
* **The run at `k = −1` is finite** (audit E).
  * In the variable `E = pe` the run map is `E ↦ (pE + p − 1)/2`. Its length is exactly `v_2((p − 2)·pe_1 + p − 1)`.
  * So it is infinite only at the fixed point `pe_1 = −(p − 1)/(p − 2)`.
  * For odd integer `e`, `pe_1 = (pe − 1)/2` is an integer, while `(p − 1)/(p − 2) = 1 + 1/(p − 2)` is not an integer for `p ≥ 5`. Hence for integer `e` the run is finite and the path returns to `(0, e*)`.
  * For non-integer `e` the run can be infinite. For example, `p = 5`, `e = −1/3` sits at `(−1, −4/15)` forever, with `f` constant.
* **The bound.** The return is `e* = (p/2)^(j+1) e_1 + A_j`, where `0 < A_j < (p/2)^(j+1)/(p − 2)`. So `|e*| ≥ (p/4)|e| − 1` in both signs:
  * `e > 0`: `e* ≥ pe/4 − 1/4`.
  * `e < 0`: `|e*| ≥ (p/2)(|e|/2 + 1/(2p) − 1/(p−2)) = (p/4)|e| + 1/4 − p/(2(p−2))`. This is `≥ (p/4)|e| − 1` because `p/(2(p−2)) ≤ 5/4` exactly when `p ≥ 10/3`; it fails at `p = 3`.
* The path never visits `(0, 0)`.
* **Iteration.** Above the fixed point `x* = 4/(p − 4)` of `x ↦ (p/4)x − 1`, the distance to `x*` grows by the factor `p/4` per lap.
  * For `p ≥ 9`, `x* < 1 ≤ |e|`.
  * For `p = 5` (`x* = 4`) and `p = 7` (`x* = 4/3`), escape_check.py (B) exhibits paths from every start with `|e| ≤ 4` to `|e'| > x*`.
* Hence from every `(0, e)` there is a positive-probability path, avoiding absorption, to a state with `|f| ≥ A_p(1/2)`. By 2, `q_p(e) ≤ 1 − 2^(−N(e))/2 < 1`.
* For `|e| ≥ A_p(η)`, 2 gives `q_p(e) ≤ η`.
* **Unequal times and unequal odd counts.** For fixed times `m, n` and fixed parity vectors with `p^o 2^(−m) ≠ p^(o') 2^(−n)`, the equation `T^m(y) = T^n(y + e)` is affine in `y` with distinct slopes, so it has at most one solution. This covers `m ≠ n` (`p` odd), and also equal times with unequal odd counts, i.e. a merge at `k_n ≠ 0`. Countably many such equations give a null set. The remaining case, equal times and equal odd counts, is the chain's absorption. ∎

**4.**
* For `p = 1`, `k` plays no role. The relation is `u = v + e`, with `e ↦ e/2` on equal parities and `e ↦ (e ± 1)/2` on unequal ones, fair.
  * So `|e|` enters `{0, ±1}` after `O(log|e|)` steps.
  * From `±1` the chain is absorbed with probability 1/2 per step.
* For `p = 3`, apply THM-4581 3 and 4. ∎

**5.**
* Absorption by time `K` depends only on `n mod 2^K`, and residues are uniform (statement 1's Terras property).
* For integers, absorption is literally a merge `T^t(n) = T^t(n + e)` with equal odd counts.
* An equal-time integer merge with unequal odd counts solves an affine equation with distinct slopes. So it occurs for at most one integer per (time, parity word); it has density zero and does not change the density (audit E found none for `n ≤ 2^16`, `t ≤ 40`). ∎

**6.**
* `T_(p,r)(rx) = r T_(p,1)(x)`, and parities agree since `r` is odd.
* The relation `u = p^k v + e` becomes `ru = p^k (rv) + re`. ∎

## Reading

* **The difference process is the size process.**
  * At zero debt, `e = u − v` is multiplied by 1/2 or `p/2` whenever the two parities agree. That is, the offset between two orbits is driven by the same multipliers as the orbits' real size.
  * Within the `px + 1` family, two orbits coalesce exactly when that common dynamics contracts. For `p = 3` the geometric mean `√3/2 < 1` contracts. For `p ≥ 5`, `√p/2 > 1` expands.
  * Outside the family, contraction alone is not enough: HYP-9244 has a contracting rank-3 map on `Z_5` that numerically does not coalesce.
  * HEURISTIC: THM-4581's Lyapunov weight `|f|^θ s^|k|`, with `ρ(θ) = 1 − √(1 − (3/4)^θ)`, is built on the factor 3/4. Here `p/4 > 1` reverses it, and `|f|^(−θ)` behaves as the escaping quantity. That is not proved: at departures it grows by `2^θ` for both coins. The proof above uses the coin walk instead.
* **A pair property that measure conjugacy does not see.**
  * Kontorovich and Lagarias note that the 2-adic `3x + 1` and `5x + 1` maps are measure-theoretically conjugate (both are the shift in parity coordinates).
  * THM-4581 and this theorem separate them by a property of pairs `(y, y + e)`. The additive structure is not preserved by the conjugacy.
  * This matches the two-point-motion picture for random iterated function systems, where synchronization goes with a negative Lyapunov exponent (e.g. Homburg–Kalle 2025).
* **Rank one.**
  * The debt group of `px + 1` is `p^Z`, of rank one. Its walk is a SRW in flip time.
  * For `p = 3` the merge tail `T^(−1/2)` is that walk's first-return tail (the upper bound at sketch level). HEURISTIC: the orphan law's measured exponents (0.42–0.63, HYP-9242) sit near the same number.
  * HYP-9244 conjectures the general picture:
    * coalescence = accessibility × contraction × recurrence of the debt walk. For `px + 1` with integer offsets there is no congruence obstruction, so accessibility is automatic here;
    * the tail is set by the rank of the debt lattice: exponential, `T^(−1/2)`, `1/log T`, then failure from rank 3 on.
  * Simulations of maps on `Z_3` and `Z_5` show all five regimes.
* **What it is not.**
  * A statement about 2-adic Haar measure and its density shadows: no individual orbit is decided.
  * For integers under `5x + 1`, merges at unequal times inside cycles are not excluded, and the theorem says nothing about divergence.
  * Novelty: no prior statement was found in brief searches (audit E). The ingredients are classical.
