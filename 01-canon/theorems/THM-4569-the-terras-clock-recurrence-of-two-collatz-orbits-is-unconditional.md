---
id: THM-4569
title: "The Terras clock: for two Collatz orbits related by x = 3^j y + c, the odd-step difference j is exactly a simple random walk run on the predictable disagreement clock N_s = #{r < s : c_r odd}; hence for Haar-random y the pair either merges or j visits every integer infinitely often (unconditionally), merging almost surely is equivalent to an archimedean box-recurrence condition, P(no merge by T) >= (1.59 + o(1)) T^(-1/2) for S19's lag-1 Mersenne pair, and computer-assisted lower bounds give P(lag-1 Mersenne merge) >= 0.3853 and P(y and y+1 merge) >= 0.5861"
status: "PROVED (elementary): the transition table, the random-walk representation, the dichotomy, the rate lower bound, the box reduction (Levy 0-1). Transitions re-verified independently here on 103,158 steps of exact orbits (j moves only at disagreements: 25,724 up, 25,907 down; disagreement density 0.5005). COMPUTER-ASSISTED: the box lower bounds (float64 value iteration from below on finite boxes, monotone; rounding far below the margins), reproduced exactly by an independent implementation on B(3,30) and B(4,100). PROVED from KNOWN facts: the orbit-equivalence form. NUMERICAL: sqrt(T) q(T) = 16.7 (lag-1) and 10.8 (y vs y+1) to T = 2e5. Found by the session's reader of openai/math #287/#248 (operator algebras and groups)."
session: mac-mini-2026-10-07-oaimath3
source: 05-knowledge/results/oai3_two_orbits_twos_and_threes_20261007.md
scripts:
  - 04-computation/experiments/oai3_20261007_readers/groups/ (tclock_chain.py, tclock_certified.c, tclock_box.py, tclock_boxbound.py, bs12_iid_returns.py, + .out)
  - 04-computation/experiments/oai3_20261007_tclock_vi.py (independent value iteration)
related:
  - THM-4564 (the same pair in the odd-step clock: alignment coupling; disagreement density 1/2 <=> variance 4 per odd step)
  - HYP-9217 (the debt recurrence law; its recurrence half is now PROVED here, its rate half remains), HYP-9218 (full-past formulation refuted; a coarse-conditioning replacement may inform the rate), HYP-9220 (y ~ y+1 a.e.)
  - THM-4556 (v) (certified Mersenne switch density 15719/131072 = 0.1199; S19: 0.1740 at K = 31)
  - Terras (1976) (the parity-vector map is a measure-preserving bijection of Z_2); Lagarias (1985, 1990); Matui (2015) (topological full group of the full shift = Thompson's V); Dougherty-Jackson-Kechris (hyperfiniteness of tail relations)
---

# THM-4569 — the Terras clock

## Setting

* `T(x) = x/2` for even `x`, and `(3x+1)/2` for odd `x`, on `Z_2`.
* Two orbits `x_s = T^s(x)`, `y_s = T^s(y)`, related by `x_s = 3^(j_s) y_s + c_s`, with `c_s ∈ Z_2` (in `Z[1/3]` when `c_0 ∈ Z[1/3]`).
* The initial relation parameters are fixed, or determined by a finite already revealed parity prefix after which the source is Haar. Arbitrary dependence of `c_0` on unrevealed source digits would not make the disagreement clock predictable.
* `p_s = y_s mod 2`, and `ε_s = c_s mod 2`.
* For Haar `y` the parities `p_s` are i.i.d. fair coins (Terras), and `ε_s` is a function of the past.

## Statements

1. **Transitions (PROVED).**
   * `x_s mod 2 = p_s ⊕ ε_s`.
   * The state `(j, c)` updates as follows:

     | `(p, ε)` | `j'` | `c'` |
     |---|---|---|
     | (0,0) | `j` | `c/2` |
     | (0,1) | `j+1` | `(3c+1)/2` |
     | (1,0) | `j` | `(3c+1−3^j)/2` |
     | (1,1) | `j−1` | `(c−3^(j−1))/2` |

2. **Random-walk representation (PROVED).**
   * `j_(s+1) − j_s = ε_s(1 − 2p_s)`, where `ε_s` is predictable and `1 − 2p_s` is a fresh fair sign.
   * By optional skipping, `j_s = j_0 + S_(N_s)`, with `S` a simple random walk and `N_s = Σ_(r<s) ε_r`.
   * So `j` is a martingale and `Var j_s = E N_s`.
3. **Dichotomy (PROVED).**
   * Up to a Haar-null set of isolated value coincidences, any merge is an equal-T-time identity hit of `(j, c) = (0, 0)`. This statement does not assert that a merge occurs almost surely; that is the box-recurrence question below, now PROVED by THM-4581.
   * If the orbits never merge, then `N_s → ∞`. Otherwise the parity vectors would eventually agree forever, and injectivity of the parity map (Terras) would force equality.
   * Hence **almost surely the pair either merges or `j` visits every integer infinitely often**. This is the recurrence half of S19's next step 2, unconditional, with no input from HYP-9218.
4. **Rate lower bound (PROVED).** `P(no merge by T) ≥ P(S avoids −j_* for T steps)`. For S19's lag-1 pair (`j_* = 2` after its forced prefix `(1,2) → (1,2) → (1,1) → (2,2) → (2,1)`) this gives `q₁(T) ≥ (2√(2/π) + o(1)) T^(−1/2) ≈ 1.596 T^(−1/2)`, so HYP-9217's exponent cannot exceed 1/2. For `x = y + 1` (`j_* = 1`) it gives `(√(2/π) + o(1)) T^(−1/2) ≈ 0.797 T^(−1/2)`. (Corrected after audit C: the earlier "0.80" rounded a lower bound up.)
5. **Box reduction (PROVED) and computer-assisted bounds.**
   * Let `B(J, C) = {|j| ≤ J, |c·3^max(0,−j)| ≤ C·3^|j|}`.
   * By Lévy's 0–1 law, the chain merges a.s. on the event that it visits `B(1,10)` infinitely often. So **almost-sure merging ⟺ box recurrence**, an archimedean condition (`|c|` does not blow up along the returns of `j`).
   * Value iteration from below on `B(9,100)` (11.8M states) gives:
     * merge probability ≥ 0.3446 from every state of `B(1,10)`;
     * `P(y and y+1 merge)` ≥ 0.5861;
     * S19's lag-1 Mersenne pair merges with probability ≥ 0.3853.
   * Each merge is decided by the exponent `a` modulo a power of 2 (S19 §1.1). So the **lower density of odd `a` with `σ(2^a − 1) = σ(2^(a−1) − 1)` is at least 0.3853**, against the previously certified 0.1740 (S19, `K = 31`, exhaustive).
6. **The real place (PROVED).** For `j ≥ 0`, `ρ = c/3^j` obeys `ρ' = ρ/2` or `ρ' = (3/2)ρ − 1/2`, up to an additive error `≤ 3^(−j)/2` (exact in rows (0,0) and (1,1)); so the statement holds as `j → +∞`. This is the real Terras iterated function system: mean multiplier 1 and tail index 1, i.e. THM-4564 (6) seen in this clock.
7. **Orbit-equivalence form (PROVED from KNOWN facts).**
   * Let `R_C` be the T-grand-orbit relation on `(Z_2, Haar)`, and `R_A` the orbit relation of `Γ_C = Z[1/6] ⋊ ⟨2, 3⟩` restricted to `Z_2`. Then `R_C ⊆ R_A`.
   * Both are ergodic, hyperfinite and of type III_(1/2). Via Lagarias's conjugacy, the Deaconu–Renault groupoid of `(Z_2, T)` is the full-2-shift groupoid. Its C*-algebra is `O_2` and its topological full group is Thompson's `V` (Matui 2015; also Nekrashevych 2004). Its orbit relation `R_C` is the ergodic hyperfinite III_(1/2) relation, with `L(R_C) ≅ R_(1/2)`, the Powers factor. (Corrected after audit C: `O_2` and `V` belong to the groupoid, not to the measured relation.)
   * The index `[R_A : R_C]` is a.e. constant. The following are equivalent:
     * the index is 1;
     * `y` and `y + 1` merge for a.e. 2-adic `y` (HYP-9220);
     * every element of `Γ_C` merges a.e.;
     * the chain started from `(0,1)` merges a.s.;
     * that chain is box-recurrent.
   * Any of these gives HYP-9217 (1)'s almost-sure part, hence `μ_2(S) = 1` and HYP-9213 (S18, Proposition 6).

## How the two orbits' exponents correlate over long gaps (with THM-4564)

The [inverse-clock bridge](../../05-knowledge/results/collatz_terras_inverse_clock_bridge_20261007.md)
supplies the separate transfer: recurrent count ties give, with conditional
probability at least one half per disjoint trial, matched odd-event times
within one step. Thus the corresponding synchronous valuation-sum debt
returns to `{-1,0,1}` infinitely often on the nonmerge alternative. This
does not bound the translation coordinate or prove box recurrence.

* In the Terras clock, every correlation between the two parity streams, at every gap, is carried by the predictable bit `ε_s = c_s mod 2`. It decides *when* the odd-step difference moves, never *which way*.
* Coupling can change the clock: the disagreement density, measured at 0.5006, is equivalent to THM-4564's variance 4 per odd step. It cannot destroy recurrence.
* The original full-past formulation of HYP-9218 is refuted: the conditioning already determines the depth. THM-4581 uses a different return-weight argument, rather than that hypothesis. The fitted constant16.7 remains distinct from the bounds below; see the [overlap integration](../../05-knowledge/results/collatz_overlap_kernel_integration_20261007.md).

---

## Update (2026-10-07, same session): box recurrence and index 1 PROVED (THM-4581)

* **Box recurrence (5) holds, so the chain merges almost surely from every admissible start** (THM-4581 (3)).
  * The weight `|c·3^(−max(j,0))|^θ s^|j|`, with `s = 2^θ(1 − √(1 − (3/4)^θ))`, has its homogeneous term balanced on flips; the additive error is bounded separately. A move of `j` toward 0 is exactly the ×3/2 branch.
  * Runs cost fresh coins per continuation beyond `v_2(3^|j| − 1) − 1` steps.
  * Together these give `E|c_return|^θ ≤ 0.634 |c|^θ + C` along the returns of `j` to 0. That is the archimedean control (5) asked for.
* **Consequences.** All equivalent statements in (7) hold, in the grand-orbit sense. THM-4581 (6b) writes out the reduction that (7) left implicit: for `g(y) = 2^a 3^b y + c_0/(2^m 3^l)`, route through `y' = 2^(M+a) y` and `3^l 2^M g(y)`. Equal-time merging fails when the multiplier has a nontrivial power of 2; `y` and `2y` a.s. never merge at equal time.
  * the index `[R_A : R_C] = 1`;
  * `y ~ y + 1` almost everywhere (HYP-9220);
  * every element of `Γ_C` merges almost everywhere.
* **The rate lower bound (4) is sharp up to logarithms:** `P(no merge by T) ≤ C T^(−1/2) (log T)^2` (THM-4581 (4), PROVED at sketch level).
* **The value-iteration lower bounds of (5) are superseded:** the merge probabilities are 1. They remain valid finite-box certificates.

**Audit (2026-10-07, independent audit C).**
* CONFIRMED: transitions (100,000 steps against direct orbits), (2), (3), box validity, and the certificates. B(3,30) and B(4,100) were reproduced exactly by independent code. B(9,100) gives 0.586175, 0.385315 and 0.344618, above the claimed 0.5861, 0.3853 and 0.3446.
* Corrected above: (4) constants, (7) groupoid wording, (6) range of validity (MISTAKE-583).
