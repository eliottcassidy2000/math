---
id: THM-4607
title: "Every expanding Matthews-Watts map fails to coalesce: for T(x) = (m_i x + r_i)/d on x = i mod d with Lambda = (1/d) sum ln(m_i/d) > 0, the pair chain u = M v + e of Haar d-adic orbits escapes from every state of large normalized offset F = |e|/(1 + M): P(ever absorbed) <= eta once F >= A(eta); hence y and y + e merge with probability q(e) -> 0 as |e| -> infinity, for every rank of the debt lattice. Proof: the step multiplier of F has conditional mean log exactly Lambda minus a convexity correction carried only by debt moves near M = 1, and the debt skeleton is a martingale with increments bounded away from 0, whose occupation of any window near log M = 0 is O(sqrt n)"
status: >
  PROVED (elementary: conditional uniformity of both digits, convexity of g(y) = ln(1 + e^y), optional sampling for the
  debt skeleton, a convex-function occupation bound, Azuma-Hoeffding, Markov's inequality on dyadic blocks, induction).
  Qualitative constants. Generalizes THM-4606 (the case d = 2, m = (1, p)) to every digit base and every rank; for d = 2 the
  departure rule is checked in z2_step_table.py (35,007 states; departures give min(m_0, m_1)/2 for both coins; exact, the
  additive term being constant: audit G); it is not needed in the proof. Independently audited 2026-10-08 (audit G): CORRECT WITH
  FIXES (scope of the reading, proof nits); MISTAKE-590. Context: the real-line synchronization dichotomy for random affine IFS
  (Homburg and Kalle, Adv. Math. 2025, arXiv:2207.09987).
session: mac-mini-2026-10-08-rank
source: 05-knowledge/results/debt_rank_least_requirements_20261008.md
scripts:
  - 04-computation/experiments/rank_20261008/z2_step_table.py (+ .out)
related:
  - THM-4606 (px+1 on Z_2; the d = 2 prototype of this argument)
  - HYP-9244 (the Polya trichotomy; this is its expanding row, proved for every map)
  - THM-4608 (the complete classification on Z_2 uses this for its expanding case)
---

# THM-4607 — every expanding Matthews–Watts map fails to coalesce

## Setting

* **The map.** `d ≥ 2`. For each `i mod d` there is a branch `T(x) = (m_i x + r_i)/d` on `x ≡ i mod d`, where:
  * `m_i ≥ 1` is prime to `d`;
  * `d | m_i i + r_i`.
* **The pair chain** is that of HYP-9244. Write `u = M v + e`.
  * `v`'s digit `j` is uniform and independent of the past `F_t`.
  * `u`'s digit is `i = π_t(j)`, where `π_t` is a permutation of `Z/d` determined by `F_t`.
  * The update is `M' = M m_i/m_j` and `e' = (m_i e + r_i − (m_i/m_j) M r_j)/d`.
  * Absorption, an equal-time merge, is `(M, e) = (1, 0)`.
* **Constants.**
  * `Λ = (1/d) Σ_i ln(m_i/d)`.
  * `R = max_i |r_i|`.
  * `Δ = max_(a,b) |ln(m_a/m_b)|`, and `δ_min` is the least nonzero value of `|ln(m_a/m_b)|`.
* **Coordinates.** `x = ln M` and `F = |e|/(1 + M)`.

## Statement

Assume `Λ > 0`.

1. **Escape lemma.**
   * For every `η ∈ (0, 1)` there is `A(η)` such that, from every state with `F_0 ≥ A(η)`, with probability at least `1 − η`, the chain is never absorbed, and `F_n ≥ F_0 e^(Λn/2 − c(η))` for all `n`.
   * This holds for every value of `M_0`.
2. **Translates.**
   * For integer offsets, `q(e) := P(y and y + e ever meet at equal time) ≤ η` once `|e| ≥ 2A(η)`. So `q(e) → 0`, and the map is not coalescent.
   * Meetings with `M ≠ 1` are null, so this is also the probability of any equal-time meeting.
3. **Small offsets.** `q(e) < 1` for every start from which some finite digit word leads, without absorption, to a state with `F ≥ A(1/2)`.
   * For `px + 1` on `Z_2` (with `r_0 = 0`, `r_1 = 1`) this holds for every nonzero integer offset (THM-4606).
   * **Rank zero, any `d`** (all `m_i = m > d`): it holds for every nonzero integer offset. The proof is greedy:
     * if `ē = 0`, then `e' = m e/d` exactly;
     * otherwise `e' = (m e + r_(j+ē) − r_j)/d`. The terms `r_(j+ē) − r_j` sum to 0 over `j`, so some `j` gives an additive term of the sign of `e`, and `|e'| ≥ (m/d)|e|`.
     * So `|e|` grows geometrically along a positive-probability path without reaching 0 (idea from audit G, for `d = 2`).
   * For other expanding maps, `q(e) < 1` at small offsets is checked by simulation only (audit G: `q(e) ≤ 0.28` in every tested case), not proved.

## Proof

**Step multiplier.**
* From `|e'| ≥ (m_i/d)|e| − R(1 + M')/d`, we get `F' ≥ μ F − R/d` with `μ = (m_i/d)(1 + M)/(1 + M')`.
* Put `g(y) = ln(1 + e^y)` and `δ = ln(m_i/m_j)`. Then

      ln μ = (1 − g'(x)) ln(m_i/d) + g'(x) ln(m_j/d) − ρ,
      ρ = g(x + δ) − g(x) − g'(x)δ ∈ [0, (δ²/2) max_[x, x+δ] g''].

* Here `[x, x + δ]` means the segment between `x` and `x + δ`, of either sign.
* Since `g'' ≤ min(1/4, e^(−|y|))`, we have `ρ ≤ (Δ²/2) φ(x)` with `φ(x) = min(1/4, e^(Δ − |x|))`. Also `ρ = 0` when `δ = 0`.
* In rank zero the debt never moves: then `ρ ≡ 0`, `δ_min` is not needed, and the skeleton step below is void.

**Exact conditional mean.**
* Given `F_t`, the digit `j` is uniform, and `i = π_t(j)` is uniform as well.
* `g'(x_t)` is `F_t`-measurable.
* So `E[ln μ_t + ρ_t | F_t] = (1 − g'(x_t))Λ + g'(x_t)Λ = Λ`, and `ψ_t := ln μ_t + ρ_t − Λ` are martingale differences bounded by `W = 2 max_i |ln(m_i/d)| + Δ`.

**The debt skeleton.**
* Let `X_k` be the value of `x` after the `k`-th time the debt moves (`δ ≠ 0`).
* Its increments lie in `[δ_min, Δ]` in absolute value.
* It is a martingale: for every `t`, `E[δ_t 1{the next move after a given move happens at t}] = E[1{no move before t} E[δ_t | F_t]] = 0`, using `δ_t 1{δ_t ≠ 0} = δ_t`.
* If the debt moves only finitely often, use the skeleton stopped at its last move. The bounds below hold for sums over the moves that actually occur.

**Occupation bound.**
* Let `H` be the even convex function with `H'' = κ min(1/4, e^(2Δ − |y|))`. Then `|H'| ≤ B_H` with `B_H = κ((2Δ + ln 4)/4 + 1/4)`.
* Taylor's formula with the integral remainder gives `E[H(X_(k+1)) − H(X_k) | past] ≥ (δ_min²/2) κ φ(X_k)`.
* Also `E H(X_N) − H(X_0) ≤ B_H Δ √N`, since the martingale increments are orthogonal.
* Hence `E Σ_(k<N) φ(X_k) ≤ C_φ √N`, with `C_φ = 2Δ((2Δ + ln 4)/4 + 1/4)/δ_min²`.
* Since there are at most `n` moves by time `n`, `E Σ_(t<n) ρ_t ≤ (Δ²/2) C_φ √n`.

**Bad events.**
* By Markov's inequality on dyadic blocks, `P(Σ_(t<2^m) ρ_t > ε 2^m) ≤ (Δ² C_φ/2) 2^(−m/2)/ε`. This is summable, so with probability `≥ 1 − η/2`, `Σ_(t<n) ρ_t ≤ 2εn + C_1(η)` for all `n`.
* By Azuma–Hoeffding, with probability `≥ 1 − η/2`, `Σ_(t<n) ψ_t ≥ −εn − C_2(η)` for all `n`.

**Induction, as in THM-4606.**
* `μ ≥ μ_min = (m_min/d)(m_min/m_max)`.
* While `F_s ≥ 2R/(d μ_min)` we have `ln F_(s+1) ≥ ln F_s + ln μ_s − 2R/(d μ_s F_s)`.
* With `ε = Λ/6`, on the good event `ln F_n ≥ ln F_0 + Λn/2 − C_1 − C_2 − Σ_(s<n) 2R/(d μ_min F_s)`. The last sum is at most 1 once `F_0 ≥ A(η)`.
* Absorption requires `e = 0`, i.e. `F = 0`. So the chain is never absorbed. ∎

**Statements 2 and 3.**
* From `(1, e_0)` we have `F_0 = |e_0|/2`, and the escape lemma applies once `|e_0| ≥ 2A(η)`.
* Meetings with `M ≠ 1` pin `v` to the fixed point of `v ↦ Mv + e`, which is null (THM-4581 3(i)).
* Statement 3 is the strong Markov property applied at the end of the finite word. ∎

## Reading

* **Least requirements.**
  * The proof needs only the uniformity of both digits given the past (the d-adic Terras property) and that the debt moves by bounded nonzero amounts.
  * It needs no rank condition, no mixing of the coupling state, and no accessibility.
* **Where the debt matters.**
  * In the normalization `F = |e|/(1 + M)`, the debt enters only through the convexity of `ln(1 + M)`, and that is felt only while `M` is near 1.
  * A martingale with increments bounded away from 0 spends only `O(√n)` moves in any window. So expansion wins whatever the debt walk does, whether recurrent (rank ≤ 2) or transient.
* **Within HYP-9244.**
  * For every expanding map, far translates are now proved to escape: `q(e) → 0`.
  * `q(e) < 1` at every offset is proved for `px + 1` and for rank zero. For other expanding maps it is proved at offsets with an escape word, and otherwise checked numerically.
  * The remaining open rows are contracting.
