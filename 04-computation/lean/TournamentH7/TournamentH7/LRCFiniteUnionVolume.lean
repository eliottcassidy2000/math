/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-10-S204)
-/
import Mathlib
import TournamentH7.LRCIntervalBridge

/-!
# Brick (iii): the finite-union volume identity for the exact `m_P` floors

mac-mini-S65 (cont.18/19) built the interval bridge — ONE anchor interval inside `goodSet ∩ safeSet`
gives `0 < witnessG2` — and honestly DEFERRED the **finite-union volume identity**, "the one remaining
measure-theory brick of the witness-floor route", needed for the FULL `m_P` floors of `hsmall3`/`hlarge`
(positivity is not enough there; one needs a quantitative floor).

This file supplies it. Given finitely many **pairwise-disjoint** anchor intervals `Ico aᵢ bᵢ ⊆ [0,1)`, all
contained in `S`, the `slowμ`-measure of `S` is at least the total length:

  `∑ᵢ (bᵢ − aᵢ) ≤ (slowμ S).toReal`.

Specialised to `S = goodSet ∩ safeSet` (with mac-mini's checkable band bounds per anchor) this is exactly
the `m_P` floor shape: `∑ᵢ (bᵢ − aᵢ) ≤ witnessG2 s`.

**Why this is on the critical path.** klein-S234's THM-685 (Kronecker transfer) proves
`|LM(q) − q·μ(S)| ≤ K(S) ≤ Σ v_l` elementarily, so a *measure floor* `μ(S) ≥ μ₀` yields live certificates
at every `q > Σv/μ₀` — collapsing the character-sum program. In klein's words, "the remaining analytic
content of the covering case = measure floors". This brick is the tool that turns exact rational anchor
tables (mac-mini's brick (ii): all 91 hk12 families, min anchor `4/637`) into those floors.

Kernel-pure: no `sorry`, no `native_decide`. Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LonelyRunner
namespace LRC14
namespace IntervalBridge

open DenseCovers MeasureTheory

/-- The `slowμ`-measure of an anchor interval inside `[0,1)` is its length. -/
theorem slowμ_Ico_eq {a b : ℝ} (h0 : 0 ≤ a) (h1 : b ≤ 1) :
    slowμ (Set.Ico a b) = ENNReal.ofReal (b - a) := by
  unfold slowμ
  rw [Measure.restrict_apply measurableSet_Ico]
  have hss : Set.Ico a b ∩ Set.Ico (0 : ℝ) 1 = Set.Ico a b := by
    apply Set.inter_eq_self_of_subset_left
    intro x hx
    exact ⟨le_trans h0 hx.1, lt_of_lt_of_le hx.2 h1⟩
  rw [hss, Real.volume_Ico]

/-- **BRICK (iii) — the finite-union volume identity.** Pairwise-disjoint anchor intervals inside `[0,1)`,
all contained in `S`, force the `slowμ`-measure of `S` to be at least their total length. This is the
quantitative upgrade of `slowμ_toReal_pos_of_Ico_subset` (which only gives positivity from one anchor). -/
theorem slowμ_toReal_ge_sum_of_disjoint_Ico {S : Set ℝ} {n : ℕ} (a b : Fin n → ℝ)
    (h0 : ∀ i, 0 ≤ a i) (h1 : ∀ i, b i ≤ 1) (hab : ∀ i, a i ≤ b i)
    (hsub : ∀ i, Set.Ico (a i) (b i) ⊆ S)
    (hdisj : Pairwise (Function.onFun Disjoint (fun i => Set.Ico (a i) (b i)))) :
    (∑ i, (b i - a i)) ≤ (slowμ S).toReal := by
  have hUmeas : slowμ (⋃ i, Set.Ico (a i) (b i)) = ∑ i, ENNReal.ofReal (b i - a i) := by
    rw [measure_iUnion hdisj (fun i => measurableSet_Ico), tsum_fintype]
    exact Finset.sum_congr rfl (fun i _ => slowμ_Ico_eq (h0 i) (h1 i))
  have hle : ENNReal.ofReal (∑ i, (b i - a i)) ≤ slowμ S := by
    rw [ENNReal.ofReal_sum_of_nonneg (fun i _ => sub_nonneg.mpr (hab i)), ← hUmeas]
    exact measure_mono (Set.iUnion_subset hsub)
  exact (ENNReal.ofReal_le_iff_le_toReal (measure_ne_top slowμ S)).mp hle

/-- **The `m_P` floor from a disjoint anchor table.** Each anchor `Ico aᵢ bᵢ` sits inside `goodSet s.2`
and satisfies mac-mini's checkable band bounds for `s.1`; the anchors are pairwise disjoint. Then the total
anchor length is a lower bound for `witnessG2 s` — the exact shape `hsmall3`/`hlarge` need. -/
theorem witnessG2_ge_sum_of_disjoint_anchors (s : Shape) {n : ℕ} (a b : Fin n → ℝ)
    (hgood : ∀ i, Set.Ico (a i) (b i) ⊆ TournamentH7.GoodSet.goodSet s.2)
    (hposP : ∀ p ∈ s.1, (0 : ℤ) < p)
    (hband : ∀ i, ∀ p ∈ s.1, ∃ j : ℤ,
      (j : ℝ) + 1 / 14 ≤ (p : ℝ) * a i ∧ (p : ℝ) * b i ≤ (j : ℝ) + 13 / 14)
    (h0 : ∀ i, 0 ≤ a i) (h1 : ∀ i, b i ≤ 1) (hab : ∀ i, a i ≤ b i)
    (hdisj : Pairwise (Function.onFun Disjoint (fun i => Set.Ico (a i) (b i)))) :
    (∑ i, (b i - a i)) ≤ witnessG2 s := by
  have hsub : ∀ i, Set.Ico (a i) (b i) ⊆ TournamentH7.GoodSet.goodSet s.2 ∩ safeSet s.1 :=
    fun i => Set.subset_inter (hgood i) (Ico_subset_safeSet_of_bounds hposP (hband i))
  exact slowμ_toReal_ge_sum_of_disjoint_Ico a b h0 h1 hab hsub hdisj

/-- Chaining form: any rational floor `m_P` dominated by the total anchor length is a floor for
`witnessG2`. (`m_P ≤ ∑ lengths ≤ witnessG2`.) -/
theorem witnessG2_ge_of_anchor_floor (s : Shape) {n : ℕ} (a b : Fin n → ℝ) (mP : ℝ)
    (hgood : ∀ i, Set.Ico (a i) (b i) ⊆ TournamentH7.GoodSet.goodSet s.2)
    (hposP : ∀ p ∈ s.1, (0 : ℤ) < p)
    (hband : ∀ i, ∀ p ∈ s.1, ∃ j : ℤ,
      (j : ℝ) + 1 / 14 ≤ (p : ℝ) * a i ∧ (p : ℝ) * b i ≤ (j : ℝ) + 13 / 14)
    (h0 : ∀ i, 0 ≤ a i) (h1 : ∀ i, b i ≤ 1) (hab : ∀ i, a i ≤ b i)
    (hdisj : Pairwise (Function.onFun Disjoint (fun i => Set.Ico (a i) (b i))))
    (hfloor : mP ≤ ∑ i, (b i - a i)) :
    mP ≤ witnessG2 s :=
  le_trans hfloor
    (witnessG2_ge_sum_of_disjoint_anchors s a b hgood hposP hband h0 h1 hab hdisj)

end IntervalBridge
end LRC14
end LonelyRunner
