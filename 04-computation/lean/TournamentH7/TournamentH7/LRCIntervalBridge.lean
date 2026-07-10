/-
  TournamentH7.LRCIntervalBridge — the interval bridge for the repaired witness-floor legs
  (mac-mini-2026-07-09-S65 cont.18).

  The engine (S65 cont.16) computes exact rational safe intervals; this file is the Lean
  consumer: an explicit interval `Ico a b` inside `goodSet ∩ safeSet` forces
  `0 < witnessG2` — the shape of the `hk12` leg of `lrc14_from_repaired_nodes`, and the
  positivity core of `hsmall3`/`hlarge` (whose full m_P floors additionally need the
  finite-union volume computation).

  Three shape-independent lemmas:
  * `slowμ_toReal_pos_of_Ico_subset` — a positive-length subinterval of `[0,1)` inside `S`
    forces `0 < (slowμ S).toReal` (monotonicity needs NO measurability; finiteness from the
    probability instance).
  * `Ico_subset_safeSet_of_bounds` — the checkable band-containment condition: if for every
    `p ∈ P` there is `j : ℤ` with `j + 1/14 ≤ p·a` and `p·b ≤ j + 13/14`, then
    `Ico a b ⊆ safeSet P` (positive `p`; `fract(px) = px − j` on the band).
  * `witnessG2_pos_of_anchor` — the composition: interval inside `goodSet` + band bounds
    ⟹ `0 < witnessG2 s`.
  Named next bricks (not claimed): `goodSet_univ_of_card_le_two` (|E| ≤ 2 ⟹ goodSet = univ:
  fract(−y) = 1 − fract(y) makes double membership in `(0, 1/7]` impossible), and the
  finite-union volume identity for the exact m_P floors.
-/
import Mathlib
import TournamentH7.LRCFourteenSkeleton

namespace LonelyRunner
namespace LRC14
namespace IntervalBridge

open DenseCovers MeasureTheory

/-- A positive-length subinterval of `[0,1)` inside `S` forces positive `slowμ`-measure
(`toReal` form, as `witnessG2` consumes it). -/
theorem slowμ_toReal_pos_of_Ico_subset {S : Set ℝ} {a b : ℝ}
    (hsub : Set.Ico a b ⊆ S) (h0 : 0 ≤ a) (h1 : b ≤ 1) (hab : a < b) :
    0 < (slowμ S).toReal := by
  have hmono : slowμ (Set.Ico a b) ≤ slowμ S := measure_mono hsub
  have hIco : slowμ (Set.Ico a b) = ENNReal.ofReal (b - a) := by
    unfold slowμ
    rw [Measure.restrict_apply' measurableSet_Ico]
    have : Set.Ico a b ∩ Set.Ico (0 : ℝ) 1 = Set.Ico a b := by
      apply Set.inter_eq_self_of_subset_left
      exact Set.Ico_subset_Ico h0 h1
    rw [this, Real.volume_Ico]
  have hpos : 0 < slowμ S :=
    lt_of_lt_of_le (by rw [hIco]; exact ENNReal.ofReal_pos.mpr (by linarith)) hmono
  have hfin : slowμ S ≠ ⊤ := measure_ne_top slowμ S
  exact ENNReal.toReal_pos hpos.ne' hfin

/-- **The checkable band-containment condition.**  If every `p ∈ P` is positive and admits an
integer `j` with `j + 1/14 ≤ p·a` and `p·b ≤ j + 13/14`, then on all of `Ico a b` every
`fract (p·x)` stays in `[1/14, 13/14]`: the interval sits inside `safeSet P`.  (On the band,
`⌊p·x⌋ = j`, so `fract (p·x) = p·x − j`.) -/
theorem Ico_subset_safeSet_of_bounds {P : List ℤ} {a b : ℝ}
    (hpos : ∀ p ∈ P, (0 : ℤ) < p)
    (hband : ∀ p ∈ P, ∃ j : ℤ, (j : ℝ) + 1 / 14 ≤ (p : ℝ) * a ∧ (p : ℝ) * b ≤ (j : ℝ) + 13 / 14) :
    Set.Ico a b ⊆ safeSet P := by
  intro x hx
  intro p hp
  obtain ⟨j, hj1, hj2⟩ := hband p hp
  have hpR : (0 : ℝ) < (p : ℝ) := by exact_mod_cast hpos p hp
  have hxa : (p : ℝ) * a ≤ (p : ℝ) * x := by
    exact mul_le_mul_of_nonneg_left hx.1 hpR.le
  have hxb : (p : ℝ) * x < (p : ℝ) * b := by
    exact mul_lt_mul_of_pos_left hx.2 hpR
  have hlo : (j : ℝ) + 1 / 14 ≤ (p : ℝ) * x := le_trans hj1 hxa
  have hhi : (p : ℝ) * x ≤ (j : ℝ) + 13 / 14 := le_of_lt (lt_of_lt_of_le hxb hj2)
  have hfloor : ⌊(p : ℝ) * x⌋ = j := by
    apply Int.floor_eq_iff.mpr
    constructor <;> [linarith; linarith]
  have hfract : Int.fract ((p : ℝ) * x) = (p : ℝ) * x - (j : ℝ) := by
    rw [Int.fract, hfloor]
  rw [hfract]
  constructor <;> linarith

/-- **The anchored positivity certificate**: an explicit interval inside `goodSet s.2` whose
band bounds hold for `s.1` forces `0 < witnessG2 s` — the exact shape of the repaired
assembly's `hk12` leg (and the positivity core of the other floor legs). -/
theorem witnessG2_pos_of_anchor (s : Shape) {a b : ℝ}
    (hgood : Set.Ico a b ⊆ TournamentH7.GoodSet.goodSet s.2)
    (hposP : ∀ p ∈ s.1, (0 : ℤ) < p)
    (hband : ∀ p ∈ s.1, ∃ j : ℤ, (j : ℝ) + 1 / 14 ≤ (p : ℝ) * a ∧ (p : ℝ) * b ≤ (j : ℝ) + 13 / 14)
    (h0 : 0 ≤ a) (h1 : b ≤ 1) (hab : a < b) :
    0 < witnessG2 s := by
  have hsafe : Set.Ico a b ⊆ safeSet s.1 := Ico_subset_safeSet_of_bounds hposP hband
  have hsub : Set.Ico a b ⊆ TournamentH7.GoodSet.goodSet s.2 ∩ safeSet s.1 :=
    Set.subset_inter hgood hsafe
  exact slowμ_toReal_pos_of_Ico_subset hsub h0 h1 hab

end IntervalBridge
end LRC14
end LonelyRunner
