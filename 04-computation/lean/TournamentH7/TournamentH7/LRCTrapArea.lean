/-
  TournamentH7.LRCTrapArea  (mac-mini-2026-07-02-S19, HYP-3876)

  THE TRAPEZOID AREA — the analytic heart of the drifting pair-floor (klein-S118 handed the
  aggregate / dense-events case to mac-mini).

  klein's `SpreadPairFloor.trap w₁ w₂ r` (LRCSpreadPairFloor Stage 1) is the EXACT overlap of two
  danger teeth at lattice residue `r`.  This file proves its integral over the support is EXACTLY
  `4·(1/14)² = 1/49`, INDEPENDENT of the speeds `w₁, w₂`:

      ∫_{-S}^{S} trap w₁ w₂ r dr = 1/49 ,   S = (1/14)(w₁+w₂).

  NOTE (mac-mini-S19): SELF-CONTAINED — re-declares `trap` verbatim from klein's
  `LonelyRunner.SpreadPairFloor` (definitionally equal), because that file does not currently
  compile against the pinned mathlib v4.30.0 (`div_le_div_of_nonneg_right`, `Int.add_mul_emod_self`
  name/signature mismatches — flagged to klein/kps).  The area transfers verbatim once it builds.
-/
import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.Tactic

namespace TournamentH7.TrapArea

/-- The trapezoid overlap of two teeth at lattice residue `r` (verbatim from
`LonelyRunner.SpreadPairFloor.trap`; `h = 1/14`).  Plateau `2h/w₂`, support `|r| < h(w₁+w₂)`. -/
noncomputable def trap (w₁ w₂ r : ℝ) : ℝ :=
  max 0 (min (2 * (1/14) / w₂) (((1/14) * (w₁ + w₂) - |r|) / (w₁ * w₂)))

theorem trap_nonneg (w₁ w₂ r : ℝ) : 0 ≤ trap w₁ w₂ r := le_max_left _ _

/-! ## The piecewise characterization of the trapezoid -/

/-- On the plateau `|r| ≤ (1/14)(w₂−w₁)` the trapezoid equals the narrower tooth width `2h/w₂`. -/
theorem trap_eq_plateau {w₁ w₂ : ℝ} (hw₁ : 0 < w₁) (hw₁₂ : w₁ ≤ w₂) {r : ℝ}
    (hr : |r| ≤ (1/14) * (w₂ - w₁)) :
    trap w₁ w₂ r = 2 * (1/14) / w₂ := by
  have hw₂ : 0 < w₂ := lt_of_lt_of_le hw₁ hw₁₂
  have hP : 0 < w₁ * w₂ := mul_pos hw₁ hw₂
  unfold trap
  have hslope_ge : 2 * (1/14) / w₂ ≤ ((1/14) * (w₁ + w₂) - |r|) / (w₁ * w₂) := by
    rw [le_div_iff₀ hP]
    have h1 : 2 * (1/14) / w₂ * (w₁ * w₂) = 2 * (1/14) * w₁ := by
      field_simp
    rw [h1]; nlinarith [hr, hw₁, hw₂]
  rw [min_eq_left hslope_ge, max_eq_right (div_nonneg (by norm_num) hw₂.le)]

/-- On the right slope `(1/14)(w₂−w₁) ≤ r ≤ (1/14)(w₁+w₂)` the trapezoid is the tent `(S−r)/(w₁w₂)`. -/
theorem trap_eq_slope_pos {w₁ w₂ : ℝ} (hw₁ : 0 < w₁) (hw₁₂ : w₁ ≤ w₂) {r : ℝ}
    (hr0 : (1/14) * (w₂ - w₁) ≤ r) (hrS : r ≤ (1/14) * (w₁ + w₂)) :
    trap w₁ w₂ r = ((1/14) * (w₁ + w₂) - r) / (w₁ * w₂) := by
  have hw₂ : 0 < w₂ := lt_of_lt_of_le hw₁ hw₁₂
  have hP : 0 < w₁ * w₂ := mul_pos hw₁ hw₂
  have hrpos : 0 ≤ r := le_trans (by positivity) hr0
  have habs : |r| = r := abs_of_nonneg hrpos
  unfold trap
  rw [habs]
  have hslope_le : ((1/14) * (w₁ + w₂) - r) / (w₁ * w₂) ≤ 2 * (1/14) / w₂ := by
    rw [div_le_iff₀ hP]
    have h1 : 2 * (1/14) / w₂ * (w₁ * w₂) = 2 * (1/14) * w₁ := by
      field_simp
    rw [h1]; nlinarith [hr0, hw₁, hw₂]
  have hslope_nonneg : 0 ≤ ((1/14) * (w₁ + w₂) - r) / (w₁ * w₂) :=
    div_nonneg (by linarith) hP.le
  rw [min_eq_right hslope_le, max_eq_right hslope_nonneg]

/-- On the left slope `−(1/14)(w₁+w₂) ≤ r ≤ −(1/14)(w₂−w₁)` the trapezoid is `(S+r)/(w₁w₂)`. -/
theorem trap_eq_slope_neg {w₁ w₂ : ℝ} (hw₁ : 0 < w₁) (hw₁₂ : w₁ ≤ w₂) {r : ℝ}
    (hrS : -((1/14) * (w₁ + w₂)) ≤ r) (hr0 : r ≤ -((1/14) * (w₂ - w₁))) :
    trap w₁ w₂ r = ((1/14) * (w₁ + w₂) + r) / (w₁ * w₂) := by
  have hw₂ : 0 < w₂ := lt_of_lt_of_le hw₁ hw₁₂
  have hP : 0 < w₁ * w₂ := mul_pos hw₁ hw₂
  have hrneg : r ≤ 0 := by nlinarith [hr0, hw₁₂]
  have habs : |r| = -r := abs_of_nonpos hrneg
  unfold trap
  rw [habs]
  have hslope_le : ((1/14) * (w₁ + w₂) - -r) / (w₁ * w₂) ≤ 2 * (1/14) / w₂ := by
    rw [div_le_iff₀ hP]
    have h1 : 2 * (1/14) / w₂ * (w₁ * w₂) = 2 * (1/14) * w₁ := by
      field_simp
    rw [h1]; nlinarith [hr0, hw₁, hw₂]
  have hslope_nonneg : 0 ≤ ((1/14) * (w₁ + w₂) - -r) / (w₁ * w₂) :=
    div_nonneg (by linarith) hP.le
  rw [min_eq_right hslope_le, max_eq_right hslope_nonneg]
  ring_nf

/-! ## Continuity and the trapezoid area = 1/49 -/

/-- The trapezoid is continuous. -/
theorem continuous_trap (w₁ w₂ : ℝ) : Continuous (trap w₁ w₂) := by
  unfold trap
  refine continuous_const.max (continuous_const.min ?_)
  exact (continuous_const.sub continuous_abs).div_const _

/-- `d/dx (x²/2) = x`. -/
private theorem hasDerivAt_sq_half (r : ℝ) : HasDerivAt (fun x : ℝ => x ^ 2 / 2) r r := by
  simpa using (hasDerivAt_pow 2 r).div_const 2

/-- A linear-piece integral: `∫_{a}^{b} (c − r) dr = c(b−a) − (b²−a²)/2`. -/
private theorem integral_sub_id (c a b : ℝ) :
    (∫ r in a..b, (c - r)) = c * (b - a) - (b ^ 2 - a ^ 2) / 2 := by
  have key : ∀ r ∈ Set.uIcc a b, HasDerivAt (fun x => c * x - x ^ 2 / 2) (c - r) r := by
    intro r _
    have h1 : HasDerivAt (fun x : ℝ => c * x) c r := by simpa using (hasDerivAt_id r).const_mul c
    exact h1.sub (hasDerivAt_sq_half r)
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt key
        ((continuous_const.sub continuous_id).intervalIntegrable _ _)]
  ring

/-- A linear-piece integral: `∫_{a}^{b} (c + r) dr = c(b−a) + (b²−a²)/2`. -/
private theorem integral_add_id (c a b : ℝ) :
    (∫ r in a..b, (c + r)) = c * (b - a) + (b ^ 2 - a ^ 2) / 2 := by
  have key : ∀ r ∈ Set.uIcc a b, HasDerivAt (fun x => c * x + x ^ 2 / 2) (c + r) r := by
    intro r _
    have h1 : HasDerivAt (fun x : ℝ => c * x) c r := by simpa using (hasDerivAt_id r).const_mul c
    exact h1.add (hasDerivAt_sq_half r)
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt key
        ((continuous_const.add continuous_id).intervalIntegrable _ _)]
  ring

/-- **THE TRAPEZOID AREA = 1/49, INDEPENDENT OF THE SPEEDS.**  The exact overlap of two danger
teeth, integrated over its support, is `4·(1/14)² = 1/49` for every pair `0 < w₁ ≤ w₂`. -/
theorem trap_integral {w₁ w₂ : ℝ} (hw₁ : 0 < w₁) (hw₁₂ : w₁ ≤ w₂) :
    ∫ r in (-((1/14) * (w₁ + w₂)))..((1/14) * (w₁ + w₂)), trap w₁ w₂ r = 1/49 := by
  have hw₂ : 0 < w₂ := lt_of_lt_of_le hw₁ hw₁₂
  have hP : 0 < w₁ * w₂ := mul_pos hw₁ hw₂
  have hPne : (w₁ * w₂) ≠ 0 := ne_of_gt hP
  have hw₂ne : w₂ ≠ 0 := ne_of_gt hw₂
  have hcont := continuous_trap w₁ w₂
  have hint : ∀ a b : ℝ, IntervalIntegrable (trap w₁ w₂) MeasureTheory.volume a b :=
    fun a b => hcont.intervalIntegrable a b
  set S : ℝ := (1/14) * (w₁ + w₂) with hSdef
  set S₀ : ℝ := (1/14) * (w₂ - w₁) with hS0def
  have hnegS0S0 : -S₀ ≤ S₀ := by rw [hS0def]; nlinarith
  have hS0S : S₀ ≤ S := by rw [hS0def, hSdef]; nlinarith
  have hnegSS0 : -S ≤ -S₀ := by linarith
  have hsplit : (∫ r in (-S)..S, trap w₁ w₂ r)
      = (∫ r in (-S)..(-S₀), trap w₁ w₂ r) + (∫ r in (-S₀)..S₀, trap w₁ w₂ r)
        + (∫ r in S₀..S, trap w₁ w₂ r) := by
    rw [intervalIntegral.integral_add_adjacent_intervals (hint (-S) (-S₀)) (hint (-S₀) S₀),
        intervalIntegral.integral_add_adjacent_intervals (hint (-S) S₀) (hint S₀ S)]
  -- plateau
  have heqonP : Set.EqOn (trap w₁ w₂) (fun _ => 2 * (1/14) / w₂) (Set.uIcc (-S₀) S₀) := by
    intro r hr
    rw [Set.uIcc_of_le hnegS0S0, Set.mem_Icc] at hr
    exact trap_eq_plateau hw₁ hw₁₂ (abs_le.mpr ⟨hr.1, hr.2⟩)
  have hplat : (∫ r in (-S₀)..S₀, trap w₁ w₂ r) = (2 * (1/14) / w₂) * (2 * S₀) := by
    rw [intervalIntegral.integral_congr heqonP, intervalIntegral.integral_const, smul_eq_mul]; ring
  -- right slope
  have heqonR : Set.EqOn (trap w₁ w₂) (fun r => (S - r) / (w₁ * w₂)) (Set.uIcc S₀ S) := by
    intro r hr
    rw [Set.uIcc_of_le hS0S, Set.mem_Icc] at hr
    exact trap_eq_slope_pos hw₁ hw₁₂ hr.1 hr.2
  have hslopeR : (∫ r in S₀..S, trap w₁ w₂ r)
      = (w₁ * w₂)⁻¹ * (S * (S - S₀) - (S ^ 2 - S₀ ^ 2) / 2) := by
    rw [intervalIntegral.integral_congr heqonR]
    have hdiv : (fun r => (S - r) / (w₁ * w₂)) = (fun r => (w₁ * w₂)⁻¹ * (S - r)) := by
      funext r; rw [div_eq_mul_inv, mul_comm]
    rw [hdiv, intervalIntegral.integral_const_mul, integral_sub_id]
  -- left slope
  have heqonL : Set.EqOn (trap w₁ w₂) (fun r => (S + r) / (w₁ * w₂)) (Set.uIcc (-S) (-S₀)) := by
    intro r hr
    rw [Set.uIcc_of_le hnegSS0, Set.mem_Icc] at hr
    exact trap_eq_slope_neg hw₁ hw₁₂ (by linarith [hr.1]) (by linarith [hr.2])
  have hslopeL : (∫ r in (-S)..(-S₀), trap w₁ w₂ r)
      = (w₁ * w₂)⁻¹ * (S * ((-S₀) - (-S)) + ((-S₀) ^ 2 - (-S) ^ 2) / 2) := by
    rw [intervalIntegral.integral_congr heqonL]
    have hdiv : (fun r => (S + r) / (w₁ * w₂)) = (fun r => (w₁ * w₂)⁻¹ * (S + r)) := by
      funext r; rw [div_eq_mul_inv, mul_comm]
    rw [hdiv, intervalIntegral.integral_const_mul, integral_add_id]
  rw [hsplit, hplat, hslopeR, hslopeL, hSdef, hS0def]
  field_simp
  ring

end TournamentH7.TrapArea
