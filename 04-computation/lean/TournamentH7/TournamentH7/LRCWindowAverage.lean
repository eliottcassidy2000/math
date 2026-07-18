/- LRCWindowAverage.lean -- opus-2026-07-17-S347 (HYP-7300 / THM-964 (M)).
   THE FUBINI POSITION STEP.  Two pieces:

   (1) window_start_measure  [PROVED, kernel-pure]: the inner integral -- the
       measure of window starts x whose window [x, x+L] contains a fixed y is
       vol (Icc (y-L) y) = L.  This is the mathematical core of the Fubini
       swap (the x-marginal of the window indicator).

   (2) window_average [PROVED]: Tonelli over the indicator product gives
       ∫⁻ x, vol (A ∩ Icc x (x+L)) = ofReal L·volA.
   (3) live_window_exists [PROVED]: vol A > 0 and L > 0 ⟹ some start x has
       positive in-window mass -- the wall's window-choice, kernel-pure. -/
import Mathlib

open MeasureTheory
open scoped ENNReal

namespace LonelyRunner.LRC14.Hunter

/-- **(1) THE INNER INTEGRAL** (kernel-pure): the starts `x` whose window
`[x, x+L]` contains a fixed `y` form `[y-L, y]`, of measure `L`. -/
theorem window_start_measure (y L : ℝ) (hL : 0 ≤ L) :
    ∫⁻ x, (Set.Icc x (x + L)).indicator (1 : ℝ → ℝ≥0∞) y = ENNReal.ofReal L := by
  have hset : (fun x => (Set.Icc x (x + L)).indicator (1 : ℝ → ℝ≥0∞) y)
      = (Set.Icc (y - L) y).indicator (1 : ℝ → ℝ≥0∞) := by
    funext x
    by_cases h : x ≤ y ∧ y ≤ x + L
    · rw [Set.indicator_of_mem (Set.mem_Icc.mpr h),
        Set.indicator_of_mem (Set.mem_Icc.mpr ⟨by linarith [h.2], h.1⟩)]
      rfl
    · have hnot : x ∉ Set.Icc (y - L) y := by
        intro hmem
        obtain ⟨h1, h2⟩ := Set.mem_Icc.mp hmem
        exact h ⟨h2, by linarith⟩
      rw [Set.indicator_of_notMem (fun hmem => h (Set.mem_Icc.mp hmem)),
        Set.indicator_of_notMem hnot]
  rw [hset, lintegral_indicator_one measurableSet_Icc, Real.volume_Icc]
  congr 1
  ring

/-- **(2) THE WINDOW-AVERAGE IDENTITY** (Tonelli): the average in-window mass
over all window starts is `ofReal L · vol A`.  Proof: the in-window mass is
the inner `y`-integral of `A.indicator 1 · (Icc x (x+L)).indicator 1`;
`lintegral_lintegral_swap` moves the `x`-integral inside, and its `x`-marginal
is `window_start_measure` (= `L`). -/
theorem window_average (A : Set ℝ) (hA : MeasurableSet A) (L : ℝ) (hL : 0 ≤ L) :
    ∫⁻ x, volume (A ∩ Set.Icc x (x + L)) = ENNReal.ofReal L * volume A := by
  -- the window slab {(x,y) : x ≤ y ≤ x+L} is measurable
  have hslab : MeasurableSet {p : ℝ × ℝ | p.1 ≤ p.2 ∧ p.2 ≤ p.1 + L} := by
    apply MeasurableSet.inter
    · exact measurableSet_le measurable_fst measurable_snd
    · exact measurableSet_le measurable_snd (measurable_fst.add measurable_const)
  -- the uncurried integrand
  have hmeas : Measurable (Function.uncurry
      fun x y => A.indicator (1 : ℝ → ℝ≥0∞) y
        * (Set.Icc x (x + L)).indicator (1 : ℝ → ℝ≥0∞) y) := by
    have h1 : Measurable (fun p : ℝ × ℝ => A.indicator (1 : ℝ → ℝ≥0∞) p.2) :=
      (measurable_one.indicator hA).comp measurable_snd
    have h2 : Measurable (fun p : ℝ × ℝ =>
        ({p : ℝ × ℝ | p.1 ≤ p.2 ∧ p.2 ≤ p.1 + L}).indicator
          (1 : ℝ × ℝ → ℝ≥0∞) p) := measurable_one.indicator hslab
    have hcongr : (Function.uncurry
        fun x y => A.indicator (1 : ℝ → ℝ≥0∞) y
          * (Set.Icc x (x + L)).indicator (1 : ℝ → ℝ≥0∞) y)
        = fun p : ℝ × ℝ => A.indicator (1 : ℝ → ℝ≥0∞) p.2
            * ({p : ℝ × ℝ | p.1 ≤ p.2 ∧ p.2 ≤ p.1 + L}).indicator
                (1 : ℝ × ℝ → ℝ≥0∞) p := by
      funext p
      simp only [Function.uncurry, Set.indicator_apply, Set.mem_Icc, Set.mem_setOf_eq,
        Pi.one_apply]
    rw [hcongr]
    exact h1.mul h2
  -- pointwise: in-window mass = ∫⁻ y, product of indicators
  have hpt : ∀ x, volume (A ∩ Set.Icc x (x + L))
      = ∫⁻ y, A.indicator (1 : ℝ → ℝ≥0∞) y
          * (Set.Icc x (x + L)).indicator (1 : ℝ → ℝ≥0∞) y := by
    intro x
    rw [← lintegral_indicator_one (hA.inter measurableSet_Icc)]
    refine lintegral_congr fun y => ?_
    rw [Set.inter_indicator_one, Pi.mul_apply]
  simp_rw [hpt]
  rw [lintegral_lintegral_swap hmeas.aemeasurable]
  have hinner : ∀ y, ∫⁻ x, A.indicator (1 : ℝ → ℝ≥0∞) y
        * (Set.Icc x (x + L)).indicator (1 : ℝ → ℝ≥0∞) y
      = A.indicator (1 : ℝ → ℝ≥0∞) y * ENNReal.ofReal L := by
    intro y
    rw [lintegral_const_mul' _ _ (by
      by_cases hy : y ∈ A
      · rw [Set.indicator_of_mem hy]; exact ENNReal.one_ne_top
      · rw [Set.indicator_of_notMem hy]; exact ENNReal.zero_ne_top),
      window_start_measure y L hL]
  simp_rw [hinner]
  rw [lintegral_mul_const' _ _ (by finiteness), lintegral_indicator_one hA, mul_comm]

/-- **THE LIVE-WINDOW COROLLARY**: if `A` has positive measure and `L > 0`,
some window start `x` gives positive in-window mass.  If every start were dead
the average `∫⁻ x, vol (A ∩ Icc x (x+L))` would be `0`, contradicting
`window_average` (= `ofReal L · vol A > 0`). -/
theorem live_window_exists (A : Set ℝ) (hA : MeasurableSet A) (L : ℝ)
    (hL : 0 < L) (hpos : 0 < volume A) :
    ∃ x : ℝ, 0 < volume (A ∩ Set.Icc x (x + L)) := by
  by_contra hcon
  push_neg at hcon
  have hzero : ∀ x, volume (A ∩ Set.Icc x (x + L)) = 0 :=
    fun x => nonpos_iff_eq_zero.mp (hcon x)
  have hint : ∫⁻ x, volume (A ∩ Set.Icc x (x + L)) = 0 := by
    simp_rw [hzero]; exact lintegral_zero
  rw [window_average A hA L (le_of_lt hL)] at hint
  exact ENNReal.mul_pos (ENNReal.ofReal_pos.mpr hL).ne' hpos.ne' |>.ne' hint

/-! ## Axiom audit (all three kernel-pure) -/
#print axioms window_start_measure
#print axioms window_average
#print axioms live_window_exists

end LonelyRunner.LRC14.Hunter
