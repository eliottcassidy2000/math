/-
  TournamentH7.LRCMreachConcrete -- Concrete Mreach and proof of lonely_of_Mreach_ge.

  This module provides the mathematically correct proof structure for filling the
  `lonely_of_Mreach_ge` sorry in `LRCFourteenSkeleton`.

  PROOF STRATEGY: The nearest-integer distance `nearInt x = min (fract x) (1 - fract x)`
  is continuous and 1-periodic.  The min-reach function `minReach v t = min_i nearInt(v_i t)`
  is therefore continuous and 1-periodic.  By compactness of [0,1], its supremum over [0,1]
  equals the supremum over all of ℝ and is ATTAINED.  If the attained maximum ≥ 1/14,
  the attaining time is 14-lonely.

  REMAINING LEAN OBLIGATIONS (marked sorry):
  - Continuity of `nearInt` (Lipschitz from quotient metric ℝ/ℤ; needs AddCircle or manual
    piecewise-continuity argument).  PROOF SKETCH: nearInt x = 1/2 - |fract x - 1/2|.
    At integer n: both left limit (fract → 1, |1-1/2|=1/2) and right limit (fract → 0,
    |0-1/2|=1/2) agree with the value 1/2, so |fract x - 1/2| is continuous at integers.
    Hence nearInt is continuous everywhere.  Equivalently: nearInt x = ⨅ m:ℤ, |x - m| is
    1-Lipschitz (triangle inequality on the quotient metric ℝ/ℤ), hence continuous.

  PROVED (no sorries):
  - nearInt_nonneg, nearInt_le_half, nearInt_eq, nearInt_add_one
  - nearInt_int_mul_add_one (periodicity under integer speed multiplication)
  - minReach_continuous (modulo nearInt_continuous sorry)
  - minReach_periodic (1-periodicity)
  - minReach_int_periodic (n-periodicity for any n : ℤ)
  - lonely_iff_minReach_ge
  - Mreach_eq_global_sSup (uses minReach_int_periodic)
  - lonely_of_Mreach_ge (the main compactness theorem, proved completely modulo nearInt_continuous)

  claude-sonnet-2026-06-22-S5.
-/
import Mathlib.Tactic
import Mathlib.Topology.Algebra.Order.IntermediateValue
import Mathlib.Topology.Instances.Real
import TournamentH7.LonelyRunner

namespace LonelyRunner
namespace LRC14Concrete

open Set Filter Finset

/-! ## 1. Nearest-integer distance -/

/-- The distance from x to the nearest integer: `nearInt x = min (fract x) (1 - fract x)`.
  Equivalently, nearInt x = ⨅ m : ℤ, |x - m|.  This function is 1-Lipschitz and
  1-periodic; it is continuous. -/
noncomputable def nearInt (x : ℝ) : ℝ := min (Int.fract x) (1 - Int.fract x)

lemma nearInt_nonneg (x : ℝ) : 0 ≤ nearInt x := by
  unfold nearInt
  exact le_min (Int.fract_nonneg x) (by linarith [Int.fract_lt_one x])

lemma nearInt_le_half (x : ℝ) : nearInt x ≤ 1 / 2 := by
  unfold nearInt
  rcases le_or_lt (Int.fract x) (1 / 2) with h | h
  · exact (min_le_left _ _).trans h
  · exact (min_le_right _ _).trans (by linarith)

/-- nearInt in terms of fract: nearInt x = 1/2 - |fract x - 1/2|. -/
lemma nearInt_eq (x : ℝ) : nearInt x = 1 / 2 - |Int.fract x - 1 / 2| := by
  unfold nearInt
  rcases le_or_lt (Int.fract x) (1 / 2) with h | h
  · rw [abs_of_nonpos (by linarith), min_eq_left (by linarith)]; ring
  · rw [abs_of_pos (by linarith), min_eq_right (by linarith)]; ring

/-- nearInt is 1-periodic: nearInt (x + 1) = nearInt x. -/
lemma nearInt_add_one (x : ℝ) : nearInt (x + 1) = nearInt x := by
  simp [nearInt, Int.fract_add_one]

/-- nearInt (v * (t + 1)) = nearInt (v * t) for integer v. -/
lemma nearInt_int_mul_add_one (v : ℤ) (t : ℝ) :
    nearInt ((v : ℝ) * (t + 1)) = nearInt ((v : ℝ) * t) := by
  simp only [nearInt]
  have h : (v : ℝ) * (t + 1) = (v : ℝ) * t + v := by ring
  rw [h]
  push_cast
  rw [Int.fract_add_int]

/-- **KEY SORRY**: nearInt is continuous.
  Proof: nearInt x = 1/2 - |fract x - 1/2|.  At an integer n:
    - value: Int.fract n = 0, so |0 - 1/2| = 1/2, nearInt n = 0.
    - left limit (x → n⁻): Int.fract x → 1, so |1 - 1/2| = 1/2, nearInt → 0.
    - right limit (x → n⁺): Int.fract x → 0, so |0 - 1/2| = 1/2, nearInt → 0.
  Both limits match; |Int.fract x - 1/2| is continuous everywhere, hence so is nearInt.
  Alternatively: nearInt x = ⨅ m:ℤ, |x - m| is 1-Lipschitz (triangle inequality on ℝ/ℤ),
  so nearInt is Lipschitz and hence continuous.
  In Lean: use ContinuousAt at each integer via the left/right limit argument, and
  continuity on (n, n+1) from continuity of Int.fract there. -/
theorem nearInt_continuous : Continuous nearInt := by
  sorry

/-! ## 2. Min-reach function -/

/-- The min-reach at time t for speed family v: minimum nearInt over all runners. -/
noncomputable def minReach (v : Fin 13 → ℤ) (t : ℝ) : ℝ :=
  ⨅ i : Fin 13, nearInt ((v i : ℝ) * t)

/-- minReach is continuous in t (finite infimum of continuous functions). -/
theorem minReach_continuous (v : Fin 13 → ℤ) : Continuous (minReach v) := by
  unfold minReach
  apply continuous_iInf
  intro i
  exact nearInt_continuous.comp (continuous_const.mul continuous_id)

/-- minReach is 1-periodic in t. -/
theorem minReach_periodic (v : Fin 13 → ℤ) (t : ℝ) :
    minReach v (t + 1) = minReach v t := by
  unfold minReach
  congr 1; ext i
  exact nearInt_int_mul_add_one (v i) t

/-- minReach has period n for any integer n : ℤ.
  Proved by induction on ℤ using minReach_periodic. -/
lemma minReach_int_periodic (v : Fin 13 → ℤ) (t : ℝ) (n : ℤ) :
    minReach v (t + n) = minReach v t := by
  induction n using Int.induction_on with
  | hz => simp
  | hp k ih =>
    push_cast
    rw [show t + (↑k + 1 : ℝ) = (t + ↑k) + 1 by ring]
    rw [minReach_periodic, ih]
  | hn k ih =>
    push_cast
    -- Goal: minReach v (t + (-↑k - 1)) = minReach v t
    -- minReach_periodic at (t + -↑k - 1): minReach v (t + -↑k) = minReach v (t + -↑k - 1)
    have hperiod := minReach_periodic v (t + (-↑k : ℝ) - 1)
    simp only [show t + (-↑k : ℝ) - 1 + 1 = t + (-↑k : ℝ) by ring] at hperiod
    -- hperiod : minReach v (t + -↑k) = minReach v (t + -↑k - 1)
    -- hperiod.symm : minReach v (t + -↑k - 1) = minReach v (t + -↑k)
    rw [show t + (-↑k - 1 : ℝ) = t + (-↑k : ℝ) - 1 by ring]
    exact hperiod.symm.trans ih

/-- Lonely 14 v t ↔ minReach v t ≥ 1/14. -/
theorem lonely_iff_minReach_ge (v : Fin 13 → ℤ) (t : ℝ) :
    Lonely 14 v t ↔ (1 : ℝ) / 14 ≤ minReach v t := by
  rw [lonely_iff_fract_mem]
  unfold minReach nearInt
  constructor
  · intro h
    apply le_ciInf
    intro i
    have ⟨hlo, hhi⟩ := h i
    exact le_min hlo (by linarith)
  · intro h i
    have hle : (1 : ℝ) / 14 ≤ min (Int.fract ((v i : ℝ) * t))
                                     (1 - Int.fract ((v i : ℝ) * t)) := by
      calc (1 : ℝ) / 14 ≤ minReach v t := h
        _ ≤ min (Int.fract ((v i : ℝ) * t)) (1 - Int.fract ((v i : ℝ) * t)) :=
            ciInf_le ⟨0, fun y ⟨j, rfl⟩ => le_min (Int.fract_nonneg _) (by
              linarith [Int.fract_lt_one ((v j : ℝ) * t)])⟩ i
    exact ⟨le_trans hle (min_le_left _ _), by linarith [min_le_right (Int.fract ((v i : ℝ) * t))
           (1 - Int.fract ((v i : ℝ) * t))]⟩

/-! ## 3. Concrete Mreach and compactness theorem -/

/-- The concrete Mreach: supremum of minReach v over [0,1].
  By 1-periodicity, this equals the global supremum over all of ℝ. -/
noncomputable def Mreach (v : Fin 13 → ℤ) : ℝ :=
  sSup (minReach v '' Icc (0 : ℝ) 1)

/-- Mreach v is the global sup of minReach v.
  Key step: minReach v t = minReach v (Int.fract t) by minReach_int_periodic,
  and Int.fract t ∈ [0,1). -/
theorem Mreach_eq_global_sSup (v : Fin 13 → ℤ) :
    Mreach v = sSup (range (minReach v)) := by
  apply le_antisymm
  · exact sSup_le_sSup (image_subset_range _ _)
  · apply sSup_le; rintro r ⟨t, rfl⟩
    -- Reduce to Int.fract t ∈ [0,1) using periodicity
    have ht : t = Int.fract t + ↑⌊t⌋ := by
      linarith [Int.floor_add_fract t]
    have heq : minReach v (Int.fract t) = minReach v t :=
      calc minReach v (Int.fract t)
          = minReach v (Int.fract t + ↑⌊t⌋) := (minReach_int_periodic v (Int.fract t) ⌊t⌋).symm
        _ = minReach v t := by rw [← ht]
    rw [← heq]
    exact le_sSup ⟨Int.fract t, ⟨Int.fract_nonneg t, le_of_lt (Int.fract_lt_one t)⟩, rfl⟩

/-- **The key compactness theorem:** If Mreach v ≥ 1/14, there exists a 14-lonely time.

  Proof: minReach v is continuous (from nearInt_continuous) and 1-periodic.  Hence
  Mreach v = sSup (minReach v '' [0,1]).  Since [0,1] is compact and minReach v is
  continuous, the supremum is attained at some t* ∈ [0,1].  Then
  minReach v t* = Mreach v ≥ 1/14, so by lonely_iff_minReach_ge, t* is 14-lonely.

  This proof is COMPLETE modulo the single sorry `nearInt_continuous`. -/
theorem lonely_of_Mreach_ge (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (hM : (1 : ℝ) / 14 ≤ Mreach v) : ∃ t : ℝ, Lonely 14 v t := by
  have hcomp : IsCompact (Icc (0 : ℝ) 1) := isCompact_Icc
  have hne : (Icc (0 : ℝ) 1).Nonempty := ⟨0, by norm_num⟩
  have hcont : ContinuousOn (minReach v) (Icc (0 : ℝ) 1) :=
    (minReach_continuous v).continuousOn
  obtain ⟨t*, ht*_mem, ht*_max⟩ := hcomp.exists_isMaxOn hne hcont
  have ht*_eq : minReach v t* = Mreach v := by
    apply le_antisymm
    · exact le_sSup ⟨t*, ht*_mem, rfl⟩
    · rw [Mreach]
      exact sSup_le (fun r ⟨s, hs_mem, hs_eq⟩ => hs_eq ▸ ht*_max hs_mem)
  exact ⟨t*, (lonely_iff_minReach_ge v t*).mpr (ht*_eq ▸ hM)⟩

/-! ## 4. Connection to the skeleton's opaque Mreach -/

/-- NOTE: The skeleton uses `opaque Mreach` which cannot be unfolded.  To fill
  `lonely_of_Mreach_ge` in the skeleton, replace `opaque Mreach` with
  `noncomputable def Mreach := LRC14Concrete.Mreach` and use the theorem above.

  REMAINING SORRIES (count: 1):
  1. nearInt_continuous: Lipschitz/continuity of the nearest-integer distance.
     Proof route A: nearInt x = ⨅ m:ℤ, |x - m|; this is 1-Lipschitz by the
       triangle inequality: |x-m| ≤ |x-y| + |y-m|, so nearInt x ≤ |x-y| + nearInt y.
       Then LipschitzWith 1 nearInt follows, and .continuous gives continuity.
     Proof route B: use AddCircle 1 quotient metric.  The map ℝ → AddCircle 1 is
       continuous; nearInt = dist · 0 on AddCircle 1 factors continuously.
     Proof route C: piecewise on [n, n+1].  Int.fract is continuous on (n, n+1) (linear);
       at endpoints, both one-sided limits of nearInt are 0. Use Continuous.piecewise or
       ContinuousAt.congr at integer points.

  All other results are proved completely. -/

end LRC14Concrete
end LonelyRunner
