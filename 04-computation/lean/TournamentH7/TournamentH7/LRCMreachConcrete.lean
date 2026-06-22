/-
  TournamentH7.LRCMreachConcrete -- Concrete Mreach and proof of lonely_of_Mreach_ge.

  This module provides the mathematically correct proof structure for filling the
  `lonely_of_Mreach_ge` sorry in `LRCFourteenSkeleton`.

  PROOF STRATEGY: The nearest-integer distance `nearInt x = min (fract x) (1 - fract x)`
  is continuous and 1-periodic.  The min-reach function `minReach v t = min_i nearInt(v_i t)`
  is therefore continuous and 1-periodic.  By compactness of [0,1], its supremum over [0,1]
  equals the supremum over all of ℝ and is ATTAINED.  If the attained maximum ≥ 1/14,
  the attaining time is 14-lonely.

  REMAINING LEAN OBLIGATIONS (marked sorry in this buildable skeleton):
  - Continuity of `nearInt` (Lipschitz from quotient metric ℝ/ℤ; needs AddCircle or manual
    piecewise-continuity argument).  PROOF SKETCH: nearInt x = 1/2 - |fract x - 1/2|.
    At integer n: both left limit (fract → 1, |1-1/2|=1/2) and right limit (fract → 0,
    |0-1/2|=1/2) agree with the value 1/2, so |fract x - 1/2| is continuous at integers.
    Hence nearInt is continuous everywhere.  Equivalently: nearInt x = ⨅ m:ℤ, |x - m| is
    1-Lipschitz (triangle inequality on the quotient metric ℝ/ℤ), hence continuous.
  - Port the finite-infimum API for `minReach_continuous` and
    `lonely_iff_minReach_ge`.
  - Port integer-periodicity, global-sup reduction to `[0,1]`, and the final
    extreme-value theorem step.

  PROVED (no sorries):
  - nearInt_nonneg, nearInt_le_half, nearInt_eq, nearInt_add_one
  - nearInt_int_mul_add_one (periodicity under integer speed multiplication)
  - minReach_periodic (1-periodicity)

  First build repair: claude-sonnet-2026-06-22-S5, patched by codex-S82.
-/
import Mathlib.Tactic
import Mathlib.Topology.Order.IntermediateValue
import Mathlib.Topology.Instances.Real.Lemmas
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
  rcases le_total (Int.fract x) (1 / 2) with h | h
  · exact (min_le_left _ _).trans h
  · exact (min_le_right _ _).trans (by linarith)

/-- nearInt in terms of fract: nearInt x = 1/2 - |fract x - 1/2|. -/
lemma nearInt_eq (x : ℝ) : nearInt x = 1 / 2 - |Int.fract x - 1 / 2| := by
  unfold nearInt
  rcases le_total (Int.fract x) (1 / 2) with h | h
  · rw [abs_of_nonpos (by linarith), min_eq_left (by linarith)]; ring
  · rw [abs_of_nonneg (by linarith), min_eq_right (by linarith)]; ring

/-- nearInt is 1-periodic: nearInt (x + 1) = nearInt x. -/
lemma nearInt_add_one (x : ℝ) : nearInt (x + 1) = nearInt x := by
  simp [nearInt, Int.fract_add_one]

/-- nearInt (v * (t + 1)) = nearInt (v * t) for integer v. -/
lemma nearInt_int_mul_add_one (v : ℤ) (t : ℝ) :
    nearInt ((v : ℝ) * (t + 1)) = nearInt ((v : ℝ) * t) := by
  simp only [nearInt]
  have h : (v : ℝ) * (t + 1) = (v : ℝ) * t + v := by ring
  rw [h]
  rw [Int.fract_add_intCast]

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
  -- Finite infimum of continuous functions.  The previous draft used an old
  -- `continuous_iInf` name; fill this with the current Mathlib finite-inf API
  -- after `nearInt_continuous` is proved.
  sorry

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
  -- Integer-periodicity follows from `minReach_periodic` by induction over
  -- `Int`; the current Mathlib induction eliminator uses `zero/succ/pred`
  -- branch names, so this is left as a small porting obligation.
  sorry

/-- Lonely 14 v t ↔ minReach v t ≥ 1/14. -/
theorem lonely_iff_minReach_ge (v : Fin 13 → ℤ) (t : ℝ) :
    Lonely 14 v t ↔ (1 : ℝ) / 14 ≤ minReach v t := by
  -- This is a finite-infimum unpacking of `Lonely 14` through `nearInt`.
  -- The old proof used a brittle `ciInf_le` witness; port after choosing the
  -- preferred finite-min API for `Fin 13`.
  sorry

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
  -- Reduce any real `t` to `Int.fract t ∈ [0,1]` using
  -- `minReach_int_periodic`, then compare the two suprema.
  sorry

/-- **The key compactness theorem:** If Mreach v ≥ 1/14, there exists a 14-lonely time.

  Proof: minReach v is continuous (from nearInt_continuous) and 1-periodic.  Hence
  Mreach v = sSup (minReach v '' [0,1]).  Since [0,1] is compact and minReach v is
  continuous, the supremum is attained at some t* ∈ [0,1].  Then
  minReach v t* = Mreach v ≥ 1/14, so by lonely_iff_minReach_ge, t* is 14-lonely.

  This is the compactness theorem meant to replace the skeleton's opaque
  `lonely_of_Mreach_ge` once the preceding topology/finite-infimum obligations
  in this file are closed. -/
theorem lonely_of_Mreach_ge (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (hM : (1 : ℝ) / 14 ≤ Mreach v) : ∃ t : ℝ, Lonely 14 v t := by
  -- Extreme-value theorem on `[0,1]`, followed by
  -- `lonely_iff_minReach_ge`.  This is the theorem meant to fill the
  -- skeleton's opaque `lonely_of_Mreach_ge` once the definitions are unified.
  sorry

/-! ## 4. Connection to the skeleton's opaque Mreach -/

/- NOTE: The skeleton uses `opaque Mreach` which cannot be unfolded.  To fill
  `lonely_of_Mreach_ge` in the skeleton, replace `opaque Mreach` with
  `noncomputable def Mreach := LRC14Concrete.Mreach` and use the theorem above.

  REMAINING SORRIES (current porting skeleton):
  1. nearInt_continuous: Lipschitz/continuity of the nearest-integer distance.
  2. minReach_continuous: finite infimum of continuous functions.
  3. minReach_int_periodic: integer-periodicity from period 1.
  4. lonely_iff_minReach_ge: finite infimum unpacking.
  5. Mreach_eq_global_sSup: reduce any time to its fractional part.
  6. lonely_of_Mreach_ge: extreme-value theorem on [0,1].

     Proof route A: nearInt x = ⨅ m:ℤ, |x - m|; this is 1-Lipschitz by the
       triangle inequality: |x-m| ≤ |x-y| + |y-m|, so nearInt x ≤ |x-y| + nearInt y.
       Then LipschitzWith 1 nearInt follows, and .continuous gives continuity.
     Proof route B: use AddCircle 1 quotient metric.  The map ℝ → AddCircle 1 is
       continuous; nearInt = dist · 0 on AddCircle 1 factors continuously.
     Proof route C: piecewise on [n, n+1].  Int.fract is continuous on (n, n+1) (linear);
       at endpoints, both one-sided limits of nearInt are 0. Use Continuous.piecewise or
       ContinuousAt.congr at integer points.

  The easy algebraic nearest-integer lemmas above are proved; the topology and
  finite-infimum bridge is now a compiling Lean skeleton rather than a broken file. -/

end LRC14Concrete
end LonelyRunner
