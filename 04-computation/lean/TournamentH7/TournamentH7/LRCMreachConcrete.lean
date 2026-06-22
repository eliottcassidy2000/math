/-
  TournamentH7.LRCMreachConcrete -- Concrete Mreach compactness bridge skeleton.

  This module provides the mathematically correct proof structure for filling the
  `lonely_of_Mreach_ge` handoff that was previously open in
  `LRCFourteenSkeleton`.

  PROOF STRATEGY: The nearest-integer distance `nearInt x = min (fract x) (1 - fract x)`
  is continuous and 1-periodic.  The min-reach function `minReach v t = min_i nearInt(v_i t)`
  is therefore continuous and 1-periodic.  By compactness of [0,1], its supremum over [0,1]
  equals the supremum over all of ℝ and is ATTAINED.  If the attained maximum ≥ 1/14,
  the attaining time is 14-lonely.

  PROVED (no sorries):
  - nearInt_nonneg, nearInt_le_half, nearInt_eq, nearInt_add_one
  - nearInt_continuous, via mathlib's endpoint-compatible `ContinuousOn.comp_fract`
  - nearInt_int_mul_add_one (periodicity under integer speed multiplication)
  - continuous_finset_inf'_real
  - minReach_continuous
  - minReach_periodic (1-periodicity)
  - minReach_int_periodic (n-periodicity for any n : ℤ)
  - minReach_le_half
  - lonely_iff_minReach_ge
  - Mreach_eq_global_sSup
  - lonely_of_Mreach_ge

  Incoming source: claude-sonnet-2026-06-22-S5.  Codex S81 repaired imports and
  made the non-typechecking proof scripts explicit theorem targets; codex-S82
  tightened the porting comments against the current Lean DAG; codex-S83 closed
  the finite-infimum/algebraic/compactness porting sorries and discharged
  nearest-integer continuity through `ContinuousOn.comp_fract`.
-/
import Mathlib.Tactic
import Mathlib.Topology.Algebra.Order.Floor
import Mathlib.Topology.Order.Compact
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
  unfold nearInt
  have h : (v : ℝ) * (t + 1) = (v : ℝ) * t + v := by ring
  rw [h]
  rw [Int.fract_add_intCast]

/-- `nearInt` is continuous.

The function `u ↦ min u (1-u)` is continuous on `[0,1]` and has matching
endpoint values, so it composes continuously with `Int.fract`. -/
theorem nearInt_continuous : Continuous nearInt := by
  have hcont : ContinuousOn (fun u : ℝ => min u (1 - u)) (Icc (0 : ℝ) 1) :=
    (continuous_id.min (continuous_const.sub continuous_id)).continuousOn
  have hend : (fun u : ℝ => min u (1 - u)) 0 = (fun u : ℝ => min u (1 - u)) 1 := by
    norm_num
  simpa [nearInt, Function.comp_def] using hcont.comp_fract'' hend

/-! ## 2. Min-reach function -/

/-- The min-reach at time t for speed family v: minimum nearInt over all runners. -/
noncomputable def minReach (v : Fin 13 → ℤ) (t : ℝ) : ℝ :=
  ⨅ i : Fin 13, nearInt ((v i : ℝ) * t)

/-- A nonempty finite infimum of continuous real-valued functions is continuous.
This avoids requiring an artificial `⊤` element for `ℝ`, so it is the right
finite-minimum API for `minReach`. -/
lemma continuous_finset_inf'_real {ι : Type*} {s : Finset ι} (hs : s.Nonempty)
    {f : ι → ℝ → ℝ} (hf : ∀ i ∈ s, Continuous (f i)) :
    Continuous (fun x => s.inf' hs (fun i => f i x)) := by
  induction hs using Finset.Nonempty.cons_induction with
  | singleton a =>
      simpa [Finset.inf'_singleton] using hf a (by simp)
  | cons a s ha hs ih =>
      have hcont : Continuous (fun x => f a x ⊓ s.inf' hs (fun i => f i x)) :=
        (hf a (by simp)).min (ih (by
          intro i hi
          exact hf i (by simp [hi])))
      convert hcont using 1
      ext x
      exact Finset.inf'_cons hs (fun i => f i x)

/-- minReach is continuous in t (finite infimum of continuous functions). -/
theorem minReach_continuous (v : Fin 13 → ℤ) : Continuous (minReach v) := by
  simpa [minReach, Finset.inf'_univ_eq_ciInf] using
    (continuous_finset_inf'_real
      (Finset.univ_nonempty : (Finset.univ : Finset (Fin 13)).Nonempty)
      (f := fun i t => nearInt ((v i : ℝ) * t)) (by
        intro i _hi
        exact nearInt_continuous.comp (continuous_const.mul continuous_id)))

/-- The concrete min-reach function is bounded above by `1/2`. -/
lemma minReach_le_half (v : Fin 13 → ℤ) (t : ℝ) : minReach v t ≤ 1 / 2 := by
  unfold minReach
  have hbdd : BddBelow (Set.range fun i : Fin 13 => nearInt ((v i : ℝ) * t)) := by
    refine ⟨0, ?_⟩
    rintro y ⟨i, rfl⟩
    exact nearInt_nonneg _
  exact (ciInf_le hbdd (0 : Fin 13)).trans (nearInt_le_half _)

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
  have hper : Function.Periodic (minReach v) (1 : ℝ) := fun x => minReach_periodic v x
  have hn : Function.Periodic (minReach v) ((n : ℝ) * (1 : ℝ)) := hper.int_mul n
  simpa using hn t

/-- Lonely 14 v t ↔ minReach v t ≥ 1/14. -/
theorem lonely_iff_minReach_ge (v : Fin 13 → ℤ) (t : ℝ) :
    Lonely 14 v t ↔ (1 : ℝ) / 14 ≤ minReach v t := by
  rw [lonely_iff_fract_mem]
  unfold minReach nearInt
  have hbdd : BddBelow
      (Set.range fun i : Fin 13 =>
        min (Int.fract ((v i : ℝ) * t)) (1 - Int.fract ((v i : ℝ) * t))) := by
    refine ⟨0, ?_⟩
    rintro y ⟨i, rfl⟩
    exact le_min (Int.fract_nonneg _)
      (by linarith [Int.fract_lt_one ((v i : ℝ) * t)])
  rw [le_ciInf_iff hbdd]
  constructor
  · intro h i
    exact le_min (h i).1 (by linarith [(h i).2])
  · intro h i
    have hle := h i
    refine ⟨?_, ?_⟩
    · exact le_trans hle (min_le_left _ _)
    · have hr : (1 : ℝ) / 14 ≤ 1 - Int.fract ((v i : ℝ) * t) :=
        le_trans hle (min_le_right _ _)
      linarith

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
  have hRangeBdd : BddAbove (range (minReach v)) := by
    refine ⟨1 / 2, ?_⟩
    rintro y ⟨t, rfl⟩
    exact minReach_le_half v t
  have hImageBdd : BddAbove (minReach v '' Icc (0 : ℝ) 1) := by
    refine ⟨1 / 2, ?_⟩
    rintro y ⟨t, _ht, rfl⟩
    exact minReach_le_half v t
  have hImageNonempty : (minReach v '' Icc (0 : ℝ) 1).Nonempty := by
    exact ⟨minReach v 0, ⟨0, by norm_num, rfl⟩⟩
  have hRangeNonempty : (range (minReach v)).Nonempty := by
    exact ⟨minReach v 0, ⟨0, rfl⟩⟩
  apply le_antisymm
  · unfold Mreach
    exact csSup_le_csSup hRangeBdd hImageNonempty (image_subset_range _ _)
  · apply csSup_le hRangeNonempty
    rintro r ⟨t, rfl⟩
    have ht : t = Int.fract t + ↑⌊t⌋ := by
      linarith [Int.floor_add_fract t]
    have heq : minReach v (Int.fract t) = minReach v t :=
      calc minReach v (Int.fract t)
          = minReach v (Int.fract t + ↑⌊t⌋) :=
              (minReach_int_periodic v (Int.fract t) ⌊t⌋).symm
        _ = minReach v t := by rw [← ht]
    rw [← heq]
    exact le_csSup hImageBdd
      ⟨Int.fract t, ⟨Int.fract_nonneg t, le_of_lt (Int.fract_lt_one t)⟩, rfl⟩

/-- **The key compactness theorem:** If Mreach v ≥ 1/14, there exists a 14-lonely time.

  Proof: minReach v is continuous and 1-periodic.  Hence
  Mreach v = sSup (minReach v '' [0,1]).  Since [0,1] is compact and minReach v is
  continuous, the supremum is attained at some t* ∈ [0,1].  Then
  minReach v t* = Mreach v ≥ 1/14, so by lonely_iff_minReach_ge, t* is 14-lonely.

  This is the compactness theorem meant to replace the skeleton's opaque
  `lonely_of_Mreach_ge` once the preceding topology/finite-infimum obligations
  in this file are closed. -/
theorem lonely_of_Mreach_ge (v : Fin 13 → ℤ) (_hv : ∀ i, v i ≠ 0)
    (hM : (1 : ℝ) / 14 ≤ Mreach v) : ∃ t : ℝ, Lonely 14 v t := by
  have hcomp : IsCompact (Icc (0 : ℝ) 1) := isCompact_Icc
  have hne : (Icc (0 : ℝ) 1).Nonempty := ⟨0, by norm_num⟩
  have hcont : ContinuousOn (minReach v) (Icc (0 : ℝ) 1) :=
    (minReach_continuous v).continuousOn
  obtain ⟨tstar, htstar_mem, htstar_max⟩ := hcomp.exists_isMaxOn hne hcont
  have hImageBdd : BddAbove (minReach v '' Icc (0 : ℝ) 1) := by
    refine ⟨1 / 2, ?_⟩
    rintro y ⟨t, _ht, rfl⟩
    exact minReach_le_half v t
  have hImageNonempty : (minReach v '' Icc (0 : ℝ) 1).Nonempty := by
    exact ⟨minReach v 0, ⟨0, by norm_num, rfl⟩⟩
  have htstar_eq : minReach v tstar = Mreach v := by
    apply le_antisymm
    · unfold Mreach
      exact le_csSup hImageBdd ⟨tstar, htstar_mem, rfl⟩
    · unfold Mreach
      exact csSup_le hImageNonempty (by
        rintro r ⟨s, hs_mem, rfl⟩
        exact htstar_max hs_mem)
  exact ⟨tstar, (lonely_iff_minReach_ge v tstar).mpr (htstar_eq ▸ hM)⟩

/-! ## 4. Connection to the skeleton's opaque Mreach -/

/- NOTE: The skeleton now defines `Mreach` as this concrete supremum and uses
  `lonely_of_Mreach_ge` directly.

  The easy algebraic nearest-integer lemmas, finite-infimum continuity bridge,
  integer-periodicity, global-sup comparison, and extreme-value handoff are now
  proved without sorries. -/

/-! ## 5. Axiom audit

The finite-infimum, periodicity, global-sup, and compactness glue below should
not inherit `sorryAx`. -/

#print axioms nearInt_continuous
#print axioms continuous_finset_inf'_real
#print axioms minReach_continuous
#print axioms minReach_int_periodic
#print axioms lonely_iff_minReach_ge
#print axioms Mreach_eq_global_sSup
#print axioms lonely_of_Mreach_ge

end LRC14Concrete
end LonelyRunner
