/-
  TournamentH7.LRCPairUpgrade  (boxeph-2026-07-09-S15, HYP-5853/5863 Lean layer 2)

  THE PAIR UPGRADE — the first μ-certificate of chain-Bonferroni in Lean:
  the joint danger set of a doubling pair `(y, 2y)` at the 1/14 band is
  contained in THREE explicit intervals of length 1/14, so its measure is
  `≤ 3/14` — strictly below the per-runner union bound `2/7 = 4/14`.
  This is the exact `+1/14`-per-doubling-pair budget upgrade (HYP-5853's
  `μ₂ = 11/14`, subset direction, which is all any Bonferroni consumer
  needs), feeding chain-Bonferroni and the chain-coarsened B5 program.

  Route (no fract case analysis): `nearInt x < c` gives an integer within
  `c` of `x` (contrapositive of `le_nearInt_of_forall_int`); for
  `y ∈ [0,1]` the witness for `2y` is forced into `{0,1,2}` by integer
  bounds, yielding the three intervals `[0,1/14) ∪ (13/28,15/28) ∪ (13/14,1]`.

  Kernel-pure target: no `sorry`, no `native_decide` (audited below).
-/
import Mathlib
import TournamentH7.LRCHembedScaleSep

namespace LonelyRunner
namespace LRC14Concrete

open Set

/-- From `nearInt x < c`, an integer within `c` of `x` (contrapositive of
`le_nearInt_of_forall_int`). -/
theorem exists_int_near_of_nearInt_lt {x c : ℝ} (h : nearInt x < c) :
    ∃ m : ℤ, |x - m| < c := by
  by_contra hcon
  push_neg at hcon
  exact absurd (le_nearInt_of_forall_int c x hcon) (not_le.mpr h)

/-- **The pair-danger inclusion.**  For `y ∈ [0,1]`, if the doubling pair
`(y, 2y)` violates the `1/14` band, then `y` lies in one of three explicit
intervals of length `1/14`. -/
theorem pair_danger_subset :
    {y : ℝ | y ∈ Icc (0:ℝ) 1 ∧ (nearInt y < 1/14 ∨ nearInt (2*y) < 1/14)}
      ⊆ Ico (0:ℝ) (1/14) ∪ Ioo (13/28 : ℝ) (15/28) ∪ Ioc (13/14 : ℝ) 1 := by
  rintro y ⟨⟨hy0, hy1⟩, hdanger⟩
  rcases hdanger with h1 | h2
  · -- nearInt y < 1/14 on [0,1]: y near 0 or near 1
    obtain ⟨m, hm⟩ := exists_int_near_of_nearInt_lt h1
    have hmlo : (-1 : ℤ) < m := by
      by_contra hc
      push_neg at hc
      have : (m : ℝ) ≤ -1 := by exact_mod_cast hc
      have := abs_lt.mp hm
      linarith
    have hmhi : m < 2 := by
      by_contra hc
      push_neg at hc
      have : (2 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hc
      have := abs_lt.mp hm
      linarith
    interval_cases m
    · push_cast at hm
      have h := abs_lt.mp hm
      exact Or.inl (Or.inl ⟨hy0, by linarith [h.2]⟩)
    · push_cast at hm
      have h := abs_lt.mp hm
      exact Or.inr ⟨by linarith [h.1], hy1⟩
  · -- nearInt (2y) < 1/14 on [0,1]: witness m ∈ {0,1,2}
    obtain ⟨m, hm⟩ := exists_int_near_of_nearInt_lt h2
    have hmlo : (-1 : ℤ) < m := by
      by_contra hc
      push_neg at hc
      have : (m : ℝ) ≤ -1 := by exact_mod_cast hc
      have := abs_lt.mp hm
      linarith
    have hmhi : m < 3 := by
      by_contra hc
      push_neg at hc
      have : (3 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hc
      have := abs_lt.mp hm
      linarith
    interval_cases m
    · push_cast at hm
      have h := abs_lt.mp hm
      exact Or.inl (Or.inl ⟨hy0, by linarith [h.2]⟩)
    · push_cast at hm
      have h := abs_lt.mp hm
      exact Or.inl (Or.inr ⟨by linarith [h.1], by linarith [h.2]⟩)
    · push_cast at hm
      have h := abs_lt.mp hm
      exact Or.inr ⟨by linarith [h.1], hy1⟩

/-- **The pair upgrade:** the joint danger measure of a doubling pair is
`≤ 3/14` — one full `1/14` below the per-runner union bound `4/14 = 2/7`.
The first exact μ-certificate of HYP-5853 in Lean. -/
theorem pair_danger_volume_le :
    MeasureTheory.volume
      {y : ℝ | y ∈ Icc (0:ℝ) 1 ∧ (nearInt y < 1/14 ∨ nearInt (2*y) < 1/14)}
      ≤ ENNReal.ofReal (3/14) := by
  refine le_trans (MeasureTheory.measure_mono pair_danger_subset) ?_
  refine le_trans (MeasureTheory.measure_union_le _ _) ?_
  have h1 : MeasureTheory.volume (Ico (0:ℝ) (1/14) ∪ Ioo (13/28 : ℝ) (15/28))
      ≤ ENNReal.ofReal (1/14) + ENNReal.ofReal (1/14) := by
    refine le_trans (MeasureTheory.measure_union_le _ _) ?_
    rw [Real.volume_Ico, Real.volume_Ioo]
    apply add_le_add <;> apply ENNReal.ofReal_le_ofReal <;> norm_num
  have h2 : MeasureTheory.volume (Ioc (13/14 : ℝ) 1) ≤ ENNReal.ofReal (1/14) := by
    rw [Real.volume_Ioc]
    apply ENNReal.ofReal_le_ofReal
    norm_num
  calc MeasureTheory.volume (Ico (0:ℝ) (1/14) ∪ Ioo (13/28 : ℝ) (15/28))
        + MeasureTheory.volume (Ioc (13/14 : ℝ) 1)
      ≤ (ENNReal.ofReal (1/14) + ENNReal.ofReal (1/14)) + ENNReal.ofReal (1/14) :=
        add_le_add h1 h2
    _ = ENNReal.ofReal (3/14) := by
        rw [← ENNReal.ofReal_add (by norm_num) (by norm_num),
          ← ENNReal.ofReal_add (by norm_num) (by norm_num)]
        norm_num

end LRC14Concrete
end LonelyRunner

-- kernel-purity audit
#print axioms LonelyRunner.LRC14Concrete.exists_int_near_of_nearInt_lt
#print axioms LonelyRunner.LRC14Concrete.pair_danger_subset
#print axioms LonelyRunner.LRC14Concrete.pair_danger_volume_le
