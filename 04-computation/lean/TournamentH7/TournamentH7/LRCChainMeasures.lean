/-
  TournamentH7.LRCChainMeasures  (boxeph-2026-07-09-S16, HYP-5853 Lean layer 3)

  The μ₃ certificate by the pair-upgrade pattern: the joint danger of a
  doubling TRIPLE `(y, 2y, 4y)` on `[0,1]` embeds in FIVE explicit
  intervals of total length `2/7 = 1 − μ₃` (`μ₃ = 5/7`, HYP-5853's exact
  value).  A triple among the accounted runners therefore costs `2/7`
  instead of the union bound `3/7` — the `+1/7`-per-triple upgrade, formal.

  Pattern (as in `LRCPairUpgrade`): a band violation at level `j` yields an
  integer witness `m` with `|2^j y − m| < 1/14`, forced into `{0,…,2^j}`;
  each witness case lands in one of the five merged intervals
    `[0,1/14) ∪ (13/56,15/56) ∪ (13/28,15/28) ∪ (41/56,43/56) ∪ (13/14,1]`
  (the even-numerator level-2 intervals are absorbed by shallower levels).

  Deeper `μ_L` (handoff, mechanical): each level `j` adds only its
  ODD-numerator intervals `(m/2^j ± 1/(14·2^j))`; the S11 script
  (lrc_chain_bonferroni_boxeph_S11.py) emits the exact merged lists —
  L=4: +4 intervals at 13/112,41/112,69/112,97/112 (total 5/14 = 1−9/14);
  L≥5: same recursion, danger 1−μ_L per the exact table.

  Kernel-pure target: no `sorry`, no `native_decide` (audited below).
-/
import Mathlib
import TournamentH7.LRCPairUpgrade

namespace LonelyRunner
namespace LRC14Concrete

open Set

/-- **The triple-danger inclusion.**  For `y ∈ [0,1]`, a `1/14`-band
violation by any of `y, 2y, 4y` forces `y` into one of five explicit
intervals. -/
theorem triple_danger_subset :
    {y : ℝ | y ∈ Icc (0:ℝ) 1 ∧
      (nearInt y < 1/14 ∨ nearInt (2*y) < 1/14 ∨ nearInt (4*y) < 1/14)}
      ⊆ Ico (0:ℝ) (1/14) ∪ Ioo (13/56 : ℝ) (15/56) ∪ Ioo (13/28 : ℝ) (15/28)
        ∪ Ioo (41/56 : ℝ) (43/56) ∪ Ioc (13/14 : ℝ) 1 := by
  rintro y ⟨⟨hy0, hy1⟩, hdanger⟩
  rcases hdanger with h0 | h1 | h2
  · obtain ⟨m, hm⟩ := exists_int_near_of_nearInt_lt h0
    have hmlo : (-1 : ℤ) < m := by
      by_contra hc; push_neg at hc
      have : (m : ℝ) ≤ -1 := by exact_mod_cast hc
      have := abs_lt.mp hm; linarith
    have hmhi : m < 2 := by
      by_contra hc; push_neg at hc
      have : (2 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hc
      have := abs_lt.mp hm; linarith
    interval_cases m
    · push_cast at hm; have h := abs_lt.mp hm
      exact Or.inl (Or.inl (Or.inl (Or.inl ⟨hy0, by linarith [h.2]⟩)))
    · push_cast at hm; have h := abs_lt.mp hm
      exact Or.inr ⟨by linarith [h.1], hy1⟩
  · obtain ⟨m, hm⟩ := exists_int_near_of_nearInt_lt h1
    have hmlo : (-1 : ℤ) < m := by
      by_contra hc; push_neg at hc
      have : (m : ℝ) ≤ -1 := by exact_mod_cast hc
      have := abs_lt.mp hm; linarith
    have hmhi : m < 3 := by
      by_contra hc; push_neg at hc
      have : (3 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hc
      have := abs_lt.mp hm; linarith
    interval_cases m
    · push_cast at hm; have h := abs_lt.mp hm
      exact Or.inl (Or.inl (Or.inl (Or.inl ⟨hy0, by linarith [h.2]⟩)))
    · push_cast at hm; have h := abs_lt.mp hm
      exact Or.inl (Or.inl (Or.inr ⟨by linarith [h.1], by linarith [h.2]⟩))
    · push_cast at hm; have h := abs_lt.mp hm
      exact Or.inr ⟨by linarith [h.1], hy1⟩
  · obtain ⟨m, hm⟩ := exists_int_near_of_nearInt_lt h2
    have hmlo : (-1 : ℤ) < m := by
      by_contra hc; push_neg at hc
      have : (m : ℝ) ≤ -1 := by exact_mod_cast hc
      have := abs_lt.mp hm; linarith
    have hmhi : m < 5 := by
      by_contra hc; push_neg at hc
      have : (5 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hc
      have := abs_lt.mp hm; linarith
    interval_cases m
    · push_cast at hm; have h := abs_lt.mp hm
      exact Or.inl (Or.inl (Or.inl (Or.inl ⟨hy0, by linarith [h.2]⟩)))
    · push_cast at hm; have h := abs_lt.mp hm
      exact Or.inl (Or.inl (Or.inl (Or.inr ⟨by linarith [h.1], by linarith [h.2]⟩)))
    · push_cast at hm; have h := abs_lt.mp hm
      exact Or.inl (Or.inl (Or.inr ⟨by linarith [h.1], by linarith [h.2]⟩))
    · push_cast at hm; have h := abs_lt.mp hm
      exact Or.inl (Or.inr ⟨by linarith [h.1], by linarith [h.2]⟩)
    · push_cast at hm; have h := abs_lt.mp hm
      exact Or.inr ⟨by linarith [h.1], hy1⟩

/-- **The triple upgrade:** the joint danger measure of a doubling triple is
`≤ 2/7 = 1 − μ₃` — a full `1/7` below the union bound `3/7`. -/
theorem triple_danger_volume_le :
    MeasureTheory.volume
      {y : ℝ | y ∈ Icc (0:ℝ) 1 ∧
        (nearInt y < 1/14 ∨ nearInt (2*y) < 1/14 ∨ nearInt (4*y) < 1/14)}
      ≤ ENNReal.ofReal (2/7) := by
  refine le_trans (MeasureTheory.measure_mono triple_danger_subset) ?_
  set A := Ico (0:ℝ) (1/14)
  set B := Ioo (13/56 : ℝ) (15/56)
  set C := Ioo (13/28 : ℝ) (15/28)
  set D := Ioo (41/56 : ℝ) (43/56)
  set E := Ioc (13/14 : ℝ) 1
  calc MeasureTheory.volume (A ∪ B ∪ C ∪ D ∪ E)
      ≤ MeasureTheory.volume (A ∪ B ∪ C ∪ D) + MeasureTheory.volume E :=
        MeasureTheory.measure_union_le _ _
    _ ≤ (MeasureTheory.volume (A ∪ B ∪ C) + MeasureTheory.volume D)
        + MeasureTheory.volume E := by
        gcongr
        exact MeasureTheory.measure_union_le _ _
    _ ≤ ((MeasureTheory.volume (A ∪ B) + MeasureTheory.volume C)
        + MeasureTheory.volume D) + MeasureTheory.volume E := by
        gcongr
        exact MeasureTheory.measure_union_le _ _
    _ ≤ (((MeasureTheory.volume A + MeasureTheory.volume B)
        + MeasureTheory.volume C) + MeasureTheory.volume D)
        + MeasureTheory.volume E := by
        gcongr
        exact MeasureTheory.measure_union_le _ _
    _ ≤ ENNReal.ofReal (2/7) := by
        rw [show A = Ico (0:ℝ) (1/14) from rfl,
          show B = Ioo (13/56 : ℝ) (15/56) from rfl,
          show C = Ioo (13/28 : ℝ) (15/28) from rfl,
          show D = Ioo (41/56 : ℝ) (43/56) from rfl,
          show E = Ioc (13/14 : ℝ) 1 from rfl,
          Real.volume_Ico, Real.volume_Ioo, Real.volume_Ioo, Real.volume_Ioo,
          Real.volume_Ioc,
          ← ENNReal.ofReal_add (by norm_num) (by norm_num),
          ← ENNReal.ofReal_add (by norm_num) (by norm_num),
          ← ENNReal.ofReal_add (by norm_num) (by norm_num),
          ← ENNReal.ofReal_add (by norm_num) (by norm_num)]
        apply ENNReal.ofReal_le_ofReal
        norm_num

end LRC14Concrete
end LonelyRunner

-- kernel-purity audit
#print axioms LonelyRunner.LRC14Concrete.triple_danger_subset
#print axioms LonelyRunner.LRC14Concrete.triple_danger_volume_le
