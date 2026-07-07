/-
  TournamentH7.LRCAP30Floor — the AP₃₀ density-floor certificate via the Farey roof
  (boxeph-2026-07-07-S2).  UNCONDITIONAL, extends death-star-S2's AP₂₀ (diameter ≤ 19)
  to diameter ≤ 29.

  death-star's `LRCAP20Certificate` proved `μ_{1/7}(AP₂₀) ≥ m_P` by exhibiting TWO empty
  arcs by hand.  Those arcs are exactly the roof-superlevel intervals of the two Farey
  cells adjacent to the nodes `0/1` and `1/1`: for the AP `{1,…,k}` the near-0 cell
  `(0/1, 1/k)` has roof `1 − (k−1)x`, giving the good interval `(0, 6/(7(k−1)))`, and the
  mirror near `1` gives `((7k−13)/(7(k−1)), 1)` — each of length `6/(7(k−1))`.

  This file DERIVES both intervals from `FareyRoofBridge.good_of_roof_gt` (opus's roof)
  instead of by hand, and pushes `k` from `20` (death-star) to `30`:

    `μ_{1/7}(AP₃₀) ≥ 6/203 + 6/203 = 12/203 ≈ 0.0591 ≥ m_P = 14249/252252 ≈ 0.0565`.

  Two endpoint intervals clear `m_P` for every `k ≤ 31` (`12/(7(k−1)) ≥ m_P`), so the
  systematic roof route reaches diameter ≤ 29 with the SAME two intervals death-star used
  for ≤ 19 — the extra reach is free once the intervals come from the roof, not by hand.

  Kernel-pure: no sorry, no native_decide.
-/
import Mathlib
import TournamentH7.LRCFareyRoofBridge
import TournamentH7.LRCTailDiameter

namespace LonelyRunner
namespace AP30Floor

open TailDiameter FareyRoofBridge MeasureTheory
open scoped ENNReal

/-- Near `x = 0`: the roof interval `(0, 6/203)` of the Farey-30 cell `(0/1, 1/30)` is good
for `AP₃₀ = {1,…,30}`.  (roof `= 1 − 29x > 1/7` iff `x < 6/203`; and `6/203 < 1/30`.) -/
theorem near0_subset_good :
    Set.Ioo (0 : ℝ) (6 / 203) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 30) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx
  obtain ⟨hx0, hxu⟩ := hx
  refine ⟨?_, le_of_lt hx0, by linarith⟩
  have h := good_of_roof_gt (p := (0 : ℤ)) (q := (1 : ℤ)) (p' := (1 : ℤ)) (q' := (30 : ℤ))
    (k := (30 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hx0])
    (by push_cast; nlinarith [hxu])
    (by push_cast; nlinarith [hxu])
  exact h

/-- Near `x = 1`: the roof interval `(197/203, 1)` of the Farey-30 cell `(29/30, 1/1)` is
good for `AP₃₀`.  (roof `= 29x − 28 > 1/7` iff `x > 197/203`.) -/
theorem near1_subset_good :
    Set.Ioo (197 / 203 : ℝ) 1 ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 30) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx
  obtain ⟨hxl, hx1⟩ := hx
  refine ⟨?_, by linarith, le_of_lt hx1⟩
  have h := good_of_roof_gt (p := (29 : ℤ)) (q := (30 : ℤ)) (p' := (1 : ℤ)) (q' := (1 : ℤ))
    (k := (30 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl])
    (by push_cast; nlinarith [hx1])
    (by push_cast; nlinarith [hxl])
  exact h

theorem near0_disjoint_near1 :
    Disjoint (Set.Ioo (0 : ℝ) (6 / 203)) (Set.Ioo (197 / 203 : ℝ) 1) := by
  apply Set.disjoint_left.mpr
  intro x hx1 hx2
  have := hx1.2; have := hx2.1; linarith

/-- **THE AP₃₀ CERTIFICATE (unconditional).**  `μ_{1/7}(AP₃₀) ≥ m_P = 14249/252252`, from the
two roof-derived endpoint intervals, each of length `6/203`. -/
theorem ap30_certificate :
    ENNReal.ofReal ((14249 : ℝ) / 252252) ≤ muGood (1 / 7) (Finset.Icc (1 : ℤ) 30) := by
  have hbound :
      ENNReal.ofReal (6 / 203 - 0) + ENNReal.ofReal (1 - 197 / 203)
        ≤ muGood (1 / 7) (Finset.Icc (1 : ℤ) 30) := by
    have hu : Set.Ioo (0 : ℝ) (6 / 203) ∪ Set.Ioo (197 / 203 : ℝ) 1
        ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 30) ∩ Set.Icc (0 : ℝ) 1 :=
      Set.union_subset near0_subset_good near1_subset_good
    calc ENNReal.ofReal (6 / 203 - 0) + ENNReal.ofReal (1 - 197 / 203)
        = volume (Set.Ioo (0 : ℝ) (6 / 203)) + volume (Set.Ioo (197 / 203 : ℝ) 1) := by
          rw [Real.volume_Ioo, Real.volume_Ioo]
      _ = volume (Set.Ioo (0 : ℝ) (6 / 203) ∪ Set.Ioo (197 / 203 : ℝ) 1) :=
          (measure_union near0_disjoint_near1 measurableSet_Ioo).symm
      _ ≤ volume (Good (1 / 7) (Finset.Icc (1 : ℤ) 30) ∩ Set.Icc (0 : ℝ) 1) := measure_mono hu
      _ = muGood (1 / 7) (Finset.Icc (1 : ℤ) 30) := rfl
  refine le_trans ?_ hbound
  rw [← ENNReal.ofReal_add (by norm_num) (by norm_num)]
  exact ENNReal.ofReal_le_ofReal (by norm_num)

/-- Rephrased on `{0,…,29}` (death-star's convention) via rotation invariance
(`good_translate`), ready to feed `TailDiameter.muGood_ge_APD` for diameter `≤ 29`. -/
theorem ap30_certificate_icc0 :
    ENNReal.ofReal ((14249 : ℝ) / 252252) ≤ muGood (1 / 7) (Finset.Icc (0 : ℤ) 29) := by
  have hE : (Finset.Icc (1 : ℤ) 30).image (fun e => e - 1) = Finset.Icc (0 : ℤ) 29 := by
    ext n
    simp only [Finset.mem_image, Finset.mem_Icc]
    constructor
    · rintro ⟨a, ⟨h1, h2⟩, rfl⟩; omega
    · intro h; exact ⟨n + 1, ⟨by omega, by omega⟩, by omega⟩
  have htr : muGood (1 / 7) (Finset.Icc (0 : ℤ) 29) = muGood (1 / 7) (Finset.Icc (1 : ℤ) 30) := by
    have := muGood_translate (1 / 7) (Finset.Icc (1 : ℤ) 30) 1
    rw [hE] at this
    exact this
  rw [htr]; exact ap30_certificate

end AP30Floor
end LonelyRunner
