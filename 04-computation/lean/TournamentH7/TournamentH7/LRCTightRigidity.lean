/-
# LRCTightRigidity — the extremal-uniqueness rigidity is a sufficient (≥ LRC-hard) route

`μ(S) = 0 ⟹ S is a dilated interval` is the **extremal-uniqueness** statement for LRC(14): the only
families whose safe set is measure-zero (tight, `Mreach = 1/14` exactly) are the dilations `c·{1,…,13}` of
the arithmetic progression.  This is **stronger than LRC(14) itself** — this file proves it *implies*
LRC(14), so it cannot be proved here.  It is left as a named hypothesis `TightRigidity`; what is proved is
the reduction

  `TightRigidity → SafeMeasureFloor → LRC14Statement`,  kernel-pure.

The mechanism is elementary: a residual family is *scale-gapped* (`GapFamily`, ratio `> 13`), whereas a
dilated interval has ratio exactly `13` (`max = 13·min`), so **a scale-gapped family is never dilated**
(`not_dilated_of_gapFamily`).  Hence if `μ = 0` forced dilation, no residual family could have `μ = 0`,
giving the measure floor.

Search evidence (`lrc14_tight_locus_search`): among all `13`-subsets of `[1,N]` for `N ≤ 18`, the ONLY
measure-zero family is `{1,…,13}` (the `c = 1` dilate; `c ≥ 2` needs `N ≥ 26`).  No non-dilated tight
family exists in range — consistent with, but not a proof of, `TightRigidity`.

kind-pasteur-2026-07-10-S127.
-/
import Mathlib
import TournamentH7.LRCResidualMeasureFloor

namespace LonelyRunner
namespace LRC14Grand

open MeasureTheory

/-- `v` is a **dilated interval**: its absolute speeds are exactly `{c, 2c, …, 13c}` for some `c > 0`. -/
def DilatedFamily (v : Fin 13 → ℤ) : Prop :=
  ∃ c : ℤ, 0 < c ∧ (Finset.univ.image fun i => |v i|) = (Finset.Icc 1 13).image (fun k => c * k)

/-- **A scale-gapped family is never a dilated interval.**  A dilated interval has every `|vᵢ| = c·k`
with `1 ≤ k ≤ 13`, so `|vᵢ| ≤ 13c ≤ 13|vⱼ|` for all `i, j` — the negation of `GapFamily`. -/
theorem not_dilated_of_gapFamily {v : Fin 13 → ℤ} (h : GapFamily v) : ¬ DilatedFamily v := by
  rintro ⟨c, hc, himg⟩
  apply h
  intro i j
  have hi : |v i| ∈ (Finset.Icc 1 13).image (fun k => c * k) := by
    rw [← himg]; exact Finset.mem_image_of_mem _ (Finset.mem_univ i)
  have hj : |v j| ∈ (Finset.Icc 1 13).image (fun k => c * k) := by
    rw [← himg]; exact Finset.mem_image_of_mem _ (Finset.mem_univ j)
  rw [Finset.mem_image] at hi hj
  obtain ⟨ki, hki, hvi⟩ := hi
  obtain ⟨kj, hkj, hvj⟩ := hj
  rw [Finset.mem_Icc] at hki hkj
  rw [← hvi, ← hvj]
  nlinarith [hki.1, hki.2, hkj.1, hkj.2, hc]

/-- **The extremal-uniqueness rigidity (open, = LRC(14) extremal uniqueness).**  Every nonzero-speed
family with a measure-zero safe set is a dilated interval `c·{1,…,13}`.  Equivalently `Mreach = 1/14` only
at the AP — stronger than LRC(14) itself. -/
def TightRigidity : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → volume (safePeriod v) = 0 → DilatedFamily v

/-- **The rigidity gives the measure floor.**  A residual family is scale-gapped, hence not dilated; so if
its safe set had measure zero the rigidity would make it dilated — contradiction.  Therefore `μ > 0`. -/
theorem safeMeasureFloor_of_tightRigidity (h : TightRigidity) : SafeMeasureFloor := by
  intro v hres
  by_cases hz : volume (safePeriod v) = 0
  · exact absurd (h v hres.1 hz) (not_dilated_of_gapFamily hres.2.2.1)
  · exact pos_iff_ne_zero.mpr hz

/-- **LRC(14) from the citation and the extremal-uniqueness rigidity.**  Kernel-pure.  This is the honest
location of `μ = 0 ⟹ dilated`: a *sufficient* route to LRC(14), hence at least as hard as the conjecture.
It is NOT proved — `TightRigidity` is the open hypothesis. -/
theorem lrc14_of_tightRigidity (cite : LRCUpTo13) (h : TightRigidity) : LRC14.LRC14Statement :=
  lrc14_of_measureFloor cite (safeMeasureFloor_of_tightRigidity h)

end LRC14Grand
end LonelyRunner
