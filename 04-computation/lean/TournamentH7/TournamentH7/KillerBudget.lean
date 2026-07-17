/- KillerBudget.lean — mac-mini-2026-07-16-S128.
   Rung three of the LRC(14) ladder: the KILLER-BUDGET COMPOSITION.
   Composing FragmentationCount.fragmentation over a finite set W of moduli
   ("killers"): the union of their bad sets meets [x, x+L] in measure at most
   Σ_{w ∈ W} (w*L + 2*lam + 1) * (2*lam/w) — the budget inequality behind the
   THM-883 multi-killer sweep and the covering-min rigidity chain, and the
   good-set floor it implies. -/
import TournamentH7.FragmentationCount

open MeasureTheory Set

namespace LRC14

/-- **The killer budget.** The union over a finite modulus set `W` of bad sets meets
    `[x, x+L]` in measure at most the sum of the per-modulus fragmentation bounds. -/
theorem killer_budget (W : Finset ℕ) (hW : ∀ w ∈ W, 0 < w) (x L lam : ℝ)
    (hlam : 0 < lam) (hL : 0 ≤ L) :
    volume ((⋃ w ∈ W, badSet w lam) ∩ Icc x (x + L))
      ≤ ∑ w ∈ W, ENNReal.ofReal (((w : ℝ) * L + 2 * lam + 1) * (2 * lam / w)) := by
  classical
  have hdistrib : (⋃ w ∈ W, badSet w lam) ∩ Icc x (x + L)
      = ⋃ w ∈ W, (badSet w lam ∩ Icc x (x + L)) := by
    ext t
    simp only [mem_inter_iff, mem_iUnion, exists_prop]
    tauto
  rw [hdistrib]
  calc volume (⋃ w ∈ W, (badSet w lam ∩ Icc x (x + L)))
      ≤ ∑ w ∈ W, volume (badSet w lam ∩ Icc x (x + L)) :=
        measure_biUnion_finset_le _ _
    _ ≤ ∑ w ∈ W, ENNReal.ofReal (((w : ℝ) * L + 2 * lam + 1) * (2 * lam / w)) := by
        refine Finset.sum_le_sum fun w hw => ?_
        exact fragmentation w (hW w hw) x L lam hlam hL

/-- **The good-set floor.** If the total budget `B` is finite, the set of points of
    `[x, x+L]` avoiding every bad arc has measure at least `L − B` (stated in the
    ENNReal subtraction form the covering chain uses). -/
theorem good_floor (W : Finset ℕ) (hW : ∀ w ∈ W, 0 < w) (x L lam : ℝ)
    (hlam : 0 < lam) (hL : 0 ≤ L) :
    ENNReal.ofReal L
        - ∑ w ∈ W, ENNReal.ofReal (((w : ℝ) * L + 2 * lam + 1) * (2 * lam / w))
      ≤ volume (Icc x (x + L) \ ⋃ w ∈ W, badSet w lam) := by
  classical
  have hvol : volume (Icc x (x + L)) = ENNReal.ofReal L := by
    rw [Real.volume_Icc]
    congr 1
    ring
  have hsplit : volume (Icc x (x + L))
      ≤ volume (Icc x (x + L) \ ⋃ w ∈ W, badSet w lam)
        + volume ((⋃ w ∈ W, badSet w lam) ∩ Icc x (x + L)) := by
    have hcover : Icc x (x + L)
        ⊆ (Icc x (x + L) \ ⋃ w ∈ W, badSet w lam)
          ∪ ((⋃ w ∈ W, badSet w lam) ∩ Icc x (x + L)) := by
      intro t ht
      by_cases hb : t ∈ ⋃ w ∈ W, badSet w lam
      · exact Or.inr ⟨hb, ht⟩
      · exact Or.inl ⟨ht, hb⟩
    exact le_trans (measure_mono hcover) (measure_union_le _ _)
  have hbudget := killer_budget W hW x L lam hlam hL
  rw [hvol] at hsplit
  exact tsub_le_iff_right.mpr
    (le_trans hsplit (add_le_add le_rfl hbudget))

end LRC14
