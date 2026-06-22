/-
  TournamentH7.LRCArcComplexity -- finite disjoint-cell core for the LRC(14)
  far-element arc-complexity bound (codex-2026-06-22-S87, HYP-2841).

  THM-546/HYP-2840 sharpen the single-far deviation by replacing the old
  `42 * sum e` arc-complexity bound with a disjoint-sector count:

      V(E') = sum_j #arcs(B_j(E')) <= #elementary cells <= 7 * sum E'.

  This file formalizes the reusable finite combinatorics behind the factor-six
  saving.  It deliberately does not define the analytic `B_j(E')` sets or their
  breakpoint partition.  Instead it proves the exact support lemma needed once
  those objects are instantiated: pairwise disjoint occupied cell sets inside a
  common finite partition cannot have total cardinality larger than the partition.
  A separate theorem packages the final `<= 7 * sumE` handoff.
-/

import Mathlib.Tactic

namespace LonelyRunner
namespace ArcComplexity

open scoped BigOperators

namespace FiniteCell

variable {ι α : Type*} [DecidableEq α]

/-- The total occupied-cell count of an indexed finite family. -/
def occupiedCount (I : Finset ι) (B : ι → Finset α) : ℕ :=
  ∑ i ∈ I, (B i).card

/-- If the occupied cell sets are pairwise disjoint on `I`, their total
cardinality is the cardinality of their finite union. -/
theorem occupiedCount_eq_card_biUnion
    (I : Finset ι) (B : ι → Finset α)
    (hdisj : (I : Set ι).PairwiseDisjoint B) :
    occupiedCount I B = (I.biUnion B).card := by
  rw [occupiedCount, Finset.card_biUnion hdisj]

/-- **Disjoint cells cost at most the ambient partition.**  This is the abstract
Lean form of the HYP-2840 factor-six improvement: the `B_j` are disjoint and all
sit inside the same elementary breakpoint-cell partition, so summing their cell
counts never multiplies the partition size by the number of sectors. -/
theorem occupiedCount_le_cells
    (I : Finset ι) (cells : Finset α) (B : ι → Finset α)
    (hsub : ∀ i ∈ I, B i ⊆ cells)
    (hdisj : (I : Set ι).PairwiseDisjoint B) :
    occupiedCount I B ≤ cells.card := by
  rw [occupiedCount_eq_card_biUnion I B hdisj]
  exact Finset.card_le_card (Finset.biUnion_subset.mpr hsub)

/-- Final arithmetic handoff for the LRC(14) arc-complexity application: once the
common breakpoint partition has at most `7 * sumE` cells, the disjoint occupied
cell count has the same bound. -/
theorem occupiedCount_le_seven_mul_sum
    (I : Finset ι) (cells : Finset α) (B : ι → Finset α) (sumE : ℕ)
    (hsub : ∀ i ∈ I, B i ⊆ cells)
    (hdisj : (I : Set ι).PairwiseDisjoint B)
    (hcells : cells.card ≤ 7 * sumE) :
    occupiedCount I B ≤ 7 * sumE :=
  (occupiedCount_le_cells I cells B hsub hdisj).trans hcells

end FiniteCell

/-! ## Axiom audit -/

#print axioms FiniteCell.occupiedCount_eq_card_biUnion
#print axioms FiniteCell.occupiedCount_le_cells
#print axioms FiniteCell.occupiedCount_le_seven_mul_sum

end ArcComplexity
end LonelyRunner
