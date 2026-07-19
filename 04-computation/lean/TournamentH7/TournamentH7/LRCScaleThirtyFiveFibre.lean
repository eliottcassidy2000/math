/-
  THM-1249: the kernel-pure consumer for the scale-35 Z/5 fibre ceiling.

  The external exact certificate proves that the independent-choice anchor
  relaxation is at most 31 on each scalar-surviving owner obligation.  This
  file checks the final logical step: the generic nested-fibre relaxation then
  bounds the literal union by 31, so it cannot cover all 35 sheets.

  The 3,002,076-row enumeration remains in the companion exact Python
  certificate.  No global n=12 sporadic-branch conclusion is asserted here.
-/
import TournamentH7.LRCNestedFibreRelaxation

namespace LonelyRunner
namespace LRC14
namespace ScaleThirtyFiveFibre

open Finset
open NestedFibreRelaxation

theorem z5_ceiling_strict : (31 : ℕ) < 35 := by
  norm_num

theorem card_le_thirty_one_not_full (U : Finset (Fin 35))
    (hU : U.card ≤ 31) : U ≠ univ := by
  intro hfull
  rw [hfull] at hU
  norm_num at hU

variable {provider choice : Type*}
variable [DecidableEq provider]

/-- The exact finite bank supplies `hceiling` at every surviving owner.  The
generic pointwise-maximum relaxation then prevents literal 35-sheet cover. -/
theorem z5_relaxation_blocks_full_cover
    (providers anchors : Finset provider)
    (options : provider → Finset choice)
    (mask : provider → choice → Finset (Fin 35))
    (chosen : provider → choice)
    (hanchors : anchors ⊆ providers)
    (hchosen : ∀ i ∈ providers \ anchors, chosen i ∈ options i)
    (hceiling :
      (anchorUnion anchors mask chosen).card +
          ∑ i ∈ providers \ anchors,
            deviationMax (anchorUnion anchors mask chosen) options mask i ≤ 31) :
    providerUnion providers mask chosen ≠ univ := by
  apply card_le_thirty_one_not_full
  exact (card_providerUnion_le_anchor_add_pointwiseMax
    providers anchors options mask chosen hanchors hchosen).trans hceiling

#print axioms z5_ceiling_strict
#print axioms card_le_thirty_one_not_full
#print axioms z5_relaxation_blocks_full_cover

end ScaleThirtyFiveFibre
end LRC14
end LonelyRunner
