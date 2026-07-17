import TournamentH7.LRCB5RelationBudget
import TournamentH7.LRCPairRatioLayerArithmetic

/-!
# The exact higher-support allowance after the weighted pair layer

The original relation-budget consumer charged the convenient pair bound
`13/30`.  The path-only THM-954 layer cake is strictly sharper.  This module
returns that exact slack to the signed support-three/four/five budget.
-/

namespace LonelyRunner
namespace LRCB5SharpPairBudget

open LRCB5RelationBudget LRCPairRatioLayerArithmetic

noncomputable section

/-- Exact signed higher-support allowance left by the path-only pair bound. -/
def pathOnlyHigherAllowance : ℝ := 541592875 / 5883525648

theorem pathOnlyHigherAllowance_eq :
    equilibrium - pairWeight *
        ((negativePairTierBoundPathOnly : ℚ) : ℝ) =
      pathOnlyHigherAllowance := by
  norm_num [equilibrium, pairWeight, negativePairTierBoundPathOnly,
    pathOnlyHigherAllowance]

/-- The weighted-ratio pair theorem returns a positive amount of budget beyond
the old `13/30` socket. -/
theorem pathOnlyHigherAllowance_sub_old_eq :
    pathOnlyHigherAllowance - 7712 / 84035 =
      8270807 / 29417628240 := by
  norm_num [pathOnlyHigherAllowance]

theorem old_allowance_lt_pathOnlyHigherAllowance :
    (7712 / 84035 : ℝ) < pathOnlyHigherAllowance := by
  rw [← sub_pos]
  rw [pathOnlyHigherAllowance_sub_old_eq]
  norm_num

/-- Sharp signed relation-budget consumer using the actual path-only
weighted-ratio bound rather than its coarse `13/30` relaxation. -/
theorem relationModel_pos_of_weightedRatioPair_split
    (mass2 mass3 mass4 mass5 : ℝ)
    (hpair :
      -((negativePairTierBoundPathOnly : ℚ) : ℝ) ≤ mass2)
    (hhigher :
      harmfulHigherContribution mass3 mass4 mass5 <
        pathOnlyHigherAllowance) :
    0 < relationModel mass2 mass3 mass4 mass5 := by
  have hweight : (0 : ℝ) ≤ pairWeight := weights_nonnegative.1
  have hpairWeighted :
      -(pairWeight * ((negativePairTierBoundPathOnly : ℚ) : ℝ)) ≤
        pairWeight * mass2 := by
    have hmul := mul_le_mul_of_nonneg_left hpair hweight
    linarith
  rw [← pathOnlyHigherAllowance_eq] at hhigher
  dsimp [relationModel, harmfulHigherContribution] at hpairWeighted hhigher ⊢
  linarith

#print axioms pathOnlyHigherAllowance_eq
#print axioms pathOnlyHigherAllowance_sub_old_eq
#print axioms relationModel_pos_of_weightedRatioPair_split

end
end LRCB5SharpPairBudget
end LonelyRunner
