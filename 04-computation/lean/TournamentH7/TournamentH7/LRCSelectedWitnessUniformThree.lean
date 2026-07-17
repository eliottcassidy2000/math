/-
  Exact parallel-class obstruction for the selected uniform `(3,3,3)` witness.

  Failure is precisely the balanced disjoint cover by the three q-three rows,
  again equality in `z(3,g;2,1)=g`.  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCSelectedWitnessCommon

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

noncomputable section

theorem not_hasThreeDetunedGoodBranch_of_uniformThreeBadPartition
    (deltaOne deltaTwo deltaThree g : ℤ) (u : ℝ)
    (hpartition :
      UniformThreeBadPartition deltaOne deltaTwo deltaThree g u) :
    ¬ HasThreeDetunedGoodBranch deltaOne deltaTwo deltaThree g u := by
  rintro ⟨c, hcIco, hcOne, hcTwo, hcThree⟩
  have hcUnion : c ∈
      detunedBadBranches deltaOne g u ∪
        detunedBadBranches deltaTwo g u ∪
        detunedBadBranches deltaThree g u := by
    rw [hpartition.cover]
    exact hcIco
  simp [hcOne, hcTwo, hcThree] at hcUnion

/-- At every modulus, uniform-q-three failure is exactly the balanced
pairwise-disjoint parallel-class partition. -/
theorem not_hasThreeDetunedGoodBranch_uniform_three_iff_badPartition
    (deltaOne deltaTwo deltaThree g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hqOne : g / (Int.gcd deltaOne g : ℤ) = 3)
    (hqTwo : g / (Int.gcd deltaTwo g : ℤ) = 3)
    (hqThree : g / (Int.gcd deltaThree g : ℤ) = 3) :
    ¬ HasThreeDetunedGoodBranch deltaOne deltaTwo deltaThree g u ↔
      UniformThreeBadPartition deltaOne deltaTwo deltaThree g u := by
  constructor
  · exact uniformThreeBadPartition_of_noGoodBranch
      deltaOne deltaTwo deltaThree g u (by omega) hqOne hqTwo hqThree
  · exact not_hasThreeDetunedGoodBranch_of_uniformThreeBadPartition
      deltaOne deltaTwo deltaThree g u

/-- Exact phase-avoidance consumer for the uniform q-three selected witness. -/
theorem uniformThree_selectedWitness_of_partition_avoidance
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g)
    (iOne iTwo iThree : Fin 13)
    (hqOne : g / (Int.gcd (v iOne) g : ℤ) = 3)
    (hqTwo : g / (Int.gcd (v iTwo) g : ℤ) = 3)
    (hqThree : g / (Int.gcd (v iThree) g : ℤ) = 3)
    (hselector : ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g iOne iTwo iThree u ∧
      ¬ UniformThreeBadPartition
        (v iOne) (v iTwo) (v iThree) g u) :
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g iOne iTwo iThree u ∧
      HasThreeDetunedGoodBranch (v iOne) (v iTwo) (v iThree) g u := by
  obtain ⟨u, hharmonic, hnotPartition⟩ := hselector
  refine ⟨u, hharmonic, ?_⟩
  by_contra hfail
  exact hnotPartition
    ((not_hasThreeDetunedGoodBranch_uniform_three_iff_badPartition
      (v iOne) (v iTwo) (v iThree) g u hg hqOne hqTwo hqThree).mp hfail)

/-! ## Axiom audit -/

#print axioms not_hasThreeDetunedGoodBranch_uniform_three_iff_badPartition
#print axioms uniformThree_selectedWitness_of_partition_avoidance

end
end LRC14Grand
end LonelyRunner
