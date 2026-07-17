import TournamentH7.LRCLeverageIdentity

/-!
# Leverage-identity consumer demo: the {1,2,3} system end to end

kind-pasteur-S128c37.  The worked template for turning a packet's depth-sweep
cell data into a kernel loneliness fact via
`LRCLeverageIdentity.goodMass_pos_of_bonferroni_pos`.

The 13 cells below are the exact depth sweep of the bad sets of speeds
`{1,2,3}` at `lam = 1/14` (weights sum to 1).  The pipeline: the level-1
truncation `B1 = 1 - S1 = 4/7` is positive, so the good mass is positive -
WITHOUT computing it - and the identity pins it exactly at `29/42` as a
cross-check.  For a real packet the same lemmas apply verbatim to its
generated cell list; only the sweep encoding remains, per the manifest.
-/

namespace LonelyRunner
namespace LRCLeverageDemo

open LRCLeverageIdentity

/-- Cell weights of the {1,2,3} sweep at lam = 1/14 (13 cells, total mass 1). -/
def demoW : Fin 13 → ℚ :=
  ![1/42, 1/84, 1/28, 5/21, 1/21, 3/28, 1/14, 3/28, 1/21, 5/21, 1/28, 1/84, 1/42]

/-- Cell depths (number of bad sets active on each cell). -/
def demoD : Fin 13 → ℕ :=
  ![3, 2, 1, 0, 1, 0, 1, 0, 1, 0, 1, 2, 3]

theorem demoW_nonneg : ∀ c, 0 ≤ demoW c := by
  intro c
  fin_cases c <;> norm_num [demoW]

theorem demoW_total : (∑ c, demoW c) = 1 := by
  norm_num [Fin.sum_univ_succ, demoW]

/-- `S1 = 3/7`: each speed's bad set has measure `1/7`. -/
theorem demo_S1 : binomMoment demoW demoD 1 = 3 / 7 := by
  norm_num [binomMoment, Fin.sum_univ_succ, demoW, demoD]

/-- The level-1 truncation `B1 = 1 - S1 = 4/7 > 0`. -/
theorem demo_B1_pos :
    0 < ∑ k ∈ Finset.range 2, ((-1 : ℚ)) ^ k * binomMoment demoW demoD k := by
  norm_num [Finset.sum_range_succ, binomMoment, Fin.sum_univ_succ, demoW, demoD]

/-- **The pipeline**: positive odd truncation implies positive good mass,
    by the certificate theorem alone. -/
theorem demo_goodMass_pos : 0 < goodMass demoW demoD :=
  goodMass_pos_of_bonferroni_pos demoW demoD 1 ⟨0, rfl⟩ demoW_nonneg demo_B1_pos

/-- Cross-check: the good mass is exactly `29/42`. -/
theorem demo_goodMass_eq : goodMass demoW demoD = 29 / 42 := by
  norm_num [goodMass, Finset.sum_filter, Fin.sum_univ_succ, demoW, demoD]

/-- Cross-check of the identity's tail at `m = 1`. -/
theorem demo_identity_check :
    ∑ k ∈ Finset.range 2, ((-1 : ℚ)) ^ k * binomMoment demoW demoD k
      = goodMass demoW demoD - leveragedTail demoW demoD 1 := by
  norm_num [Finset.sum_range_succ, binomMoment, goodMass, leveragedTail,
    Finset.sum_filter, Fin.sum_univ_succ, demoW, demoD]

/-! ## Axiom audit -/

#print axioms demo_goodMass_pos
#print axioms demo_identity_check

end LRCLeverageDemo
end LonelyRunner
