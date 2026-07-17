/-
  Exact obstruction-avoidance interfaces for all three selected witnesses.

  Fixed-phase Zarankiewicz counting is now complete: every failure is a
  saturated parallel-class partition.  The remaining analytic task is to make
  one harmonic-good phase avoid that partition.  This module packages those
  three minimal selectors and proves that they imply the original supplier
  interfaces.  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCSelectedWitnessTwoTwo
import TournamentH7.LRCSelectedWitnessTwoFourFour
import TournamentH7.LRCSelectedWitnessUniformThree

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner

/-- Harmonic-good phase avoidance of the exhaustive `(2,2,q)` obstruction. -/
def TwoTwoParallelAvoidanceSupply : Prop :=
  ∀ (v : Fin 13 → ℤ), (∀ i, v i ≠ 0) → ∀ (g : ℤ), 2 ≤ g →
    ∀ first second third : Fin 13,
    first ≠ second → first ≠ third → second ≠ third →
    (∀ j, j ≠ first → j ≠ second → j ≠ third → g ∣ v j) →
    ¬ g ∣ v third →
    g / (Int.gcd (v first) g : ℤ) = 2 →
    g / (Int.gcd (v second) g : ℤ) = 2 →
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g first second third u ∧
      ¬ TwoTwoParallelObstruction
        (v first) (v second) (v third) g u

/-- Harmonic-good phase avoidance of the complete `(2,4,4)` prefix code. -/
def TwoFourFourParallelAvoidanceSupply : Prop :=
  ∀ (v : Fin 13 → ℤ), (∀ i, v i ≠ 0) → ∀ (g : ℤ), 2 ≤ g →
    ∀ two fourA fourB : Fin 13,
    two ≠ fourA → two ≠ fourB → fourA ≠ fourB →
    (∀ j, j ≠ two → j ≠ fourA → j ≠ fourB → g ∣ v j) →
    g / (Int.gcd (v two) g : ℤ) = 2 →
    g / (Int.gcd (v fourA) g : ℤ) = 4 →
    g / (Int.gcd (v fourB) g : ℤ) = 4 →
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g two fourA fourB u ∧
      ¬ TwoFourFourParallelPartition (v two) (v fourA) (v fourB) g u

/-- Harmonic-good phase avoidance of the balanced uniform-q-three partition. -/
def UniformThreeParallelAvoidanceSupply : Prop :=
  ∀ (v : Fin 13 → ℤ), (∀ i, v i ≠ 0) → ∀ (g : ℤ), 2 ≤ g →
    ∀ first second third : Fin 13,
    first ≠ second → first ≠ third → second ≠ third →
    (∀ j, j ≠ first → j ≠ second → j ≠ third → g ∣ v j) →
    g / (Int.gcd (v first) g : ℤ) = 3 →
    g / (Int.gcd (v second) g : ℤ) = 3 →
    g / (Int.gcd (v third) g : ℤ) = 3 →
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g first second third u ∧
      ¬ UniformThreeBadPartition (v first) (v second) (v third) g u

theorem twoTwoSelectedWitnessSupply_of_parallelAvoidance
    (hsupply : TwoTwoParallelAvoidanceSupply) :
    TwoTwoSelectedWitnessSupply := by
  intro v hv g hg first second third h12 h13 h23 hdvd hthird hqFirst hqSecond
  apply twoTwo_selectedWitness_of_parallelObstruction_avoidance
    v g hg first second third hthird hqFirst hqSecond
  exact hsupply v hv g hg first second third h12 h13 h23 hdvd hthird
    hqFirst hqSecond

theorem twoFourFourSelectedWitnessSupply_of_parallelAvoidance
    (hsupply : TwoFourFourParallelAvoidanceSupply) :
    TwoFourFourSelectedWitnessSupply := by
  intro v hv g hg two fourA fourB hTwoA hTwoB hAB hdvd hqTwo hqFourA hqFourB
  apply twoFourFour_selectedWitness_of_partition_avoidance
    v g hg two fourA fourB hqTwo hqFourA hqFourB
  exact hsupply v hv g hg two fourA fourB hTwoA hTwoB hAB hdvd
    hqTwo hqFourA hqFourB

theorem uniformThreeSelectedWitnessSupply_of_parallelAvoidance
    (hsupply : UniformThreeParallelAvoidanceSupply) :
    UniformThreeSelectedWitnessSupply := by
  intro v hv g hg first second third h12 h13 h23 hdvd hqFirst hqSecond hqThird
  apply uniformThree_selectedWitness_of_partition_avoidance
    v g hg first second third hqFirst hqSecond hqThird
  exact hsupply v hv g hg first second third h12 h13 h23 hdvd
    hqFirst hqSecond hqThird

/-! ## Axiom audit -/

#print axioms twoTwoSelectedWitnessSupply_of_parallelAvoidance
#print axioms twoFourFourSelectedWitnessSupply_of_parallelAvoidance
#print axioms uniformThreeSelectedWitnessSupply_of_parallelAvoidance

end LRC14Grand
end LonelyRunner
