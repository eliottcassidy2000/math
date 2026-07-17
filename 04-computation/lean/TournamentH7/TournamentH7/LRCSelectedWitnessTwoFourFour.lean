/-
  Exact parallel-class obstruction for the selected `(2,4,4)` witness.

  The three bad rows fail precisely when they are the disjoint complete binary
  prefix code `{1/2,1/4,1/4}`.  This is equality in the fixed-phase incidence
  value `z(3,g;2,1)=g`; static Zarankiewicz counting has no remaining slack.
  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCSelectedWitnessCommon

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

noncomputable section

/-- The exact pair-codegree-zero extremizer on the parallel-class circle for
a q-two/q-four/q-four triple. -/
def TwoFourFourParallelPartition
    (deltaTwo deltaFourA deltaFourB g : ℤ) (u : ℝ) : Prop :=
  ∃ classTwo classFourA classFourB,
    classTwo ∈ detunedBadBranches deltaTwo g u ∧
    classFourA ∈ detunedBadBranches deltaFourA g u ∧
    classFourB ∈ detunedBadBranches deltaFourB g u ∧
    detunedBadBranches deltaTwo g u =
      {x ∈ Finset.Ico (0 : ℤ) g | x ≡ classTwo [ZMOD 2]} ∧
    detunedBadBranches deltaFourA g u =
      {x ∈ Finset.Ico (0 : ℤ) g | x ≡ classFourA [ZMOD 4]} ∧
    detunedBadBranches deltaFourB g u =
      {x ∈ Finset.Ico (0 : ℤ) g | x ≡ classFourB [ZMOD 4]} ∧
    ¬ classTwo ≡ classFourA [ZMOD 2] ∧
    ¬ classTwo ≡ classFourB [ZMOD 2] ∧
    ¬ classFourA ≡ classFourB [ZMOD 4]

/-- A parity fiber modulo four consists of exactly two residue classes. -/
theorem modEq_four_or_modEq_four_of_modEq_two
    {a b c : ℤ} (hca : c ≡ a [ZMOD 2]) (hba : b ≡ a [ZMOD 2])
    (hab : ¬ a ≡ b [ZMOD 4]) :
    c ≡ a [ZMOD 4] ∨ c ≡ b [ZMOD 4] := by
  change c % 4 = a % 4 ∨ c % 4 = b % 4
  change c % 2 = a % 2 at hca
  change b % 2 = a % 2 at hba
  change ¬ a % 4 = b % 4 at hab
  omega

theorem not_hasThreeDetunedGoodBranch_of_twoFourFourParallelPartition
    (deltaTwo deltaFourA deltaFourB g : ℤ) (u : ℝ)
    (hpartition :
      TwoFourFourParallelPartition deltaTwo deltaFourA deltaFourB g u) :
    ¬ HasThreeDetunedGoodBranch deltaTwo deltaFourA deltaFourB g u := by
  obtain ⟨classTwo, classFourA, classFourB, -, -, -, hrowTwo,
    hrowFourA, hrowFourB, hparityTwoFourA, hparityTwoFourB,
    hdistinctFour⟩ := hpartition
  rintro ⟨c, hcIco, hcTwo, hcFourA, hcFourB⟩
  by_cases hparTwo : c ≡ classTwo [ZMOD 2]
  · apply hcTwo
    rw [hrowTwo, Finset.mem_filter]
    exact ⟨hcIco, hparTwo⟩
  · have hnotFourA : ¬ classFourA ≡ classTwo [ZMOD 2] :=
      fun h => hparityTwoFourA h.symm
    have hnotFourB : ¬ classFourB ≡ classTwo [ZMOD 2] :=
      fun h => hparityTwoFourB h.symm
    have hparFourA : c ≡ classFourA [ZMOD 2] :=
      modEq_two_of_not_modEq hparTwo hnotFourA
    have hparFourBA : classFourB ≡ classFourA [ZMOD 2] :=
      modEq_two_of_not_modEq hnotFourB hnotFourA
    rcases modEq_four_or_modEq_four_of_modEq_two
        hparFourA hparFourBA hdistinctFour with hmodFourA | hmodFourB
    · apply hcFourA
      rw [hrowFourA, Finset.mem_filter]
      exact ⟨hcIco, hmodFourA⟩
    · apply hcFourB
      rw [hrowFourB, Finset.mem_filter]
      exact ⟨hcIco, hmodFourB⟩

/-- Failure is equivalent to the exact three-class partition. -/
theorem not_hasThreeDetunedGoodBranch_two_four_four_iff_parallelPartition
    (deltaTwo deltaFourA deltaFourB g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hqTwo : g / (Int.gcd deltaTwo g : ℤ) = 2)
    (hqFourA : g / (Int.gcd deltaFourA g : ℤ) = 4)
    (hqFourB : g / (Int.gcd deltaFourB g : ℤ) = 4) :
    ¬ HasThreeDetunedGoodBranch deltaTwo deltaFourA deltaFourB g u ↔
      TwoFourFourParallelPartition deltaTwo deltaFourA deltaFourB g u := by
  constructor
  · intro hfail
    obtain ⟨classTwo, classFourA, classFourB, hcTwo, hcFourA, hcFourB,
      hrowTwo, hrowFourA, hrowFourB, hparityTwoFourA, hparityTwoFourB,
      hdistinctFour⟩ :=
        qTwo_four_four_failure_normal_form deltaTwo deltaFourA deltaFourB
          g u hg hqTwo hqFourA hqFourB hfail
    exact ⟨classTwo, classFourA, classFourB, hcTwo, hcFourA, hcFourB,
      hrowTwo, hrowFourA, hrowFourB, hparityTwoFourA, hparityTwoFourB,
      hdistinctFour⟩
  · exact not_hasThreeDetunedGoodBranch_of_twoFourFourParallelPartition
      deltaTwo deltaFourA deltaFourB g u

/-- Exact phase-avoidance consumer for the `(2,4,4)` selected witness. -/
theorem twoFourFour_selectedWitness_of_partition_avoidance
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g)
    (iTwo iFourA iFourB : Fin 13)
    (hqTwo : g / (Int.gcd (v iTwo) g : ℤ) = 2)
    (hqFourA : g / (Int.gcd (v iFourA) g : ℤ) = 4)
    (hqFourB : g / (Int.gcd (v iFourB) g : ℤ) = 4)
    (hselector : ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g iTwo iFourA iFourB u ∧
      ¬ TwoFourFourParallelPartition
        (v iTwo) (v iFourA) (v iFourB) g u) :
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g iTwo iFourA iFourB u ∧
      HasThreeDetunedGoodBranch (v iTwo) (v iFourA) (v iFourB) g u := by
  obtain ⟨u, hharmonic, hnotPartition⟩ := hselector
  refine ⟨u, hharmonic, ?_⟩
  by_contra hfail
  exact hnotPartition
    ((not_hasThreeDetunedGoodBranch_two_four_four_iff_parallelPartition
      (v iTwo) (v iFourA) (v iFourB) g u hg hqTwo hqFourA hqFourB).mp
      hfail)

/-! ## Axiom audit -/

#print axioms not_hasThreeDetunedGoodBranch_two_four_four_iff_parallelPartition
#print axioms twoFourFour_selectedWitness_of_partition_avoidance

end
end LRC14Grand
end LonelyRunner
