/-
  Exact compatibility among the two-fresh q22 and q244 obstruction walls.
  Two q244 failures sharing the old q4 pair force equal fresh half-rows, so
  they cannot coexist with q22 opposition. The wall-type nerve has no triple
  face; separate wall widths or a Helly/common-wall quotient lose this fact.
  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCPairTowerReduction

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

/-- Two failing q244 partitions with the same old q4 pair have the same
fresh q2 half-row. -/
theorem qTwo_badRows_eq_of_two_q244_failures
    (δ₂a δ₂b δ₄a δ₄b g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hq2a : g / (Int.gcd δ₂a g : ℤ) = 2)
    (hq2b : g / (Int.gcd δ₂b g : ℤ) = 2)
    (hq4a : g / (Int.gcd δ₄a g : ℤ) = 4)
    (hq4b : g / (Int.gcd δ₄b g : ℤ) = 4)
    (hfailA : ¬ HasThreeDetunedGoodBranch δ₂a δ₄a δ₄b g u)
    (hfailB : ¬ HasThreeDetunedGoodBranch δ₂b δ₄a δ₄b g u) :
    detunedBadBranches δ₂a g u = detunedBadBranches δ₂b g u := by
  have hg0 : (0 : ℤ) < g := by omega
  have hg1 : (1 : ℤ) ≤ g := by omega
  have hbad2a := two_mul_badCount_eq δ₂a g hg0 hq2a
  have hbad2b := two_mul_badCount_eq δ₂b g hg0 hq2b
  have hbad4a :=
    four_mul_badCount_eq_of_reducedDenominator_four δ₄a g hg0 hq4a
  have hbad4b :=
    four_mul_badCount_eq_of_reducedDenominator_four δ₄b g hg0 hq4b
  have hbudgetA :
      DetunedD3.badCount δ₂a g + DetunedD3.badCount δ₄a g +
          DetunedD3.badCount δ₄b g ≤ g.toNat := by
    omega
  have hbudgetB :
      DetunedD3.badCount δ₂b g + DetunedD3.badCount δ₄a g +
          DetunedD3.badCount δ₄b g ≤ g.toNat := by
    omega
  have hAOld : Disjoint (detunedBadBranches δ₂a g u)
      (detunedBadBranches δ₄a g u ∪ detunedBadBranches δ₄b g u) := by
    rw [Finset.disjoint_left]
    intro branch hbranchA hbranchOld
    rcases Finset.mem_union.mp hbranchOld with hbranch4a | hbranch4b
    · have hinter :
          (detunedBadBranches δ₂a g u ∩
            detunedBadBranches δ₄a g u).Nonempty :=
        ⟨branch, Finset.mem_inter.mpr ⟨hbranchA, hbranch4a⟩⟩
      exact hfailA (hasThreeDetunedGoodBranch_of_pairOverlap
        δ₂a δ₄a δ₄b g u hg1 hbudgetA (Or.inl hinter))
    · have hinter :
          (detunedBadBranches δ₂a g u ∩
            detunedBadBranches δ₄b g u).Nonempty :=
        ⟨branch, Finset.mem_inter.mpr ⟨hbranchA, hbranch4b⟩⟩
      exact hfailA (hasThreeDetunedGoodBranch_of_pairOverlap
        δ₂a δ₄a δ₄b g u hg1 hbudgetA (Or.inr (Or.inl hinter)))
  have hBOld : Disjoint (detunedBadBranches δ₂b g u)
      (detunedBadBranches δ₄a g u ∪ detunedBadBranches δ₄b g u) := by
    rw [Finset.disjoint_left]
    intro branch hbranchB hbranchOld
    rcases Finset.mem_union.mp hbranchOld with hbranch4a | hbranch4b
    · have hinter :
          (detunedBadBranches δ₂b g u ∩
            detunedBadBranches δ₄a g u).Nonempty :=
        ⟨branch, Finset.mem_inter.mpr ⟨hbranchB, hbranch4a⟩⟩
      exact hfailB (hasThreeDetunedGoodBranch_of_pairOverlap
        δ₂b δ₄a δ₄b g u hg1 hbudgetB (Or.inl hinter))
    · have hinter :
          (detunedBadBranches δ₂b g u ∩
            detunedBadBranches δ₄b g u).Nonempty :=
        ⟨branch, Finset.mem_inter.mpr ⟨hbranchB, hbranch4b⟩⟩
      exact hfailB (hasThreeDetunedGoodBranch_of_pairOverlap
        δ₂b δ₄a δ₄b g u hg1 hbudgetB (Or.inr (Or.inl hinter)))
  have hcoverA : Finset.Ico (0 : ℤ) g ⊆
      detunedBadBranches δ₂a g u ∪
        (detunedBadBranches δ₄a g u ∪ detunedBadBranches δ₄b g u) := by
    intro branch hbranch
    by_contra hnot
    simp only [Finset.mem_union, not_or] at hnot
    exact hfailA ⟨branch, hbranch, hnot.1, hnot.2.1, hnot.2.2⟩
  have hcoverB : Finset.Ico (0 : ℤ) g ⊆
      detunedBadBranches δ₂b g u ∪
        (detunedBadBranches δ₄a g u ∪ detunedBadBranches δ₄b g u) := by
    intro branch hbranch
    by_contra hnot
    simp only [Finset.mem_union, not_or] at hnot
    exact hfailB ⟨branch, hbranch, hnot.1, hnot.2.1, hnot.2.2⟩
  ext branch
  constructor
  · intro hbranchA
    have hbranchIco := detunedBadBranches_subset_Ico δ₂a g u hbranchA
    rcases Finset.mem_union.mp (hcoverB hbranchIco) with hbranchB | hbranchOld
    · exact hbranchB
    · exact (Finset.disjoint_left.mp hAOld hbranchA hbranchOld).elim
  · intro hbranchB
    have hbranchIco := detunedBadBranches_subset_Ico δ₂b g u hbranchB
    rcases Finset.mem_union.mp (hcoverA hbranchIco) with hbranchA | hbranchOld
    · exact hbranchA
    · exact (Finset.disjoint_left.mp hBOld hbranchB hbranchOld).elim

/-- Consequently simultaneous q244 failures exclude q22 opposition between
the two fresh rows. -/
theorem not_twoTwoOpposition_of_two_q244_failures
    (δ₂a δ₂b δ₄a δ₄b g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hq2a : g / (Int.gcd δ₂a g : ℤ) = 2)
    (hq2b : g / (Int.gcd δ₂b g : ℤ) = 2)
    (hq4a : g / (Int.gcd δ₄a g : ℤ) = 4)
    (hq4b : g / (Int.gcd δ₄b g : ℤ) = 4)
    (hfailA : ¬ HasThreeDetunedGoodBranch δ₂a δ₄a δ₄b g u)
    (hfailB : ¬ HasThreeDetunedGoodBranch δ₂b δ₄a δ₄b g u) :
    ¬ TwoTwoThreePhaseOpposition δ₂a δ₂b g u := by
  intro hopposition
  have hrows := qTwo_badRows_eq_of_two_q244_failures
    δ₂a δ₂b δ₄a δ₄b g u hg hq2a hq2b hq4a hq4b hfailA hfailB
  obtain ⟨branch, hbranchA⟩ := hopposition.1
  have hbranchB : branch ∈ detunedBadBranches δ₂b g u := by
    rw [← hrows]
    exact hbranchA
  exact Finset.disjoint_left.mp hopposition.2.2 hbranchA hbranchB

/-- On a q22 opposition edge, failure of one endpoint against a fixed old q4
pair forces the other endpoint to have a three-row good branch.  Equivalently,
the q244-failure vertices form an independent set in the q22 opposition graph
at every fixed phase. -/
theorem hasThreeDetunedGoodBranch_of_twoTwoOpposition_and_q244_failure
    (δ₂a δ₂b δ₄a δ₄b g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hq2a : g / (Int.gcd δ₂a g : ℤ) = 2)
    (hq2b : g / (Int.gcd δ₂b g : ℤ) = 2)
    (hq4a : g / (Int.gcd δ₄a g : ℤ) = 4)
    (hq4b : g / (Int.gcd δ₄b g : ℤ) = 4)
    (hopposition : TwoTwoThreePhaseOpposition δ₂a δ₂b g u)
    (hfailA : ¬ HasThreeDetunedGoodBranch δ₂a δ₄a δ₄b g u) :
    HasThreeDetunedGoodBranch δ₂b δ₄a δ₄b g u := by
  by_contra hfailB
  exact (not_twoTwoOpposition_of_two_q244_failures
    δ₂a δ₂b δ₄a δ₄b g u hg hq2a hq2b hq4a hq4b hfailA hfailB)
      hopposition

/-- Pair-tower specialization of the independent-failure law.  Fresh rows
crossing the first divisibility wall have reduced denominator two at `2g`,
while the two inherited base-q2 rows become the fixed q4 pair. -/
theorem pairTower_q244_escape_on_opposite_fresh_pair
    (v : Fin 13 → ℤ) (g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (oldA oldB freshA freshB : Fin 13)
    (hqOldA : g / (Int.gcd (v oldA) g : ℤ) = 2)
    (hqOldB : g / (Int.gcd (v oldB) g : ℤ) = 2)
    (hfreshA : g ∣ v freshA ∧ ¬ 2 * g ∣ v freshA)
    (hfreshB : g ∣ v freshB ∧ ¬ 2 * g ∣ v freshB)
    (hopposition :
      TwoTwoThreePhaseOpposition (v freshA) (v freshB) (2 * g) u)
    (hfailA :
      ¬ HasThreeDetunedGoodBranch
        (v freshA) (v oldA) (v oldB) (2 * g) u) :
    HasThreeDetunedGoodBranch
      (v freshB) (v oldA) (v oldB) (2 * g) u := by
  have hqFreshA := odd_harmonic_lifts_to_q_two
    (v freshA) g hg hfreshA.1 hfreshA.2
  have hqFreshB := odd_harmonic_lifts_to_q_two
    (v freshB) g hg hfreshB.1 hfreshB.2
  have hqOldA' := reducedDenominator_double_eq_four_of_eq_two
    (v oldA) g (by omega) hqOldA
  have hqOldB' := reducedDenominator_double_eq_four_of_eq_two
    (v oldB) g (by omega) hqOldB
  exact hasThreeDetunedGoodBranch_of_twoTwoOpposition_and_q244_failure
    (v freshA) (v freshB) (v oldA) (v oldB) (2 * g) u (by omega)
      hqFreshA hqFreshB hqOldA' hqOldB' hopposition hfailA

#print axioms qTwo_badRows_eq_of_two_q244_failures
#print axioms not_twoTwoOpposition_of_two_q244_failures
#print axioms hasThreeDetunedGoodBranch_of_twoTwoOpposition_and_q244_failure
#print axioms pairTower_q244_escape_on_opposite_fresh_pair

end LRC14Grand
end LonelyRunner
