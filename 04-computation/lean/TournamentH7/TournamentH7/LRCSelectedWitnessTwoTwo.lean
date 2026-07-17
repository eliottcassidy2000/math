/-
  Exact parallel-class obstruction for the selected `(2,2,q)` witness.

  At qx >= 3 only opposition of the named q-two rows is fatal; at qx = 2 any
  opposite pair is fatal.  These are saturated `{1/2,1/2}` prefix codes and
  equality cases of the fixed-phase incidence bound.  Static row incidence
  preserves coverage but loses the cyclic phase moment, so a manufactured
  tie-tournament adds no information.  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCSelectedWitnessCommon

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

noncomputable section

/-- If the third row has reduced denominator at least three, its universal
bad degree is strictly below one half of the branch circle. -/
theorem two_mul_badCount_lt_of_reducedDenominator_ge_three
    (δ g : ℤ) (hg : 2 ≤ g)
    (hq3 : 3 ≤ g / (Int.gcd δ g : ℤ)) :
    2 * DetunedD3.badCount δ g < g.toNat := by
  have hg0 : (0 : ℤ) < g := by omega
  by_cases heq3 : g / (Int.gcd δ g : ℤ) = 3
  · have hbad := three_mul_badCount_eq δ g hg0 heq3
    omega
  · have hq4 : 4 ≤ g / (Int.gcd δ g : ℤ) := by omega
    have hbad := seven_mul_badCount_le_two_mul δ g hg0 hq4
    omega

/-- Exact fixed-phase TwoTwo closure when the third denominator is at least
three. Away from opposite q-two parity rows, their overlap/emptiness plus the
strict half-bound for the third row leaves a common branch. -/
theorem hasThreeDetunedGoodBranch_two_two_qgeThree_of_not_opposition
    (δ₂a δ₂b δₓ g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hq2a : g / (Int.gcd δ₂a g : ℤ) = 2)
    (hq2b : g / (Int.gcd δ₂b g : ℤ) = 2)
    (hqx3 : 3 ≤ g / (Int.gcd δₓ g : ℤ))
    (hnot : ¬ TwoTwoThreePhaseOpposition δ₂a δ₂b g u) :
    HasThreeDetunedGoodBranch δ₂a δ₂b δₓ g u := by
  have hg0 : (0 : ℤ) < g := by omega
  have hg1 : (1 : ℤ) ≤ g := by omega
  have hbad2a := two_mul_badCount_eq δ₂a g hg0 hq2a
  have hbad2b := two_mul_badCount_eq δ₂b g hg0 hq2b
  have hsmall := two_mul_badCount_lt_of_reducedDenominator_ge_three
    δₓ g hg hqx3
  have hcard2a := detunedBadBranches_card_le δ₂a g hg1 u
  have hcard2b := detunedBadBranches_card_le δ₂b g hg1 u
  have hcardx := detunedBadBranches_card_le δₓ g hg1 u
  by_cases hrow2a : (detunedBadBranches δ₂a g u).Nonempty
  · by_cases hrow2b : (detunedBadBranches δ₂b g u).Nonempty
    · by_cases hoverlap :
        (detunedBadBranches δ₂a g u ∩
          detunedBadBranches δ₂b g u).Nonempty
      · have hroweq := detunedBadBranches_eq_of_overlap_same_reducedDenominator
          δ₂a δ₂b g 2 u hg1 (by norm_num) hq2a hq2b hoverlap
        have hinter : detunedBadBranches δ₂a g u ∩
            detunedBadBranches δ₂b g u = detunedBadBranches δ₂a g u := by
          rw [hroweq]
          simp
        have hintercard := two_mul_detunedBadBranches_card_eq_of_nonempty
          δ₂a g u hg1 hq2a hrow2a
        apply hasThreeDetunedGoodBranch_of_overlapDebt
          δ₂a δ₂b δₓ g u hg1
        unfold ThreeDetunedOverlapDebtPaid
        left
        rw [hinter]
        omega
      · exfalso
        apply hnot
        refine ⟨hrow2a, hrow2b, Finset.disjoint_left.mpr ?_⟩
        intro c hc2a hc2b
        exact hoverlap ⟨c, Finset.mem_inter.mpr ⟨hc2a, hc2b⟩⟩
    · apply hasThreeDetunedGoodBranch_of_card_sum_lt δ₂a δ₂b δₓ g u
      have hzero : (detunedBadBranches δ₂b g u).card = 0 := by
        rw [Finset.card_eq_zero]
        exact Finset.not_nonempty_iff_eq_empty.mp hrow2b
      omega
  · apply hasThreeDetunedGoodBranch_of_card_sum_lt δ₂a δ₂b δₓ g u
    have hzero : (detunedBadBranches δ₂a g u).card = 0 := by
      rw [Finset.card_eq_zero]
      exact Finset.not_nonempty_iff_eq_empty.mp hrow2a
    omega

/-- For qx at least three, the elementary Zarankiewicz extremizer is exact:
failure occurs iff the two q-two rows are the two disjoint parity classes. -/
theorem not_hasThreeDetunedGoodBranch_two_two_qgeThree_iff_opposition
    (δ₂a δ₂b δₓ g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hq2a : g / (Int.gcd δ₂a g : ℤ) = 2)
    (hq2b : g / (Int.gcd δ₂b g : ℤ) = 2)
    (hqx3 : 3 ≤ g / (Int.gcd δₓ g : ℤ)) :
    ¬ HasThreeDetunedGoodBranch δ₂a δ₂b δₓ g u ↔
      TwoTwoThreePhaseOpposition δ₂a δ₂b g u := by
  constructor
  · intro hfail
    by_contra hnot
    exact hfail
      (hasThreeDetunedGoodBranch_two_two_qgeThree_of_not_opposition
        δ₂a δ₂b δₓ g u hg hq2a hq2b hqx3 hnot)
  · exact not_hasThreeDetunedGoodBranch_two_two_three_of_opposition
      δ₂a δ₂b δₓ g u (by omega) hq2a hq2b

/-- The exact extra selector for the qx>=3 TwoTwo residual. The citation
supplies harmonic-good phases, but one of them must avoid parity opposition. -/
theorem twoTwo_selectedWitness_of_nonopposition_selector
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g)
    (i₂a i₂b iₓ : Fin 13)
    (hq2a : g / (Int.gcd (v i₂a) g : ℤ) = 2)
    (hq2b : g / (Int.gcd (v i₂b) g : ℤ) = 2)
    (hqx3 : 3 ≤ g / (Int.gcd (v iₓ) g : ℤ))
    (hselector : ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₂a i₂b iₓ u ∧
      ¬ TwoTwoThreePhaseOpposition (v i₂a) (v i₂b) g u) :
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₂a i₂b iₓ u ∧
      HasThreeDetunedGoodBranch (v i₂a) (v i₂b) (v iₓ) g u := by
  obtain ⟨u, hharm, hnot⟩ := hselector
  exact ⟨u, hharm,
    hasThreeDetunedGoodBranch_two_two_qgeThree_of_not_opposition
      (v i₂a) (v i₂b) (v iₓ) g u hg hq2a hq2b hqx3 hnot⟩

/-- In the all-q-two regime, obstruction means that some pair of nonempty
rows occupies the two opposite parity classes. -/
def ThreeQTwoPairOpposition
    (δ₂a δ₂b δ₂c g : ℤ) (u : ℝ) : Prop :=
  TwoTwoThreePhaseOpposition δ₂a δ₂b g u ∨
  TwoTwoThreePhaseOpposition δ₂a δ₂c g u ∨
  TwoTwoThreePhaseOpposition δ₂b δ₂c g u

theorem qTwo_inter_nonempty_of_nonempty_of_not_opposition
    (δ₂a δ₂b g : ℤ) (u : ℝ)
    (hrowa : (detunedBadBranches δ₂a g u).Nonempty)
    (hrowb : (detunedBadBranches δ₂b g u).Nonempty)
    (hnot : ¬ TwoTwoThreePhaseOpposition δ₂a δ₂b g u) :
    (detunedBadBranches δ₂a g u ∩
      detunedBadBranches δ₂b g u).Nonempty := by
  by_contra hnointer
  apply hnot
  refine ⟨hrowa, hrowb, Finset.disjoint_left.mpr ?_⟩
  intro c hca hcb
  exact hnointer ⟨c, Finset.mem_inter.mpr ⟨hca, hcb⟩⟩

theorem two_mul_qTwo_inter_card_eq_of_overlap
    (δ₂a δ₂b g : ℤ) (u : ℝ) (hg : 1 ≤ g)
    (hq2a : g / (Int.gcd δ₂a g : ℤ) = 2)
    (hq2b : g / (Int.gcd δ₂b g : ℤ) = 2)
    (hrowa : (detunedBadBranches δ₂a g u).Nonempty)
    (hoverlap : (detunedBadBranches δ₂a g u ∩
      detunedBadBranches δ₂b g u).Nonempty) :
    2 * (detunedBadBranches δ₂a g u ∩
      detunedBadBranches δ₂b g u).card = g.toNat := by
  have hroweq := detunedBadBranches_eq_of_overlap_same_reducedDenominator
    δ₂a δ₂b g 2 u hg (by norm_num) hq2a hq2b hoverlap
  have hinter : detunedBadBranches δ₂a g u ∩
      detunedBadBranches δ₂b g u = detunedBadBranches δ₂a g u := by
    rw [hroweq]
    simp
  rw [hinter]
  exact two_mul_detunedBadBranches_card_eq_of_nonempty
    δ₂a g u hg hq2a hrowa

theorem hasThreeDetunedGoodBranch_three_qTwo_of_no_pair_opposition
    (δ₂a δ₂b δ₂c g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hq2a : g / (Int.gcd δ₂a g : ℤ) = 2)
    (hq2b : g / (Int.gcd δ₂b g : ℤ) = 2)
    (hq2c : g / (Int.gcd δ₂c g : ℤ) = 2)
    (hnot : ¬ ThreeQTwoPairOpposition δ₂a δ₂b δ₂c g u) :
    HasThreeDetunedGoodBranch δ₂a δ₂b δ₂c g u := by
  have hg0 : (0 : ℤ) < g := by omega
  have hg1 : (1 : ℤ) ≤ g := by omega
  have hbad2a := two_mul_badCount_eq δ₂a g hg0 hq2a
  have hbad2b := two_mul_badCount_eq δ₂b g hg0 hq2b
  have hbad2c := two_mul_badCount_eq δ₂c g hg0 hq2c
  have hcarda := detunedBadBranches_card_le δ₂a g hg1 u
  have hcardb := detunedBadBranches_card_le δ₂b g hg1 u
  have hcardc := detunedBadBranches_card_le δ₂c g hg1 u
  have hbranches : (Finset.Ico (0 : ℤ) g).card = g.toNat := by
    rw [Int.card_Ico]
    congr 1
    omega
  simp only [ThreeQTwoPairOpposition, not_or] at hnot
  rcases hnot with ⟨hnotab, hnotac, hnotbc⟩
  by_cases hrowa : (detunedBadBranches δ₂a g u).Nonempty
  · by_cases hrowb : (detunedBadBranches δ₂b g u).Nonempty
    · by_cases hrowc : (detunedBadBranches δ₂c g u).Nonempty
      · have hab := qTwo_inter_nonempty_of_nonempty_of_not_opposition
          δ₂a δ₂b g u hrowa hrowb hnotab
        have hac := qTwo_inter_nonempty_of_nonempty_of_not_opposition
          δ₂a δ₂c g u hrowa hrowc hnotac
        have hcardab := two_mul_qTwo_inter_card_eq_of_overlap
          δ₂a δ₂b g u hg1 hq2a hq2b hrowa hab
        have hcardac := two_mul_qTwo_inter_card_eq_of_overlap
          δ₂a δ₂c g u hg1 hq2a hq2c hrowa hac
        apply exists_outside_three_of_two_first_overlapDebt
          (Finset.Ico (0 : ℤ) g)
          (detunedBadBranches δ₂a g u)
          (detunedBadBranches δ₂b g u)
          (detunedBadBranches δ₂c g u)
        rw [hbranches]
        omega
      · apply hasThreeDetunedGoodBranch_of_exactOverlapDebt
          δ₂a δ₂b δ₂c g u
        unfold ThreeDetunedExactOverlapDebtPaid
        left
        have hab := qTwo_inter_nonempty_of_nonempty_of_not_opposition
          δ₂a δ₂b g u hrowa hrowb hnotab
        have hcardab := two_mul_qTwo_inter_card_eq_of_overlap
          δ₂a δ₂b g u hg1 hq2a hq2b hrowa hab
        have hzero : (detunedBadBranches δ₂c g u).card = 0 := by
          rw [Finset.card_eq_zero]
          exact Finset.not_nonempty_iff_eq_empty.mp hrowc
        omega
    · by_cases hrowc : (detunedBadBranches δ₂c g u).Nonempty
      · apply hasThreeDetunedGoodBranch_of_exactOverlapDebt
          δ₂a δ₂b δ₂c g u
        unfold ThreeDetunedExactOverlapDebtPaid
        right
        left
        have hac := qTwo_inter_nonempty_of_nonempty_of_not_opposition
          δ₂a δ₂c g u hrowa hrowc hnotac
        have hcardac := two_mul_qTwo_inter_card_eq_of_overlap
          δ₂a δ₂c g u hg1 hq2a hq2c hrowa hac
        have hzero : (detunedBadBranches δ₂b g u).card = 0 := by
          rw [Finset.card_eq_zero]
          exact Finset.not_nonempty_iff_eq_empty.mp hrowb
        omega
      · apply hasThreeDetunedGoodBranch_of_card_sum_lt
          δ₂a δ₂b δ₂c g u
        have hzerob : (detunedBadBranches δ₂b g u).card = 0 := by
          rw [Finset.card_eq_zero]
          exact Finset.not_nonempty_iff_eq_empty.mp hrowb
        have hzeroc : (detunedBadBranches δ₂c g u).card = 0 := by
          rw [Finset.card_eq_zero]
          exact Finset.not_nonempty_iff_eq_empty.mp hrowc
        omega
  · by_cases hrowb : (detunedBadBranches δ₂b g u).Nonempty
    · by_cases hrowc : (detunedBadBranches δ₂c g u).Nonempty
      · apply hasThreeDetunedGoodBranch_of_exactOverlapDebt
          δ₂a δ₂b δ₂c g u
        unfold ThreeDetunedExactOverlapDebtPaid
        right
        right
        have hbc := qTwo_inter_nonempty_of_nonempty_of_not_opposition
          δ₂b δ₂c g u hrowb hrowc hnotbc
        have hcardbc := two_mul_qTwo_inter_card_eq_of_overlap
          δ₂b δ₂c g u hg1 hq2b hq2c hrowb hbc
        have hzero : (detunedBadBranches δ₂a g u).card = 0 := by
          rw [Finset.card_eq_zero]
          exact Finset.not_nonempty_iff_eq_empty.mp hrowa
        omega
      · apply hasThreeDetunedGoodBranch_of_card_sum_lt
          δ₂a δ₂b δ₂c g u
        have hzeroa : (detunedBadBranches δ₂a g u).card = 0 := by
          rw [Finset.card_eq_zero]
          exact Finset.not_nonempty_iff_eq_empty.mp hrowa
        have hzeroc : (detunedBadBranches δ₂c g u).card = 0 := by
          rw [Finset.card_eq_zero]
          exact Finset.not_nonempty_iff_eq_empty.mp hrowc
        omega
    · apply hasThreeDetunedGoodBranch_of_card_sum_lt
        δ₂a δ₂b δ₂c g u
      have hzeroa : (detunedBadBranches δ₂a g u).card = 0 := by
        rw [Finset.card_eq_zero]
        exact Finset.not_nonempty_iff_eq_empty.mp hrowa
      have hzerob : (detunedBadBranches δ₂b g u).card = 0 := by
        rw [Finset.card_eq_zero]
        exact Finset.not_nonempty_iff_eq_empty.mp hrowb
      omega

/-- Exact all-q-two obstruction: failure iff at least one pair occupies
opposite parity classes and already covers the circle. -/
theorem not_hasThreeDetunedGoodBranch_three_qTwo_iff_pair_opposition
    (δ₂a δ₂b δ₂c g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hq2a : g / (Int.gcd δ₂a g : ℤ) = 2)
    (hq2b : g / (Int.gcd δ₂b g : ℤ) = 2)
    (hq2c : g / (Int.gcd δ₂c g : ℤ) = 2) :
    ¬ HasThreeDetunedGoodBranch δ₂a δ₂b δ₂c g u ↔
      ThreeQTwoPairOpposition δ₂a δ₂b δ₂c g u := by
  constructor
  · intro hfail
    by_contra hnot
    exact hfail
      (hasThreeDetunedGoodBranch_three_qTwo_of_no_pair_opposition
        δ₂a δ₂b δ₂c g u hg hq2a hq2b hq2c hnot)
  · rintro (hab | hac | hbc) hgood
    · exact (not_hasThreeDetunedGoodBranch_two_two_three_of_opposition
        δ₂a δ₂b δ₂c g u (by omega) hq2a hq2b hab) hgood
    · apply (not_hasThreeDetunedGoodBranch_two_two_three_of_opposition
        δ₂a δ₂c δ₂b g u (by omega) hq2a hq2c hac)
      obtain ⟨c, hcIco, hca, hcb, hcc⟩ := hgood
      exact ⟨c, hcIco, hca, hcc, hcb⟩
    · apply (not_hasThreeDetunedGoodBranch_two_two_three_of_opposition
        δ₂b δ₂c δ₂a g u (by omega) hq2b hq2c hbc)
      obtain ⟨c, hcIco, hca, hcb, hcc⟩ := hgood
      exact ⟨c, hcIco, hcb, hcc, hca⟩

/-- Exhaustive parallel-circle obstruction for the exact `(2,2,q)` pattern.
If qx=2, any opposite q-two pair is fatal; if qx>=3, only opposition of the
two selected q-two rows can saturate. -/
def TwoTwoParallelObstruction
    (δ₂a δ₂b δₓ g : ℤ) (u : ℝ) : Prop :=
  (g / (Int.gcd δₓ g : ℤ) = 2 ∧
    ThreeQTwoPairOpposition δ₂a δ₂b δₓ g u) ∨
  (3 ≤ g / (Int.gcd δₓ g : ℤ) ∧
    TwoTwoThreePhaseOpposition δ₂a δ₂b g u)

theorem not_hasThreeDetunedGoodBranch_two_two_iff_parallelObstruction
    (δ₂a δ₂b δₓ g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hδx : ¬ g ∣ δₓ)
    (hq2a : g / (Int.gcd δ₂a g : ℤ) = 2)
    (hq2b : g / (Int.gcd δ₂b g : ℤ) = 2) :
    ¬ HasThreeDetunedGoodBranch δ₂a δ₂b δₓ g u ↔
      TwoTwoParallelObstruction δ₂a δ₂b δₓ g u := by
  have hqx2 := reducedDetuningDenominator_ge_two hg hδx
  by_cases hqx : g / (Int.gcd δₓ g : ℤ) = 2
  · constructor
    · intro hfail
      exact Or.inl ⟨hqx,
        (not_hasThreeDetunedGoodBranch_three_qTwo_iff_pair_opposition
          δ₂a δ₂b δₓ g u hg hq2a hq2b hqx).mp hfail⟩
    · rintro (⟨-, hopposition⟩ | ⟨hqx3, -⟩)
      · exact (not_hasThreeDetunedGoodBranch_three_qTwo_iff_pair_opposition
          δ₂a δ₂b δₓ g u hg hq2a hq2b hqx).mpr hopposition
      · omega
  · have hqx3 : 3 ≤ g / (Int.gcd δₓ g : ℤ) := by omega
    constructor
    · intro hfail
      exact Or.inr ⟨hqx3,
        (not_hasThreeDetunedGoodBranch_two_two_qgeThree_iff_opposition
          δ₂a δ₂b δₓ g u hg hq2a hq2b hqx3).mp hfail⟩
    · rintro (⟨hqx', -⟩ | ⟨-, hopposition⟩)
      · exact (hqx hqx').elim
      · exact (not_hasThreeDetunedGoodBranch_two_two_qgeThree_iff_opposition
          δ₂a δ₂b δₓ g u hg hq2a hq2b hqx3).mpr hopposition

/-- With the exact nonmultiplicity hypothesis on the third selected row, no
further denominator assumption is needed. The sole remaining supplier is a
harmonic-good phase avoiding `TwoTwoParallelObstruction`. -/
theorem twoTwo_selectedWitness_of_parallelObstruction_avoidance
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g)
    (i₂a i₂b iₓ : Fin 13) (hδx : ¬ g ∣ v iₓ)
    (hq2a : g / (Int.gcd (v i₂a) g : ℤ) = 2)
    (hq2b : g / (Int.gcd (v i₂b) g : ℤ) = 2)
    (hselector : ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₂a i₂b iₓ u ∧
      ¬ TwoTwoParallelObstruction (v i₂a) (v i₂b) (v iₓ) g u) :
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₂a i₂b iₓ u ∧
      HasThreeDetunedGoodBranch (v i₂a) (v i₂b) (v iₓ) g u := by
  obtain ⟨u, hharm, hnotObstruction⟩ := hselector
  refine ⟨u, hharm, ?_⟩
  by_contra hfail
  exact hnotObstruction
    ((not_hasThreeDetunedGoodBranch_two_two_iff_parallelObstruction
      (v i₂a) (v i₂b) (v iₓ) g u hg hδx hq2a hq2b).mp hfail)

/-! ## Axiom audit -/

#print axioms not_hasThreeDetunedGoodBranch_two_two_qgeThree_iff_opposition
#print axioms not_hasThreeDetunedGoodBranch_three_qTwo_iff_pair_opposition
#print axioms not_hasThreeDetunedGoodBranch_two_two_iff_parallelObstruction
#print axioms twoTwo_selectedWitness_of_parallelObstruction_avoidance

end

end LRC14Grand
end LonelyRunner
