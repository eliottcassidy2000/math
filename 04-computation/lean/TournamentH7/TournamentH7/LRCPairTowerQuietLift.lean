/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: codex (LRC multi-agent project, 2026-07-17-S47)
-/
import TournamentH7.LRCPairTowerValuation

/-!
# The quiet two-row lift is the proved valuation-gap-three leaf

Suppose a nongeneric pair at modulus `g` exposes exactly two fresh rows at
the first doubling.  If the divisible core is quiet at the next two walls,
then the four exceptional rows have reduced denominators `(8,8,16,16)` at
modulus `8g`.  This module proves that assertion from the actual pair-tower
predicates and feeds it to `lonely14_of_pairTower_two_min_gap_three`.

Thus the remaining many-lift supplier can exclude the concrete condition

`liftFailureCard v g = 2`, `TwoAdicLiftTerminates v (2*g)`, and
`TwoAdicLiftTerminates v (4*g)`.

Tournament-analysis audit: vertices here are divisibility-wall events.  The
observable is the first wall crossed by a speed; it gives a valuation-layer
partition, not a binary orientation.  A tournament gauge would forget the
common modulus needed by the four-row debt estimate.  The quotient preserves
the exact two fresh rows and the three-wall quiet tail, while discarding their
phase positions.  The challenged assumption is that two first-wall events
always remain a hard pair-tower branch: two subsequent empty walls force the
already closed `(8,8,16,16)` leaf.

No `sorry`; no `native_decide`.
-/

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

noncomputable section

/-! ## Odd quotients are coprime to the required dyadic powers -/

private theorem gcd_eight_eq_one_of_not_two_dvd (a : ℤ)
    (ha : ¬ (2 : ℤ) ∣ a) :
    Int.gcd a 8 = 1 := by
  have haNat : ¬ (2 : ℕ) ∣ a.natAbs := by
    intro hdiv
    apply ha
    exact_mod_cast hdiv
  have hcop := Nat.prime_two.coprime_pow_of_not_dvd (m := 3) haNat
  rw [Nat.coprime_iff_gcd_eq_one] at hcop
  simpa [Int.gcd] using hcop

private theorem gcd_sixteen_eq_one_of_not_two_dvd (a : ℤ)
    (ha : ¬ (2 : ℤ) ∣ a) :
    Int.gcd a 16 = 1 := by
  have haNat : ¬ (2 : ℕ) ∣ a.natAbs := by
    intro hdiv
    apply ha
    exact_mod_cast hdiv
  have hcop := Nat.prime_two.coprime_pow_of_not_dvd (m := 4) haNat
  rw [Nat.coprime_iff_gcd_eq_one] at hcop
  simpa [Int.gcd] using hcop

/-! ## Extraction of the genuine quiet leaf -/

/-- The genuine strict pair-tower leaf.  Two original q-two rows and exactly
two first-wall rows are extracted from the cardinality predicates.  Quietness
at the following two walls forces all remaining rows to be divisible by
`8g`; after writing `g = 2d`, this is exactly the factorized
valuation-gap-three theorem. -/
theorem lonely14_of_pairTower_quietLift
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (hg : 2 ≤ g)
    (hcard : nonMultCard v g = 2)
    (hnongeneric : ¬ genericCount v g)
    (hfreshCard : liftFailureCard v g = 2)
    (hquietTwo : TwoAdicLiftTerminates v (2 * g))
    (hquietFour : TwoAdicLiftTerminates v (4 * g)) :
    ∃ t : ℝ, Lonely 14 v t := by
  obtain ⟨j₁, j₂, hj₁j₂, hdvdBase, hj₁Not, hj₂Not, hqj₁, hqj₂⟩ :=
    nonGeneric_pair_q_two v g hg hcard hnongeneric

  rw [liftFailureCard] at hfreshCard
  obtain ⟨i₁, i₂, hi₁i₂, hfresh⟩ := Finset.card_eq_two.mp hfreshCard
  have hi₁Wall : g ∣ v i₁ ∧ ¬ 2 * g ∣ v i₁ := by
    have : i₁ ∈ Finset.univ.filter (fun i => g ∣ v i ∧ ¬ 2 * g ∣ v i) := by
      rw [hfresh]
      simp
    simpa using this
  have hi₂Wall : g ∣ v i₂ ∧ ¬ 2 * g ∣ v i₂ := by
    have : i₂ ∈ Finset.univ.filter (fun i => g ∣ v i ∧ ¬ 2 * g ∣ v i) := by
      rw [hfresh]
      simp
    simpa using this

  have hi₁j₁ : i₁ ≠ j₁ := by
    intro h
    subst i₁
    exact hj₁Not hi₁Wall.1
  have hi₁j₂ : i₁ ≠ j₂ := by
    intro h
    subst i₁
    exact hj₂Not hi₁Wall.1
  have hi₂j₁ : i₂ ≠ j₁ := by
    intro h
    subst i₂
    exact hj₁Not hi₂Wall.1
  have hi₂j₂ : i₂ ≠ j₂ := by
    intro h
    subst i₂
    exact hj₂Not hi₂Wall.1

  let d : ℤ := (Int.gcd (v j₁) g : ℤ)
  have hgEq : g = 2 * d := by
    have hfactor := Int.mul_ediv_cancel' (Int.gcd_dvd_right (v j₁) g)
    rw [hqj₁] at hfactor
    simpa [d, mul_comm] using hfactor.symm
  have hd : 1 ≤ d := by omega

  have hgcdj₂Eq : (Int.gcd (v j₂) g : ℤ) = d := by
    have hfactor := Int.mul_ediv_cancel' (Int.gcd_dvd_right (v j₂) g)
    rw [hqj₂] at hfactor
    have hgEq' : g = (Int.gcd (v j₂) g : ℤ) * 2 := hfactor.symm
    rw [hgEq] at hgEq'
    omega

  let m₁ : ℤ := v i₁ / g
  let m₂ : ℤ := v i₂ / g
  let a₁ : ℤ := v j₁ / d
  let a₂ : ℤ := v j₂ / d

  have hi₁FactorG : v i₁ = g * m₁ :=
    (Int.mul_ediv_cancel' hi₁Wall.1).symm
  have hi₂FactorG : v i₂ = g * m₂ :=
    (Int.mul_ediv_cancel' hi₂Wall.1).symm
  have hi₁Factor : v i₁ = (2 * d) * m₁ := by
    simpa [hgEq] using hi₁FactorG
  have hi₂Factor : v i₂ = (2 * d) * m₂ := by
    simpa [hgEq] using hi₂FactorG

  have hdvdj₁ : d ∣ v j₁ := by
    simpa [d] using (Int.gcd_dvd_left (v j₁) g)
  have hdvdj₂' : d ∣ v j₂ := by
    rw [← hgcdj₂Eq]
    exact_mod_cast Int.gcd_dvd_left (v j₂) g
  have hj₁Factor : v j₁ = d * a₁ :=
    (Int.mul_ediv_cancel' hdvdj₁).symm
  have hj₂Factor : v j₂ = d * a₂ :=
    (Int.mul_ediv_cancel' hdvdj₂').symm

  have hm₁Odd : ¬ (2 : ℤ) ∣ m₁ := by
    intro htwo
    obtain ⟨k, hk⟩ := htwo
    apply hi₁Wall.2
    refine ⟨k, ?_⟩
    rw [hi₁FactorG, hk]
    ring
  have hm₂Odd : ¬ (2 : ℤ) ∣ m₂ := by
    intro htwo
    obtain ⟨k, hk⟩ := htwo
    apply hi₂Wall.2
    refine ⟨k, ?_⟩
    rw [hi₂FactorG, hk]
    ring
  have ha₁Odd : ¬ (2 : ℤ) ∣ a₁ := by
    intro htwo
    obtain ⟨k, hk⟩ := htwo
    apply hj₁Not
    refine ⟨k, ?_⟩
    rw [hj₁Factor, hk, hgEq]
    ring
  have ha₂Odd : ¬ (2 : ℤ) ∣ a₂ := by
    intro htwo
    obtain ⟨k, hk⟩ := htwo
    apply hj₂Not
    refine ⟨k, ?_⟩
    rw [hj₂Factor, hk, hgEq]
    ring

  have hdvdTail : ∀ k, k ≠ i₁ → k ≠ i₂ → k ≠ j₁ → k ≠ j₂ →
      16 * d ∣ v k := by
    intro k hki₁ hki₂ hkj₁ hkj₂
    have hkg : g ∣ v k := hdvdBase k hkj₁ hkj₂
    have hk2g : 2 * g ∣ v k := by
      by_contra hkNot
      have hmem : k ∈ Finset.univ.filter
          (fun i => g ∣ v i ∧ ¬ 2 * g ∣ v i) := by
        simp [hkg, hkNot]
      rw [hfresh] at hmem
      simp only [Finset.mem_insert, Finset.mem_singleton] at hmem
      rcases hmem with h | h
      · exact hki₁ h
      · exact hki₂ h
    have hk4g : 4 * g ∣ v k := by
      convert hquietTwo k hk2g using 1 <;> ring
    have hk8g : 8 * g ∣ v k := by
      convert hquietFour k hk4g using 1 <;> ring
    rw [hgEq] at hk8g
    convert hk8g using 1 <;> ring

  exact lonely14_of_pairTower_two_min_gap_three
    cite v hv d hd i₁ i₂ j₁ j₂ hi₁i₂ hi₁j₁ hi₁j₂
    hi₂j₁ hi₂j₂ hj₁j₂ m₁ m₂ a₁ a₂ hi₁Factor hi₂Factor
    hj₁Factor hj₂Factor (gcd_eight_eq_one_of_not_two_dvd m₁ hm₁Odd)
    (gcd_eight_eq_one_of_not_two_dvd m₂ hm₂Odd)
    (gcd_sixteen_eq_one_of_not_two_dvd a₁ ha₁Odd)
    (gcd_sixteen_eq_one_of_not_two_dvd a₂ ha₂Odd) hdvdTail

/-! ## Remove the leaf from the remaining supplier -/

/-- The concrete quiet-lift condition, separated from base pair and
nongenericity hypotheses so it can be negated by the remaining supplier. -/
def PairTowerQuietLiftCondition (v : Fin 13 → ℤ) (g : ℤ) : Prop :=
  liftFailureCard v g = 2 ∧
    TwoAdicLiftTerminates v (2 * g) ∧
    TwoAdicLiftTerminates v (4 * g)

/-- Remaining many-failure supplier after deleting the proved genuine quiet
two-row leaf. -/
def ManyLiftFailureBeyondQuietLiftSupply : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → ∀ g : ℤ, 2 ≤ g →
    nonMultCard v g = 2 → 2 ≤ liftFailureCard v g →
    ¬ genericCount v g → ¬ PairTowerQuietLiftCondition v g →
    ∃ t : ℝ, Lonely 14 v t

/-- Reinserting the proved quiet leaf recovers the full first-wall
many-failure supplier. -/
theorem manyLiftFailurePairTowerSupply_of_beyondQuietLift
    (cite : LRCUpTo13)
    (hremaining : ManyLiftFailureBeyondQuietLiftSupply) :
    ManyLiftFailurePairTowerSupply := by
  intro v hv g hg hcard hmany hnongeneric
  by_cases hquiet : PairTowerQuietLiftCondition v g
  · exact lonely14_of_pairTower_quietLift cite v hv g hg hcard hnongeneric
      hquiet.1 hquiet.2.1 hquiet.2.2
  · exact hremaining v hv g hg hcard hmany hnongeneric hquiet

/-- Capstone with the genuine two-row quiet leaf removed from the pair-tower
hypothesis and the concrete signed-deviation dense-core socket retained. -/
theorem lrc14_from_pairTowerBeyondQuietLift_and_selectedWitnessSupplies_and_deviationBudget
    (cite : LRCUpTo13)
    (hremaining : ManyLiftFailureBeyondQuietLiftSupply)
    (h22 : TwoTwoSelectedWitnessSupply)
    (h244 : TwoFourFourSelectedWitnessSupply)
    (h333 : UniformThreeSelectedWitnessSupply)
    (hsupply : DenseCoreDeviationBudgetSupply) :
    LRC14.LRC14Statement :=
  lrc14_from_manyLiftFailure_and_selectedWitnessSupplies_and_deviationBudget
    cite (manyLiftFailurePairTowerSupply_of_beyondQuietLift cite hremaining)
    h22 h244 h333 hsupply

/-! ## Axiom audit -/

#print axioms lonely14_of_pairTower_quietLift
#print axioms manyLiftFailurePairTowerSupply_of_beyondQuietLift
#print axioms lrc14_from_pairTowerBeyondQuietLift_and_selectedWitnessSupplies_and_deviationBudget

end
end LRC14Grand
end LonelyRunner
