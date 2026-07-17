/-
  TournamentH7.LRCSelectedWitnessResidual

  Restates the three selected-phase sockets in the form consumed by the
  endgame: with the nongeneric branch-level hypothesis visible.  This is a
  proof-local normalization, not automatically a weaker mathematical theorem:
  the exceptional denominator patterns themselves force nongenericity.  The
  converse adapters below make that equivalence explicit rather than claiming
  an artificial reduction of the frontier.

  Tournament-analysis audit: the useful carrier is the bipartite incidence
  hypergraph between harmonic-safe phase cells and detuned branch obligations.
  A pairwise orientation on runners, residues, or cells loses the existential
  common-incidence predicate, so no arbitrary tournament gauge or tie
  Hamiltonian path is promoted.  The challenged assumption is that a supplier
  should cover every denominator pattern: the proof quotient records the
  nongeneric hypothesis at its point of use, while the pattern arithmetic
  explains why it is redundant in the supplier statements.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCEndgameTwoEight

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

noncomputable section

/-- A coordinate whose reduced denominator is bigger than one is genuinely
detuned.  This elementary bridge is shared by the three automatic
nongenericity calculations below and the pair-tower valuation leaves. -/
theorem not_dvd_of_one_lt_reducedDenominator
    (x g q : ℤ) (hg : 0 < g)
    (hq : g / (Int.gcd x g : ℤ) = q) (hq1 : 1 < q) :
    ¬ g ∣ x := by
  intro hdvd
  have habs : g.natAbs ∣ x.natAbs :=
    Int.natAbs_dvd_natAbs.mpr hdvd
  have hgcd : Int.gcd x g = g.natAbs := by
    rw [Int.gcd_def]
    exact Nat.gcd_eq_right_iff_dvd.mpr habs
  have hgabs : (g.natAbs : ℤ) = g := by
    rw [Int.natAbs_of_nonneg (le_of_lt hg)]
  have hqone : g / (Int.gcd x g : ℤ) = 1 := by
    rw [hgcd, hgabs]
    exact Int.ediv_self (ne_of_gt hg)
  omega

/-- If the three selected bad rows already carry total degree at least `g`,
then their exact three-row decomposition cannot satisfy the generic strict
union bound. -/
theorem not_genericCount_of_three_of_badCount_sum_ge
    (v : Fin 13 → ℤ) (g : ℤ) (i₁ i₂ i₃ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃)
    (hδ1 : ¬ g ∣ v i₁) (hδ2 : ¬ g ∣ v i₂) (hδ3 : ¬ g ∣ v i₃)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j)
    (hge : g.toNat ≤ DetunedD3.badCount (v i₁) g +
      DetunedD3.badCount (v i₂) g + DetunedD3.badCount (v i₃) g) :
    ¬ genericCount v g := by
  intro hgeneric
  have hfilter :
      Finset.univ.filter (fun j : Fin 13 => ¬ g ∣ v j) = {i₁, i₂, i₃} := by
    ext j
    simp only [Finset.mem_filter, Finset.mem_univ, true_and,
      Finset.mem_insert, Finset.mem_singleton]
    constructor
    · intro hj
      by_cases hj1 : j = i₁
      · exact Or.inl hj1
      by_cases hj2 : j = i₂
      · exact Or.inr (Or.inl hj2)
      by_cases hj3 : j = i₃
      · exact Or.inr (Or.inr hj3)
      exact False.elim (hj (hdvd j hj1 hj2 hj3))
    · rintro (rfl | rfl | rfl)
      · exact hδ1
      · exact hδ2
      · exact hδ3
  rw [genericCount, hfilter,
    Finset.sum_insert (by simp [h12, h13]),
    Finset.sum_insert (by simp [h23]), Finset.sum_singleton] at hgeneric
  omega

/-- The two q-two rows alone saturate the generic branch budget. -/
theorem not_genericCount_two_two_of_reducedDenominators
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g)
    (i₂a i₂b iₓ : Fin 13)
    (h2ab : i₂a ≠ i₂b) (h2ax : i₂a ≠ iₓ) (h2bx : i₂b ≠ iₓ)
    (hdvd : ∀ j, j ≠ i₂a → j ≠ i₂b → j ≠ iₓ → g ∣ v j)
    (hδx : ¬ g ∣ v iₓ)
    (hq2a : g / (Int.gcd (v i₂a) g : ℤ) = 2)
    (hq2b : g / (Int.gcd (v i₂b) g : ℤ) = 2) :
    ¬ genericCount v g := by
  have hg0 : (0 : ℤ) < g := by omega
  have hδa := not_dvd_of_one_lt_reducedDenominator
    (v i₂a) g 2 hg0 hq2a (by omega)
  have hδb := not_dvd_of_one_lt_reducedDenominator
    (v i₂b) g 2 hg0 hq2b (by omega)
  apply not_genericCount_of_three_of_badCount_sum_ge
    v g i₂a i₂b iₓ h2ab h2ax h2bx hδa hδb hδx hdvd
  have hbadA := two_mul_badCount_eq (v i₂a) g hg0 hq2a
  have hbadB := two_mul_badCount_eq (v i₂b) g hg0 hq2b
  omega

/-- The q-two and two q-four rows exactly saturate the generic branch budget. -/
theorem not_genericCount_two_four_four_of_reducedDenominators
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g)
    (i₂ i₄a i₄b : Fin 13)
    (h24a : i₂ ≠ i₄a) (h24b : i₂ ≠ i₄b) (h4ab : i₄a ≠ i₄b)
    (hdvd : ∀ j, j ≠ i₂ → j ≠ i₄a → j ≠ i₄b → g ∣ v j)
    (hq2 : g / (Int.gcd (v i₂) g : ℤ) = 2)
    (hq4a : g / (Int.gcd (v i₄a) g : ℤ) = 4)
    (hq4b : g / (Int.gcd (v i₄b) g : ℤ) = 4) :
    ¬ genericCount v g := by
  have hg0 : (0 : ℤ) < g := by omega
  have hδ2 := not_dvd_of_one_lt_reducedDenominator
    (v i₂) g 2 hg0 hq2 (by omega)
  have hδ4a := not_dvd_of_one_lt_reducedDenominator
    (v i₄a) g 4 hg0 hq4a (by omega)
  have hδ4b := not_dvd_of_one_lt_reducedDenominator
    (v i₄b) g 4 hg0 hq4b (by omega)
  apply not_genericCount_of_three_of_badCount_sum_ge
    v g i₂ i₄a i₄b h24a h24b h4ab hδ2 hδ4a hδ4b hdvd
  have hbad2 := two_mul_badCount_eq (v i₂) g hg0 hq2
  have hbad4a := four_mul_badCount_eq_of_reducedDenominator_four
    (v i₄a) g hg0 hq4a
  have hbad4b := four_mul_badCount_eq_of_reducedDenominator_four
    (v i₄b) g hg0 hq4b
  omega

/-- Three q-three rows exactly saturate the generic branch budget. -/
theorem not_genericCount_uniform_three_of_reducedDenominators
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g)
    (i₁ i₂ i₃ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j)
    (hq1 : g / (Int.gcd (v i₁) g : ℤ) = 3)
    (hq2 : g / (Int.gcd (v i₂) g : ℤ) = 3)
    (hq3 : g / (Int.gcd (v i₃) g : ℤ) = 3) :
    ¬ genericCount v g := by
  have hg0 : (0 : ℤ) < g := by omega
  have hδ1 := not_dvd_of_one_lt_reducedDenominator
    (v i₁) g 3 hg0 hq1 (by omega)
  have hδ2 := not_dvd_of_one_lt_reducedDenominator
    (v i₂) g 3 hg0 hq2 (by omega)
  have hδ3 := not_dvd_of_one_lt_reducedDenominator
    (v i₃) g 3 hg0 hq3 (by omega)
  apply not_genericCount_of_three_of_badCount_sum_ge
    v g i₁ i₂ i₃ h12 h13 h23 hδ1 hδ2 hδ3 hdvd
  have hbad1 := three_mul_badCount_eq (v i₁) g hg0 hq1
  have hbad2 := three_mul_badCount_eq (v i₂) g hg0 hq2
  have hbad3 := three_mul_badCount_eq (v i₃) g hg0 hq3
  omega

/-- The `(2,2,q)` selected-phase theorem only on the nongeneric residual. -/
def ResidualTwoTwoSelectedWitnessSupply : Prop :=
  ∀ (v : Fin 13 → ℤ), (∀ i, v i ≠ 0) → ∀ (g : ℤ), 2 ≤ g →
    ¬ genericCount v g →
    ∀ i₂a i₂b iₓ : Fin 13,
    i₂a ≠ i₂b → i₂a ≠ iₓ → i₂b ≠ iₓ →
    (∀ j, j ≠ i₂a → j ≠ i₂b → j ≠ iₓ → g ∣ v j) →
    ¬ g ∣ v iₓ →
    g / (Int.gcd (v i₂a) g : ℤ) = 2 →
    g / (Int.gcd (v i₂b) g : ℤ) = 2 →
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₂a i₂b iₓ u ∧
      HasThreeDetunedGoodBranch (v i₂a) (v i₂b) (v iₓ) g u

/-- The `(2,4,4)` selected-phase theorem only on the nongeneric residual. -/
def ResidualTwoFourFourSelectedWitnessSupply : Prop :=
  ∀ (v : Fin 13 → ℤ), (∀ i, v i ≠ 0) → ∀ (g : ℤ), 2 ≤ g →
    ¬ genericCount v g →
    ∀ i₂ i₄a i₄b : Fin 13,
    i₂ ≠ i₄a → i₂ ≠ i₄b → i₄a ≠ i₄b →
    (∀ j, j ≠ i₂ → j ≠ i₄a → j ≠ i₄b → g ∣ v j) →
    g / (Int.gcd (v i₂) g : ℤ) = 2 →
    g / (Int.gcd (v i₄a) g : ℤ) = 4 →
    g / (Int.gcd (v i₄b) g : ℤ) = 4 →
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₂ i₄a i₄b u ∧
      HasThreeDetunedGoodBranch (v i₂) (v i₄a) (v i₄b) g u

/-- The uniform-q-three selected-phase theorem only on the nongeneric residual. -/
def ResidualUniformThreeSelectedWitnessSupply : Prop :=
  ∀ (v : Fin 13 → ℤ), (∀ i, v i ≠ 0) → ∀ (g : ℤ), 2 ≤ g →
    ¬ genericCount v g →
    ∀ i₁ i₂ i₃ : Fin 13,
    i₁ ≠ i₂ → i₁ ≠ i₃ → i₂ ≠ i₃ →
    (∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j) →
    g / (Int.gcd (v i₁) g : ℤ) = 3 →
    g / (Int.gcd (v i₂) g : ℤ) = 3 →
    g / (Int.gcd (v i₃) g : ℤ) = 3 →
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₁ i₂ i₃ u ∧
      HasThreeDetunedGoodBranch (v i₁) (v i₂) (v i₃) g u

/-- Strong all-level suppliers restrict to the residual suppliers. -/
theorem residualSelectedWitnessSupplies_of_selectedWitnessSupplies
    (h22 : TwoTwoSelectedWitnessSupply)
    (h244 : TwoFourFourSelectedWitnessSupply)
    (h333 : UniformThreeSelectedWitnessSupply) :
    ResidualTwoTwoSelectedWitnessSupply ∧
      ResidualTwoFourFourSelectedWitnessSupply ∧
      ResidualUniformThreeSelectedWitnessSupply := by
  refine ⟨?_, ?_, ?_⟩
  · intro v hv g hg _
    exact h22 v hv g hg
  · intro v hv g hg _
    exact h244 v hv g hg
  · intro v hv g hg _
    exact h333 v hv g hg

/-- In the `(2,2,q)` pattern the residual hypothesis is automatic, so the
residual supplier recovers the original supplier. -/
theorem twoTwoSelectedWitnessSupply_of_residual
    (hres : ResidualTwoTwoSelectedWitnessSupply) :
    TwoTwoSelectedWitnessSupply := by
  intro v hv g hg i₂a i₂b iₓ h2ab h2ax h2bx hdvd hδx hq2a hq2b
  exact hres v hv g hg
    (not_genericCount_two_two_of_reducedDenominators
      v g hg i₂a i₂b iₓ h2ab h2ax h2bx hdvd hδx hq2a hq2b)
    i₂a i₂b iₓ h2ab h2ax h2bx hdvd hδx hq2a hq2b

/-- In the `(2,4,4)` pattern the residual hypothesis is automatic. -/
theorem twoFourFourSelectedWitnessSupply_of_residual
    (hres : ResidualTwoFourFourSelectedWitnessSupply) :
    TwoFourFourSelectedWitnessSupply := by
  intro v hv g hg i₂ i₄a i₄b h24a h24b h4ab hdvd hq2 hq4a hq4b
  exact hres v hv g hg
    (not_genericCount_two_four_four_of_reducedDenominators
      v g hg i₂ i₄a i₄b h24a h24b h4ab hdvd hq2 hq4a hq4b)
    i₂ i₄a i₄b h24a h24b h4ab hdvd hq2 hq4a hq4b

/-- In the uniform `(3,3,3)` pattern the residual hypothesis is automatic. -/
theorem uniformThreeSelectedWitnessSupply_of_residual
    (hres : ResidualUniformThreeSelectedWitnessSupply) :
    UniformThreeSelectedWitnessSupply := by
  intro v hv g hg i₁ i₂ i₃ h12 h13 h23 hdvd hq1 hq2 hq3
  exact hres v hv g hg
    (not_genericCount_uniform_three_of_reducedDenominators
      v g hg i₁ i₂ i₃ h12 h13 h23 hdvd hq1 hq2 hq3)
    i₁ i₂ i₃ h12 h13 h23 hdvd hq1 hq2 hq3

/-- The residual and original selected-witness supplier triples are exactly
equivalent; the residual presentation exposes a proof-local hypothesis but
does not shrink the mathematical obligation. -/
theorem selectedWitnessSupplies_of_residualSelectedWitnessSupplies
    (h22 : ResidualTwoTwoSelectedWitnessSupply)
    (h244 : ResidualTwoFourFourSelectedWitnessSupply)
    (h333 : ResidualUniformThreeSelectedWitnessSupply) :
    TwoTwoSelectedWitnessSupply ∧ TwoFourFourSelectedWitnessSupply ∧
      UniformThreeSelectedWitnessSupply :=
  ⟨twoTwoSelectedWitnessSupply_of_residual h22,
    twoFourFourSelectedWitnessSupply_of_residual h244,
    uniformThreeSelectedWitnessSupply_of_residual h333⟩

/-- The three residual suppliers, together with the pair tower, discharge the
exact final-residue dispatch without asking for selected witnesses on generic
levels. -/
theorem deepExceptionalDetunedDispatchFinalResidues_of_residualSelectedWitnessSupplies
    (hpairs : NonterminatingPairTowerSupply)
    (h22 : ResidualTwoTwoSelectedWitnessSupply)
    (h244 : ResidualTwoFourFourSelectedWitnessSupply)
    (h333 : ResidualUniformThreeSelectedWitnessSupply) :
    DeepExceptionalDetunedDispatchFinalResidues := by
  intro v hv g hg hcase hnongeneric
  rcases hcase with ⟨hcard, hnonterm⟩ | ⟨hcard, hpattern⟩
  · exact hpairs v hv g hg hcard hnonterm hnongeneric
  · rcases hpattern with h22pattern | h244pattern | hallThree
    · obtain ⟨i₂a, i₂b, iₓ, h2ab, h2ax, h2bx, hδ2a, hδ2b, hδx,
        hq2a, hq2b⟩ := h22pattern
      have hdvd := dvd_of_nonMultCard_three_of_selected
        v g hcard i₂a i₂b iₓ h2ab h2ax h2bx hδ2a hδ2b hδx
      apply lonely14_of_three_detuned_selectedWitness v g hg i₂a i₂b iₓ hdvd
      exact h22 v hv g hg hnongeneric i₂a i₂b iₓ h2ab h2ax h2bx hdvd
        hδx hq2a hq2b
    · obtain ⟨i₂, i₄a, i₄b, h24a, h24b, h4ab, hδ2, hδ4a, hδ4b,
        hq2, hq4a, hq4b⟩ := h244pattern
      have hdvd := dvd_of_nonMultCard_three_of_selected
        v g hcard i₂ i₄a i₄b h24a h24b h4ab hδ2 hδ4a hδ4b
      apply lonely14_of_three_detuned_selectedWitness v g hg i₂ i₄a i₄b hdvd
      exact h244 v hv g hg hnongeneric i₂ i₄a i₄b h24a h24b h4ab hdvd
        hq2 hq4a hq4b
    · have hcard' := hcard
      rw [nonMultCard] at hcard'
      obtain ⟨i₁, i₂, i₃, h12, h13, h23, hfilter⟩ :=
        Finset.card_eq_three.mp hcard'
      have hδ1 : ¬ g ∣ v i₁ := by
        have : i₁ ∈ Finset.univ.filter (fun i => ¬ g ∣ v i) := by
          rw [hfilter]
          simp
        simpa using this
      have hδ2 : ¬ g ∣ v i₂ := by
        have : i₂ ∈ Finset.univ.filter (fun i => ¬ g ∣ v i) := by
          rw [hfilter]
          simp
        simpa using this
      have hδ3 : ¬ g ∣ v i₃ := by
        have : i₃ ∈ Finset.univ.filter (fun i => ¬ g ∣ v i) := by
          rw [hfilter]
          simp
        simpa using this
      have hdvd := dvd_of_nonMultCard_three_of_selected
        v g hcard i₁ i₂ i₃ h12 h13 h23 hδ1 hδ2 hδ3
      apply lonely14_of_three_detuned_selectedWitness v g hg i₁ i₂ i₃ hdvd
      exact h333 v hv g hg hnongeneric i₁ i₂ i₃ h12 h13 h23 hdvd
        (hallThree i₁ hδ1) (hallThree i₂ hδ2) (hallThree i₃ hδ3)

/-- Sharp residual-selected capstone with the analytic relation-budget socket. -/
theorem lrc14_from_residualSelectedWitnessSupplies_and_relationBudget
    (cite : LRCUpTo13)
    (hpairs : NonterminatingPairTowerSupply)
    (h22 : ResidualTwoTwoSelectedWitnessSupply)
    (h244 : ResidualTwoFourFourSelectedWitnessSupply)
    (h333 : ResidualUniformThreeSelectedWitnessSupply)
    (hsupply : DenseCoreRelationBudgetSupply) :
    LRC14.LRC14Statement :=
  lrc14_from_finalResidues_and_relationBudget cite
    (deepExceptionalDetunedDispatchFinalResidues_of_residualSelectedWitnessSupplies
      hpairs h22 h244 h333) hsupply

/-- Sharp residual-selected capstone with THM-940's concrete deviation socket. -/
theorem lrc14_from_residualSelectedWitnessSupplies_and_deviationBudget
    (cite : LRCUpTo13)
    (hpairs : NonterminatingPairTowerSupply)
    (h22 : ResidualTwoTwoSelectedWitnessSupply)
    (h244 : ResidualTwoFourFourSelectedWitnessSupply)
    (h333 : ResidualUniformThreeSelectedWitnessSupply)
    (hsupply : DenseCoreDeviationBudgetSupply) :
    LRC14.LRC14Statement :=
  lrc14_from_twoThree_detuned_and_denseCore_dissociated_B5 cite
    (deepExceptionalDetunedDispatchTwoThree_of_finalResidues cite
      (deepExceptionalDetunedDispatchFinalResidues_of_residualSelectedWitnessSupplies
        hpairs h22 h244 h333))
    (denseCoreDissociatedB5Supply_of_deviationBudget hsupply)

/-- Residual-selected capstone with the incoming clean-modulus normalized
coverage/depth certificate. -/
theorem lrc14_from_residualSelectedWitnessSupplies_and_normalizedRelationBudget
    (cite : LRCUpTo13)
    (hpairs : NonterminatingPairTowerSupply)
    (h22 : ResidualTwoTwoSelectedWitnessSupply)
    (h244 : ResidualTwoFourFourSelectedWitnessSupply)
    (h333 : ResidualUniformThreeSelectedWitnessSupply)
    (hsupply : DenseCoreNormalizedRelationBudgetSupply) :
    LRC14.LRC14Statement :=
  lrc14_from_twoThree_detuned_and_denseCore_dissociated_B5 cite
    (deepExceptionalDetunedDispatchTwoThree_of_finalResidues cite
      (deepExceptionalDetunedDispatchFinalResidues_of_residualSelectedWitnessSupplies
        hpairs h22 h244 h333))
    (denseCoreDissociatedB5Supply_of_normalizedRelationBudget hsupply)

/-- Residual-selected capstone with THM-950's unconditional uniform census.
The pair and phase selectors are needed only on their final residue families,
while the dense core needs one usable modulus whose live count beats `792`
times its full depth-at-least-six count. -/
theorem lrc14_from_residualSelectedWitnessSupplies_and_censusB5
    (cite : LRCUpTo13)
    (hpairs : NonterminatingPairTowerSupply)
    (h22 : ResidualTwoTwoSelectedWitnessSupply)
    (h244 : ResidualTwoFourFourSelectedWitnessSupply)
    (h333 : ResidualUniformThreeSelectedWitnessSupply)
    (hsupply : DenseCoreCensusB5Supply) :
    LRC14.LRC14Statement :=
  lrc14_from_twoThree_detuned_and_denseCore_dissociated_B5 cite
    (deepExceptionalDetunedDispatchTwoThree_of_finalResidues cite
      (deepExceptionalDetunedDispatchFinalResidues_of_residualSelectedWitnessSupplies
        hpairs h22 h244 h333))
    (denseCoreDissociatedB5Supply_of_census hsupply)

/-- Sharp residual-selected composition with the exact weighted census.  The
dense core now pays `choose (bandCount - 1) 5` at each multiplier rather than
the uniform depth-thirteen envelope. -/
theorem lrc14_from_residualSelectedWitnessSupplies_and_weightedCensusB5
    (cite : LRCUpTo13)
    (hpairs : NonterminatingPairTowerSupply)
    (h22 : ResidualTwoTwoSelectedWitnessSupply)
    (h244 : ResidualTwoFourFourSelectedWitnessSupply)
    (h333 : ResidualUniformThreeSelectedWitnessSupply)
    (hsupply : DenseCoreWeightedCensusB5Supply) :
    LRC14.LRC14Statement :=
  lrc14_from_twoThree_detuned_and_denseCore_dissociated_B5 cite
    (deepExceptionalDetunedDispatchTwoThree_of_finalResidues cite
      (deepExceptionalDetunedDispatchFinalResidues_of_residualSelectedWitnessSupplies
        hpairs h22 h244 h333))
    (denseCoreDissociatedB5Supply_of_weightedCensus hsupply)

/-! ## Axiom audit -/

#print axioms residualSelectedWitnessSupplies_of_selectedWitnessSupplies
#print axioms not_genericCount_of_three_of_badCount_sum_ge
#print axioms not_genericCount_two_two_of_reducedDenominators
#print axioms not_genericCount_two_four_four_of_reducedDenominators
#print axioms not_genericCount_uniform_three_of_reducedDenominators
#print axioms selectedWitnessSupplies_of_residualSelectedWitnessSupplies
#print axioms deepExceptionalDetunedDispatchFinalResidues_of_residualSelectedWitnessSupplies
#print axioms lrc14_from_residualSelectedWitnessSupplies_and_relationBudget
#print axioms lrc14_from_residualSelectedWitnessSupplies_and_deviationBudget
#print axioms lrc14_from_residualSelectedWitnessSupplies_and_normalizedRelationBudget
#print axioms lrc14_from_residualSelectedWitnessSupplies_and_censusB5
#print axioms lrc14_from_residualSelectedWitnessSupplies_and_weightedCensusB5

end
end LRC14Grand
end LonelyRunner
