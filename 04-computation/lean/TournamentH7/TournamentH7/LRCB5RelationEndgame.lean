import TournamentH7.LRCEndgameParameterDischargeTwoThree
import TournamentH7.LRCB5NormalizedBridge

/-!
# Relation-budget certificates feed the chain-dense LRC(14) endgame

This module is the formal consumer for THM-935.  A certificate supplies a
positive modulus, exact-support relation masses, the exact identity between
their signed model and the concrete integer `B5`, and the sharp one-sided
horizon-thirty tail split.  The checked budget then gives `B5 > 0`, which feeds
the chain-dense endgame without any additional analytic assumption.
-/

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open LRCB5RelationBudget
open LRCB5NormalizedBridge
open LRCB5CleanModulus
open scoped Classical

/-- A proof-producing THM-935 certificate for one concrete speed tuple. -/
structure B5RelationBudgetCertificate (v : Fin 13 → ℤ) where
  q : ℕ
  one_lt_q : 1 < q
  mass2 : ℝ
  mass3 : ℝ
  mass4 : ℝ
  mass5 : ℝ
  b5_eq_scaled_model :
    (LRC14Concrete.B5 v q : ℝ) =
      ((q : ℝ) - 1) * relationModel mass2 mass3 mass4 mass5
  pair_mass_lower_bound : -(13 / 30 : ℝ) ≤ mass2
  harmful_higher_budget :
    harmfulHigherContribution mass3 mass4 mass5 < 7712 / 84035

/-- Every relation-budget certificate produces the integer positivity witness
consumed by the discrete Bonferroni endgame. -/
theorem B5RelationBudgetCertificate.b5_pos {v : Fin 13 → ℤ}
    (certificate : B5RelationBudgetCertificate v) :
    0 < LRC14Concrete.B5 v certificate.q := by
  have hmodel : 0 < relationModel certificate.mass2 certificate.mass3
      certificate.mass4 certificate.mass5 :=
    relationModel_pos_of_signed_horizon_thirty_split _ _ _ _
      certificate.pair_mass_lower_bound certificate.harmful_higher_budget
  have hqR : (1 : ℝ) < certificate.q := by
    exact_mod_cast certificate.one_lt_q
  have hscale : (0 : ℝ) < (certificate.q : ℝ) - 1 := by
    linarith
  have hreal : (0 : ℝ) < (LRC14Concrete.B5 v certificate.q : ℝ) := by
    rw [certificate.b5_eq_scaled_model]
    positivity
  exact_mod_cast hreal

/-- Concrete THM-940 certificate at the canonical clean modulus.  The modulus,
coprimality, singleton sign, and finite-height modular-to-exact relation bridge
are automatic; only the second factorial-moment floor and harmful depth budget
remain. -/
structure NormalizedB5RelationBudgetCertificate (v : Fin 13 → ℤ) where
  height : ℕ
  pair_depth_budget : 1703 / 1470 ≤
    normalizedPairDepthMoment v (cleanModulus v height)
  harmful_depth_budget :
    harmfulDepthMoment v (cleanModulus v height) /
      ((cleanModulus v height : ℚ) - 1) < -(65218 / 84035)

/-- The normalized certificate feeds the concrete integer B5 directly; no
separate relation-mass identification premise remains. -/
theorem NormalizedB5RelationBudgetCertificate.b5_pos {v : Fin 13 → ℤ}
    (certificate : NormalizedB5RelationBudgetCertificate v)
    (hv : ∀ i, v i ≠ 0) :
    0 < LRC14Concrete.B5 v (cleanModulus v certificate.height) :=
  B5_pos_at_cleanModulus_of_depthMoment_budgets v hv certificate.height
    certificate.pair_depth_budget certificate.harmful_depth_budget

/-- THM-935 certificate supply only on the primitive, dissociated,
chain-dense core isolated by the current endgame. -/
def DenseCoreRelationBudgetSupply : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.tupleGcd v = 1 →
    LRC14.CoveringFamily v → GapFamily v →
    (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) →
    (∀ i j, i ≠ j → |v i| ≠ |v j|) →
    (∃ i, 23 ≤ |v i|) →
    (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13,
      (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
    (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ),
      (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
      (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1 / 13 - 1 / 14 ∧
      (∀ i, k i ≠ 0) ∧ (Finset.univ.image k).card ≤ 12) →
    (∀ g : ℤ, 2 ≤ g →
      nonMultCard v g ≠ 2 ∧ nonMultCard v g ≠ 3) →
    (∃ σ : Equiv.Perm (Fin 13), Monotone (fun i => |v (σ i)|) ∧
      ChainDenseCore (fun i => |v (σ i)|)) →
    Nonempty (B5RelationBudgetCertificate v)

/-- Stronger concrete supplier on the same primitive, dissociated,
chain-dense core.  Its open mathematics is exactly a coprime modulus together
with the normalized pair and signed higher-support estimates. -/
def DenseCoreNormalizedRelationBudgetSupply : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → LRC14.tupleGcd v = 1 →
    LRC14.CoveringFamily v → GapFamily v →
    (∀ i, ∃ j, j ≠ i ∧ |v i| ≤ 13 * |v j|) →
    (∀ i j, i ≠ j → |v i| ≠ |v j|) →
    (∃ i, 23 ≤ |v i|) →
    (∀ g : ℤ, 2 ≤ g → ∀ i₀ : Fin 13,
      (∀ j, j ≠ i₀ → g ∣ v j) → g ∣ v i₀) →
    (¬ ∃ (L : ℤ) (k a : Fin 13 → ℤ) (A : ℝ),
      (∀ i, v i = a i + L * k i) ∧ 0 < (L : ℝ) ∧
      (∀ i, |(a i : ℝ)| ≤ A) ∧ A / (L : ℝ) ≤ 1 / 13 - 1 / 14 ∧
      (∀ i, k i ≠ 0) ∧ (Finset.univ.image k).card ≤ 12) →
    (∀ g : ℤ, 2 ≤ g →
      nonMultCard v g ≠ 2 ∧ nonMultCard v g ≠ 3) →
    (∃ σ : Equiv.Perm (Fin 13), Monotone (fun i => |v (σ i)|) ∧
      ChainDenseCore (fun i => |v (σ i)|)) →
    Nonempty (NormalizedB5RelationBudgetCertificate v)

/-- The structured relation-budget supplier discharges the raw positive-B5
supplier used by `LRCDenseCoreEndgame`. -/
theorem denseCoreDissociatedB5Supply_of_relationBudget
    (hsupply : DenseCoreRelationBudgetSupply) :
    DenseCoreDissociatedB5Supply := by
  intro v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc hcore
  obtain ⟨certificate⟩ :=
    hsupply v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc hcore
  exact ⟨certificate.q, lt_trans Nat.zero_lt_one certificate.one_lt_q,
    certificate.b5_pos⟩

/-- The normalized THM-940 supplier discharges the raw positive-B5 supplier
without routing through the older abstract equality certificate. -/
theorem denseCoreDissociatedB5Supply_of_normalizedRelationBudget
    (hsupply : DenseCoreNormalizedRelationBudgetSupply) :
    DenseCoreDissociatedB5Supply := by
  intro v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc hcore
  obtain ⟨certificate⟩ :=
    hsupply v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc hcore
  exact ⟨cleanModulus v certificate.height,
    lt_of_lt_of_le (by omega) (fourteen_le_cleanModulus v certificate.height hv),
    certificate.b5_pos hv⟩

/-- **THM-935-shaped current capstone.**  The sanctioned LRC(≤13) citation,
the sharply reduced `q=2,3`/two-adic exceptional dispatch, and relation-budget
certificates only on the primitive dissociated chain-dense core imply LRC(14). -/
theorem lrc14_from_four_detuned_and_relationBudget
    (cite : LRCUpTo13) (hdeep : DeepExceptionalDetunedDispatchFour)
    (hsupply : DenseCoreRelationBudgetSupply) : LRC14.LRC14Statement :=
  lrc14_from_four_detuned_and_denseCore_dissociated_B5 cite hdeep
    (denseCoreDissociatedB5Supply_of_relationBudget hsupply)

/-- **Sharpest THM-935-shaped capstone.**  Degree arithmetic reduces the
three-detuned residual to a `q = 2` row with a distinct `q ≤ 8` companion or
the uniform `(3,3,3)` pattern; relation-budget certificates are required only
on the primitive dissociated chain-dense core. -/
theorem lrc14_from_twoThree_detuned_and_relationBudget
    (cite : LRCUpTo13) (hdeep : DeepExceptionalDetunedDispatchTwoThree)
    (hsupply : DenseCoreRelationBudgetSupply) : LRC14.LRC14Statement :=
  lrc14_from_twoThree_detuned_and_denseCore_dissociated_B5 cite hdeep
    (denseCoreDissociatedB5Supply_of_relationBudget hsupply)

/-- Concrete normalized-THM-940 capstone.  The B5 side now exposes only the
second factorial-moment floor and signed higher-depth budget at a coprime
modulus. -/
theorem lrc14_from_twoThree_detuned_and_normalizedRelationBudget
    (cite : LRCUpTo13) (hdeep : DeepExceptionalDetunedDispatchTwoThree)
    (hsupply : DenseCoreNormalizedRelationBudgetSupply) :
    LRC14.LRC14Statement :=
  lrc14_from_twoThree_detuned_and_denseCore_dissociated_B5 cite hdeep
    (denseCoreDissociatedB5Supply_of_normalizedRelationBudget hsupply)

/-! ## Axiom audit -/

#print axioms B5RelationBudgetCertificate.b5_pos
#print axioms NormalizedB5RelationBudgetCertificate.b5_pos
#print axioms denseCoreDissociatedB5Supply_of_relationBudget
#print axioms denseCoreDissociatedB5Supply_of_normalizedRelationBudget
#print axioms lrc14_from_four_detuned_and_relationBudget
#print axioms lrc14_from_twoThree_detuned_and_relationBudget
#print axioms lrc14_from_twoThree_detuned_and_normalizedRelationBudget

end LRC14Grand
end LonelyRunner
