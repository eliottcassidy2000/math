import TournamentH7.LRCEndgameParameterDischargeTwoThree
import TournamentH7.LRCB5RelationBudget

/-!
# Relation-budget certificates feed the chain-dense LRC(14) endgame

This module is the formal consumer for THM-935.  A certificate supplies a
positive modulus, exact-support relation masses, the exact identity between
their signed model and the concrete integer `B5`, and the proved quarter / open
three-quarter tail split.  The checked budget then gives `B5 > 0`, which feeds
the chain-dense endgame without any additional analytic assumption.
-/

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open LRCB5RelationBudget
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
  pair_budget : pairWeight * |mass2| ≤ equilibrium / 4
  higher_budget : higherRelationDebt mass3 mass4 mass5 < 3 * equilibrium / 4

/-- Every relation-budget certificate produces the integer positivity witness
consumed by the discrete Bonferroni endgame. -/
theorem B5RelationBudgetCertificate.b5_pos {v : Fin 13 → ℤ}
    (certificate : B5RelationBudgetCertificate v) :
    0 < LRC14Concrete.B5 v certificate.q := by
  have hmodel : 0 < relationModel certificate.mass2 certificate.mass3
      certificate.mass4 certificate.mass5 :=
    relationModel_pos_of_quarter_threeQuarter_split _ _ _ _
      certificate.pair_budget certificate.higher_budget
  have hqR : (1 : ℝ) < certificate.q := by
    exact_mod_cast certificate.one_lt_q
  have hscale : (0 : ℝ) < (certificate.q : ℝ) - 1 := by
    linarith
  have hreal : (0 : ℝ) < (LRC14Concrete.B5 v certificate.q : ℝ) := by
    rw [certificate.b5_eq_scaled_model]
    positivity
  exact_mod_cast hreal

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

/-! ## Axiom audit -/

#print axioms B5RelationBudgetCertificate.b5_pos
#print axioms denseCoreDissociatedB5Supply_of_relationBudget
#print axioms lrc14_from_four_detuned_and_relationBudget
#print axioms lrc14_from_twoThree_detuned_and_relationBudget

end LRC14Grand
end LonelyRunner
