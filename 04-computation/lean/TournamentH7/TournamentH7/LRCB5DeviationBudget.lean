/-
  TournamentH7.LRCB5DeviationBudget

  Proof-facing names for THM-940's concrete deviation layers, the sharp
  one-sided positivity criterion, and an honest layer-by-layer adapter to the
  signed relation-shaped model.  The adapter keeps the finite-grid singleton
  defect explicit instead of silently absorbing it into an unconstrained
  mass.  Its support masses are normalized *discrete deviation layers*; a
  future comparison with continuous Fourier masses must add quadrature-error
  fields separately.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCB5SubsetExpansion
import TournamentH7.LRCB5RelationBudget

namespace LonelyRunner
namespace LRC14Concrete

open Finset
open LRCB5RelationBudget

/-- The total THM-940 deviation on subsets of one fixed cardinality. -/
def deviationLayer (v : Fin 13 → ℤ) (q k : ℕ) : ℚ :=
  ∑ T ∈ (Finset.univ : Finset (Fin 13)).powersetCard k, deviation v q T

/-- The signed depth-five deviation ledger, before adding the equilibrium. -/
def signedDeviationLedger (v : Fin 13 → ℤ) (q : ℕ) : ℚ :=
  ∑ k ∈ Finset.range 6, (-1 : ℚ) ^ k * deviationLayer v q k

/-- THM-940 in the named layer coordinates. -/
theorem B5_eq_equilibrium_add_signedDeviationLedger
    (v : Fin 13 → ℤ) (q : ℕ) :
    (B5 v q : ℚ) =
      ((q : ℚ) - 1) * (2052 / 16807) + signedDeviationLedger v q := by
  simpa [signedDeviationLedger, deviationLayer] using
    B5_eq_equilibrium_add_deviation v q

/-- The sharp sufficient condition is one-sided: favorable negative deviation
need not be charged through an absolute value. -/
theorem B5_pos_of_signedDeviationLedger_lower
    (v : Fin 13 → ℤ) (q : ℕ)
    (hlower : -(((q : ℚ) - 1) * (2052 / 16807)) <
      signedDeviationLedger v q) :
    0 < B5 v q := by
  have hid := B5_eq_equilibrium_add_signedDeviationLedger v q
  have hQ : (0 : ℚ) < (B5 v q : ℚ) := by
    rw [hid]
    linarith
  exact_mod_cast hQ

/-- The empty-support layer vanishes as soon as the multiplier grid is
nonempty. -/
theorem deviationLayer_zero (v : Fin 13 → ℤ) (q : ℕ) (hq : 1 ≤ q) :
    deviationLayer v q 0 = 0 := by
  rw [deviationLayer, Finset.powersetCard_zero, Finset.sum_singleton,
    deviation_empty v q hq]

/-- Explicit six-layer expansion, used by the bridge below. -/
theorem signedDeviationLedger_expand (v : Fin 13 → ℤ) (q : ℕ) :
    signedDeviationLedger v q =
      deviationLayer v q 0 - deviationLayer v q 1 +
      deviationLayer v q 2 - deviationLayer v q 3 +
      deviationLayer v q 4 - deviationLayer v q 5 := by
  rw [signedDeviationLedger]
  norm_num [Finset.sum_range_succ]
  ring

/-- Layer-by-layer finite-grid identification with the signed relation-shaped
model.  `singletonDefect` is exactly normalized support-one deviation.  It is
not a quadrature error, and the fields do not yet identify continuous
Fourier/singular-series masses with these discrete support coordinates. -/
structure B5NormalizedRelationDeviationBridge (v : Fin 13 → ℤ) where
  q : ℕ
  one_lt_q : 1 < q
  singletonDefect : ℝ
  mass2 : ℝ
  mass3 : ℝ
  mass4 : ℝ
  mass5 : ℝ
  singleton_id :
    -((deviationLayer v q 1 : ℚ) : ℝ) =
      ((q : ℝ) - 1) * singletonDefect
  pair_id :
    ((deviationLayer v q 2 : ℚ) : ℝ) =
      ((q : ℝ) - 1) * pairWeight * mass2
  triple_id :
    ((deviationLayer v q 3 : ℚ) : ℝ) =
      ((q : ℝ) - 1) * tripleWeight * mass3
  quad_id :
    ((deviationLayer v q 4 : ℚ) : ℝ) =
      -((q : ℝ) - 1) * quadWeight * mass4
  quint_id :
    ((deviationLayer v q 5 : ℚ) : ℝ) =
      ((q : ℝ) - 1) * mass5
  pair_mass_lower_bound : -(13 / 30 : ℝ) ≤ mass2
  harmful_with_singleton_budget :
    harmfulHigherContribution mass3 mass4 mass5 - singletonDefect < 7712 / 84035

/-- The five support identities recover the concrete `B5`, including the
finite-grid singleton defect. -/
theorem B5NormalizedRelationDeviationBridge.b5_eq_scaled_model_add_singleton
    {v : Fin 13 → ℤ}
    (certificate : B5NormalizedRelationDeviationBridge v) :
    (B5 v certificate.q : ℝ) =
      ((certificate.q : ℝ) - 1) *
        (relationModel certificate.mass2 certificate.mass3
          certificate.mass4 certificate.mass5 + certificate.singletonDefect) := by
  have hidQ := B5_eq_equilibrium_add_signedDeviationLedger v certificate.q
  have hidCast := congrArg (fun x : ℚ => (x : ℝ)) hidQ
  have hidR : (B5 v certificate.q : ℝ) =
      ((certificate.q : ℝ) - 1) * (2052 / 16807 : ℝ) +
        (signedDeviationLedger v certificate.q : ℝ) := by
    norm_num at hidCast ⊢
    exact hidCast
  have hexpandQ := signedDeviationLedger_expand v certificate.q
  have hexpandR : (signedDeviationLedger v certificate.q : ℝ) =
      (deviationLayer v certificate.q 0 : ℝ) -
      (deviationLayer v certificate.q 1 : ℝ) +
      (deviationLayer v certificate.q 2 : ℝ) -
      (deviationLayer v certificate.q 3 : ℝ) +
      (deviationLayer v certificate.q 4 : ℝ) -
      (deviationLayer v certificate.q 5 : ℝ) := by
    exact_mod_cast hexpandQ
  have hzeroQ := deviationLayer_zero v certificate.q
    (le_of_lt certificate.one_lt_q)
  have hzeroR : (deviationLayer v certificate.q 0 : ℝ) = 0 := by
    exact_mod_cast hzeroQ
  rw [hidR, hexpandR, hzeroR, zero_sub, certificate.singleton_id,
    certificate.pair_id, certificate.triple_id, certificate.quad_id,
    certificate.quint_id]
  dsimp [relationModel, equilibrium]
  ring

/-- The signed horizon-thirty arithmetic remains valid after explicitly
charging the grid defect. -/
theorem B5NormalizedRelationDeviationBridge.b5_pos
    {v : Fin 13 → ℤ}
    (certificate : B5NormalizedRelationDeviationBridge v) :
    0 < B5 v certificate.q := by
  have hpairWeight : (0 : ℝ) ≤ pairWeight := weights_nonnegative.1
  have hpairWeighted : -(52 / 1715 : ℝ) ≤
      pairWeight * certificate.mass2 := by
    have hmul := mul_le_mul_of_nonneg_left
      certificate.pair_mass_lower_bound hpairWeight
    norm_num [pairWeight] at hmul ⊢
    exact hmul
  have hmodelGrid : 0 <
      relationModel certificate.mass2 certificate.mass3
        certificate.mass4 certificate.mass5 + certificate.singletonDefect := by
    have hbudget := certificate.harmful_with_singleton_budget
    dsimp [relationModel, harmfulHigherContribution, equilibrium, pairWeight,
      tripleWeight, quadWeight] at hpairWeighted hbudget ⊢
    linarith
  have hscale : (0 : ℝ) < (certificate.q : ℝ) - 1 := by
    have hqR : (1 : ℝ) < certificate.q := by
      exact_mod_cast certificate.one_lt_q
    linarith
  have hreal : (0 : ℝ) < (B5 v certificate.q : ℝ) := by
    rw [certificate.b5_eq_scaled_model_add_singleton]
    exact mul_pos hscale hmodelGrid
  exact_mod_cast hreal

/-! ## Axiom audit -/

#print axioms B5_eq_equilibrium_add_signedDeviationLedger
#print axioms B5_pos_of_signedDeviationLedger_lower
#print axioms deviationLayer_zero
#print axioms signedDeviationLedger_expand
#print axioms B5NormalizedRelationDeviationBridge.b5_eq_scaled_model_add_singleton
#print axioms B5NormalizedRelationDeviationBridge.b5_pos

end LRC14Concrete
end LonelyRunner
