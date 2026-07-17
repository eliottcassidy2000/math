/-
  TournamentH7.LRCB5CertificateAudit

  Semantic audit of `B5RelationBudgetCertificate`.

  The four fields named `mass2`, ..., `mass5` are currently unconstrained real
  numbers: the structure does not identify them with the exact-support relation
  masses or with the concrete subset-deviation ledger.  Consequently the signed
  budget is not an additional certificate condition.  Given any positive `B5` at
  a modulus `q > 1`, put the pair mass at its permitted lower endpoint and use the
  unconstrained support-five field to absorb all remaining slack.  This produces
  a `B5RelationBudgetCertificate`; conversely the existing consumer recovers
  positive `B5`.  Thus the current structure is logically equivalent to the raw
  positive-`B5` witness it was intended to refine.

  This audit applies only to the older abstract structure.  The subsequently
  introduced `NormalizedB5RelationBudgetCertificate` is the honest replacement:
  its fields are inequalities on the concretely defined pair-depth and harmful-
  depth moments at `cleanModulus`, so there is no free mass coordinate with which
  to manufacture a witness.

  Assumption/tournament audit: the natural vertices here are proof obligations
  (scaled identity, pair floor, signed higher budget), not runners.  The pairwise
  observable "obligation A logically constrains obligation B" is not a canonical
  orientation: the free `mass5` coordinate makes the three obligations mutually
  satisfiable as soon as `B5 > 0`.  An index-order gauge would manufacture a
  tournament while destroying exactly the missing semantic datum, namely whether
  a mass comes from an actual relation-support stratum.  No tournament fingerprint
  is therefore asserted.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCB5RelationEndgame

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open LRCB5RelationBudget

noncomputable section

/-- Any positive concrete `B5` value at `q > 1` manufactures the current
`B5RelationBudgetCertificate`.  The construction uses no relation-mass data:
`mass2` is fixed at `-13/30`, `mass3 = mass4 = 0`, and the free `mass5`
coordinate absorbs the remaining difference. -/
def b5RelationBudgetCertificate_of_b5_pos
    {v : Fin 13 → ℤ} {q : ℕ} (hq : 1 < q)
    (hb5 : 0 < LRC14Concrete.B5 v q) :
    B5RelationBudgetCertificate v := by
  let scale : ℝ := (q : ℝ) - 1
  let target : ℝ := (LRC14Concrete.B5 v q : ℝ) / scale
  let artificialMass5 : ℝ :=
    equilibrium - pairWeight * (13 / 30 : ℝ) - target
  have hscale : 0 < scale := by
    have hqReal : (1 : ℝ) < q := by exact_mod_cast hq
    dsimp [scale]
    linarith
  have htarget : 0 < target := by
    dsimp [target]
    positivity
  have hscaledTarget :
      (LRC14Concrete.B5 v q : ℝ) = scale * target := by
    dsimp [target]
    field_simp [ne_of_gt hscale]
  have hmodel :
      relationModel (-(13 / 30 : ℝ)) 0 0 artificialMass5 = target := by
    dsimp [relationModel, artificialMass5]
    ring
  have hmass5 : artificialMass5 = 7712 / 84035 - target := by
    dsimp [artificialMass5]
    rw [remaining_after_pair_tail_eq]
  have hharmful :
      harmfulHigherContribution 0 0 artificialMass5 < 7712 / 84035 := by
    rw [show harmfulHigherContribution 0 0 artificialMass5 = artificialMass5 by
      simp [harmfulHigherContribution]]
    rw [hmass5]
    linarith
  refine
    { q := q
      one_lt_q := hq
      mass2 := -(13 / 30 : ℝ)
      mass3 := 0
      mass4 := 0
      mass5 := artificialMass5
      b5_eq_scaled_model := ?_
      pair_mass_lower_bound := le_rfl
      harmful_higher_budget := hharmful }
  rw [hmodel]
  simpa [scale] using hscaledTarget

/-- Exact semantic audit: because the mass fields have no realization
predicate, inhabiting the current relation-budget certificate is equivalent to
exhibiting a positive `B5` value at some modulus greater than one. -/
theorem nonempty_b5RelationBudgetCertificate_iff
    (v : Fin 13 → ℤ) :
    Nonempty (B5RelationBudgetCertificate v) ↔
      ∃ q : ℕ, 1 < q ∧ 0 < LRC14Concrete.B5 v q := by
  constructor
  · rintro ⟨certificate⟩
    exact ⟨certificate.q, certificate.one_lt_q, certificate.b5_pos⟩
  · rintro ⟨q, hq, hb5⟩
    exact ⟨b5RelationBudgetCertificate_of_b5_pos hq hb5⟩

/-! ## Axiom audit -/

#print axioms b5RelationBudgetCertificate_of_b5_pos
#print axioms nonempty_b5RelationBudgetCertificate_iff

end
end LRC14Grand
end LonelyRunner
