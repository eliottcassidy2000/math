/-
  TournamentH7.LRCDenseCoreEndgame

  Orthogonal composition of the two sharpest reductions on the primitive
  LRC(14) residual:

  * the detuned peel removes every family with two or three nonmultiples;
  * the chain-split dichotomy removes every family with a citable ratio-3 tail.

  Consequently positive B5 is needed only on primitive, dissociated families
  whose sorted absolute speeds lie in `ChainDenseCore`.  Loneliness from the
  chain branch is used directly; it is not incorrectly converted into a B5
  certificate.
-/

import TournamentH7.LRCEndgameParameterDischargeFour
import TournamentH7.LRCChainDichotomy

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

/-- Positive depth-five Bonferroni supply on the intersection of the primitive
dissociated residual and the explicit chain-dense core.  This is strictly
weaker than `DissociatedB5Supply`: the supplier may use the sorted dense-core
certificate. -/
def DenseCoreDissociatedB5Supply : Prop :=
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
    ∃ q : ℕ, 0 < q ∧ 0 < LRC14Concrete.B5 v q

/-- The older dissociated B5 supplier implies the chain-dense supplier by
forgetting the extra certificate. -/
theorem denseCoreDissociatedB5Supply_of_dissociated
    (hB5 : DissociatedB5Supply) : DenseCoreDissociatedB5Supply := by
  intro v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc _hcore
  exact hB5 v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdissoc

/-- The primitive residual follows from the multi-detuned dispatch and B5 only
on the chain-dense dissociated core.  The order matters: detuned families are
peeled first, then the sorted chain dichotomy either supplies loneliness
directly or hands its dense-core certificate to the B5 supplier. -/
theorem residualObligationPrimitive_of_denseCore_dissoc_B5
    (cite : LRCUpTo13) (hMD : MultiDetunedDispatch)
    (hB5 : DenseCoreDissociatedB5Supply) : ResidualObligationPrimitive := by
  intro v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse _hcres
  by_cases hdet : ∃ g : ℤ, 2 ≤ g ∧
      (nonMultCard v g = 2 ∨ nonMultCard v g = 3)
  · exact hMD v hv hdet
  · push_neg at hdet
    set va : Fin 13 → ℤ := fun i => |v i| with hva
    have hva_pos : ∀ i, 0 < va i := fun i => abs_pos.mpr (hv i)
    set σ : Equiv.Perm (Fin 13) := Tuple.sort va with hσ
    set w : Fin 13 → ℤ := va ∘ σ with hw
    have hw_mono : Monotone w := Tuple.monotone_sort va
    have hw_pos : ∀ i, 0 < w i := fun i => hva_pos (σ i)
    rcases lonely_or_denseCore cite w hw_pos hw_mono with ⟨t, ht⟩ | hcore
    · refine ⟨t, (LRC14.lonely_abs_iff 14 v t).mp ?_⟩
      exact (LRC14.lonely_comp_equiv 14 va t σ).mp ht
    · obtain ⟨q, hq, hB5pos⟩ :=
        hB5 v hv hgcd hcov hgap hcomp hdist hlarge hdiv hcoarse hdet
          ⟨σ, hw_mono, hcore⟩
      exact LRC14Concrete.lonely_of_Mreach_ge v hv
        (LRC14Concrete.mreach_ge_of_B5_pos v q hq hB5pos)

/-- **Current narrowest machine-checked endgame surface.**  The sanctioned
LRC(≤13) citation, the deep exceptional detuned dispatch, and positive B5 only
on the primitive dissociated chain-dense core imply LRC(14). -/
theorem lrc14_from_deep_detuned_and_denseCore_dissociated_B5
    (cite : LRCUpTo13) (hdeep : DeepExceptionalDetunedDispatch)
    (hB5 : DenseCoreDissociatedB5Supply) : LRC14.LRC14Statement :=
  lrc14_grand_assembly_primitive cite
    (residualObligationPrimitive_of_denseCore_dissoc_B5 cite
      (multiDetunedDispatch_of_exceptional cite
        (exceptionalDetunedDispatch_of_deep cite hdeep))
      hB5)

/-- **Sharpened narrowest machine-checked endgame surface.**  After the exact
`q ≥ 4` bad-count arithmetic, exceptional triples need dispatch only when one
reduced denominator is `2` or `3`; positive B5 is still required only on the
primitive dissociated chain-dense core. -/
theorem lrc14_from_four_detuned_and_denseCore_dissociated_B5
    (cite : LRCUpTo13) (hdeep : DeepExceptionalDetunedDispatchFour)
    (hB5 : DenseCoreDissociatedB5Supply) : LRC14.LRC14Statement :=
  lrc14_from_deep_detuned_and_denseCore_dissociated_B5 cite
    (deepExceptionalDetunedDispatch_of_four cite hdeep) hB5

/-! ## Axiom audit -/

#print axioms denseCoreDissociatedB5Supply_of_dissociated
#print axioms residualObligationPrimitive_of_denseCore_dissoc_B5
#print axioms lrc14_from_deep_detuned_and_denseCore_dissociated_B5
#print axioms lrc14_from_four_detuned_and_denseCore_dissociated_B5

end LRC14Grand
end LonelyRunner
