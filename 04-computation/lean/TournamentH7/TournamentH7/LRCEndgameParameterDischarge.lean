/-
  TournamentH7.LRCEndgameParameterDischarge

  Maximal currently proved composition of the primitive/dissociated LRC(14)
  endgame.  The generic two- and three-detuned branches were already closed;
  this module also removes two concrete pieces from their exceptional residue:

  * a two-detuned branch whose 2-adic lift terminates immediately;
  * a three-detuned branch whose three reduced denominators are all at least 8.

  What remains is named exactly: the nonterminating two-adic tower or a
  small-denominator exceptional triple, together with positive B5 supply on
  the dissociated primitive residual.  No mathematical residual is hidden in
  assembly glue.
-/

import TournamentH7.LRCDetunedDispatchReduce
import TournamentH7.LRCTwoDetunedLift
import TournamentH7.LRCThreeDetunedCoarse

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

/-- The first 2-adic lift terminates: every speed divisible by `g` is already
divisible by `2 * g`. -/
def TwoAdicLiftTerminates (v : Fin 13 → ℤ) (g : ℤ) : Prop :=
  ∀ i, g ∣ v i → 2 * g ∣ v i

/-- A detuned coordinate has reduced denominator below the all-fine threshold
`8` used by `lonely14_of_three_detuned_coarse`. -/
def HasSmallDetuningDenominator (v : Fin 13 → ℤ) (g : ℤ) : Prop :=
  ∃ i, ¬ g ∣ v i ∧ g / (Int.gcd (v i) g : ℤ) < 8

/-- The genuinely deep part of `ExceptionalDetunedDispatch`: either a
nonterminating two-coordinate 2-adic tower, or a three-coordinate exceptional
branch containing a reduced denominator below `8`. -/
def DeepExceptionalDetunedDispatch : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → ∀ g : ℤ, 2 ≤ g →
    ((nonMultCard v g = 2 ∧ ¬ TwoAdicLiftTerminates v g) ∨
      (nonMultCard v g = 3 ∧ HasSmallDetuningDenominator v g)) →
    ¬ genericCount v g →
    ∃ t : ℝ, Lonely 14 v t

/-- Positive depth-five Bonferroni supply only on the primitive dissociated
residual, the sharpest B5 surface currently consumed by the assembly. -/
def DissociatedB5Supply : Prop :=
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
    ∃ q : ℕ, 0 < q ∧ 0 < LRC14Concrete.B5 v q

/-- Exactly two nonmultiples, with no odd harmonic multipliers, are handled by
the proved mod-`2g` lift. -/
theorem lonely14_of_nonMultCard_two_liftTerminates (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (hcard : nonMultCard v g = 2) (hterm : TwoAdicLiftTerminates v g) :
    ∃ t : ℝ, Lonely 14 v t := by
  rw [nonMultCard] at hcard
  obtain ⟨i₁, i₂, h12, hfilter⟩ := Finset.card_eq_two.mp hcard
  have hmem : ∀ j,
      j ∈ Finset.univ.filter (fun i => ¬ g ∣ v i) ↔ ¬ g ∣ v j := by
    intro j
    simp
  have hδ1 : ¬ g ∣ v i₁ := (hmem i₁).mp (by rw [hfilter]; simp)
  have hδ2 : ¬ g ∣ v i₂ := (hmem i₂).mp (by rw [hfilter]; simp)
  have hlift : ∀ j, j ≠ i₁ → j ≠ i₂ → 2 * g ∣ v j := by
    intro j hj1 hj2
    apply hterm
    by_contra hnot
    have hj : j ∈ Finset.univ.filter (fun i => ¬ g ∣ v i) := (hmem j).mpr hnot
    rw [hfilter] at hj
    simp only [Finset.mem_insert, Finset.mem_singleton] at hj
    rcases hj with h | h
    · exact hj1 h
    · exact hj2 h
  exact DetunedD2.lonely14_of_two_detuned_lift2g
    cite v hv g hg i₁ i₂ h12 hlift hδ1 hδ2

/-- Exactly three nonmultiples whose reduced denominators are all at least `8`
are handled by the proved all-fine coarse dispatch. -/
theorem lonely14_of_nonMultCard_three_noSmallDenominator (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (hcard : nonMultCard v g = 3)
    (hsmall : ¬ HasSmallDetuningDenominator v g) :
    ∃ t : ℝ, Lonely 14 v t := by
  rw [nonMultCard] at hcard
  obtain ⟨i₁, i₂, i₃, h12, h13, h23, hfilter⟩ :=
    Finset.card_eq_three.mp hcard
  have hmem : ∀ j,
      j ∈ Finset.univ.filter (fun i => ¬ g ∣ v i) ↔ ¬ g ∣ v j := by
    intro j
    simp
  have hδ1 : ¬ g ∣ v i₁ := (hmem i₁).mp (by rw [hfilter]; simp)
  have hδ2 : ¬ g ∣ v i₂ := (hmem i₂).mp (by rw [hfilter]; simp)
  have hδ3 : ¬ g ∣ v i₃ := (hmem i₃).mp (by rw [hfilter]; simp)
  have hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j := by
    intro j hj1 hj2 hj3
    by_contra hnot
    have hj : j ∈ Finset.univ.filter (fun i => ¬ g ∣ v i) := (hmem j).mpr hnot
    rw [hfilter] at hj
    simp only [Finset.mem_insert, Finset.mem_singleton] at hj
    rcases hj with h | h | h
    · exact hj1 h
    · exact hj2 h
    · exact hj3 h
  have hq1 : 8 ≤ g / (Int.gcd (v i₁) g : ℤ) := by
    by_contra hnot
    exact hsmall ⟨i₁, hδ1, lt_of_not_ge hnot⟩
  have hq2 : 8 ≤ g / (Int.gcd (v i₂) g : ℤ) := by
    by_contra hnot
    exact hsmall ⟨i₂, hδ2, lt_of_not_ge hnot⟩
  have hq3 : 8 ≤ g / (Int.gcd (v i₃) g : ℤ) := by
    by_contra hnot
    exact hsmall ⟨i₃, hδ3, lt_of_not_ge hnot⟩
  exact DetunedD3.lonely14_of_three_detuned_coarse cite v hv g hg
    i₁ i₂ i₃ h12 h13 h23 hdvd hδ1 hδ2 hδ3 hq1 hq2 hq3

/-- The full exceptional detuned dispatch follows from only its deep residual:
the terminating pair and all-fine triple branches are theorems. -/
theorem exceptionalDetunedDispatch_of_deep (cite : LRCUpTo13)
    (hdeep : DeepExceptionalDetunedDispatch) : ExceptionalDetunedDispatch := by
  intro v hv g hg hcard hnongeneric
  rcases hcard with htwo | hthree
  · by_cases hterm : TwoAdicLiftTerminates v g
    · exact lonely14_of_nonMultCard_two_liftTerminates cite v hv g hg htwo hterm
    · exact hdeep v hv g hg (Or.inl ⟨htwo, hterm⟩) hnongeneric
  · by_cases hsmall : HasSmallDetuningDenominator v g
    · exact hdeep v hv g hg (Or.inr ⟨hthree, hsmall⟩) hnongeneric
    · exact lonely14_of_nonMultCard_three_noSmallDenominator
        cite v hv g hg hthree hsmall

/-- **Current sharp endgame surface.**  The sanctioned LRC(≤13) citation, deep
exceptional detuned dispatch, and dissociated positive-B5 supply imply LRC(14).
Every intervening primitive peel, generic detuned branch, terminating lift,
all-fine triple, B5 consumer, and witness-attainment step is proved. -/
theorem lrc14_from_deep_detuned_and_dissociated_B5 (cite : LRCUpTo13)
    (hdeep : DeepExceptionalDetunedDispatch) (hB5 : DissociatedB5Supply) :
    LRC14.LRC14Statement :=
  lrc14_from_B5_dissoc cite
    (multiDetunedDispatch_of_exceptional cite
      (exceptionalDetunedDispatch_of_deep cite hdeep))
    hB5

/-! ## Axiom audit -/

#print axioms lonely14_of_nonMultCard_two_liftTerminates
#print axioms lonely14_of_nonMultCard_three_noSmallDenominator
#print axioms exceptionalDetunedDispatch_of_deep
#print axioms lrc14_from_deep_detuned_and_dissociated_B5

end LRC14Grand
end LonelyRunner
