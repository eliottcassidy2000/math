/-
  TournamentH7.LRCEndgameParameterDischargeFour

  Sharpens the small-denominator branch of
  `DeepExceptionalDetunedDispatch`.  The earlier coarse corollary used the
  sufficient threshold `q = g / gcd(delta,g) >= 8`.  The exact integer count

      badCount(delta,g) = gcd(delta,g) * (floor(q / 7) + 1)

  satisfies the stronger uniform estimate

      7 * badCount(delta,g) <= 2 * g

  already for every `q >= 4`.  Hence three detuned coordinates with all
  reduced denominators at least four have total bad count at most `6g/7 < g`
  and are handled by the proved generic three-detuned dispatch.

  The resulting named residual contains only:

  * the nonterminating two-adic pair tower; or
  * a non-generic triple with at least one reduced denominator in `{2,3}`.

  Assumption challenge / quotient audit: the finite vertices used here are the
  three detuned coordinates, not all runners or circle arcs.  Quotienting a
  coordinate to its reduced denominator preserves the normalized bad-count
  upper bound needed by `genericCount`, but destroys its phase and identity;
  therefore the quotient is used only to certify the counting inequality.
  The challenged assumption is precisely that the prior threshold `8` was a
  structural boundary: exact floor arithmetic lowers it to `4`.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCEndgameParameterDischarge

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

/-- A detuned coordinate has reduced denominator strictly below four. -/
def HasSubFourDetuningDenominator (v : Fin 13 → ℤ) (g : ℤ) : Prop :=
  ∃ i, ¬ g ∣ v i ∧ g / (Int.gcd (v i) g : ℤ) < 4

/-- A nonmultiple has reduced denominator at least two. -/
theorem reducedDetuningDenominator_ge_two {δ g : ℤ} (hg : 2 ≤ g)
    (hδ : ¬ g ∣ δ) :
    2 ≤ g / (Int.gcd δ g : ℤ) := by
  have hg0 : (0 : ℤ) < g := by omega
  have hdvdg : ((Int.gcd δ g : ℤ)) ∣ g := Int.gcd_dvd_right δ g
  have hdpos : (0 : ℤ) < (Int.gcd δ g : ℤ) := by
    have : 0 < Int.gcd δ g := by
      rw [Int.gcd_pos_iff]
      right
      omega
    exact_mod_cast this
  have hqnonneg : (0 : ℤ) ≤ g / (Int.gcd δ g : ℤ) :=
    Int.ediv_nonneg (le_of_lt hg0) (le_of_lt hdpos)
  have hfactor : (Int.gcd δ g : ℤ) * (g / (Int.gcd δ g : ℤ)) = g :=
    Int.mul_ediv_cancel' hdvdg
  have hqzero : g / (Int.gcd δ g : ℤ) ≠ 0 := by
    intro hzero
    rw [hzero, mul_zero] at hfactor
    omega
  have hqone : g / (Int.gcd δ g : ℤ) ≠ 1 := by
    intro hone
    have hgcd : (Int.gcd δ g : ℤ) = g := by
      rw [hone, mul_one] at hfactor
      exact hfactor
    have hdvdδ : ((Int.gcd δ g : ℤ)) ∣ δ := Int.gcd_dvd_left δ g
    exact hδ (by simpa [hgcd] using hdvdδ)
  omega

/-- Thus the sub-four condition is exactly the finite denominator alphabet
`{2,3}`. -/
theorem hasSubFourDetuningDenominator_iff_two_or_three
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g) :
    HasSubFourDetuningDenominator v g ↔
      ∃ i, ¬ g ∣ v i ∧
        (g / (Int.gcd (v i) g : ℤ) = 2 ∨
          g / (Int.gcd (v i) g : ℤ) = 3) := by
  constructor
  · rintro ⟨i, hδ, hlt⟩
    have hge := reducedDetuningDenominator_ge_two hg hδ
    exact ⟨i, hδ, by omega⟩
  · rintro ⟨i, hδ, hq⟩
    exact ⟨i, hδ, by omega⟩

private theorem gcd_toNat_mul_reducedDenominator (δ g : ℤ) (hg : 0 < g) :
    g.toNat = (Int.gcd δ g) * (g / (Int.gcd δ g : ℤ)).toNat := by
  have hdvd : ((Int.gcd δ g : ℤ)) ∣ g := Int.gcd_dvd_right δ g
  have hdz : (0 : ℤ) < (Int.gcd δ g : ℤ) := by
    have : 0 < Int.gcd δ g := by
      rw [Int.gcd_pos_iff]
      right
      omega
    exact_mod_cast this
  have hqnn : (0 : ℤ) ≤ g / (Int.gcd δ g : ℤ) :=
    Int.ediv_nonneg (le_of_lt hg) (le_of_lt hdz)
  have heq : g = (Int.gcd δ g : ℤ) * (g / (Int.gcd δ g : ℤ)) :=
    (Int.mul_ediv_cancel' hdvd).symm
  have hcast : (g.toNat : ℤ) =
      ((Int.gcd δ g) * (g / (Int.gcd δ g : ℤ)).toNat : ℕ) := by
    rw [Int.toNat_of_nonneg (le_of_lt hg)]
    push_cast
    rw [Int.toNat_of_nonneg hqnn]
    exact heq
  exact_mod_cast hcast

/-- Exact-floor improvement of the old all-fine estimate: reduced denominator
`q >= 4` already gives `badCount / g <= 2/7`. -/
theorem seven_mul_badCount_le_two_mul (δ g : ℤ) (hg : 0 < g)
    (hq : 4 ≤ g / (Int.gcd δ g : ℤ)) :
    7 * DetunedD3.badCount δ g ≤ 2 * g.toNat := by
  have hQ4 : 4 ≤ (g / (Int.gcd δ g : ℤ)).toNat := by omega
  have hG : g.toNat =
      (Int.gcd δ g) * (g / (Int.gcd δ g : ℤ)).toNat :=
    gcd_toNat_mul_reducedDenominator δ g hg
  have hseven :
      7 * ((g / (Int.gcd δ g : ℤ)).toNat / 7 + 1) ≤
        2 * (g / (Int.gcd δ g : ℤ)).toNat := by
    omega
  calc
    7 * DetunedD3.badCount δ g =
        (Int.gcd δ g) *
          (7 * ((g / (Int.gcd δ g : ℤ)).toNat / 7 + 1)) := by
            rw [DetunedD3.badCount]
            ring
    _ ≤ (Int.gcd δ g) *
        (2 * (g / (Int.gcd δ g : ℤ)).toNat) :=
      Nat.mul_le_mul_left _ hseven
    _ = 2 * ((Int.gcd δ g) *
        (g / (Int.gcd δ g : ℤ)).toNat) := by ring
    _ = 2 * g.toNat := by rw [hG]

/-- With exactly three detuned coordinates, absence of a sub-four denominator
forces the generic bad-count inequality. -/
theorem genericCount_of_three_noSubFour
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g)
    (hcard : nonMultCard v g = 3)
    (hsub : ¬ HasSubFourDetuningDenominator v g) :
    genericCount v g := by
  have hg0 : (0 : ℤ) < g := by omega
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
  have hq1 : 4 ≤ g / (Int.gcd (v i₁) g : ℤ) := by
    by_contra hnot
    exact hsub ⟨i₁, hδ1, lt_of_not_ge hnot⟩
  have hq2 : 4 ≤ g / (Int.gcd (v i₂) g : ℤ) := by
    by_contra hnot
    exact hsub ⟨i₂, hδ2, lt_of_not_ge hnot⟩
  have hq3 : 4 ≤ g / (Int.gcd (v i₃) g : ℤ) := by
    by_contra hnot
    exact hsub ⟨i₃, hδ3, lt_of_not_ge hnot⟩
  have hb1 := seven_mul_badCount_le_two_mul (v i₁) g hg0 hq1
  have hb2 := seven_mul_badCount_le_two_mul (v i₂) g hg0 hq2
  have hb3 := seven_mul_badCount_le_two_mul (v i₃) g hg0 hq3
  rw [genericCount, hfilter,
    Finset.sum_insert (by simp [h12, h13]),
    Finset.sum_insert (by simp [h23]), Finset.sum_singleton]
  omega

/-- The newly closed honest subcase: every exactly-three-detuned family whose
three reduced denominators are at least four is lonely. -/
theorem lonely14_of_nonMultCard_three_noSubFour (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (hcard : nonMultCard v g = 3)
    (hsub : ¬ HasSubFourDetuningDenominator v g) :
    ∃ t : ℝ, Lonely 14 v t :=
  lonely14_of_nonMultCard_three cite v hv g hg hcard
    (genericCount_of_three_noSubFour v g hg hcard hsub)

/-- The sharpened deep residual: a nonterminating pair tower, or an exceptional
triple carrying a reduced denominator in the exact alphabet `{2,3}`. -/
def DeepExceptionalDetunedDispatchFour : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → ∀ g : ℤ, 2 ≤ g →
    ((nonMultCard v g = 2 ∧ ¬ TwoAdicLiftTerminates v g) ∨
      (nonMultCard v g = 3 ∧ HasSubFourDetuningDenominator v g)) →
    ¬ genericCount v g →
    ∃ t : ℝ, Lonely 14 v t

/-- The sub-four residual supplies the previous sub-eight residual because all
triples with denominators at least four are now theorems. -/
theorem deepExceptionalDetunedDispatch_of_four (cite : LRCUpTo13)
    (hdeep : DeepExceptionalDetunedDispatchFour) :
    DeepExceptionalDetunedDispatch := by
  intro v hv g hg hcase hnongeneric
  rcases hcase with htwo | ⟨hthree, _⟩
  · exact hdeep v hv g hg (Or.inl htwo) hnongeneric
  · by_cases hsub : HasSubFourDetuningDenominator v g
    · exact hdeep v hv g hg (Or.inr ⟨hthree, hsub⟩) hnongeneric
    · exact lonely14_of_nonMultCard_three_noSubFour
        cite v hv g hg hthree hsub

/-- Sharpened current endgame surface: only the nonterminating pair tower, a
triple with a `q=2` or `q=3` coordinate, and dissociated positive-B5 supply
remain as named mathematical parameters. -/
theorem lrc14_from_four_detuned_and_dissociated_B5 (cite : LRCUpTo13)
    (hdeep : DeepExceptionalDetunedDispatchFour)
    (hB5 : DissociatedB5Supply) :
    LRC14.LRC14Statement :=
  lrc14_from_deep_detuned_and_dissociated_B5 cite
    (deepExceptionalDetunedDispatch_of_four cite hdeep) hB5

/-! ## Axiom audit -/

#print axioms reducedDetuningDenominator_ge_two
#print axioms hasSubFourDetuningDenominator_iff_two_or_three
#print axioms seven_mul_badCount_le_two_mul
#print axioms genericCount_of_three_noSubFour
#print axioms lonely14_of_nonMultCard_three_noSubFour
#print axioms deepExceptionalDetunedDispatch_of_four
#print axioms lrc14_from_four_detuned_and_dissociated_B5

end LRC14Grand
end LonelyRunner
