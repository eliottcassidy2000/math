/-
  Shared harmonic-phase facts for the three selected-witness residues.

  The sanctioned LRC(<=13) citation supplies a phase clearing the ten harmonic
  quotients by `1/11`, hence by `1/14`.  The missing content of each supplier is
  exactly that this nonempty harmonic-good set is not contained in its explicit
  parallel-class failure locus.  The *strict* periodic unit budget does not
  close ten distinct frequencies: at the target gap its inequality handles at
  most six.  Compactness at the exact endpoint does handle seven distinct
  absolute frequency classes; that sharpened quotient statement is proved in
  `LRCSelectedWitnessFrequencyEndpoint`.

  Assumption challenge: harmonic quotient frequencies, rather than runners or
  branch arcs, are the vertices preserved by this reduction.  Their danger-set
  intersection preserves harmonic clearance but forgets which detuned branch
  survives.  A tournament gauge on frequency labels would add no information,
  so no artificial tie path is introduced.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCEndgameTwoEight

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

noncomputable section

/-- The citation supplies a harmonic-good phase for every nonzero harmonic
quotient family.  It does not select that phase relative to the detuned
parallel-class obstruction. -/
theorem exists_threeDetunedHarmonicGoodAt_of_cite
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (i₁ i₂ i₃ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j) :
    ∃ u : ℝ, ThreeDetunedHarmonicGoodAt v g i₁ i₂ i₃ u := by
  have hselectedCard : ({i₁, i₂, i₃} : Finset (Fin 13)).card = 3 :=
    Finset.card_eq_three.mpr ⟨i₁, i₂, i₃, h12, h13, h23, rfl⟩
  have hcard :
      (Finset.univ \ ({i₁, i₂, i₃} : Finset (Fin 13))).card = 10 := by
    rw [Finset.card_sdiff, Finset.inter_univ, hselectedCard, Finset.card_univ,
      Fintype.card_fin]
  set embedding : Fin 10 ↪o Fin 13 :=
    (Finset.univ \ ({i₁, i₂, i₃} : Finset (Fin 13))).orderEmbOfFin hcard
      with hembedding
  have hembeddingMem : ∀ index,
      embedding index ∈ Finset.univ \ ({i₁, i₂, i₃} : Finset (Fin 13)) :=
    fun index =>
      (Finset.univ \ ({i₁, i₂, i₃} : Finset (Fin 13))).orderEmbOfFin_mem
        hcard index
  have hembeddingNe : ∀ index,
      embedding index ≠ i₁ ∧ embedding index ≠ i₂ ∧ embedding index ≠ i₃ := by
    intro index
    have hmem := hembeddingMem index
    rw [Finset.mem_sdiff, Finset.mem_insert, Finset.mem_insert,
      Finset.mem_singleton] at hmem
    push Not at hmem
    exact hmem.2
  set quotient : Fin 10 → ℤ := fun index => v (embedding index) / g
    with hquotient
  have hquotientNonzero : ∀ index, quotient index ≠ 0 := by
    intro index hzero
    exact hv (embedding index) (by
      have hdivides := hdvd (embedding index) (hembeddingNe index).1
        (hembeddingNe index).2.1 (hembeddingNe index).2.2
      have hfactor : v (embedding index) = g * quotient index :=
        (Int.mul_ediv_cancel' hdivides).symm
      rw [hfactor, hzero, mul_zero])
  obtain ⟨u, hu⟩ := cite 10 (by norm_num) quotient hquotientNonzero
  refine ⟨u, ?_⟩
  intro index hindex1 hindex2 hindex3 integer
  have hindexMem :
      index ∈ Finset.univ \ ({i₁, i₂, i₃} : Finset (Fin 13)) := by
    rw [Finset.mem_sdiff, Finset.mem_insert, Finset.mem_insert,
      Finset.mem_singleton]
    exact ⟨Finset.mem_univ index, by
      push Not
      exact ⟨hindex1, hindex2, hindex3⟩⟩
  obtain ⟨source, hsource⟩ : ∃ source, embedding source = index := by
    have hrange : index ∈ Set.range embedding := by
      rw [hembedding, Finset.range_orderEmbOfFin]
      exact hindexMem
    exact hrange
  have hquotientEq : quotient source = v index / g := by
    show v (embedding source) / g = v index / g
    rw [hsource]
  calc
    (1 : ℝ) / 14 ≤ 1 / 11 := by norm_num
    _ ≤ |(quotient source : ℝ) * u - integer| := hu source integer
    _ = |((v index / g : ℤ) : ℝ) * u - integer| := by rw [hquotientEq]

/-- The exact extra selector is that the harmonic-good phase set is not
contained in the detuned failure set. -/
theorem selectedWitness_iff_harmonicGood_not_subset_failure
    (v : Fin 13 → ℤ) (g : ℤ) (i₁ i₂ i₃ : Fin 13) :
    (∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₁ i₂ i₃ u ∧
      HasThreeDetunedGoodBranch (v i₁) (v i₂) (v i₃) g u) ↔
    ¬ (∀ u : ℝ, ThreeDetunedHarmonicGoodAt v g i₁ i₂ i₃ u →
      ¬ HasThreeDetunedGoodBranch (v i₁) (v i₂) (v i₃) g u) := by
  constructor
  · rintro ⟨u, hharmonic, hgood⟩ hall
    exact hall u hharmonic hgood
  · intro hescape
    push Not at hescape
    exact hescape

/-- Any actual LRC(14) witness canonically produces a selected harmonic phase:
take `u = g*t` and the zero branch. -/
theorem selectedWitness_of_lonely
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g) (i₁ i₂ i₃ : Fin 13)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j)
    (hlonely : ∃ t : ℝ, Lonely 14 v t) :
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₁ i₂ i₃ u ∧
      HasThreeDetunedGoodBranch (v i₁) (v i₂) (v i₃) g u := by
  obtain ⟨t, ht⟩ := hlonely
  let u : ℝ := (g : ℝ) * t
  have hgReal : (g : ℝ) ≠ 0 := by
    exact_mod_cast (show g ≠ 0 by omega)
  have hharmonic : ThreeDetunedHarmonicGoodAt v g i₁ i₂ i₃ u := by
    intro index hindex1 hindex2 hindex3 integer
    have hfactorInt : v index = g * (v index / g) :=
      (Int.mul_ediv_cancel' (hdvd index hindex1 hindex2 hindex3)).symm
    have hfactorReal : (v index : ℝ) =
        (g : ℝ) * ((v index / g : ℤ) : ℝ) := by
      exact_mod_cast hfactorInt
    have heq : ((v index / g : ℤ) : ℝ) * u - integer =
        (v index : ℝ) * t - integer := by
      dsimp [u]
      rw [hfactorReal]
      ring
    rw [heq]
    exact ht index integer
  have hphase (index : Fin 13) :
      (v index : ℝ) * ((u + ((0 : ℤ) : ℝ)) / (g : ℝ)) =
        (v index : ℝ) * t := by
    dsimp [u]
    field_simp
    ring
  have hgood : HasThreeDetunedGoodBranch
      (v i₁) (v i₂) (v i₃) g u := by
    refine ⟨0, Finset.mem_Ico.mpr ⟨by omega, by omega⟩, ?_, ?_, ?_⟩
    · intro hbad
      obtain ⟨integer, hnear⟩ := (Finset.mem_filter.mp hbad).2
      apply (not_lt_of_ge (ht i₁ integer))
      rwa [hphase i₁] at hnear
    · intro hbad
      obtain ⟨integer, hnear⟩ := (Finset.mem_filter.mp hbad).2
      apply (not_lt_of_ge (ht i₂ integer))
      rwa [hphase i₂] at hnear
    · intro hbad
      obtain ⟨integer, hnear⟩ := (Finset.mem_filter.mp hbad).2
      apply (not_lt_of_ge (ht i₃ integer))
      rwa [hphase i₃] at hnear
  exact ⟨u, hharmonic, hgood⟩

/-- For a fixed three-detuned decomposition, selected-witness existence is
exactly equivalent to the LRC(14) conclusion.  It is therefore an exact
residual restatement, not a weaker consequence of the current unit budget. -/
theorem selectedWitness_iff_exists_lonely
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g) (i₁ i₂ i₃ : Fin 13)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j) :
    (∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₁ i₂ i₃ u ∧
      HasThreeDetunedGoodBranch (v i₁) (v i₂) (v i₃) g u) ↔
    ∃ t : ℝ, Lonely 14 v t := by
  constructor
  · exact lonely14_of_three_detuned_selectedWitness
      v g hg i₁ i₂ i₃ hdvd
  · exact selectedWitness_of_lonely v g hg i₁ i₂ i₃ hdvd

/-- The strict periodic unit-budget inequality cannot directly supply ten
*distinct* harmonic frequencies at the target gap: its raw ten-label budget is
already supercritical.  Repeated or sign-opposite quotients must first be
collapsed to absolute frequency classes. -/
theorem unitBudget_target_gap_ten_arithmetic :
    ¬ (2 * ((1 : ℝ) / 14) * 10 < 1) := by
  norm_num

/-- At most six frequencies fit the *strict* unit-budget inequality at gap
`1/14`.  This does not contradict the compactness endpoint for seven distinct
absolute frequencies, where equality `2 * (1/14) * 7 = 1` is retained as a
limit. -/
theorem unitBudget_target_gap_card_bound {cardinality : ℕ}
    (hbudget : 2 * ((1 : ℝ) / 14) * cardinality < 1) :
    cardinality ≤ 6 := by
  norm_num at hbudget
  have hreal : (cardinality : ℝ) < 7 := by
    nlinarith [hbudget]
  have hlt : cardinality < 7 := by
    exact_mod_cast hreal
  omega

/-! ## Axiom audit -/

#print axioms exists_threeDetunedHarmonicGoodAt_of_cite
#print axioms selectedWitness_iff_harmonicGood_not_subset_failure
#print axioms selectedWitness_iff_exists_lonely
#print axioms unitBudget_target_gap_card_bound

end

end LRC14Grand
end LonelyRunner
