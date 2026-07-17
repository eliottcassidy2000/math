/-
  TournamentH7.LRCEndgameParameterDischargeTwoThree

  Sharpens the exactly-three-detuned part of
  `DeepExceptionalDetunedDispatchFour`.  View the `g` branch choices as the
  right vertices of a bipartite incidence graph and the three detuned
  coordinates as the left vertices.  A left vertex of reduced denominator
  `q = 3` has exact bad degree `g / 3`; one with `q >= 4` has bad degree at
  most `2g / 7`.  Consequently the mixed extremal pattern

      (q₁, q₂, q₃) = (3, 3, q),  q >= 4,

  has total bad degree at most `20g / 21 < g`, so its three bad neighborhoods
  cannot cover the parallel-class circle.  The same argument handles every
  pattern without a `q = 2` coordinate except the uniform `(3,3,3)` pattern.

  Thus a non-generic triple has an exact structural normalization:

  * one reduced denominator is two and a distinct one is at most eight; or
  * all three reduced denominators are three.

  The companion cutoff is sharp for this degree-only method: `(2,8,8)` has
  total bad degree exactly `g`, whereas `(2,q₂,q₃)` is generic when both
  companion denominators are at least nine.

  Zarankiewicz audit: a stronger `K_{s,t}` incidence estimate would require
  pairwise bad-neighborhood intersections.  Reduced denominators preserve the
  row degrees used here but destroy the phases controlling those intersections,
  so no unsupported `K_{2,t}`-free claim is made.  The remaining q=2-with-small-
  companion and uniform-q=3 patterns are precisely where row-degree arithmetic
  can saturate and a parity/two-adic or phase-intersection argument is still
  required.

  Assumption challenge: the finite vertices are detuned coordinates and branch
  classes, not runners or fixed arcs.  This quotient preserves the generic
  branch-cover predicate through its degree sum, but not the identity or cyclic
  placement of bad branches.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCDenseCoreEndgame

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

/-- At least one detuned coordinate has reduced denominator exactly two. -/
def HasTwoDetuningDenominator (v : Fin 13 → ℤ) (g : ℤ) : Prop :=
  ∃ i, ¬ g ∣ v i ∧ g / (Int.gcd (v i) g : ℤ) = 2

/-- A `q = 2` detuned coordinate has a distinct detuned companion with reduced
denominator below nine. -/
def HasTwoWithSubNineCompanion (v : Fin 13 → ℤ) (g : ℤ) : Prop :=
  ∃ i j, i ≠ j ∧ ¬ g ∣ v i ∧ ¬ g ∣ v j ∧
    g / (Int.gcd (v i) g : ℤ) = 2 ∧
    g / (Int.gcd (v j) g : ℤ) < 9

/-- Every detuned coordinate has reduced denominator exactly three. -/
def AllDetuningDenominatorsThree (v : Fin 13 → ℤ) (g : ℤ) : Prop :=
  ∀ i, ¬ g ∣ v i → g / (Int.gcd (v i) g : ℤ) = 3

/-- A `q = 2` bad row occupies exactly one half of the branch classes. -/
theorem two_mul_badCount_eq (δ g : ℤ) (hg : 0 < g)
    (hq : g / (Int.gcd δ g : ℤ) = 2) :
    2 * DetunedD3.badCount δ g = g.toNat := by
  obtain ⟨hbad, hfactor⟩ := badCount_of_q_two (δ := δ) (g := g) hg hq
  omega

/-- A `q = 3` bad row occupies exactly one third of the branch classes. -/
theorem three_mul_badCount_eq (δ g : ℤ) (hg : 0 < g)
    (hq : g / (Int.gcd δ g : ℤ) = 3) :
    3 * DetunedD3.badCount δ g = g.toNat := by
  have hdvd : ((Int.gcd δ g : ℤ)) ∣ g := Int.gcd_dvd_right δ g
  have htoNat : (g / (Int.gcd δ g : ℤ)).toNat = 3 := by
    rw [hq]
    rfl
  have hbad : DetunedD3.badCount δ g = Int.gcd δ g := by
    rw [DetunedD3.badCount, htoNat]
    norm_num
  have hfactor := Int.mul_ediv_cancel' hdvd
  rw [hq] at hfactor
  rw [hbad]
  omega

private theorem gcd_toNat_mul_reducedDenominator_twoThree
    (δ g : ℤ) (hg : 0 < g) :
    g.toNat = (Int.gcd δ g) *
      (g / (Int.gcd δ g : ℤ)).toNat := by
  have hdvd : ((Int.gcd δ g : ℤ)) ∣ g := Int.gcd_dvd_right δ g
  have hdz : (0 : ℤ) < (Int.gcd δ g : ℤ) := by
    have : 0 < Int.gcd δ g := by
      rw [Int.gcd_pos_iff]
      right
      omega
    exact_mod_cast this
  have hqnn : (0 : ℤ) ≤ g / (Int.gcd δ g : ℤ) :=
    Int.ediv_nonneg (le_of_lt hg) (le_of_lt hdz)
  have heq : g = (Int.gcd δ g : ℤ) *
      (g / (Int.gcd δ g : ℤ)) :=
    (Int.mul_ediv_cancel' hdvd).symm
  have hcast : (g.toNat : ℤ) =
      ((Int.gcd δ g) * (g / (Int.gcd δ g : ℤ)).toNat : ℕ) := by
    rw [Int.toNat_of_nonneg (le_of_lt hg)]
    push_cast
    rw [Int.toNat_of_nonneg hqnn]
    exact heq
  exact_mod_cast hcast

/-- Reduced denominator nine is the sharp threshold above which a row has bad
degree at most `2g/9`. -/
theorem nine_mul_badCount_le_two_mul (δ g : ℤ) (hg : 0 < g)
    (hq : 9 ≤ g / (Int.gcd δ g : ℤ)) :
    9 * DetunedD3.badCount δ g ≤ 2 * g.toNat := by
  have hQ9 : 9 ≤ (g / (Int.gcd δ g : ℤ)).toNat := by omega
  have hG : g.toNat = (Int.gcd δ g) *
      (g / (Int.gcd δ g : ℤ)).toNat :=
    gcd_toNat_mul_reducedDenominator_twoThree δ g hg
  have hnine :
      9 * ((g / (Int.gcd δ g : ℤ)).toNat / 7 + 1) ≤
        2 * (g / (Int.gcd δ g : ℤ)).toNat := by
    omega
  calc
    9 * DetunedD3.badCount δ g =
        (Int.gcd δ g) *
          (9 * ((g / (Int.gcd δ g : ℤ)).toNat / 7 + 1)) := by
            rw [DetunedD3.badCount]
            ring
    _ ≤ (Int.gcd δ g) *
        (2 * (g / (Int.gcd δ g : ℤ)).toNat) :=
      Nat.mul_le_mul_left _ hnine
    _ = 2 * ((Int.gcd δ g) *
        (g / (Int.gcd δ g : ℤ)).toNat) := by ring
    _ = 2 * g.toNat := by rw [hG]

/-- The common scaled row-degree bound for every reduced denominator at least
three.  It is exact at `q = 3`. -/
private theorem twentyone_mul_badCount_le_seven_mul (δ g : ℤ) (hg : 0 < g)
    (hq : 3 ≤ g / (Int.gcd δ g : ℤ)) :
    21 * DetunedD3.badCount δ g ≤ 7 * g.toNat := by
  by_cases hthree : g / (Int.gcd δ g : ℤ) = 3
  · have h := three_mul_badCount_eq δ g hg hthree
    omega
  · have hfour : 4 ≤ g / (Int.gcd δ g : ℤ) := by omega
    have h := seven_mul_badCount_le_two_mul δ g hg hfour
    omega

/-- The strict scaled row-degree bound once the reduced denominator is at
least four. -/
private theorem twentyone_mul_badCount_le_six_mul (δ g : ℤ) (hg : 0 < g)
    (hq : 4 ≤ g / (Int.gcd δ g : ℤ)) :
    21 * DetunedD3.badCount δ g ≤ 6 * g.toNat := by
  have h := seven_mul_badCount_le_two_mul δ g hg hq
  omega

/-- With exactly three detuned coordinates, excluding `q = 2` and the uniform
`(3,3,3)` pattern forces the generic bad-count inequality. -/
theorem genericCount_of_three_noTwo_notAllThree
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g)
    (hcard : nonMultCard v g = 3)
    (hnoTwo : ¬ HasTwoDetuningDenominator v g)
    (hnotThree : ¬ AllDetuningDenominatorsThree v g) :
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
  have hq1 : 3 ≤ g / (Int.gcd (v i₁) g : ℤ) := by
    have hge := reducedDetuningDenominator_ge_two hg hδ1
    have hne : g / (Int.gcd (v i₁) g : ℤ) ≠ 2 := by
      intro htwo
      exact hnoTwo ⟨i₁, hδ1, htwo⟩
    omega
  have hq2 : 3 ≤ g / (Int.gcd (v i₂) g : ℤ) := by
    have hge := reducedDetuningDenominator_ge_two hg hδ2
    have hne : g / (Int.gcd (v i₂) g : ℤ) ≠ 2 := by
      intro htwo
      exact hnoTwo ⟨i₂, hδ2, htwo⟩
    omega
  have hq3 : 3 ≤ g / (Int.gcd (v i₃) g : ℤ) := by
    have hge := reducedDetuningDenominator_ge_two hg hδ3
    have hne : g / (Int.gcd (v i₃) g : ℤ) ≠ 2 := by
      intro htwo
      exact hnoTwo ⟨i₃, hδ3, htwo⟩
    omega
  have hnotAll : ¬
      (g / (Int.gcd (v i₁) g : ℤ) = 3 ∧
        g / (Int.gcd (v i₂) g : ℤ) = 3 ∧
        g / (Int.gcd (v i₃) g : ℤ) = 3) := by
    rintro ⟨hq1eq, hq2eq, hq3eq⟩
    apply hnotThree
    intro i hδ
    have hi : i ∈ Finset.univ.filter (fun j => ¬ g ∣ v j) :=
      (hmem i).mpr hδ
    rw [hfilter] at hi
    simp only [Finset.mem_insert, Finset.mem_singleton] at hi
    rcases hi with rfl | rfl | rfl
    · exact hq1eq
    · exact hq2eq
    · exact hq3eq
  have hstrict :
      4 ≤ g / (Int.gcd (v i₁) g : ℤ) ∨
        4 ≤ g / (Int.gcd (v i₂) g : ℤ) ∨
        4 ≤ g / (Int.gcd (v i₃) g : ℤ) := by
    by_contra hnone
    push Not at hnone
    exact hnotAll ⟨by omega, by omega, by omega⟩
  have hb1 := twentyone_mul_badCount_le_seven_mul (v i₁) g hg0 hq1
  have hb2 := twentyone_mul_badCount_le_seven_mul (v i₂) g hg0 hq2
  have hb3 := twentyone_mul_badCount_le_seven_mul (v i₃) g hg0 hq3
  rw [genericCount, hfilter,
    Finset.sum_insert (by simp [h12, h13]),
    Finset.sum_insert (by simp [h23]), Finset.sum_singleton]
  rcases hstrict with hq1strict | hq2strict | hq3strict
  · have hs1 := twentyone_mul_badCount_le_six_mul (v i₁) g hg0 hq1strict
    omega
  · have hs2 := twentyone_mul_badCount_le_six_mul (v i₂) g hg0 hq2strict
    omega
  · have hs3 := twentyone_mul_badCount_le_six_mul (v i₃) g hg0 hq3strict
    omega

/-- A triple with a `q = 2` row is generic if both other detuned rows have
reduced denominator at least nine. -/
theorem genericCount_of_three_two_noSubNineCompanion
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g)
    (hcard : nonMultCard v g = 3)
    (htwo : HasTwoDetuningDenominator v g)
    (hnoCompanion : ¬ HasTwoWithSubNineCompanion v g) :
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
  obtain ⟨i, hδi, hqi⟩ := htwo
  have hi : i ∈ Finset.univ.filter (fun j => ¬ g ∣ v j) :=
    (hmem i).mpr hδi
  rw [hfilter] at hi
  simp only [Finset.mem_insert, Finset.mem_singleton] at hi
  have hcompanion : ∀ j, ¬ g ∣ v j → j ≠ i →
      9 ≤ g / (Int.gcd (v j) g : ℤ) := by
    intro j hδj hji
    by_contra hnot
    exact hnoCompanion
      ⟨i, j, hji.symm, hδi, hδj, hqi, lt_of_not_ge hnot⟩
  rw [genericCount, hfilter,
    Finset.sum_insert (by simp [h12, h13]),
    Finset.sum_insert (by simp [h23]), Finset.sum_singleton]
  rcases hi with hi | hi | hi
  · subst i
    have hb1 := two_mul_badCount_eq (v i₁) g hg0 hqi
    have hb2 := nine_mul_badCount_le_two_mul (v i₂) g hg0
      (hcompanion i₂ hδ2 h12.symm)
    have hb3 := nine_mul_badCount_le_two_mul (v i₃) g hg0
      (hcompanion i₃ hδ3 h13.symm)
    omega
  · subst i
    have hb1 := nine_mul_badCount_le_two_mul (v i₁) g hg0
      (hcompanion i₁ hδ1 h12)
    have hb2 := two_mul_badCount_eq (v i₂) g hg0 hqi
    have hb3 := nine_mul_badCount_le_two_mul (v i₃) g hg0
      (hcompanion i₃ hδ3 h23.symm)
    omega
  · subst i
    have hb1 := nine_mul_badCount_le_two_mul (v i₁) g hg0
      (hcompanion i₁ hδ1 h13)
    have hb2 := nine_mul_badCount_le_two_mul (v i₂) g hg0
      (hcompanion i₂ hδ2 h23)
    have hb3 := two_mul_badCount_eq (v i₃) g hg0 hqi
    omega

/-- Exact denominator-pattern classification of a non-generic triple. -/
theorem nonGeneric_three_hasTwo_or_allThree
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g)
    (hcard : nonMultCard v g = 3) (hnongeneric : ¬ genericCount v g) :
    HasTwoDetuningDenominator v g ∨
      AllDetuningDenominatorsThree v g := by
  by_cases htwo : HasTwoDetuningDenominator v g
  · exact Or.inl htwo
  · right
    by_contra hthree
    exact hnongeneric
      (genericCount_of_three_noTwo_notAllThree v g hg hcard htwo hthree)

/-- Sharp one-coordinate refinement of the denominator classification: if the
exceptional triple has a `q = 2` row, it also has a distinct row with `q ≤ 8`.
The cutoff is sharp for degree arithmetic because `(2,8,8)` has total bad
degree exactly `g`. -/
theorem nonGeneric_three_hasTwoCompanion_or_allThree
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g)
    (hcard : nonMultCard v g = 3) (hnongeneric : ¬ genericCount v g) :
    HasTwoWithSubNineCompanion v g ∨
      AllDetuningDenominatorsThree v g := by
  by_cases htwo : HasTwoDetuningDenominator v g
  · by_cases hcompanion : HasTwoWithSubNineCompanion v g
    · exact Or.inl hcompanion
    · exact False.elim (hnongeneric
        (genericCount_of_three_two_noSubNineCompanion
          v g hg hcard htwo hcompanion))
  · exact Or.inr ((nonGeneric_three_hasTwo_or_allThree
      v g hg hcard hnongeneric).resolve_left htwo)

/-- The newly closed honest subcase of the exactly-three-detuned dispatch. -/
theorem lonely14_of_nonMultCard_three_noTwo_notAllThree (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (hcard : nonMultCard v g = 3)
    (hnoTwo : ¬ HasTwoDetuningDenominator v g)
    (hnotThree : ¬ AllDetuningDenominatorsThree v g) :
    ∃ t : ℝ, Lonely 14 v t :=
  lonely14_of_nonMultCard_three cite v hv g hg hcard
    (genericCount_of_three_noTwo_notAllThree
      v g hg hcard hnoTwo hnotThree)

/-- The exact `q = 2/3` deep residual.  The triple branch now consists only of
a `q = 2` coordinate with a distinct `q ≤ 8` companion, or the uniform
`(3,3,3)` pattern. -/
def DeepExceptionalDetunedDispatchTwoThree : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → ∀ g : ℤ, 2 ≤ g →
    ((nonMultCard v g = 2 ∧ ¬ TwoAdicLiftTerminates v g) ∨
      (nonMultCard v g = 3 ∧
        (HasTwoWithSubNineCompanion v g ∨
          AllDetuningDenominatorsThree v g))) →
    ¬ genericCount v g →
    ∃ t : ℝ, Lonely 14 v t

/-- The exact `q = 2/3` residual supplies the prior sub-four residual because
every other triple denominator pattern is now a theorem. -/
theorem deepExceptionalDetunedDispatchFour_of_twoThree
    (hdeep : DeepExceptionalDetunedDispatchTwoThree) :
    DeepExceptionalDetunedDispatchFour := by
  intro v hv g hg hcase hnongeneric
  rcases hcase with htwo | ⟨hthree, _⟩
  · exact hdeep v hv g hg (Or.inl htwo) hnongeneric
  · have hpattern :=
      nonGeneric_three_hasTwoCompanion_or_allThree
        v g hg hthree hnongeneric
    exact hdeep v hv g hg (Or.inr ⟨hthree, hpattern⟩) hnongeneric

/-- Current narrowest machine-checked endgame surface: the nonterminating pair
tower, the `q = 2`-with-`q ≤ 8`-companion/uniform-`q = 3` triple residual, and
positive B5 only on the primitive dissociated chain-dense core. -/
theorem lrc14_from_twoThree_detuned_and_denseCore_dissociated_B5
    (cite : LRCUpTo13)
    (hdeep : DeepExceptionalDetunedDispatchTwoThree)
    (hB5 : DenseCoreDissociatedB5Supply) :
    LRC14.LRC14Statement :=
  lrc14_from_four_detuned_and_denseCore_dissociated_B5 cite
    (deepExceptionalDetunedDispatchFour_of_twoThree hdeep) hB5

/-! ## Axiom audit -/

#print axioms two_mul_badCount_eq
#print axioms three_mul_badCount_eq
#print axioms nine_mul_badCount_le_two_mul
#print axioms genericCount_of_three_noTwo_notAllThree
#print axioms genericCount_of_three_two_noSubNineCompanion
#print axioms nonGeneric_three_hasTwo_or_allThree
#print axioms nonGeneric_three_hasTwoCompanion_or_allThree
#print axioms lonely14_of_nonMultCard_three_noTwo_notAllThree
#print axioms deepExceptionalDetunedDispatchFour_of_twoThree
#print axioms lrc14_from_twoThree_detuned_and_denseCore_dissociated_B5

end LRC14Grand
end LonelyRunner
