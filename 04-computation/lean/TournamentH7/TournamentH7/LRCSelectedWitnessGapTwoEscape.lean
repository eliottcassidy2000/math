import TournamentH7.LRCPairTowerGapTwo
import TournamentH7.LRCPairTowerGapTwoFrequency
import TournamentH7.LRCLadderFattening

/-!
# Dynamic escape for the `(4,4,8,8)` pair-tower wall

Four exceptional rows leave nine harmonic rows.  The `LRC(9)` citation
therefore supplies a closed harmonic-good interval of radius `1 / (35 * B)`.
Any nonzero scalar failure frequency whose `3/14` danger chambers fit inside
that interval image must miss one chamber and hence select a common good
branch.  This is the dynamic consumer needed after the static parallel-class
classification in `LRCPairTowerGapTwo`.

Tournament-analysis audit: the vertices are phase chambers, not runners.  The
binary relation records whether two chambers are paired by integer translation;
the scalar-frequency ordering is a Hamiltonian path through these chambers.
Quotienting to runners loses the phase width and the integer chamber labels.
-/

namespace LonelyRunner
namespace LRCPairTowerGapTwo

open LonelyRunner LRCPairTowerGapTwoFrequency
open scoped Classical

noncomputable section

private theorem exists_integer_clear_in_interval_four
    (lower upper radius : ℝ)
    (hradius0 : 0 ≤ radius) (hradiusHalf : 2 * radius ≤ 1)
    (hlowerUpper : lower ≤ upper)
    (hwidth : 2 * radius ≤ upper - lower) :
    ∃ point ∈ Set.Icc lower upper,
      ∀ integer : ℤ, radius ≤ |point - integer| := by
  by_cases hlowerClear : ∀ integer : ℤ, radius ≤ |lower - integer|
  · exact ⟨lower, ⟨le_rfl, hlowerUpper⟩, hlowerClear⟩
  · push Not at hlowerClear
    obtain ⟨nearest, hnearest⟩ := hlowerClear
    refine ⟨(nearest : ℝ) + radius, ?_, ?_⟩
    · rw [abs_lt] at hnearest
      constructor <;> linarith
    · intro integer
      rcases lt_trichotomy integer nearest with hlt | heq | hgt
      · have hstep : (integer : ℝ) + 1 ≤ nearest := by exact_mod_cast hlt
        rw [abs_of_nonneg]
        · linarith
        · linarith
      · subst integer
        simp [abs_of_nonneg hradius0]
      · have hstep : (nearest : ℝ) + 1 ≤ integer := by exact_mod_cast hgt
        rw [abs_of_nonpos]
        · linarith
        · linarith

/-- A scalar failure wall cannot cover a harmonic-good interval once its image
contains a complete closed complement of one integer danger chamber. -/
theorem selectedWitness_of_four_harmonicInterval_and_scalarFailureWall
    (v : Fin 13 → ℤ) (g : ℤ) (i₁ i₂ i₃ i₄ : Fin 13)
    (center intervalRadius dangerRadius : ℝ) (frequency : ℤ)
    (hintervalRadius : 0 ≤ intervalRadius)
    (hdangerRadius : 0 ≤ dangerRadius)
    (hdangerHalf : 2 * dangerRadius ≤ 1)
    (hfrequency : frequency ≠ 0)
    (hinterval : ∀ u ∈ Set.Icc (center - intervalRadius)
        (center + intervalRadius),
      FourDetunedHarmonicGoodAt v g i₁ i₂ i₃ i₄ u)
    (hfailureWall : ∀ u : ℝ,
      ¬ HasFourDetunedGoodBranch (v i₁) (v i₂) (v i₃) (v i₄) g u →
      ∃ integer : ℤ,
        |(frequency : ℝ) * u - integer| < dangerRadius)
    (hwidth : 2 * dangerRadius ≤
      2 * intervalRadius * |(frequency : ℝ)|) :
    ∃ u : ℝ,
      FourDetunedHarmonicGoodAt v g i₁ i₂ i₃ i₄ u ∧
      HasFourDetunedGoodBranch (v i₁) (v i₂) (v i₃) (v i₄) g u := by
  have hfrequencyReal : (frequency : ℝ) ≠ 0 := by exact_mod_cast hfrequency
  rcases lt_or_gt_of_ne hfrequencyReal with hnegative | hpositive
  · let lower : ℝ := (frequency : ℝ) * (center + intervalRadius)
    let upper : ℝ := (frequency : ℝ) * (center - intervalRadius)
    have hlowerUpper : lower ≤ upper := by
      dsimp [lower, upper]
      nlinarith
    have hspan : upper - lower =
        2 * intervalRadius * |(frequency : ℝ)| := by
      rw [abs_of_neg hnegative]
      dsimp [lower, upper]
      ring
    obtain ⟨point, hpoint, hclear⟩ :=
      exists_integer_clear_in_interval_four lower upper dangerRadius
        hdangerRadius hdangerHalf hlowerUpper (by rwa [hspan])
    let u : ℝ := point / (frequency : ℝ)
    have hphase : (frequency : ℝ) * u = point := by
      dsimp [u]
      field_simp
    have hu : u ∈ Set.Icc (center - intervalRadius)
        (center + intervalRadius) := by
      dsimp [lower, upper] at hpoint
      constructor <;> nlinarith [hpoint.1, hpoint.2]
    refine ⟨u, hinterval u hu, ?_⟩
    by_contra hfail
    obtain ⟨integer, hinteger⟩ := hfailureWall u hfail
    rw [hphase] at hinteger
    exact (not_lt_of_ge (hclear integer)) hinteger
  · let lower : ℝ := (frequency : ℝ) * (center - intervalRadius)
    let upper : ℝ := (frequency : ℝ) * (center + intervalRadius)
    have hlowerUpper : lower ≤ upper := by
      dsimp [lower, upper]
      nlinarith
    have hspan : upper - lower =
        2 * intervalRadius * |(frequency : ℝ)| := by
      rw [abs_of_pos hpositive]
      dsimp [lower, upper]
      ring
    obtain ⟨point, hpoint, hclear⟩ :=
      exists_integer_clear_in_interval_four lower upper dangerRadius
        hdangerRadius hdangerHalf hlowerUpper (by rwa [hspan])
    let u : ℝ := point / (frequency : ℝ)
    have hphase : (frequency : ℝ) * u = point := by
      dsimp [u]
      field_simp
    have hu : u ∈ Set.Icc (center - intervalRadius)
        (center + intervalRadius) := by
      dsimp [lower, upper] at hpoint
      constructor <;> nlinarith [hpoint.1, hpoint.2]
    refine ⟨u, hinterval u hu, ?_⟩
    by_contra hfail
    obtain ⟨integer, hinteger⟩ := hfailureWall u hfail
    rw [hphase] at hinteger
    exact (not_lt_of_ge (hclear integer)) hinteger

/-- Four selected rows leave nine quotient speeds.  An `LRC(9)` witness has
margin `1/10`, so it remains `1/14`-good throughout the exact radius
`(1/10 - 1/14) / B = 1/(35B)`. -/
theorem exists_fourDetunedHarmonicGoodAt_sharpInterval_of_cite
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (i₁ i₂ i₃ i₄ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h14 : i₁ ≠ i₄)
    (h23 : i₂ ≠ i₃) (h24 : i₂ ≠ i₄) (h34 : i₃ ≠ i₄)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → j ≠ i₄ → g ∣ v j)
    (B : ℝ) (hB0 : 0 < B)
    (hB : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → j ≠ i₄ →
      |((v j / g : ℤ) : ℝ)| ≤ B) :
    ∃ center : ℝ,
      ∀ u ∈ Set.Icc (center - 1 / (35 * B))
          (center + 1 / (35 * B)),
        FourDetunedHarmonicGoodAt v g i₁ i₂ i₃ i₄ u := by
  have hselectedCard :
      ({i₁, i₂, i₃, i₄} : Finset (Fin 13)).card = 4 := by
    simp [h12, h13, h14, h23, h24, h34]
  have hcard :
      (Finset.univ \ ({i₁, i₂, i₃, i₄} : Finset (Fin 13))).card = 9 := by
    rw [Finset.card_sdiff, Finset.inter_univ, hselectedCard,
      Finset.card_univ, Fintype.card_fin]
  set embedding : Fin 9 ↪o Fin 13 :=
    (Finset.univ \ ({i₁, i₂, i₃, i₄} : Finset (Fin 13))).orderEmbOfFin hcard
      with hembedding
  have hembeddingMem : ∀ index,
      embedding index ∈ Finset.univ \
        ({i₁, i₂, i₃, i₄} : Finset (Fin 13)) :=
    fun index =>
      (Finset.univ \ ({i₁, i₂, i₃, i₄} : Finset (Fin 13))).orderEmbOfFin_mem
        hcard index
  have hembeddingNe : ∀ index,
      embedding index ≠ i₁ ∧ embedding index ≠ i₂ ∧
        embedding index ≠ i₃ ∧ embedding index ≠ i₄ := by
    intro index
    have hmem := hembeddingMem index
    rw [Finset.mem_sdiff, Finset.mem_insert, Finset.mem_insert,
      Finset.mem_insert, Finset.mem_singleton] at hmem
    push Not at hmem
    exact hmem.2
  set quotient : Fin 9 → ℤ := fun index => v (embedding index) / g
    with hquotient
  have hquotientNonzero : ∀ index, quotient index ≠ 0 := by
    intro index hzero
    exact hv (embedding index) (by
      have hdivides := hdvd (embedding index) (hembeddingNe index).1
        (hembeddingNe index).2.1 (hembeddingNe index).2.2.1
        (hembeddingNe index).2.2.2
      have hfactor : v (embedding index) = g * quotient index :=
        (Int.mul_ediv_cancel' hdivides).symm
      rw [hfactor, hzero, mul_zero])
  have hquotientBound : ∀ index, |(quotient index : ℝ)| ≤ B := by
    intro index
    exact hB (embedding index) (hembeddingNe index).1
      (hembeddingNe index).2.1 (hembeddingNe index).2.2.1
      (hembeddingNe index).2.2.2
  obtain ⟨center, hcenter⟩ :=
    cite 9 (by norm_num) quotient hquotientNonzero
  refine ⟨center, ?_⟩
  intro u hu
  have habs : |u - center| ≤ 1 / (35 * B) := by
    rw [abs_le]
    constructor <;> linarith [hu.1, hu.2]
  have hdrift : B * |u - center| ≤ 1 / 35 := by
    calc
      B * |u - center| ≤ B * (1 / (35 * B)) :=
        mul_le_mul_of_nonneg_left habs (le_of_lt hB0)
      _ = 1 / 35 := by field_simp
  have htransport : Lonely 14 quotient u := by
    intro index integer
    apply lonely_band_transport 10 quotient center u B (1 / 14)
      hcenter hquotientBound
    norm_num
    linarith
  intro index hindex1 hindex2 hindex3 hindex4 integer
  have hindexMem :
      index ∈ Finset.univ \
        ({i₁, i₂, i₃, i₄} : Finset (Fin 13)) := by
    rw [Finset.mem_sdiff, Finset.mem_insert, Finset.mem_insert,
      Finset.mem_insert, Finset.mem_singleton]
    exact ⟨Finset.mem_univ index, by
      push Not
      exact ⟨hindex1, hindex2, hindex3, hindex4⟩⟩
  obtain ⟨source, hsource⟩ : ∃ source, embedding source = index := by
    have hrange : index ∈ Set.range embedding := by
      rw [hembedding, Finset.range_orderEmbOfFin]
      exact hindexMem
    exact hrange
  have hquotientEq : quotient source = v index / g := by
    show v (embedding source) / g = v index / g
    rw [hsource]
  simpa [hquotientEq] using htransport source integer

/-- Exact large-frequency socket for the gap-two wall.  The sharp width
condition is `15B ≤ 2|F|`; the simpler sufficient condition `8B ≤ |F|` may be
derived by arithmetic by downstream callers. -/
theorem selectedWitness_of_cite_and_four_threeFourteenths_scalarWall
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (i₁ i₂ i₃ i₄ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h14 : i₁ ≠ i₄)
    (h23 : i₂ ≠ i₃) (h24 : i₂ ≠ i₄) (h34 : i₃ ≠ i₄)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → j ≠ i₄ → g ∣ v j)
    (B : ℝ) (hB0 : 0 < B)
    (hB : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → j ≠ i₄ →
      |((v j / g : ℤ) : ℝ)| ≤ B)
    (frequency : ℤ) (hfrequency : frequency ≠ 0)
    (hfrequencyLarge : 15 * B ≤ 2 * |(frequency : ℝ)|)
    (hfailureWall : ∀ u : ℝ,
      ¬ HasFourDetunedGoodBranch (v i₁) (v i₂) (v i₃) (v i₄) g u →
      ∃ integer : ℤ,
        |(frequency : ℝ) * u - integer| < 3 / 14) :
    ∃ u : ℝ,
      FourDetunedHarmonicGoodAt v g i₁ i₂ i₃ i₄ u ∧
      HasFourDetunedGoodBranch (v i₁) (v i₂) (v i₃) (v i₄) g u := by
  obtain ⟨center, hinterval⟩ :=
    exists_fourDetunedHarmonicGoodAt_sharpInterval_of_cite
      cite v hv g i₁ i₂ i₃ i₄ h12 h13 h14 h23 h24 h34 hdvd B hB0 hB
  apply selectedWitness_of_four_harmonicInterval_and_scalarFailureWall
    v g i₁ i₂ i₃ i₄ center (1 / (35 * B)) (3 / 14) frequency
  · positivity
  · norm_num
  · norm_num
  · exact hfrequency
  · exact hinterval
  · exact hfailureWall
  · have hBne : B ≠ 0 := ne_of_gt hB0
    have hscaled : (15 / 2 : ℝ) * B ≤ |(frequency : ℝ)| := by
      linarith
    calc
      2 * (3 / 14 : ℝ) =
          2 * (1 / (35 * B)) * ((15 / 2 : ℝ) * B) := by
            field_simp [hBne]
            ring
      _ ≤ 2 * (1 / (35 * B)) * |(frequency : ℝ)| :=
        mul_le_mul_of_nonneg_left hscaled (by positivity)

/-- Composition with the static q4488 matching theorem.  All normalized
numerators are fixed; only the three row offsets and their bad witnesses vary
with phase.  The remaining dynamic producer must establish `hmatching` and
the nonzero/large-frequency alternatives. -/
theorem selectedWitness_of_cite_and_fourEightEight_matchingWall
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (i₁ i₂ i₃ i₄ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h14 : i₁ ≠ i₄)
    (h23 : i₂ ≠ i₃) (h24 : i₂ ≠ i₄) (h34 : i₃ ≠ i₄)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → j ≠ i₄ → g ∣ v j)
    (B : ℝ) (hB0 : 0 < B)
    (hB : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → j ≠ i₄ →
      |((v j / g : ℤ) : ℝ)| ≤ B)
    (a₄ a₈a a₈b residue : ℤ)
    (hres₄ : (4 : ℤ) ∣ a₄ - residue)
    (hres₈a : (8 : ℤ) ∣ a₈a - residue)
    (hres₈b : (8 : ℤ) ∣ a₈b - residue)
    (hfrequency : fourEightEightPhaseFrequency a₄ a₈a a₈b ≠ 0)
    (hfrequencyLarge :
      15 * B ≤ 2 * |(fourEightEightPhaseFrequency a₄ a₈a a₈b : ℝ)|)
    (hmatching : ∀ u : ℝ,
      ¬ HasFourDetunedGoodBranch (v i₁) (v i₂) (v i₃) (v i₄) g u →
      ∃ c₄ c₈a c₈b : ℤ,
        (8 : ℤ) ∣ -2 * c₄ + c₈a + c₈b ∧
        (∃ n : ℤ,
          |(a₄ : ℝ) * ((u + (c₄ : ℝ)) / 4) - n| < 1 / 14) ∧
        (∃ n : ℤ,
          |(a₈a : ℝ) * ((u + (c₈a : ℝ)) / 8) - n| < 1 / 14) ∧
        (∃ n : ℤ,
          |(a₈b : ℝ) * ((u + (c₈b : ℝ)) / 8) - n| < 1 / 14)) :
    ∃ u : ℝ,
      FourDetunedHarmonicGoodAt v g i₁ i₂ i₃ i₄ u ∧
      HasFourDetunedGoodBranch (v i₁) (v i₂) (v i₃) (v i₄) g u := by
  apply selectedWitness_of_cite_and_four_threeFourteenths_scalarWall
    cite v hv g i₁ i₂ i₃ i₄ h12 h13 h14 h23 h24 h34 hdvd B hB0 hB
      (fourEightEightPhaseFrequency a₄ a₈a a₈b)
      hfrequency hfrequencyLarge
  intro u hfail
  obtain ⟨c₄, c₈a, c₈b, hcode, hbad₄, hbad₈a, hbad₈b⟩ :=
    hmatching u hfail
  exact frequency_bad_of_fourEightEight_matching
    a₄ a₈a a₈b residue c₄ c₈a c₈b u
      hres₄ hres₈a hres₈b hcode hbad₄ hbad₈a hbad₈b

/-- Direct pair-tower consumer: once the four-row failure wall has a nonzero
large scalar frequency, the selected phase produces an actual LRC(14) time. -/
theorem lonely14_of_gapTwo_cite_and_threeFourteenths_scalarWall
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (hg : 2 ≤ g) (i₁ i₂ i₃ i₄ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h14 : i₁ ≠ i₄)
    (h23 : i₂ ≠ i₃) (h24 : i₂ ≠ i₄) (h34 : i₃ ≠ i₄)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → j ≠ i₄ → g ∣ v j)
    (B : ℝ) (hB0 : 0 < B)
    (hB : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → j ≠ i₄ →
      |((v j / g : ℤ) : ℝ)| ≤ B)
    (frequency : ℤ) (hfrequency : frequency ≠ 0)
    (hfrequencyLarge : 15 * B ≤ 2 * |(frequency : ℝ)|)
    (hfailureWall : ∀ u : ℝ,
      ¬ HasFourDetunedGoodBranch (v i₁) (v i₂) (v i₃) (v i₄) g u →
      ∃ integer : ℤ,
        |(frequency : ℝ) * u - integer| < 3 / 14) :
    ∃ t : ℝ, Lonely 14 v t := by
  apply lonely14_of_four_detuned_selectedWitness v g hg i₁ i₂ i₃ i₄ hdvd
  exact selectedWitness_of_cite_and_four_threeFourteenths_scalarWall
    cite v hv g i₁ i₂ i₃ i₄ h12 h13 h14 h23 h24 h34 hdvd B hB0 hB
      frequency hfrequency hfrequencyLarge hfailureWall

#print axioms selectedWitness_of_four_harmonicInterval_and_scalarFailureWall
#print axioms exists_fourDetunedHarmonicGoodAt_sharpInterval_of_cite
#print axioms selectedWitness_of_cite_and_four_threeFourteenths_scalarWall
#print axioms selectedWitness_of_cite_and_fourEightEight_matchingWall
#print axioms lonely14_of_gapTwo_cite_and_threeFourteenths_scalarWall

end
end LRCPairTowerGapTwo
end LonelyRunner
