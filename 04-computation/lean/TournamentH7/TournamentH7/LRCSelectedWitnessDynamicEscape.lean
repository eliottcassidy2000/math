import TournamentH7.LRCEndgameUniformThreePhase
import TournamentH7.LRCSelectedWitnessTwoFourFourFrequency

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

noncomputable section

/-- An interval at least as long as one open integer-danger chamber contains
a point on or outside that chamber. -/
theorem exists_integer_clear_in_interval
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

/-- Abstract dynamic escape lemma.  If every failed selected branch forces an
integer scalar frequency into an open chamber of radius `dangerRadius`, then
a harmonic-good interval whose scalar image is at least one full chamber long
contains a successful phase. -/
theorem selectedWitness_of_harmonicInterval_and_scalarFailureWall
    (v : Fin 13 → ℤ) (g : ℤ) (i₁ i₂ i₃ : Fin 13)
    (center intervalRadius dangerRadius : ℝ) (frequency : ℤ)
    (hintervalRadius : 0 ≤ intervalRadius)
    (hdangerRadius : 0 ≤ dangerRadius)
    (hdangerHalf : 2 * dangerRadius ≤ 1)
    (hfrequency : frequency ≠ 0)
    (hinterval : ∀ u ∈ Set.Icc (center - intervalRadius)
        (center + intervalRadius),
      ThreeDetunedHarmonicGoodAt v g i₁ i₂ i₃ u)
    (hfailureWall : ∀ u : ℝ,
      ¬ HasThreeDetunedGoodBranch (v i₁) (v i₂) (v i₃) g u →
      ∃ integer : ℤ,
        |(frequency : ℝ) * u - integer| < dangerRadius)
    (hwidth : 2 * dangerRadius ≤
      2 * intervalRadius * |(frequency : ℝ)|) :
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₁ i₂ i₃ u ∧
      HasThreeDetunedGoodBranch (v i₁) (v i₂) (v i₃) g u := by
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
      exists_integer_clear_in_interval lower upper dangerRadius
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
      exists_integer_clear_in_interval lower upper dangerRadius
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

/-- Local Lipschitz transport, specialized only in the next theorem. -/
private theorem lonely_band_transport_dynamic
    {ι : Type*} (n : ℕ) (w : ι → ℤ) (u₀ u B clearance : ℝ)
    (hlonely : Lonely n w u₀) (hB : ∀ index, |(w index : ℝ)| ≤ B)
    (hclearance : clearance ≤ 1 / n - B * |u - u₀|) :
    ∀ index, ∀ integer : ℤ,
      clearance ≤ |(w index : ℝ) * u - integer| := by
  intro index integer
  have hbase : (1 : ℝ) / n ≤
      |(w index : ℝ) * u₀ - integer| := hlonely index integer
  have htriangle :
      |(w index : ℝ) * u₀ - integer| -
          |(w index : ℝ) * u₀ - (w index : ℝ) * u| ≤
        |(w index : ℝ) * u - integer| := by
    have h := abs_sub_abs_le_abs_sub
      ((w index : ℝ) * u₀ - integer) ((w index : ℝ) * u - integer)
    have hrearrange :
        ((w index : ℝ) * u₀ - integer) -
            ((w index : ℝ) * u - integer) =
          (w index : ℝ) * u₀ - (w index : ℝ) * u := by ring
    rw [hrearrange] at h
    linarith
  have hdrift :
      |(w index : ℝ) * u₀ - (w index : ℝ) * u| ≤
        B * |u - u₀| := by
    have hfactor :
        (w index : ℝ) * u₀ - (w index : ℝ) * u =
          (w index : ℝ) * (u₀ - u) := by ring
    rw [hfactor, abs_mul, abs_sub_comm]
    exact mul_le_mul_of_nonneg_right (hB index) (abs_nonneg _)
  linarith

/-- Sharp ten-frequency citation fattening.  The `LRC(10)` witness has margin
`1/11`; therefore the harmonic-good set contains a closed interval of radius
`(1/11-1/14)/B = 3/(154B)`. -/
theorem exists_threeDetunedHarmonicGoodAt_sharpInterval_of_cite
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (i₁ i₂ i₃ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j)
    (B : ℝ) (hB0 : 0 < B)
    (hB : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ →
      |((v j / g : ℤ) : ℝ)| ≤ B) :
    ∃ center : ℝ,
      ∀ u ∈ Set.Icc (center - 3 / (154 * B))
          (center + 3 / (154 * B)),
        ThreeDetunedHarmonicGoodAt v g i₁ i₂ i₃ u := by
  have hselectedCard : ({i₁, i₂, i₃} : Finset (Fin 13)).card = 3 :=
    Finset.card_eq_three.mpr ⟨i₁, i₂, i₃, h12, h13, h23, rfl⟩
  have hcard :
      (Finset.univ \ ({i₁, i₂, i₃} : Finset (Fin 13))).card = 10 := by
    rw [Finset.card_sdiff, Finset.inter_univ, hselectedCard,
      Finset.card_univ, Fintype.card_fin]
  set embedding : Fin 10 ↪o Fin 13 :=
    (Finset.univ \ ({i₁, i₂, i₃} : Finset (Fin 13))).orderEmbOfFin hcard
      with hembedding
  have hembeddingMem : ∀ index,
      embedding index ∈ Finset.univ \
        ({i₁, i₂, i₃} : Finset (Fin 13)) :=
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
  have hquotientBound : ∀ index, |(quotient index : ℝ)| ≤ B := by
    intro index
    exact hB (embedding index) (hembeddingNe index).1
      (hembeddingNe index).2.1 (hembeddingNe index).2.2
  obtain ⟨center, hcenter⟩ :=
    cite 10 (by norm_num) quotient hquotientNonzero
  refine ⟨center, ?_⟩
  intro u hu
  have habs : |u - center| ≤ 3 / (154 * B) := by
    rw [abs_le]
    constructor <;> linarith [hu.1, hu.2]
  have hdrift : B * |u - center| ≤ 3 / 154 := by
    calc
      B * |u - center| ≤ B * (3 / (154 * B)) :=
        mul_le_mul_of_nonneg_left habs (le_of_lt hB0)
      _ = 3 / 154 := by field_simp
  have htransport : Lonely 14 quotient u := by
    intro index integer
    apply lonely_band_transport_dynamic 11 quotient center u B (1 / 14)
      hcenter hquotientBound
    norm_num
    linarith
  intro index hindex1 hindex2 hindex3 integer
  have hindexMem :
      index ∈ Finset.univ \
        ({i₁, i₂, i₃} : Finset (Fin 13)) := by
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
  simpa [hquotientEq] using htransport source integer

/-- Repo-facing composition: the sharp citation interval escapes any scalar
failure wall of radius `3/14` once the wall frequency is at least eleven times
the largest harmonic quotient.  This is the common q244/q333 large-frequency
socket. -/
theorem selectedWitness_of_cite_and_threeFourteenths_scalarWall
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (i₁ i₂ i₃ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j)
    (B : ℝ) (hB0 : 0 < B)
    (hB : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ →
      |((v j / g : ℤ) : ℝ)| ≤ B)
    (frequency : ℤ) (hfrequency : frequency ≠ 0)
    (hfrequencyLarge : 11 * B ≤ |(frequency : ℝ)|)
    (hfailureWall : ∀ u : ℝ,
      ¬ HasThreeDetunedGoodBranch (v i₁) (v i₂) (v i₃) g u →
      ∃ integer : ℤ,
        |(frequency : ℝ) * u - integer| < 3 / 14) :
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₁ i₂ i₃ u ∧
      HasThreeDetunedGoodBranch (v i₁) (v i₂) (v i₃) g u := by
  obtain ⟨center, hinterval⟩ :=
    exists_threeDetunedHarmonicGoodAt_sharpInterval_of_cite
      cite v hv g i₁ i₂ i₃ h12 h13 h23 hdvd B hB0 hB
  apply selectedWitness_of_harmonicInterval_and_scalarFailureWall
    v g i₁ i₂ i₃ center (3 / (154 * B)) (3 / 14) frequency
  · positivity
  · norm_num
  · norm_num
  · exact hfrequency
  · exact hinterval
  · exact hfailureWall
  · have hBne : B ≠ 0 := ne_of_gt hB0
    calc
      2 * (3 / 14 : ℝ) =
          2 * (3 / (154 * B)) * (11 * B) := by
            field_simp [hBne]
            ring
      _ ≤ 2 * (3 / (154 * B)) * |(frequency : ℝ)| :=
        mul_le_mul_of_nonneg_left hfrequencyLarge (by positivity)

/-- Direct primitive q333 corollary.  The existing cyclic-obstruction theorem
supplies the scalar wall; only the quantitative nonzero-frequency separation
`11B ≤ |F|` is new. -/
theorem uniformThree_selectedWitness_of_cite_and_frequency_large
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (i₁ i₂ i₃ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → (3 : ℤ) ∣ v j)
    (a : Fin 3 → ℤ)
    (hv₁ : v i₁ = a 0) (hv₂ : v i₂ = a 1) (hv₃ : v i₃ = a 2)
    (hunit : ∀ i, ¬ (3 : ℤ) ∣ a i)
    (hresidue : ∀ i, (3 : ℤ) ∣ a i - 1)
    (B : ℝ) (hB0 : 0 < B)
    (hB : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ →
      |((v j / 3 : ℤ) : ℝ)| ≤ B)
    (hfrequency : threePhaseFrequency a ≠ 0)
    (hfrequencyLarge :
      11 * B ≤ |(threePhaseFrequency a : ℝ)|) :
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v 3 i₁ i₂ i₃ u ∧
      HasThreeDetunedGoodBranch (v i₁) (v i₂) (v i₃) 3 u := by
  apply selectedWitness_of_cite_and_threeFourteenths_scalarWall
    cite v hv 3 i₁ i₂ i₃ h12 h13 h23 hdvd B hB0 hB
      (threePhaseFrequency a) hfrequency hfrequencyLarge
  intro u hfail
  have hfailNormalized :
      ¬ HasThreeDetunedGoodBranch (a 0) (a 1) (a 2) 3 u := by
    simpa [hv₁, hv₂, hv₃] using hfail
  have hobstruction : ThreeClassCyclicObstruction a u :=
    (noThreeDetunedGoodBranch_three_iff_cyclicObstruction
      a u hunit).mp hfailNormalized
  exact frequency_bad_of_cyclicObstruction a u hresidue hobstruction

/-- Direct q244 socket in existing branch predicates.  The companion module's
parallel-partition scalar lemma supplies `hfailureWall`; this theorem performs
the sharp citation/interval escape. -/
theorem twoFourFour_selectedWitness_of_cite_and_frequency_large
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (i₂ i₄a i₄b : Fin 13)
    (h24a : i₂ ≠ i₄a) (h24b : i₂ ≠ i₄b) (h4ab : i₄a ≠ i₄b)
    (hdvd : ∀ j, j ≠ i₂ → j ≠ i₄a → j ≠ i₄b → g ∣ v j)
    (_hq2 : g / (Int.gcd (v i₂) g : ℤ) = 2)
    (_hq4a : g / (Int.gcd (v i₄a) g : ℤ) = 4)
    (_hq4b : g / (Int.gcd (v i₄b) g : ℤ) = 4)
    (B : ℝ) (hB0 : 0 < B)
    (hB : ∀ j, j ≠ i₂ → j ≠ i₄a → j ≠ i₄b →
      |((v j / g : ℤ) : ℝ)| ≤ B)
    (frequency : ℤ) (hfrequency : frequency ≠ 0)
    (hfrequencyLarge : 11 * B ≤ |(frequency : ℝ)|)
    (hfailureWall : ∀ u : ℝ,
      ¬ HasThreeDetunedGoodBranch (v i₂) (v i₄a) (v i₄b) g u →
      ∃ integer : ℤ,
        |(frequency : ℝ) * u - integer| < 3 / 14) :
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₂ i₄a i₄b u ∧
      HasThreeDetunedGoodBranch (v i₂) (v i₄a) (v i₄b) g u :=
  selectedWitness_of_cite_and_threeFourteenths_scalarWall
    cite v hv g i₂ i₄a i₄b h24a h24b h4ab hdvd B hB0 hB
      frequency hfrequency hfrequencyLarge hfailureWall

/-- Complete q244 frequency split.  Denominator data automatically produces a
normalized scalar wall.  A nonzero wall frequency of size at least `11B`
selects the witness; the only residual is an explicit zero/small support-three
frequency. -/
theorem twoFourFour_selectedWitness_or_small_normalized_frequency
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (hg : 2 ≤ g) (i₂ i₄a i₄b : Fin 13)
    (h24a : i₂ ≠ i₄a) (h24b : i₂ ≠ i₄b) (h4ab : i₄a ≠ i₄b)
    (hdvd : ∀ j, j ≠ i₂ → j ≠ i₄a → j ≠ i₄b → g ∣ v j)
    (hq2 : g / (Int.gcd (v i₂) g : ℤ) = 2)
    (hq4a : g / (Int.gcd (v i₄a) g : ℤ) = 4)
    (hq4b : g / (Int.gcd (v i₄b) g : ℤ) = 4)
    (B : ℝ) (hB0 : 0 < B)
    (hB : ∀ j, j ≠ i₂ → j ≠ i₄a → j ≠ i₄b →
      |((v j / g : ℤ) : ℝ)| ≤ B) :
    (∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₂ i₄a i₄b u ∧
      HasThreeDetunedGoodBranch (v i₂) (v i₄a) (v i₄b) g u) ∨
    ∃ a₂ a₄a a₄b : ℤ,
      (4 : ℤ) ∣ a₂ - 1 ∧ (4 : ℤ) ∣ a₄a - 1 ∧
      (4 : ℤ) ∣ a₄b - 1 ∧
      (∀ u : ℝ,
        ¬ HasThreeDetunedGoodBranch (v i₂) (v i₄a) (v i₄b) g u →
        ∃ n : ℤ,
          |(twoFourFourPhaseFrequency a₂ a₄a a₄b : ℝ) * u - n| < 3 / 14) ∧
      (twoFourFourPhaseFrequency a₂ a₄a a₄b = 0 ∨
        |(twoFourFourPhaseFrequency a₂ a₄a a₄b : ℝ)| < 11 * B) := by
  obtain ⟨a₂, a₄a, a₄b, hres₂, hres₄a, hres₄b, hwall⟩ :=
    exists_twoFourFour_normalized_failureWall
      (v i₂) (v i₄a) (v i₄b) g hg hq2 hq4a hq4b
  by_cases hzero : twoFourFourPhaseFrequency a₂ a₄a a₄b = 0
  · exact Or.inr ⟨a₂, a₄a, a₄b, hres₂, hres₄a, hres₄b,
      hwall, Or.inl hzero⟩
  by_cases hlarge : 11 * B ≤
      |(twoFourFourPhaseFrequency a₂ a₄a a₄b : ℝ)|
  · exact Or.inl (twoFourFour_selectedWitness_of_cite_and_frequency_large
      cite v hv g i₂ i₄a i₄b h24a h24b h4ab hdvd hq2 hq4a hq4b
        B hB0 hB (twoFourFourPhaseFrequency a₂ a₄a a₄b)
        hzero hlarge hwall)
  · exact Or.inr ⟨a₂, a₄a, a₄b, hres₂, hres₄a, hres₄b,
      hwall, Or.inr (lt_of_not_ge hlarge)⟩

end

end LRC14Grand
end LonelyRunner
