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

/-- Cancellation-robust companion socket.  A q-two bad row gives a scalar wall
of radius only `1/7`; the sharp citation interval escapes it as soon as
`3|frequency| ≥ 22B`.  This threshold is strictly weaker than the combined
q244 frequency threshold and remains useful when the three-row first moment
cancels. -/
theorem selectedWitness_of_cite_and_oneSeventh_scalarWall
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (i₁ i₂ i₃ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j)
    (B : ℝ) (hB0 : 0 < B)
    (hB : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ →
      |((v j / g : ℤ) : ℝ)| ≤ B)
    (frequency : ℤ) (hfrequency : frequency ≠ 0)
    (hfrequencyLarge : 22 * B ≤ 3 * |(frequency : ℝ)|)
    (hfailureWall : ∀ u : ℝ,
      ¬ HasThreeDetunedGoodBranch (v i₁) (v i₂) (v i₃) g u →
      ∃ integer : ℤ,
        |(frequency : ℝ) * u - integer| < 1 / 7) :
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₁ i₂ i₃ u ∧
      HasThreeDetunedGoodBranch (v i₁) (v i₂) (v i₃) g u := by
  obtain ⟨center, hinterval⟩ :=
    exists_threeDetunedHarmonicGoodAt_sharpInterval_of_cite
      cite v hv g i₁ i₂ i₃ h12 h13 h23 hdvd B hB0 hB
  apply selectedWitness_of_harmonicInterval_and_scalarFailureWall
    v g i₁ i₂ i₃ center (3 / (154 * B)) (1 / 7) frequency
  · positivity
  · norm_num
  · norm_num
  · exact hfrequency
  · exact hinterval
  · exact hfailureWall
  · have hBne : B ≠ 0 := ne_of_gt hB0
    have hscaled : (22 / 3 : ℝ) * B ≤ |(frequency : ℝ)| := by
      linarith
    calc
      2 * (1 / 7 : ℝ) =
          2 * (3 / (154 * B)) * ((22 / 3) * B) := by
            field_simp [hBne]
            ring
      _ ≤ 2 * (3 / (154 * B)) * |(frequency : ℝ)| :=
        mul_le_mul_of_nonneg_left hscaled (by positivity)

/-- Q-four analogue of the preceding cancellation-robust socket.  A normalized
q-four bad row has danger radius `2/7` and is escapable when
`3|frequency| ≥ 44B`. -/
theorem selectedWitness_of_cite_and_twoSevenths_scalarWall
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (i₁ i₂ i₃ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → g ∣ v j)
    (B : ℝ) (hB0 : 0 < B)
    (hB : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ →
      |((v j / g : ℤ) : ℝ)| ≤ B)
    (frequency : ℤ) (hfrequency : frequency ≠ 0)
    (hfrequencyLarge : 44 * B ≤ 3 * |(frequency : ℝ)|)
    (hfailureWall : ∀ u : ℝ,
      ¬ HasThreeDetunedGoodBranch (v i₁) (v i₂) (v i₃) g u →
      ∃ integer : ℤ,
        |(frequency : ℝ) * u - integer| < 2 / 7) :
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₁ i₂ i₃ u ∧
      HasThreeDetunedGoodBranch (v i₁) (v i₂) (v i₃) g u := by
  obtain ⟨center, hinterval⟩ :=
    exists_threeDetunedHarmonicGoodAt_sharpInterval_of_cite
      cite v hv g i₁ i₂ i₃ h12 h13 h23 hdvd B hB0 hB
  apply selectedWitness_of_harmonicInterval_and_scalarFailureWall
    v g i₁ i₂ i₃ center (3 / (154 * B)) (2 / 7) frequency
  · positivity
  · norm_num
  · norm_num
  · exact hfrequency
  · exact hinterval
  · exact hfailureWall
  · have hBne : B ≠ 0 := ne_of_gt hB0
    have hscaled : (44 / 3 : ℝ) * B ≤ |(frequency : ℝ)| := by
      linarith
    calc
      2 * (2 / 7 : ℝ) =
          2 * (3 / (154 * B)) * ((44 / 3) * B) := by
            field_simp [hBne]
            ring
      _ ≤ 2 * (3 / (154 * B)) * |(frequency : ℝ)| :=
        mul_le_mul_of_nonneg_left hscaled (by positivity)

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

/-- Sign-robust primitive q333 socket.  The selected speeds need only be units
modulo three; `a` is their coordinatewise sign normalization to residue one.
This removes the artificial requirement in the preceding theorem that the
original tuple already use those three signs. -/
theorem uniformThree_selectedWitness_of_cite_and_normalized_frequency_large
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (i₁ i₂ i₃ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → (3 : ℤ) ∣ v j)
    (hunit₁ : ¬ (3 : ℤ) ∣ v i₁)
    (hunit₂ : ¬ (3 : ℤ) ∣ v i₂)
    (hunit₃ : ¬ (3 : ℤ) ∣ v i₃)
    (ε a : Fin 3 → ℤ)
    (hε : ∀ i, ε i = 1 ∨ ε i = -1)
    (ha₁ : a 0 = ε 0 * v i₁)
    (ha₂ : a 1 = ε 1 * v i₂)
    (ha₃ : a 2 = ε 2 * v i₃)
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
  let p : Fin 3 → ℤ := ![v i₁, v i₂, v i₃]
  have hpunit : ∀ i, ¬ (3 : ℤ) ∣ p i := by
    intro i
    fin_cases i
    · simpa [p] using hunit₁
    · simpa [p] using hunit₂
    · simpa [p] using hunit₃
  have ha : ∀ i, a i = ε i * p i := by
    intro i
    fin_cases i
    · simpa [p] using ha₁
    · simpa [p] using ha₂
    · simpa [p] using ha₃
  have haEq : (fun i => ε i * p i) = a := by
    funext i
    exact (ha i).symm
  apply selectedWitness_of_cite_and_threeFourteenths_scalarWall
    cite v hv 3 i₁ i₂ i₃ h12 h13 h23 hdvd B hB0 hB
      (threePhaseFrequency a) hfrequency hfrequencyLarge
  intro u hfail
  have hpFail :
      ¬ HasThreeDetunedGoodBranch (p 0) (p 1) (p 2) 3 u := by
    simpa [p] using hfail
  have hpObstruction : ThreeClassCyclicObstruction p u :=
    (noThreeDetunedGoodBranch_three_iff_cyclicObstruction
      p u hpunit).mp hpFail
  have haObstruction : ThreeClassCyclicObstruction a u := by
    have hsigned :=
      (cyclicObstruction_sign_iff p ε u hε).mpr hpObstruction
    rwa [haEq] at hsigned
  exact frequency_bad_of_cyclicObstruction a u hresidue haObstruction

/-- Complete primitive q333 frequency split.  The normalization and scalar
wall are now produced from the three nondivisibility hypotheses themselves.
The zero branch exposes its exact unit-coefficient support-three relation;
otherwise the only residue is a nonzero integer frequency smaller than `11B`.
-/
theorem uniformThree_selectedWitness_or_small_normalized_frequency
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (i₁ i₂ i₃ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → (3 : ℤ) ∣ v j)
    (hunit₁ : ¬ (3 : ℤ) ∣ v i₁)
    (hunit₂ : ¬ (3 : ℤ) ∣ v i₂)
    (hunit₃ : ¬ (3 : ℤ) ∣ v i₃)
    (B : ℝ) (hB0 : 0 < B)
    (hB : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ →
      |((v j / 3 : ℤ) : ℝ)| ≤ B) :
    (∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v 3 i₁ i₂ i₃ u ∧
      HasThreeDetunedGoodBranch (v i₁) (v i₂) (v i₃) 3 u) ∨
    ∃ ε a : Fin 3 → ℤ,
      (∀ i, ε i = 1 ∨ ε i = -1) ∧
      a 0 = ε 0 * v i₁ ∧ a 1 = ε 1 * v i₂ ∧
      a 2 = ε 2 * v i₃ ∧
      (∀ i, (3 : ℤ) ∣ a i - 1) ∧
      (∀ u : ℝ,
        ¬ HasThreeDetunedGoodBranch (v i₁) (v i₂) (v i₃) 3 u →
        ∃ n : ℤ,
          |(threePhaseFrequency a : ℝ) * u - n| < 3 / 14) ∧
      ((threePhaseFrequency a = 0 ∧
          ε 0 * v i₁ + ε 1 * v i₂ + ε 2 * v i₃ = 0) ∨
        (threePhaseFrequency a ≠ 0 ∧
          |(threePhaseFrequency a : ℝ)| < 11 * B)) := by
  let p : Fin 3 → ℤ := ![v i₁, v i₂, v i₃]
  have hpunit : ∀ i, ¬ (3 : ℤ) ∣ p i := by
    intro i
    fin_cases i
    · simpa [p] using hunit₁
    · simpa [p] using hunit₂
    · simpa [p] using hunit₃
  obtain ⟨ε, a, hε, ha, hresidue, hpwall⟩ :=
    exists_uniformThree_normalized_failureWall p hpunit
  have ha₁ : a 0 = ε 0 * v i₁ := by simpa [p] using ha 0
  have ha₂ : a 1 = ε 1 * v i₂ := by simpa [p] using ha 1
  have ha₃ : a 2 = ε 2 * v i₃ := by simpa [p] using ha 2
  have hwall : ∀ u : ℝ,
      ¬ HasThreeDetunedGoodBranch (v i₁) (v i₂) (v i₃) 3 u →
      ∃ n : ℤ,
        |(threePhaseFrequency a : ℝ) * u - n| < 3 / 14 := by
    intro u hfail
    apply hpwall u
    simpa [p] using hfail
  by_cases hzero : threePhaseFrequency a = 0
  · have hrelationA : a 0 + a 1 + a 2 = 0 :=
      (threePhaseFrequency_eq_zero_iff_three_term_relation
        a hresidue).mp hzero
    have hrelation :
        ε 0 * v i₁ + ε 1 * v i₂ + ε 2 * v i₃ = 0 := by
      rw [← ha₁, ← ha₂, ← ha₃]
      exact hrelationA
    exact Or.inr ⟨ε, a, hε, ha₁, ha₂, ha₃, hresidue, hwall,
      Or.inl ⟨hzero, hrelation⟩⟩
  by_cases hlarge : 11 * B ≤ |(threePhaseFrequency a : ℝ)|
  · exact Or.inl
      (uniformThree_selectedWitness_of_cite_and_normalized_frequency_large
        cite v hv i₁ i₂ i₃ h12 h13 h23 hdvd hunit₁ hunit₂ hunit₃
        ε a hε ha₁ ha₂ ha₃ hresidue B hB0 hB hzero hlarge)
  · exact Or.inr ⟨ε, a, hε, ha₁, ha₂, ha₃, hresidue, hwall,
      Or.inr ⟨hzero, lt_of_not_ge hlarge⟩⟩

/-- Honest closed q333 subcase: after normalization, a coordinate whose
absolute-value excess over the other two is at least `33B` forces the large
frequency gate and hence selects a joint harmonic/detuned witness.  This is
the quantitative form naturally supplied by a sufficiently separated
chain-dense ladder top. -/
theorem uniformThree_selectedWitness_of_cite_and_normalized_third_gap
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (i₁ i₂ i₃ : Fin 13)
    (h12 : i₁ ≠ i₂) (h13 : i₁ ≠ i₃) (h23 : i₂ ≠ i₃)
    (hdvd : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ → (3 : ℤ) ∣ v j)
    (hunit₁ : ¬ (3 : ℤ) ∣ v i₁)
    (hunit₂ : ¬ (3 : ℤ) ∣ v i₂)
    (hunit₃ : ¬ (3 : ℤ) ∣ v i₃)
    (ε a : Fin 3 → ℤ)
    (hε : ∀ i, ε i = 1 ∨ ε i = -1)
    (ha₁ : a 0 = ε 0 * v i₁)
    (ha₂ : a 1 = ε 1 * v i₂)
    (ha₃ : a 2 = ε 2 * v i₃)
    (hresidue : ∀ i, (3 : ℤ) ∣ a i - 1)
    (B : ℝ) (hB0 : 0 < B)
    (hB : ∀ j, j ≠ i₁ → j ≠ i₂ → j ≠ i₃ →
      |((v j / 3 : ℤ) : ℝ)| ≤ B)
    (hgap : 33 * B ≤
      |(a 2 : ℝ)| - |(a 0 : ℝ)| - |(a 1 : ℝ)|) :
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v 3 i₁ i₂ i₃ u ∧
      HasThreeDetunedGoodBranch (v i₁) (v i₂) (v i₃) 3 u := by
  have hlarge : 11 * B ≤ |(threePhaseFrequency a : ℝ)| :=
    threePhaseFrequency_large_of_third_gap a hresidue B hgap
  have hnonzero : threePhaseFrequency a ≠ 0 := by
    intro hzero
    rw [hzero] at hlarge
    norm_num at hlarge
    linarith
  exact uniformThree_selectedWitness_of_cite_and_normalized_frequency_large
    cite v hv i₁ i₂ i₃ h12 h13 h23 hdvd hunit₁ hunit₂ hunit₃
      ε a hε ha₁ ha₂ ha₃ hresidue B hB0 hB hnonzero hlarge

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

/-- A relation-lattice formulation of the q244 large-frequency exit.  The
normalized frequency is canonically one of the eight signed sums of the three
selected speeds divided by `g`.  Thus a gap of `11*g*B` for those signed sums
closes the branch without mentioning the auxiliary normalized numerators.

The universal quantifier only ranges over signs in `{+1,-1}`; it is a finite
four-check condition up to simultaneous sign reversal. -/
theorem twoFourFour_selectedWitness_of_cite_and_signed_sum_gap
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (hg : 2 ≤ g) (i₂ i₄a i₄b : Fin 13)
    (h24a : i₂ ≠ i₄a) (h24b : i₂ ≠ i₄b) (h4ab : i₄a ≠ i₄b)
    (hdvd : ∀ j, j ≠ i₂ → j ≠ i₄a → j ≠ i₄b → g ∣ v j)
    (hq2 : g / (Int.gcd (v i₂) g : ℤ) = 2)
    (hq4a : g / (Int.gcd (v i₄a) g : ℤ) = 4)
    (hq4b : g / (Int.gcd (v i₄b) g : ℤ) = 4)
    (B : ℝ) (hB0 : 0 < B)
    (hB : ∀ j, j ≠ i₂ → j ≠ i₄a → j ≠ i₄b →
      |((v j / g : ℤ) : ℝ)| ≤ B)
    (hsignedGap : ∀ signTwo signFourA signFourB : ℤ,
      (signTwo = 1 ∨ signTwo = -1) →
      (signFourA = 1 ∨ signFourA = -1) →
      (signFourB = 1 ∨ signFourB = -1) →
      (g : ℝ) * (11 * B) ≤
        |((signTwo * v i₂ + signFourA * v i₄a +
          signFourB * v i₄b : ℤ) : ℝ)|) :
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₂ i₄a i₄b u ∧
      HasThreeDetunedGoodBranch (v i₂) (v i₄a) (v i₄b) g u := by
  obtain ⟨signTwo, signFourA, signFourB,
      numeratorTwo, numeratorFourA, numeratorFourB,
      hsignTwo, hsignFourA, hsignFourB,
      -, -, -, -, -, -, hidentity, -, -, -, hwall⟩ :=
    exists_twoFourFour_normalized_signed_failureWall
      (v i₂) (v i₄a) (v i₄b) g hg hq2 hq4a hq4b
  let frequency : ℤ := twoFourFourPhaseFrequency
    numeratorTwo numeratorFourA numeratorFourB
  let signedSum : ℤ :=
    signTwo * v i₂ + signFourA * v i₄a + signFourB * v i₄b
  have hidentity' : signedSum = g * frequency := hidentity
  have hgap : (g : ℝ) * (11 * B) ≤ |(signedSum : ℝ)| :=
    hsignedGap signTwo signFourA signFourB
      hsignTwo hsignFourA hsignFourB
  have hgReal : (0 : ℝ) < g := by exact_mod_cast (show (0 : ℤ) < g by omega)
  have hfrequency : frequency ≠ 0 := by
    intro hzero
    have hsumZero : signedSum = 0 := by rw [hidentity', hzero, mul_zero]
    rw [hsumZero, Int.cast_zero, abs_zero] at hgap
    have hgapPositive : (0 : ℝ) < (g : ℝ) * (11 * B) := by positivity
    linarith
  have hidentityReal : (signedSum : ℝ) = (g : ℝ) * (frequency : ℝ) := by
    exact_mod_cast hidentity'
  have hfrequencyLarge : 11 * B ≤ |(frequency : ℝ)| := by
    rw [hidentityReal, abs_mul, abs_of_pos hgReal] at hgap
    exact le_of_mul_le_mul_left hgap hgReal
  apply twoFourFour_selectedWitness_of_cite_and_frequency_large
    cite v hv g i₂ i₄a i₄b h24a h24b h4ab hdvd
      hq2 hq4a hq4b B hB0 hB frequency hfrequency hfrequencyLarge
  exact hwall

/-- Cancellation-robust q244 exit from the q-two row alone.  Since
`2|v(i₂)| = g|a₂|`, the concrete speed gap

`11*g*B ≤ 3|v(i₂)|`

is exactly the threshold `22B ≤ 3|a₂|` needed for the q-two wall of radius
`1/7` to cross the citation interval.  This closes cases in which the combined
frequency is small or zero because of cancellation but the q-two numerator is
large. -/
theorem twoFourFour_selectedWitness_of_cite_and_qTwo_speed_gap
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (g : ℤ) (hg : 2 ≤ g) (i₂ i₄a i₄b : Fin 13)
    (h24a : i₂ ≠ i₄a) (h24b : i₂ ≠ i₄b) (h4ab : i₄a ≠ i₄b)
    (hdvd : ∀ j, j ≠ i₂ → j ≠ i₄a → j ≠ i₄b → g ∣ v j)
    (hq2 : g / (Int.gcd (v i₂) g : ℤ) = 2)
    (hq4a : g / (Int.gcd (v i₄a) g : ℤ) = 4)
    (hq4b : g / (Int.gcd (v i₄b) g : ℤ) = 4)
    (B : ℝ) (hB0 : 0 < B)
    (hB : ∀ j, j ≠ i₂ → j ≠ i₄a → j ≠ i₄b →
      |((v j / g : ℤ) : ℝ)| ≤ B)
    (hspeedGap : 11 * (g : ℝ) * B ≤ 3 * |(v i₂ : ℝ)|) :
    ∃ u : ℝ,
      ThreeDetunedHarmonicGoodAt v g i₂ i₄a i₄b u ∧
      HasThreeDetunedGoodBranch (v i₂) (v i₄a) (v i₄b) g u := by
  obtain ⟨signTwo, signFourA, signFourB,
      numeratorTwo, numeratorFourA, numeratorFourB,
      hsignTwo, -, -, -, -, -, hscaleTwo, -, -, -, hqTwoWall, -, -, -⟩ :=
    exists_twoFourFour_normalized_signed_failureWall
      (v i₂) (v i₄a) (v i₄b) g hg hq2 hq4a hq4b
  have habsScale :
      (2 : ℝ) * |(v i₂ : ℝ)| = (g : ℝ) * |(numeratorTwo : ℝ)| :=
    normalizedQTwo_abs_scale (v i₂) g signTwo numeratorTwo
      (by omega) hsignTwo hscaleTwo
  have hgReal : (0 : ℝ) < g := by exact_mod_cast (show (0 : ℤ) < g by omega)
  have hfrequency : numeratorTwo ≠ 0 := by
    intro hzero
    rw [hzero, Int.cast_zero, abs_zero, mul_zero] at habsScale
    have hspeedReal : (v i₂ : ℝ) ≠ 0 := by exact_mod_cast hv i₂
    have hspeedAbs : (0 : ℝ) < |(v i₂ : ℝ)| := abs_pos.mpr hspeedReal
    linarith
  have hfrequencyLarge : 22 * B ≤ 3 * |(numeratorTwo : ℝ)| := by
    have hscaled :
        (g : ℝ) * (22 * B) ≤
          (g : ℝ) * (3 * |(numeratorTwo : ℝ)|) := by
      nlinarith [hspeedGap, habsScale]
    exact le_of_mul_le_mul_left hscaled hgReal
  exact selectedWitness_of_cite_and_oneSeventh_scalarWall
    cite v hv g i₂ i₄a i₄b h24a h24b h4ab hdvd B hB0 hB
      numeratorTwo hfrequency hfrequencyLarge hqTwoWall

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

/-- Arithmetic sharpening of the q244 residual after both scalar exits.  The
q-two numerator is forced below `22B/3`.  The zero combined-frequency branch
is an exact signed support-three relation among the selected speeds; otherwise
the same signed sum is a nonzero multiple of `g`, bounded between `g` and
`11*g*B`.  The normalized scalars are linked back to all three original speeds
by exact scale identities. -/
theorem twoFourFour_selectedWitness_or_signed_small_relation
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
    ∃ signTwo signFourA signFourB
        numeratorTwo numeratorFourA numeratorFourB : ℤ,
      (signTwo = 1 ∨ signTwo = -1) ∧
      (signFourA = 1 ∨ signFourA = -1) ∧
      (signFourB = 1 ∨ signFourB = -1) ∧
      (4 : ℤ) ∣ numeratorTwo - 1 ∧
      (4 : ℤ) ∣ numeratorFourA - 1 ∧
      (4 : ℤ) ∣ numeratorFourB - 1 ∧
      2 * (signTwo * v i₂) = g * numeratorTwo ∧
      4 * (signFourA * v i₄a) = g * numeratorFourA ∧
      4 * (signFourB * v i₄b) = g * numeratorFourB ∧
      signTwo * v i₂ + signFourA * v i₄a + signFourB * v i₄b =
        g * twoFourFourPhaseFrequency
          numeratorTwo numeratorFourA numeratorFourB ∧
      (∀ u : ℝ,
        ¬ HasThreeDetunedGoodBranch (v i₂) (v i₄a) (v i₄b) g u →
        ∃ integer : ℤ,
          |(twoFourFourPhaseFrequency
              numeratorTwo numeratorFourA numeratorFourB : ℝ) * u - integer| <
            3 / 14) ∧
      3 * |(numeratorTwo : ℝ)| < 22 * B ∧
      (signTwo * v i₂ + signFourA * v i₄a + signFourB * v i₄b = 0 ∨
        (signTwo * v i₂ + signFourA * v i₄a +
            signFourB * v i₄b ≠ 0 ∧
          (g : ℝ) ≤
            |((signTwo * v i₂ + signFourA * v i₄a +
              signFourB * v i₄b : ℤ) : ℝ)| ∧
          |((signTwo * v i₂ + signFourA * v i₄a +
            signFourB * v i₄b : ℤ) : ℝ)| <
              (g : ℝ) * (11 * B))) := by
  obtain ⟨signTwo, signFourA, signFourB,
      numeratorTwo, numeratorFourA, numeratorFourB,
      hsignTwo, hsignFourA, hsignFourB,
      hresTwo, hresFourA, hresFourB,
      hscaleTwo, hscaleFourA, hscaleFourB, hidentity,
      hqTwoWall, -, -, hwall⟩ :=
    exists_twoFourFour_normalized_signed_failureWall
      (v i₂) (v i₄a) (v i₄b) g hg hq2 hq4a hq4b
  let frequency : ℤ := twoFourFourPhaseFrequency
    numeratorTwo numeratorFourA numeratorFourB
  let signedSum : ℤ :=
    signTwo * v i₂ + signFourA * v i₄a + signFourB * v i₄b
  have hidentity' : signedSum = g * frequency := hidentity
  have habsScale :
      (2 : ℝ) * |(v i₂ : ℝ)| = (g : ℝ) * |(numeratorTwo : ℝ)| :=
    normalizedQTwo_abs_scale (v i₂) g signTwo numeratorTwo
      (by omega) hsignTwo hscaleTwo
  have hnumeratorTwo : numeratorTwo ≠ 0 := by
    intro hzero
    rw [hzero, Int.cast_zero, abs_zero, mul_zero] at habsScale
    have hspeedReal : (v i₂ : ℝ) ≠ 0 := by exact_mod_cast hv i₂
    linarith [abs_pos.mpr hspeedReal]
  by_cases hqTwoLarge : 22 * B ≤ 3 * |(numeratorTwo : ℝ)|
  · exact Or.inl (selectedWitness_of_cite_and_oneSeventh_scalarWall
      cite v hv g i₂ i₄a i₄b h24a h24b h4ab hdvd B hB0 hB
        numeratorTwo hnumeratorTwo hqTwoLarge hqTwoWall)
  have hqTwoSmall : 3 * |(numeratorTwo : ℝ)| < 22 * B :=
    lt_of_not_ge hqTwoLarge
  by_cases hzero : frequency = 0
  · right
    refine ⟨signTwo, signFourA, signFourB,
      numeratorTwo, numeratorFourA, numeratorFourB,
      hsignTwo, hsignFourA, hsignFourB,
      hresTwo, hresFourA, hresFourB,
      hscaleTwo, hscaleFourA, hscaleFourB, hidentity, hwall,
      hqTwoSmall, Or.inl ?_⟩
    exact show signedSum = 0 by rw [hidentity', hzero, mul_zero]
  by_cases hlarge : 11 * B ≤ |(frequency : ℝ)|
  · left
    exact twoFourFour_selectedWitness_of_cite_and_frequency_large
      cite v hv g i₂ i₄a i₄b h24a h24b h4ab hdvd
        hq2 hq4a hq4b B hB0 hB frequency hzero hlarge hwall
  · right
    have hsmall : |(frequency : ℝ)| < 11 * B := lt_of_not_ge hlarge
    have hsumNonzero : signedSum ≠ 0 := by
      intro hsumZero
      have hproductZero : g * frequency = 0 := by rw [← hidentity', hsumZero]
      exact hzero ((mul_eq_zero.mp hproductZero).resolve_left (by omega))
    have hgReal : (0 : ℝ) < g := by
      exact_mod_cast (show (0 : ℤ) < g by omega)
    have hidentityReal : (signedSum : ℝ) = (g : ℝ) * (frequency : ℝ) := by
      exact_mod_cast hidentity'
    have hfrequencyAbsOne : (1 : ℝ) ≤ |(frequency : ℝ)| := by
      rcases lt_or_gt_of_ne hzero with hnegative | hpositive
      · rw [abs_of_neg]
        · exact_mod_cast (show (1 : ℤ) ≤ -frequency by omega)
        · exact_mod_cast hnegative
      · rw [abs_of_pos]
        · exact_mod_cast (show (1 : ℤ) ≤ frequency by omega)
        · exact_mod_cast hpositive
    have hsumLower : (g : ℝ) ≤ |(signedSum : ℝ)| := by
      rw [hidentityReal, abs_mul, abs_of_pos hgReal]
      simpa only [mul_one] using
        mul_le_mul_of_nonneg_left hfrequencyAbsOne (le_of_lt hgReal)
    have hsumSmall : |(signedSum : ℝ)| < (g : ℝ) * (11 * B) := by
      rw [hidentityReal, abs_mul, abs_of_pos hgReal]
      exact mul_lt_mul_of_pos_left hsmall hgReal
    refine ⟨signTwo, signFourA, signFourB,
      numeratorTwo, numeratorFourA, numeratorFourB,
      hsignTwo, hsignFourA, hsignFourB,
      hresTwo, hresFourA, hresFourB,
      hscaleTwo, hscaleFourA, hscaleFourB, hidentity, hwall,
      hqTwoSmall, Or.inr ?_⟩
    exact ⟨hsumNonzero, hsumLower, hsumSmall⟩

end

end LRC14Grand
end LonelyRunner
