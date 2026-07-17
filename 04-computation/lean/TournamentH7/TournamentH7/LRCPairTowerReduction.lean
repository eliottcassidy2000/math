/-
  TournamentH7.LRCPairTowerReduction

  The first exact finite-layer reduction of the residual q-two pair tower.
  A nongeneric two-row level is necessarily `(2,2)`.  At the doubled modulus,
  every old q-two row becomes q-four, while every newly exposed odd harmonic
  row is q-two.  Consequently a lift exposing exactly one new row is the
  already isolated `(2,4,4)` selected-witness problem.  The genuinely new pair
  residual starts only when the first lift has at least four nonmultiples.

  Tournament-analysis audit: vertices are valuation layers and speed indices,
  not runners-as-points or circle arcs.  The observable is whether an index
  crosses the divisibility wall between `g` and `2g`; thresholding gives a
  layer-membership relation, not an orientation.  An index-order gauge could
  manufacture a tournament and a tie Hamiltonian path, but would discard the
  q-two/q-four transition used by the proof.  The valuation quotient preserves
  exact nonmultiple counts and reduced-denominator transitions; it destroys
  phase placement on the parallel-class circle.  The challenged assumption is
  that a failed first lift creates an infinite tower: finite tuples instead
  create a strictly larger exposed layer, and the singleton layer closes here.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCSelectedWitnessObstructions

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open scoped Classical

/-- Coordinates which leave the divisible core when the modulus doubles. -/
def liftFailureCard (v : Fin 13 → ℤ) (g : ℤ) : ℕ :=
  (Finset.univ.filter fun i => g ∣ v i ∧ ¬ 2 * g ∣ v i).card

/-- Exact partition of the doubled-modulus nonmultiples into old
nonmultiples and newly exposed coordinates. -/
theorem nonMultCard_double_eq_add_liftFailureCard
    (v : Fin 13 → ℤ) (g : ℤ) :
    nonMultCard v (2 * g) = nonMultCard v g + liftFailureCard v g := by
  let old := Finset.univ.filter fun i => ¬ g ∣ v i
  let fresh := Finset.univ.filter fun i => g ∣ v i ∧ ¬ 2 * g ∣ v i
  have hdisjoint : Disjoint old fresh := by
    rw [Finset.disjoint_left]
    intro i hiOld hiFresh
    have hiOld' : ¬ g ∣ v i := by simpa [old] using hiOld
    have hiFresh' : g ∣ v i ∧ ¬ 2 * g ∣ v i := by
      simpa [fresh] using hiFresh
    exact hiOld' hiFresh'.1
  have hunion :
      Finset.univ.filter (fun i => ¬ 2 * g ∣ v i) = old ∪ fresh := by
    ext i
    simp only [Finset.mem_filter, Finset.mem_univ, true_and,
      Finset.mem_union, old, fresh]
    constructor
    · intro hnot
      by_cases hdiv : g ∣ v i
      · exact Or.inr ⟨hdiv, hnot⟩
      · exact Or.inl hdiv
    · rintro (hnot | ⟨-, hnot⟩)
      · intro hdouble
        exact hnot (dvd_trans ⟨2, by ring⟩ hdouble)
      · exact hnot
  unfold nonMultCard liftFailureCard
  rw [hunion, Finset.card_union_of_disjoint hdisjoint]

/-- Nontermination exposes at least one coordinate on the next 2-adic sheet. -/
theorem liftFailureCard_pos_of_not_twoAdicLiftTerminates
    (v : Fin 13 → ℤ) (g : ℤ)
    (hnonterm : ¬ TwoAdicLiftTerminates v g) :
    0 < liftFailureCard v g := by
  unfold TwoAdicLiftTerminates at hnonterm
  push Not at hnonterm
  obtain ⟨i, hdiv, hnot⟩ := hnonterm
  rw [liftFailureCard, Finset.card_pos]
  exact ⟨i, by simp [hdiv, hnot]⟩

/-- Exact first-step pair-tower dichotomy. -/
theorem pairTower_firstStep_dichotomy
    (v : Fin 13 → ℤ) (g : ℤ)
    (hcard : nonMultCard v g = 2)
    (hnonterm : ¬ TwoAdicLiftTerminates v g) :
    (liftFailureCard v g = 1 ∧ nonMultCard v (2 * g) = 3) ∨
      (2 ≤ liftFailureCard v g ∧ 4 ≤ nonMultCard v (2 * g)) := by
  have hpos := liftFailureCard_pos_of_not_twoAdicLiftTerminates v g hnonterm
  have heq := nonMultCard_double_eq_add_liftFailureCard v g
  rw [hcard] at heq
  by_cases hone : liftFailureCard v g = 1
  · exact Or.inl ⟨hone, by omega⟩
  · exact Or.inr ⟨by omega, by omega⟩

/-- A q-two row at modulus `g` becomes exactly q-four at modulus `2g`. -/
theorem q_two_lifts_to_q_four (delta g : ℤ) (hg : 2 ≤ g)
    (hq : g / (Int.gcd delta g : ℤ) = 2) (hdelta : ¬ g ∣ delta) :
    2 * g / (Int.gcd delta (2 * g) : ℤ) = 4 := by
  let divisor : ℤ := (Int.gcd delta g : ℤ)
  let liftedDivisor : ℤ := (Int.gcd delta (2 * g) : ℤ)
  have hdivisorPos : 0 < divisor := by
    dsimp [divisor]
    exact_mod_cast (show 0 < Int.gcd delta g by
      rw [Int.gcd_pos_iff]
      right
      omega)
  have hliftedPos : 0 < liftedDivisor := by
    dsimp [liftedDivisor]
    exact_mod_cast (show 0 < Int.gcd delta (2 * g) by
      rw [Int.gcd_pos_iff]
      right
      omega)
  have hgFactor : g = divisor * 2 := by
    have hfactor := Int.mul_ediv_cancel' (Int.gcd_dvd_right delta g)
    rw [hq] at hfactor
    simpa [divisor, mul_comm] using hfactor.symm
  have hdividesNat : Int.gcd delta g ∣ Int.gcd delta (2 * g) := by
    rw [Int.dvd_gcd_iff]
    exact ⟨Int.gcd_dvd_left delta g,
      dvd_trans (Int.gcd_dvd_right delta g) ⟨2, by ring⟩⟩
  have hdivides : divisor ∣ liftedDivisor := by
    obtain ⟨factor, hfactor⟩ := hdividesNat
    refine ⟨factor, ?_⟩
    dsimp [divisor, liftedDivisor]
    exact_mod_cast hfactor
  obtain ⟨factor, hfactor⟩ := hdivides
  have hfactorPos : 0 < factor := by nlinarith
  have hliftedDvd : liftedDivisor ∣ 2 * g :=
    Int.gcd_dvd_right delta (2 * g)
  obtain ⟨cofactor, hcofactor⟩ := hliftedDvd
  have hcofactorPos : 0 < cofactor := by nlinarith
  have hcancel : divisor * 4 = divisor * (factor * cofactor) := by
    calc
      divisor * 4 = 2 * g := by rw [hgFactor]; ring
      _ = liftedDivisor * cofactor := hcofactor
      _ = divisor * (factor * cofactor) := by rw [hfactor]; ring
  have hproduct : 4 = factor * cofactor :=
    mul_left_cancel₀ (ne_of_gt hdivisorPos) hcancel
  have hfactorLe : factor ≤ 4 := by
    have hnonneg : 0 ≤ (cofactor - 1) * factor :=
      mul_nonneg (by omega) (by omega)
    nlinarith
  have hfactorCases : factor = 1 ∨ factor = 2 ∨ factor = 4 := by
    interval_cases factor <;> omega
  have hfactorOne : factor = 1 := by
    rcases hfactorCases with rfl | rfl | rfl
    · rfl
    · exfalso
      apply hdelta
      have hliftedDelta : liftedDivisor ∣ delta :=
        Int.gcd_dvd_left delta (2 * g)
      have hliftedEq : liftedDivisor = g := by rw [hfactor, hgFactor]
      rwa [hliftedEq] at hliftedDelta
    · exfalso
      apply hdelta
      have hliftedDelta : liftedDivisor ∣ delta :=
        Int.gcd_dvd_left delta (2 * g)
      exact dvd_trans ⟨2, by rw [hfactor, hgFactor]; ring⟩ hliftedDelta
  have hliftedEq : liftedDivisor = divisor := by
    rw [hfactor, hfactorOne, mul_one]
  apply Int.ediv_eq_of_eq_mul_left (ne_of_gt hliftedPos)
  rw [hliftedEq, hgFactor]
  ring

/-- A harmonic row which crosses the `g`/`2g` wall is exactly q-two at the
doubled modulus. -/
theorem odd_harmonic_lifts_to_q_two (delta g : ℤ) (hg : 2 ≤ g)
    (hdiv : g ∣ delta) (hnot : ¬ 2 * g ∣ delta) :
    2 * g / (Int.gcd delta (2 * g) : ℤ) = 2 := by
  let liftedDivisor : ℤ := (Int.gcd delta (2 * g) : ℤ)
  have hliftedPos : 0 < liftedDivisor := by
    dsimp [liftedDivisor]
    exact_mod_cast (show 0 < Int.gcd delta (2 * g) by
      rw [Int.gcd_pos_iff]
      right
      omega)
  have hgNonneg : 0 ≤ g := by omega
  have hgCast : (g.toNat : ℤ) = g := by rw [Int.toNat_of_nonneg hgNonneg]
  have hgDvdNat : g.toNat ∣ Int.gcd delta (2 * g) := by
    rw [Int.dvd_gcd_iff]
    constructor
    · simpa [hgCast] using hdiv
    · rw [hgCast]
      exact ⟨2, by ring⟩
  have hgDvd : g ∣ liftedDivisor := by
    obtain ⟨factor, hfactor⟩ := hgDvdNat
    refine ⟨factor, ?_⟩
    dsimp [liftedDivisor]
    have hfactorCast : (Int.gcd delta (2 * g) : ℤ) =
        (g.toNat : ℤ) * (factor : ℤ) := by
      exact_mod_cast hfactor
    simpa [hgCast] using hfactorCast
  obtain ⟨factor, hfactor⟩ := hgDvd
  have hfactorPos : 0 < factor := by nlinarith
  have hliftedDvd : liftedDivisor ∣ 2 * g :=
    Int.gcd_dvd_right delta (2 * g)
  obtain ⟨cofactor, hcofactor⟩ := hliftedDvd
  have hcofactorPos : 0 < cofactor := by nlinarith
  have hcancel : g * 2 = g * (factor * cofactor) := by
    calc
      g * 2 = 2 * g := by ring
      _ = liftedDivisor * cofactor := hcofactor
      _ = g * (factor * cofactor) := by rw [hfactor]; ring
  have hproduct : 2 = factor * cofactor :=
    mul_left_cancel₀ (by omega : g ≠ 0) hcancel
  have hfactorLe : factor ≤ 2 := by
    have hnonneg : 0 ≤ (cofactor - 1) * factor :=
      mul_nonneg (by omega) (by omega)
    nlinarith
  have hfactorCases : factor = 1 ∨ factor = 2 := by omega
  have hfactorOne : factor = 1 := by
    rcases hfactorCases with rfl | rfl
    · rfl
    · exfalso
      apply hnot
      have hliftedDelta : liftedDivisor ∣ delta :=
        Int.gcd_dvd_left delta (2 * g)
      have hliftedEq : liftedDivisor = 2 * g := by rw [hfactor]; ring
      rwa [hliftedEq] at hliftedDelta
  have hliftedEq : liftedDivisor = g := by
    rw [hfactor, hfactorOne, mul_one]
  apply Int.ediv_eq_of_eq_mul_left (ne_of_gt hliftedPos)
  rw [hliftedEq]

private theorem two_mul_badCount_lt_of_three_le
    (delta g : ℤ) (hg : 0 < g)
    (hq : 3 ≤ g / (Int.gcd delta g : ℤ)) :
    2 * DetunedD3.badCount delta g < g.toNat := by
  by_cases hthree : g / (Int.gcd delta g : ℤ) = 3
  · have hcount := three_mul_badCount_eq delta g hg hthree
    omega
  · have hfour : 4 ≤ g / (Int.gcd delta g : ℤ) := by omega
    have hcount := seven_mul_badCount_le_two_mul delta g hg hfour
    have hgNat : 0 < g.toNat := by omega
    omega

private theorem two_mul_badCount_le_of_two_le
    (delta g : ℤ) (hg : 0 < g)
    (hq : 2 ≤ g / (Int.gcd delta g : ℤ)) :
    2 * DetunedD3.badCount delta g ≤ g.toNat := by
  by_cases htwo : g / (Int.gcd delta g : ℤ) = 2
  · exact (two_mul_badCount_eq delta g hg htwo).le
  · exact (two_mul_badCount_lt_of_three_le delta g hg (by omega)).le

/-- Exact denominator normal form for a nongeneric two-row level. -/
theorem nonGeneric_pair_q_two
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g)
    (hcard : nonMultCard v g = 2) (hnongeneric : ¬ genericCount v g) :
    ∃ first second : Fin 13, first ≠ second ∧
      (∀ index, index ≠ first → index ≠ second → g ∣ v index) ∧
      ¬ g ∣ v first ∧ ¬ g ∣ v second ∧
      g / (Int.gcd (v first) g : ℤ) = 2 ∧
      g / (Int.gcd (v second) g : ℤ) = 2 := by
  have hgPos : 0 < g := by omega
  rw [nonMultCard] at hcard
  obtain ⟨first, second, hne, hfilter⟩ := Finset.card_eq_two.mp hcard
  have hmem : ∀ index,
      index ∈ Finset.univ.filter (fun i => ¬ g ∣ v i) ↔ ¬ g ∣ v index := by
    intro index
    simp
  have hfirst : ¬ g ∣ v first := (hmem first).mp (by rw [hfilter]; simp)
  have hsecond : ¬ g ∣ v second := (hmem second).mp (by rw [hfilter]; simp)
  have hdvd : ∀ index, index ≠ first → index ≠ second → g ∣ v index := by
    intro index hindexFirst hindexSecond
    by_contra hnot
    have hindex := (hmem index).mpr hnot
    rw [hfilter] at hindex
    simp only [Finset.mem_insert, Finset.mem_singleton] at hindex
    rcases hindex with h | h
    · exact hindexFirst h
    · exact hindexSecond h
  have hqFirst : 2 ≤ g / (Int.gcd (v first) g : ℤ) :=
    reducedDetuningDenominator_ge_two hg hfirst
  have hqSecond : 2 ≤ g / (Int.gcd (v second) g : ℤ) :=
    reducedDetuningDenominator_ge_two hg hsecond
  have hsum : DetunedD3.badCount (v first) g +
      DetunedD3.badCount (v second) g ≥ g.toNat := by
    rw [genericCount, hfilter,
      Finset.sum_insert (by simp [hne]), Finset.sum_singleton] at hnongeneric
    omega
  have hqFirstEq : g / (Int.gcd (v first) g : ℤ) = 2 := by
    by_contra hnot
    have hboundFirst := two_mul_badCount_lt_of_three_le
      (v first) g hgPos (by omega)
    have hboundSecond := two_mul_badCount_le_of_two_le
      (v second) g hgPos hqSecond
    omega
  have hqSecondEq : g / (Int.gcd (v second) g : ℤ) = 2 := by
    by_contra hnot
    have hboundFirst := two_mul_badCount_le_of_two_le
      (v first) g hgPos hqFirst
    have hboundSecond := two_mul_badCount_lt_of_three_le
      (v second) g hgPos (by omega)
    omega
  exact ⟨first, second, hne, hdvd, hfirst, hsecond, hqFirstEq, hqSecondEq⟩

/-- The remaining pair supplier after the first lift exposes at least two new
odd harmonic rows, hence at least four nonmultiples in total. -/
def FourOrMorePairLiftSupply : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → ∀ g : ℤ, 2 ≤ g →
    nonMultCard v g = 2 → ¬ TwoAdicLiftTerminates v g →
    ¬ genericCount v g → 4 ≤ nonMultCard v (2 * g) →
    ∃ t : ℝ, Lonely 14 v t

/-- Sharp residual supplier: at least two coordinates cross the first
divisibility wall.  The doubled cardinality bound is then automatic. -/
def ManyLiftFailurePairTowerSupply : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → ∀ g : ℤ, 2 ≤ g →
    nonMultCard v g = 2 → 2 ≤ liftFailureCard v g →
    ¬ genericCount v g → ∃ t : ℝ, Lonely 14 v t

/-- Harmonic safety on coordinates which remain in the divisible core after
doubling. -/
def PairTowerCoreGoodAt (v : Fin 13 → ℤ) (g : ℤ) (u : ℝ) : Prop :=
  ∀ j, 2 * g ∣ v j → ∀ n : ℤ,
    (1 : ℝ) / 14 ≤ |((v j / (2 * g) : ℤ) : ℝ) * u - n|

/-- Exact fixed-phase obstruction on the first doubled sheet.  Either two
active fresh q-two rows occupy opposite parity sheets, or one fresh q-two row
and the inherited q-four pair form the saturated `(2,4,4)` prefix code. -/
def PairTowerParallelObstruction
    (v : Fin 13 → ℤ) (g : ℤ) (first second : Fin 13) (u : ℝ) : Prop :=
  (∃ i j : Fin 13,
      (g ∣ v i ∧ ¬ 2 * g ∣ v i) ∧
      (g ∣ v j ∧ ¬ 2 * g ∣ v j) ∧ i ≠ j ∧
      TwoTwoThreePhaseOpposition (v i) (v j) (2 * g) u) ∨
  (∃ i : Fin 13,
      (g ∣ v i ∧ ¬ 2 * g ∣ v i) ∧
      ¬ HasThreeDetunedGoodBranch
        (v i) (v first) (v second) (2 * g) u)

/-- Minimal dynamic pair-tower crux.  Vertices are harmonic-safe components
and divisibility-wall events; the selector must avoid both complete binary
prefix codes on one component. -/
def ManyLiftFailurePhaseSelector : Prop :=
  ∀ (v : Fin 13 → ℤ) (g : ℤ), 2 ≤ g →
    ∀ first second : Fin 13,
    first ≠ second →
    ¬ g ∣ v first → ¬ g ∣ v second →
    (∀ j, j ≠ first → j ≠ second → g ∣ v j) →
    g / (Int.gcd (v first) g : ℤ) = 2 →
    g / (Int.gcd (v second) g : ℤ) = 2 →
    2 ≤ liftFailureCard v g →
    ∃ u : ℝ,
      PairTowerCoreGoodAt v g u ∧
      ¬ PairTowerParallelObstruction v g first second u

/-- The wall-event phase selector implies the many-crossing loneliness
supplier.  This is the exact fixed-phase block-gluing consumer. -/
theorem manyLiftFailurePairTowerSupply_of_phaseSelector
    (hselector : ManyLiftFailurePhaseSelector) :
    ManyLiftFailurePairTowerSupply := by
  intro v _ g hg hcard hmany hnongeneric
  obtain ⟨first, second, hne, hdvd, hfirst, hsecond,
      hqFirst, hqSecond⟩ :=
    nonGeneric_pair_q_two v g hg hcard hnongeneric
  obtain ⟨u, hcore, hnoObstruction⟩ :=
    hselector v g hg first second hne hfirst hsecond hdvd
      hqFirst hqSecond hmany
  have hfreshNonempty :
      (Finset.univ.filter fun i => g ∣ v i ∧ ¬ 2 * g ∣ v i).Nonempty := by
    rw [← Finset.card_pos]
    simpa [liftFailureCard] using (show 0 < liftFailureCard v g by omega)
  obtain ⟨defaultIndex, hdefaultMem⟩ := hfreshNonempty
  have hdefault : g ∣ v defaultIndex ∧ ¬ 2 * g ∣ v defaultIndex := by
    simpa using hdefaultMem
  obtain ⟨fresh, hfresh, hmode⟩ :
      ∃ i : Fin 13, (g ∣ v i ∧ ¬ 2 * g ∣ v i) ∧
        ((detunedBadBranches (v i) (2 * g) u).Nonempty ∨
          ∀ j : Fin 13, (g ∣ v j ∧ ¬ 2 * g ∣ v j) →
            ¬ (detunedBadBranches (v j) (2 * g) u).Nonempty) := by
    by_cases hactive : ∃ i : Fin 13,
        (g ∣ v i ∧ ¬ 2 * g ∣ v i) ∧
          (detunedBadBranches (v i) (2 * g) u).Nonempty
    · obtain ⟨i, hi, hrow⟩ := hactive
      exact ⟨i, hi, Or.inl hrow⟩
    · refine ⟨defaultIndex, hdefault, Or.inr ?_⟩
      intro j hj hrow
      exact hactive ⟨j, hj, hrow⟩
  have hgood :
      HasThreeDetunedGoodBranch
        (v fresh) (v first) (v second) (2 * g) u := by
    by_contra hfail
    exact hnoObstruction (Or.inr ⟨fresh, hfresh, hfail⟩)
  obtain ⟨branch, hbranch, hfreshGood, hfirstGood, hsecondGood⟩ := hgood
  have hallFreshGood : ∀ j : Fin 13,
      (g ∣ v j ∧ ¬ 2 * g ∣ v j) →
      branch ∉ detunedBadBranches (v j) (2 * g) u := by
    intro j hj
    rcases hmode with hfreshActive | hallEmpty
    · by_cases hjActive : (detunedBadBranches (v j) (2 * g) u).Nonempty
      · by_cases hjFresh : j = fresh
        · subst j
          exact hfreshGood
        · have hnotOpposition :
              ¬ TwoTwoThreePhaseOpposition
                (v fresh) (v j) (2 * g) u := by
            intro hopposition
            exact hnoObstruction (Or.inl
              ⟨fresh, j, hfresh, hj, Ne.symm hjFresh, hopposition⟩)
          have hnotDisjoint :
              ¬ Disjoint
                (detunedBadBranches (v fresh) (2 * g) u)
                (detunedBadBranches (v j) (2 * g) u) := by
            intro hdisjoint
            exact hnotOpposition ⟨hfreshActive, hjActive, hdisjoint⟩
          have hinter :
              (detunedBadBranches (v fresh) (2 * g) u ∩
                detunedBadBranches (v j) (2 * g) u).Nonempty := by
            obtain ⟨c, hleft, hright⟩ :=
              Finset.not_disjoint_iff.mp hnotDisjoint
            exact ⟨c, Finset.mem_inter.mpr ⟨hleft, hright⟩⟩
          have hqFresh := odd_harmonic_lifts_to_q_two
            (v fresh) g hg hfresh.1 hfresh.2
          have hqJ := odd_harmonic_lifts_to_q_two
            (v j) g hg hj.1 hj.2
          have hrowsEq :=
            detunedBadBranches_eq_of_overlap_same_reducedDenominator
              (v fresh) (v j) (2 * g) 2 u (by omega) (by norm_num)
              hqFresh hqJ hinter
          rw [hrowsEq] at hfreshGood
          exact hfreshGood
      · intro hjBranch
        exact hjActive ⟨branch, hjBranch⟩
    · intro hjBranch
      exact hallEmpty j hj ⟨branch, hjBranch⟩
  refine ⟨(u + (branch : ℝ)) / ((2 * g : ℤ) : ℝ), ?_⟩
  intro i n
  by_cases hiFirst : i = first
  · subst i
    exact not_lt.mp fun hlt => hfirstGood (by
      rw [detunedBadBranches, Finset.mem_filter]
      exact ⟨hbranch, n, hlt⟩)
  by_cases hiSecond : i = second
  · subst i
    exact not_lt.mp fun hlt => hsecondGood (by
      rw [detunedBadBranches, Finset.mem_filter]
      exact ⟨hbranch, n, hlt⟩)
  have hdivG : g ∣ v i := hdvd i hiFirst hiSecond
  by_cases hdivDouble : 2 * g ∣ v i
  · have hvalue :
        (v i : ℝ) * ((u + (branch : ℝ)) / ((2 * g : ℤ) : ℝ)) - n =
          ((v i / (2 * g) : ℤ) : ℝ) * u -
            (((n - (v i / (2 * g)) * branch : ℤ)) : ℝ) := by
      have hspeed : (v i : ℝ) =
          (2 * g : ℝ) * ((v i / (2 * g) : ℤ) : ℝ) := by
        have : v i = 2 * g * (v i / (2 * g)) :=
          (Int.mul_ediv_cancel' hdivDouble).symm
        exact_mod_cast this
      rw [hspeed]
      push_cast
      field_simp
      ring
    rw [hvalue]
    exact hcore i hdivDouble (n - (v i / (2 * g)) * branch)
  · have hbranchGood := hallFreshGood i ⟨hdivG, hdivDouble⟩
    exact not_lt.mp fun hlt => hbranchGood (by
      rw [detunedBadBranches, Finset.mem_filter]
      exact ⟨hbranch, n, hlt⟩)

/-- The exact cardinality identity turns the many-wall-event supplier into
the earlier four-or-more doubled-level interface. -/
theorem fourOrMorePairLiftSupply_of_manyLiftFailure
    (hmany : ManyLiftFailurePairTowerSupply) :
    FourOrMorePairLiftSupply := by
  intro v hv g hg hcard _ hnongeneric hfour
  apply hmany v hv g hg hcard
  · have heq := nonMultCard_double_eq_add_liftFailureCard v g
    rw [hcard] at heq
    omega
  · exact hnongeneric

/-- The singleton exposed layer is exactly `(2,4,4)` and is consumed by its
selected-witness supply.  Thus only the four-or-more first lift remains. -/
theorem nonterminatingPairTowerSupply_of_fourOrMore_and_twoFourFour
    (hmany : FourOrMorePairLiftSupply)
    (h244 : TwoFourFourSelectedWitnessSupply) :
    NonterminatingPairTowerSupply := by
  intro v hv g hg hcard hnonterm hnongeneric
  obtain ⟨first, second, hne, hdvd, hfirst, hsecond, hqFirst, hqSecond⟩ :=
    nonGeneric_pair_q_two v g hg hcard hnongeneric
  unfold TwoAdicLiftTerminates at hnonterm
  push Not at hnonterm
  obtain ⟨newIndex, hdivNew, hnotNew⟩ := hnonterm
  have hnewFirst : newIndex ≠ first := by
    intro h
    subst newIndex
    exact hfirst hdivNew
  have hnewSecond : newIndex ≠ second := by
    intro h
    subst newIndex
    exact hsecond hdivNew
  have hfirstLift : ¬ 2 * g ∣ v first := by
    intro h
    exact hfirst (dvd_trans ⟨2, by ring⟩ h)
  have hsecondLift : ¬ 2 * g ∣ v second := by
    intro h
    exact hsecond (dvd_trans ⟨2, by ring⟩ h)
  have hsubset : ({newIndex, first, second} : Finset (Fin 13)) ⊆
      Finset.univ.filter (fun i => ¬ 2 * g ∣ v i) := by
    intro index hindex
    simp only [Finset.mem_insert, Finset.mem_singleton] at hindex
    rw [Finset.mem_filter]
    refine ⟨Finset.mem_univ index, ?_⟩
    rcases hindex with rfl | rfl | rfl
    · exact hnotNew
    · exact hfirstLift
    · exact hsecondLift
  have hsetCard : ({newIndex, first, second} : Finset (Fin 13)).card = 3 := by
    simp [hnewFirst, hnewSecond, hne]
  have hthree : 3 ≤ nonMultCard v (2 * g) := by
    rw [nonMultCard]
    rw [← hsetCard]
    exact Finset.card_le_card hsubset
  by_cases hcardThree : nonMultCard v (2 * g) = 3
  · have hgLift : 2 ≤ 2 * g := by omega
    have hdvdLift := dvd_of_nonMultCard_three_of_selected
      v (2 * g) hcardThree newIndex first second hnewFirst hnewSecond hne
        hnotNew hfirstLift hsecondLift
    have hqNew := odd_harmonic_lifts_to_q_two
      (v newIndex) g hg hdivNew hnotNew
    have hqFirstLift := q_two_lifts_to_q_four
      (v first) g hg hqFirst hfirst
    have hqSecondLift := q_two_lifts_to_q_four
      (v second) g hg hqSecond hsecond
    apply lonely14_of_three_detuned_selectedWitness
      v (2 * g) hgLift newIndex first second hdvdLift
    exact h244 v hv (2 * g) hgLift newIndex first second
      hnewFirst hnewSecond hne hdvdLift hqNew hqFirstLift hqSecondLift
  · exact hmany v hv g hg hcard (by
      intro hterm
      exact hnotNew (hterm newIndex hdivNew)) hnongeneric (by omega)

/-- Exact pair-tower reduction: the singleton wall event is `(2,4,4)`;
only two-or-more simultaneous wall events remain. -/
theorem nonterminatingPairTowerSupply_of_manyLiftFailure_and_twoFourFour
    (hmany : ManyLiftFailurePairTowerSupply)
    (h244 : TwoFourFourSelectedWitnessSupply) :
    NonterminatingPairTowerSupply :=
  nonterminatingPairTowerSupply_of_fourOrMore_and_twoFourFour
    (fourOrMorePairLiftSupply_of_manyLiftFailure hmany) h244

/-- Final pair-tower reduction to the dynamic wall-event selector and the
already named singleton `(2,4,4)` supplier. -/
theorem nonterminatingPairTowerSupply_of_phaseSelector_and_twoFourFour
    (hselector : ManyLiftFailurePhaseSelector)
    (h244 : TwoFourFourSelectedWitnessSupply) :
    NonterminatingPairTowerSupply :=
  nonterminatingPairTowerSupply_of_manyLiftFailure_and_twoFourFour
    (manyLiftFailurePairTowerSupply_of_phaseSelector hselector) h244

/-- Refined supplier-level capstone: the pair obligation begins only at a
first doubled level with at least four nonmultiples. -/
theorem lrc14_from_fourOrMorePairLift_and_selectedWitnessSupplies_and_relationBudget
    (cite : LRCUpTo13)
    (hmany : FourOrMorePairLiftSupply)
    (h22 : TwoTwoSelectedWitnessSupply)
    (h244 : TwoFourFourSelectedWitnessSupply)
    (h333 : UniformThreeSelectedWitnessSupply)
    (hsupply : DenseCoreRelationBudgetSupply) :
    LRC14.LRC14Statement :=
  lrc14_from_selectedWitnessSupplies_and_relationBudget cite
    (nonterminatingPairTowerSupply_of_fourOrMore_and_twoFourFour hmany h244)
    h22 h244 h333 hsupply

/-- Sharpest current pair-tower capstone, stated directly in terms of the
number of coordinates crossing the first doubled-modulus wall. -/
theorem lrc14_from_manyLiftFailure_and_selectedWitnessSupplies_and_relationBudget
    (cite : LRCUpTo13)
    (hmany : ManyLiftFailurePairTowerSupply)
    (h22 : TwoTwoSelectedWitnessSupply)
    (h244 : TwoFourFourSelectedWitnessSupply)
    (h333 : UniformThreeSelectedWitnessSupply)
    (hsupply : DenseCoreRelationBudgetSupply) :
    LRC14.LRC14Statement :=
  lrc14_from_selectedWitnessSupplies_and_relationBudget cite
    (nonterminatingPairTowerSupply_of_manyLiftFailure_and_twoFourFour hmany h244)
    h22 h244 h333 hsupply

/-- Sharpest current checked capstone.  The remaining mathematical inputs are
the dynamic pair wall-event selector, the three phase-selection suppliers, and
the normalized pair/higher-depth deviation budget on the dense core. -/
theorem lrc14_from_pairPhaseSelector_and_selectedWitnessSupplies_and_normalizedBudget
    (cite : LRCUpTo13)
    (hselector : ManyLiftFailurePhaseSelector)
    (h22 : TwoTwoSelectedWitnessSupply)
    (h244 : TwoFourFourSelectedWitnessSupply)
    (h333 : UniformThreeSelectedWitnessSupply)
    (hsupply : DenseCoreNormalizedRelationBudgetSupply) :
    LRC14.LRC14Statement :=
  lrc14_from_twoThree_detuned_and_denseCore_dissociated_B5 cite
    (deepExceptionalDetunedDispatchTwoThree_of_finalResidues cite
      (deepExceptionalDetunedDispatchFinalResidues_of_selectedWitnessSupplies
        (nonterminatingPairTowerSupply_of_phaseSelector_and_twoFourFour
          hselector h244)
        h22 h244 h333))
    (denseCoreDissociatedB5Supply_of_normalizedRelationBudget hsupply)

/-- Fully obstruction-normalized capstone.  Static parallel-class counting has
been discharged exactly; the three phase inputs ask only for a harmonic-good
phase outside the corresponding saturated partition. -/
theorem lrc14_from_pairPhaseSelector_and_parallelAvoidanceSupplies_and_normalizedBudget
    (cite : LRCUpTo13)
    (hselector : ManyLiftFailurePhaseSelector)
    (h22 : TwoTwoParallelAvoidanceSupply)
    (h244 : TwoFourFourParallelAvoidanceSupply)
    (h333 : UniformThreeParallelAvoidanceSupply)
    (hsupply : DenseCoreNormalizedRelationBudgetSupply) :
    LRC14.LRC14Statement :=
  lrc14_from_pairPhaseSelector_and_selectedWitnessSupplies_and_normalizedBudget
    cite hselector
    (twoTwoSelectedWitnessSupply_of_parallelAvoidance h22)
    (twoFourFourSelectedWitnessSupply_of_parallelAvoidance h244)
    (uniformThreeSelectedWitnessSupply_of_parallelAvoidance h333)
    hsupply

/-! ## Axiom audit -/

#print axioms q_two_lifts_to_q_four
#print axioms odd_harmonic_lifts_to_q_two
#print axioms nonMultCard_double_eq_add_liftFailureCard
#print axioms pairTower_firstStep_dichotomy
#print axioms nonGeneric_pair_q_two
#print axioms fourOrMorePairLiftSupply_of_manyLiftFailure
#print axioms manyLiftFailurePairTowerSupply_of_phaseSelector
#print axioms nonterminatingPairTowerSupply_of_fourOrMore_and_twoFourFour
#print axioms nonterminatingPairTowerSupply_of_manyLiftFailure_and_twoFourFour
#print axioms nonterminatingPairTowerSupply_of_phaseSelector_and_twoFourFour
#print axioms lrc14_from_fourOrMorePairLift_and_selectedWitnessSupplies_and_relationBudget
#print axioms lrc14_from_manyLiftFailure_and_selectedWitnessSupplies_and_relationBudget
#print axioms lrc14_from_pairPhaseSelector_and_selectedWitnessSupplies_and_normalizedBudget
#print axioms lrc14_from_pairPhaseSelector_and_parallelAvoidanceSupplies_and_normalizedBudget

end LRC14Grand
end LonelyRunner
