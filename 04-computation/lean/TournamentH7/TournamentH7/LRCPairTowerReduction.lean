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
import TournamentH7.LRCLadderFattening

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


/-- A half-harmonic row at `g` becomes a quarter-harmonic row at `2g`. -/
theorem reducedDenominator_double_eq_four_of_eq_two
    (δ g : ℤ) (hg : 0 < g)
    (hq : g / (Int.gcd δ g : ℤ) = 2) :
    2 * g / (Int.gcd δ (2 * g) : ℤ) = 4 := by
  set d : ℤ := (Int.gcd δ g : ℤ) with hd
  set p : ℤ := δ / d with hp
  set q : ℤ := g / d with hqDef
  have hdpos : 0 < d := by
    rw [hd]
    have : Int.gcd δ g ≠ 0 := by
      intro hzero
      rw [Int.gcd_eq_zero_iff] at hzero
      omega
    exact_mod_cast Nat.pos_of_ne_zero this
  have hd0 : d ≠ 0 := ne_of_gt hdpos
  have hddδ : d ∣ δ := by
    rw [hd]
    exact Int.gcd_dvd_left δ g
  have hddg : d ∣ g := by
    rw [hd]
    exact Int.gcd_dvd_right δ g
  have hδdp : δ = d * p := (Int.mul_ediv_cancel' hddδ).symm
  have hgdq : g = d * q := (Int.mul_ediv_cancel' hddg).symm
  have hq2 : q = 2 := by simpa [q, d, hd] using hq
  have hcop : IsCoprime p q := by
    set A := Int.gcdA δ g
    set B := Int.gcdB δ g
    have hbez : d = δ * A + g * B := by
      rw [hd]
      exact Int.gcd_eq_gcd_ab δ g
    refine ⟨A, B, ?_⟩
    have hfac : d * (A * p + B * q) = d * 1 := by
      have h1 : d * (A * p + B * q) = (d * p) * A + (d * q) * B := by ring
      rw [h1, ← hδdp, ← hgdq, ← hbez, mul_one]
    exact mul_left_cancel₀ hd0 hfac
  have hcop2 : IsCoprime p (2 : ℤ) := by simpa [hq2] using hcop
  have hcop4 : IsCoprime p (4 : ℤ) := by
    simpa using hcop2.pow_right (n := 2)
  have hgcdp4 : Int.gcd p 4 = 1 :=
    Int.isCoprime_iff_gcd_eq_one.mp hcop4
  have hgcdp2 : Int.gcd p 2 = 1 :=
    Int.isCoprime_iff_gcd_eq_one.mp hcop2
  have hgcdNat : Int.gcd δ (2 * g) = Int.gcd δ g := by
    rw [hδdp, hgdq]
    have hdouble : 2 * (d * q) = d * (2 * q) := by ring
    rw [hdouble, hq2]
    simp only [Int.gcd_mul_left]
    norm_num [hgcdp4, hgcdp2]
  have hgcdZ : (Int.gcd δ (2 * g) : ℤ) = d := by
    rw [hgcdNat, hd]
  rw [hgcdZ]
  apply Int.ediv_eq_of_eq_mul_right hd0
  rw [hgdq, hq2]
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

/-- At the zero phase every detuned bad row contains branch zero. -/
theorem zero_mem_detunedBadBranches
    (δ modulus : ℤ) (hmodulus : 1 ≤ modulus) :
    (0 : ℤ) ∈ detunedBadBranches δ modulus 0 := by
  rw [detunedBadBranches, Finset.mem_filter]
  refine ⟨Finset.mem_Ico.mpr ⟨by omega, by omega⟩, 0, ?_⟩
  norm_num

/-- The zero phase cannot carry either complete binary prefix obstruction.
All q-two rows meet at branch zero, and in q244 the same overlap pays the
exact saturated degree budget. -/
theorem not_pairTowerParallelObstruction_zero
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g)
    (i₄a i₄b : Fin 13)
    (hq4aBase : g / (Int.gcd (v i₄a) g : ℤ) = 2)
    (hq4bBase : g / (Int.gcd (v i₄b) g : ℤ) = 2) :
    ¬ PairTowerParallelObstruction v g i₄a i₄b 0 := by
  have hdoublePos : (0 : ℤ) < 2 * g := by omega
  have hdoubleOne : (1 : ℤ) ≤ 2 * g := by omega
  have hq4a := reducedDenominator_double_eq_four_of_eq_two
    (v i₄a) g (by omega) hq4aBase
  have hq4b := reducedDenominator_double_eq_four_of_eq_two
    (v i₄b) g (by omega) hq4bBase
  rintro (h22 | h244)
  · obtain ⟨i, j, hi, hj, -, hOpposition⟩ := h22
    have hzeroI := zero_mem_detunedBadBranches (v i) (2 * g) hdoubleOne
    have hzeroJ := zero_mem_detunedBadBranches (v j) (2 * g) hdoubleOne
    exact Finset.disjoint_left.mp hOpposition.2.2 hzeroI hzeroJ
  · obtain ⟨i, hi, hfail⟩ := h244
    have hq2 := odd_harmonic_lifts_to_q_two (v i) g hg hi.1 hi.2
    have hbad2 := two_mul_badCount_eq (v i) (2 * g) hdoublePos hq2
    have hbad4a := four_mul_badCount_eq_of_reducedDenominator_four
      (v i₄a) (2 * g) hdoublePos hq4a
    have hbad4b := four_mul_badCount_eq_of_reducedDenominator_four
      (v i₄b) (2 * g) hdoublePos hq4b
    have hbudget :
        DetunedD3.badCount (v i) (2 * g) +
            DetunedD3.badCount (v i₄a) (2 * g) +
            DetunedD3.badCount (v i₄b) (2 * g) ≤ (2 * g).toNat := by
      omega
    have hinter :
        (detunedBadBranches (v i) (2 * g) 0 ∩
          detunedBadBranches (v i₄a) (2 * g) 0).Nonempty := by
      refine ⟨0, Finset.mem_inter.mpr ⟨?_, ?_⟩⟩
      · exact zero_mem_detunedBadBranches (v i) (2 * g) hdoubleOne
      · exact zero_mem_detunedBadBranches (v i₄a) (2 * g) hdoubleOne
    exact hfail (hasThreeDetunedGoodBranch_of_pairOverlap
      (v i) (v i₄a) (v i₄b) (2 * g) 0 hdoubleOne hbudget
      (Or.inl hinter))

/-- A genuine sufficient condition for the many-failure selector: if the
doubled harmonic core is empty, the zero phase is core-good and simultaneously
escapes every q22 and q244 parallel-class obstruction. -/
theorem pairTowerPhaseSelector_of_emptyCore
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g)
    (i₄a i₄b : Fin 13)
    (hq4aBase : g / (Int.gcd (v i₄a) g : ℤ) = 2)
    (hq4bBase : g / (Int.gcd (v i₄b) g : ℤ) = 2)
    (hemptyCore : ∀ j, ¬ 2 * g ∣ v j) :
    ∃ u : ℝ,
      PairTowerCoreGoodAt v g u ∧
      ¬ PairTowerParallelObstruction v g i₄a i₄b u := by
  refine ⟨0, ?_, not_pairTowerParallelObstruction_zero
    v g hg i₄a i₄b hq4aBase hq4bBase⟩
  intro j hj
  exact (hemptyCore j hj).elim

/-- If the two old q-two rows are the only base nonmultiples and all eleven
remaining coordinates fail the first lift, then the doubled harmonic core is
empty.  The cardinal statement is the intrinsic maximal-failure form of the
previous sufficient condition. -/
theorem emptyCore_of_liftFailureCard_eleven
    (v : Fin 13 → ℤ) (g : ℤ) (i₄a i₄b : Fin 13)
    (h4ab : i₄a ≠ i₄b)
    (hδ4a : ¬ g ∣ v i₄a) (hδ4b : ¬ g ∣ v i₄b)
    (hlift : liftFailureCard v g = 11) :
    ∀ j, ¬ 2 * g ∣ v j := by
  let fresh : Finset (Fin 13) :=
    Finset.univ.filter fun i => g ∣ v i ∧ ¬ 2 * g ∣ v i
  let rest : Finset (Fin 13) :=
    Finset.univ \ ({i₄a, i₄b} : Finset (Fin 13))
  have hfreshCard : fresh.card = 11 := by
    simpa [fresh, liftFailureCard] using hlift
  have hrestCard : rest.card = 11 := by
    dsimp [rest]
    rw [Finset.card_sdiff, Finset.inter_univ, Finset.card_pair h4ab,
      Finset.card_univ, Fintype.card_fin]
  have hfreshSubset : fresh ⊆ rest := by
    intro i hi
    have hi' : g ∣ v i ∧ ¬ 2 * g ∣ v i := by
      simpa [fresh] using hi
    simp only [rest, Finset.mem_sdiff, Finset.mem_insert, Finset.mem_singleton]
    refine ⟨Finset.mem_univ i, ?_⟩
    push Not
    constructor
    · intro hieq
      subst i
      exact hδ4a hi'.1
    · intro hieq
      subst i
      exact hδ4b hi'.1
  have hfreshEq : fresh = rest := by
    apply Finset.eq_of_subset_of_card_le hfreshSubset
    omega
  intro j hjDouble
  have hjBase : g ∣ v j := by
    obtain ⟨factor, hfactor⟩ := hjDouble
    refine ⟨2 * factor, ?_⟩
    rw [hfactor]
    ring
  have hj4a : j ≠ i₄a := by
    intro hj
    subst j
    exact hδ4a hjBase
  have hj4b : j ≠ i₄b := by
    intro hj
    subst j
    exact hδ4b hjBase
  have hjRest : j ∈ rest := by
    simp only [rest, Finset.mem_sdiff, Finset.mem_insert, Finset.mem_singleton]
    exact ⟨Finset.mem_univ j, by push Not; exact ⟨hj4a, hj4b⟩⟩
  have hjFresh : j ∈ fresh := by rw [hfreshEq]; exact hjRest
  have hjFresh' : g ∣ v j ∧ ¬ 2 * g ∣ v j := by
    simpa [fresh] using hjFresh
  exact hjFresh'.2 hjDouble

/-- Direct maximal-failure selector socket.  In the pair-tower premises,
`liftFailureCard = 11` forces the empty-core situation and is therefore
selected at phase zero. -/
theorem pairTowerPhaseSelector_of_liftFailureCard_eleven
    (v : Fin 13 → ℤ) (g : ℤ) (hg : 2 ≤ g)
    (i₄a i₄b : Fin 13) (h4ab : i₄a ≠ i₄b)
    (hδ4a : ¬ g ∣ v i₄a) (hδ4b : ¬ g ∣ v i₄b)
    (hq4aBase : g / (Int.gcd (v i₄a) g : ℤ) = 2)
    (hq4bBase : g / (Int.gcd (v i₄b) g : ℤ) = 2)
    (hlift : liftFailureCard v g = 11) :
    ∃ u : ℝ,
      PairTowerCoreGoodAt v g u ∧
      ¬ PairTowerParallelObstruction v g i₄a i₄b u :=
  pairTowerPhaseSelector_of_emptyCore v g hg i₄a i₄b hq4aBase hq4bBase
    (emptyCore_of_liftFailureCard_eleven
      v g i₄a i₄b h4ab hδ4a hδ4b hlift)

/-- Indices of the doubled harmonic core. -/
def pairTowerCoreIndices (v : Fin 13 → ℤ) (g : ℤ) : Finset (Fin 13) :=
  Finset.univ.filter fun j => 2 * g ∣ v j

/-- Exact citation fattening radius for a core of the displayed cardinality. -/
noncomputable def pairTowerCitationRadius (coreCard : ℕ) (B : ℝ) : ℝ :=
  (1 / ((coreCard + 1 : ℕ) : ℝ) - 1 / 14) / B

/-- The citation gives a whole harmonic-good interval for the doubled core,
with the sharp margin `1/(coreCard+1)-1/14` rather than the generic ten-row
margin. -/
theorem exists_pairTowerCoreGoodAt_sharpInterval_of_cite
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ j, v j ≠ 0)
    (g : ℤ) (coreCard : ℕ)
    (hcoreCard : (pairTowerCoreIndices v g).card = coreCard)
    (hcard : coreCard ≤ 12)
    (B : ℝ) (hB0 : 0 < B)
    (hB : ∀ j, 2 * g ∣ v j →
      |((v j / (2 * g) : ℤ) : ℝ)| ≤ B) :
    ∃ center : ℝ,
      ∀ u ∈ Set.Icc
          (center - pairTowerCitationRadius coreCard B)
          (center + pairTowerCitationRadius coreCard B),
        PairTowerCoreGoodAt v g u := by
  let core := pairTowerCoreIndices v g
  set embedding : Fin coreCard ↪o Fin 13 :=
    core.orderEmbOfFin hcoreCard with hembedding
  have hembeddingMem : ∀ index, embedding index ∈ core :=
    fun index => core.orderEmbOfFin_mem hcoreCard index
  have hembeddingDivides : ∀ index, 2 * g ∣ v (embedding index) := by
    intro index
    have hmem := hembeddingMem index
    simpa [core, pairTowerCoreIndices] using hmem
  set quotient : Fin coreCard → ℤ :=
    fun index => v (embedding index) / (2 * g) with hquotient
  have hquotientNonzero : ∀ index, quotient index ≠ 0 := by
    intro index hzero
    have hfactor :
        v (embedding index) = (2 * g) * quotient index :=
      (Int.mul_ediv_cancel' (hembeddingDivides index)).symm
    apply hv (embedding index)
    rw [hfactor, hzero, mul_zero]
  have hquotientBound : ∀ index, |(quotient index : ℝ)| ≤ B := by
    intro index
    exact hB (embedding index) (hembeddingDivides index)
  obtain ⟨center, hcenter⟩ :=
    cite coreCard hcard quotient hquotientNonzero
  refine ⟨center, ?_⟩
  intro u hu
  have habs : |u - center| ≤ pairTowerCitationRadius coreCard B := by
    rw [abs_le]
    constructor <;> linarith [hu.1, hu.2]
  have hdrift :
      B * |u - center| ≤
        1 / ((coreCard + 1 : ℕ) : ℝ) - 1 / 14 := by
    calc
      B * |u - center| ≤
          B * pairTowerCitationRadius coreCard B :=
        mul_le_mul_of_nonneg_left habs (le_of_lt hB0)
      _ = 1 / ((coreCard + 1 : ℕ) : ℝ) - 1 / 14 := by
        simp only [pairTowerCitationRadius]
        field_simp [ne_of_gt hB0]
  have htransport : Lonely 14 quotient u := by
    intro index integer
    apply lonely_band_transport (coreCard + 1) quotient center u B (1 / 14)
      hcenter hquotientBound
    linarith
  intro j hjDouble integer
  have hjCore : j ∈ core := by
    simpa [core, pairTowerCoreIndices] using hjDouble
  obtain ⟨source, hsource⟩ : ∃ source, embedding source = j := by
    have hrange : j ∈ Set.range embedding := by
      rw [hembedding, Finset.range_orderEmbOfFin]
      exact hjCore
    exact hrange
  have hquotientEq : quotient source = v j / (2 * g) := by
    show v (embedding source) / (2 * g) = v j / (2 * g)
    rw [hsource]
  simpa [hquotientEq] using htransport source integer

/-- An interval at least as wide as one open integer-danger chamber contains
a point outside that chamber. -/
theorem exists_pairTower_integer_clear_in_interval
    (lower upper dangerRadius : ℝ)
    (hdanger0 : 0 ≤ dangerRadius) (hdangerHalf : 2 * dangerRadius ≤ 1)
    (hlowerUpper : lower ≤ upper)
    (hwidth : 2 * dangerRadius ≤ upper - lower) :
    ∃ point ∈ Set.Icc lower upper,
      ∀ integer : ℤ, dangerRadius ≤ |point - integer| := by
  by_cases hlowerClear :
      ∀ integer : ℤ, dangerRadius ≤ |lower - integer|
  · exact ⟨lower, ⟨le_rfl, hlowerUpper⟩, hlowerClear⟩
  · push Not at hlowerClear
    obtain ⟨nearest, hnearest⟩ := hlowerClear
    refine ⟨(nearest : ℝ) + dangerRadius, ?_, ?_⟩
    · rw [abs_lt] at hnearest
      constructor <;> linarith
    · intro integer
      rcases lt_trichotomy integer nearest with hlt | heq | hgt
      · have hstep : (integer : ℝ) + 1 ≤ nearest := by exact_mod_cast hlt
        rw [abs_of_nonneg]
        · linarith
        · linarith
      · subst integer
        simp [abs_of_nonneg hdanger0]
      · have hstep : (nearest : ℝ) + 1 ≤ integer := by exact_mod_cast hgt
        rw [abs_of_nonpos]
        · linarith
        · linarith

/-- Dynamic pair-tower escape.  A harmonic-good interval crosses any *single
common affine obstruction wall* whose scalar image spans one whole danger
chamber.  This is the exact nonchaining hypothesis missing from independent
q22/q244 moment bounds. -/
theorem pairTowerPhaseSelector_of_coreInterval_and_commonWall
    (v : Fin 13 → ℤ) (g : ℤ) (i₄a i₄b : Fin 13)
    (center intervalRadius dangerRadius moment : ℝ) (frequency : ℤ)
    (hintervalRadius : 0 ≤ intervalRadius)
    (hdanger0 : 0 ≤ dangerRadius)
    (hdangerHalf : 2 * dangerRadius ≤ 1)
    (hfrequency : frequency ≠ 0)
    (hinterval : ∀ u ∈ Set.Icc
        (center - intervalRadius) (center + intervalRadius),
      PairTowerCoreGoodAt v g u)
    (hfailureWall : ∀ u : ℝ,
      PairTowerParallelObstruction v g i₄a i₄b u →
      ∃ integer : ℤ,
        |(frequency : ℝ) * u + moment - integer| < dangerRadius)
    (hwidth : 2 * dangerRadius ≤
      2 * intervalRadius * |(frequency : ℝ)|) :
    ∃ u : ℝ,
      PairTowerCoreGoodAt v g u ∧
      ¬ PairTowerParallelObstruction v g i₄a i₄b u := by
  have hfrequencyReal : (frequency : ℝ) ≠ 0 := by exact_mod_cast hfrequency
  rcases lt_or_gt_of_ne hfrequencyReal with hnegative | hpositive
  · let lower : ℝ :=
      (frequency : ℝ) * (center + intervalRadius) + moment
    let upper : ℝ :=
      (frequency : ℝ) * (center - intervalRadius) + moment
    have hlowerUpper : lower ≤ upper := by
      dsimp [lower, upper]
      nlinarith
    have hspan : upper - lower =
        2 * intervalRadius * |(frequency : ℝ)| := by
      rw [abs_of_neg hnegative]
      dsimp [lower, upper]
      ring
    obtain ⟨point, hpoint, hclear⟩ :=
      exists_pairTower_integer_clear_in_interval
        lower upper dangerRadius hdanger0 hdangerHalf hlowerUpper
        (by rwa [hspan])
    let u : ℝ := (point - moment) / (frequency : ℝ)
    have hphase : (frequency : ℝ) * u + moment = point := by
      dsimp [u]
      field_simp
      ring
    have hu : u ∈ Set.Icc
        (center - intervalRadius) (center + intervalRadius) := by
      dsimp [lower, upper] at hpoint
      constructor <;> nlinarith [hpoint.1, hpoint.2]
    refine ⟨u, hinterval u hu, ?_⟩
    intro hfail
    obtain ⟨integer, hinteger⟩ := hfailureWall u hfail
    rw [hphase] at hinteger
    exact (not_lt_of_ge (hclear integer)) hinteger

  · let lower : ℝ :=
      (frequency : ℝ) * (center - intervalRadius) + moment
    let upper : ℝ :=
      (frequency : ℝ) * (center + intervalRadius) + moment
    have hlowerUpper : lower ≤ upper := by
      dsimp [lower, upper]
      nlinarith
    have hspan : upper - lower =
        2 * intervalRadius * |(frequency : ℝ)| := by
      rw [abs_of_pos hpositive]
      dsimp [lower, upper]
      ring
    obtain ⟨point, hpoint, hclear⟩ :=
      exists_pairTower_integer_clear_in_interval
        lower upper dangerRadius hdanger0 hdangerHalf hlowerUpper
        (by rwa [hspan])
    let u : ℝ := (point - moment) / (frequency : ℝ)
    have hphase : (frequency : ℝ) * u + moment = point := by
      dsimp [u]
      field_simp
      ring
    have hu : u ∈ Set.Icc
        (center - intervalRadius) (center + intervalRadius) := by
      dsimp [lower, upper] at hpoint
      constructor <;> nlinarith [hpoint.1, hpoint.2]
    refine ⟨u, hinterval u hu, ?_⟩
    intro hfail
    obtain ⟨integer, hinteger⟩ := hfailureWall u hfail
    rw [hphase] at hinteger
    exact (not_lt_of_ge (hclear integer)) hinteger

/-- Citation-facing common-wall socket.  The exact core-cardinality margin
funds the wall crossing; the displayed width inequality is the only numerical
condition. -/
theorem pairTowerPhaseSelector_of_cite_and_commonWall
    (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ j, v j ≠ 0)
    (g : ℤ) (i₄a i₄b : Fin 13)
    (coreCard : ℕ)
    (hcoreCard : (pairTowerCoreIndices v g).card = coreCard)
    (hcard : coreCard ≤ 12)
    (B : ℝ) (hB0 : 0 < B)
    (hB : ∀ j, 2 * g ∣ v j →
      |((v j / (2 * g) : ℤ) : ℝ)| ≤ B)
    (dangerRadius moment : ℝ) (frequency : ℤ)
    (hdanger0 : 0 ≤ dangerRadius)
    (hdangerHalf : 2 * dangerRadius ≤ 1)
    (hfrequency : frequency ≠ 0)
    (hfailureWall : ∀ u : ℝ,
      PairTowerParallelObstruction v g i₄a i₄b u →
      ∃ integer : ℤ,
        |(frequency : ℝ) * u + moment - integer| < dangerRadius)
    (hwidth : 2 * dangerRadius ≤
      2 * pairTowerCitationRadius coreCard B *
        |(frequency : ℝ)|) :
    ∃ u : ℝ,
      PairTowerCoreGoodAt v g u ∧
      ¬ PairTowerParallelObstruction v g i₄a i₄b u := by
  obtain ⟨center, hinterval⟩ :=
    exists_pairTowerCoreGoodAt_sharpInterval_of_cite
      cite v hv g coreCard hcoreCard hcard B hB0 hB
  have hcardReal : ((coreCard : ℝ) + 1) ≤ 13 := by
    have : (coreCard : ℝ) ≤ 12 := by exact_mod_cast hcard
    linarith
  have hcardPositive : (0 : ℝ) < (coreCard : ℝ) + 1 := by positivity
  have hgap : (1 : ℝ) / 14 ≤ 1 / ((coreCard + 1 : ℕ) : ℝ) := by
    have hcast : ((coreCard + 1 : ℕ) : ℝ) = (coreCard : ℝ) + 1 := by
      push_cast
      ring
    rw [hcast]
    calc
      (1 : ℝ) / 14 ≤ 1 / 13 := by norm_num
      _ ≤ 1 / ((coreCard : ℝ) + 1) :=
        one_div_le_one_div_of_le hcardPositive hcardReal
  have hintervalRadius : 0 ≤ pairTowerCitationRadius coreCard B := by
    rw [pairTowerCitationRadius]
    exact div_nonneg (sub_nonneg.mpr hgap) (le_of_lt hB0)
  exact pairTowerPhaseSelector_of_coreInterval_and_commonWall
    v g i₄a i₄b center (pairTowerCitationRadius coreCard B)
      dangerRadius moment frequency hintervalRadius hdanger0 hdangerHalf
      hfrequency hinterval hfailureWall hwidth

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
#print axioms not_pairTowerParallelObstruction_zero
#print axioms pairTowerPhaseSelector_of_emptyCore
#print axioms emptyCore_of_liftFailureCard_eleven
#print axioms pairTowerPhaseSelector_of_liftFailureCard_eleven
#print axioms exists_pairTowerCoreGoodAt_sharpInterval_of_cite
#print axioms pairTowerPhaseSelector_of_coreInterval_and_commonWall
#print axioms pairTowerPhaseSelector_of_cite_and_commonWall

end LRC14Grand
end LonelyRunner
