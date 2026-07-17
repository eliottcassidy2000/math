import TournamentH7.LRCSevenOverlapRelations
import TournamentH7.LRCAlignedResonance

/-!
# Zero-color witness gluing

At fixed `(q,p)`, a zero determinant says that two witness-speed vectors have
the same rational slope.  For nonzero speeds this relation is reflexive,
symmetric, and transitive.  Hence every connected zero-color stalk is already
a zero-color clique and carries one common rational witness parameter.

This closes the local gluing primitive left open by `LRCAlignedResonance`.
It does not bound how many multiplier events share the same parameter.  The
faithful vertices are witness-speed vectors; the zero/nonzero determinant is
the binary switch and connected components are the tie classes.  A tournament
orientation of runners preserves only an order between classes and destroys
the common rational parameter inside each class.
-/

namespace LonelyRunner
namespace LRC14Concrete

open scoped Classical

/-- Rational nearest-integer witness slope at one multiplier event. -/
def witnessSlope (v : Fin 13 → ℤ) (q p : ℕ) (index : Fin 13) : ℚ :=
  (failWitness v q p index : ℚ) / (v index : ℚ)

/-- Zero determinant is exactly equality of rational witness slopes. -/
theorem overlapDet_eq_zero_iff_witnessSlope_eq
    (v : Fin 13 → ℤ) (q p : ℕ) (first second : Fin 13)
    (hfirst : v first ≠ 0) (hsecond : v second ≠ 0) :
    overlapDet v q p first second = 0 ↔
      witnessSlope v q p first = witnessSlope v q p second := by
  have hfirstQ : (v first : ℚ) ≠ 0 := by exact_mod_cast hfirst
  have hsecondQ : (v second : ℚ) ≠ 0 := by exact_mod_cast hsecond
  constructor
  · intro hzero
    rw [witnessSlope, witnessSlope,
      div_eq_div_iff hfirstQ hsecondQ]
    have hcross :
        v first * failWitness v q p second =
          v second * failWitness v q p first := by
      exact sub_eq_zero.mp hzero
    exact_mod_cast (show
      failWitness v q p first * v second =
        failWitness v q p second * v first by
      simpa [mul_comm] using hcross.symm)
  · intro hslope
    rw [witnessSlope, witnessSlope,
      div_eq_div_iff hfirstQ hsecondQ] at hslope
    have hcross :
        failWitness v q p first * v second =
          failWitness v q p second * v first := by
      exact_mod_cast hslope
    exact sub_eq_zero.mpr (by simpa [mul_comm] using hcross.symm)

theorem overlapDet_zero_refl
    (v : Fin 13 → ℤ) (q p : ℕ) (index : Fin 13) :
    overlapDet v q p index index = 0 := by
  simp [overlapDet]

theorem overlapDet_zero_symm
    (v : Fin 13 → ℤ) (q p : ℕ) (first second : Fin 13)
    (hzero : overlapDet v q p first second = 0) :
    overlapDet v q p second first = 0 := by
  rw [overlapDet_skew, hzero, neg_zero]

/-- Transitivity is the Plücker identity with the middle speed cancelled. -/
theorem overlapDet_zero_trans
    (v : Fin 13 → ℤ) (q p : ℕ) (first middle last : Fin 13)
    (hmiddle : v middle ≠ 0)
    (hfirstMiddle : overlapDet v q p first middle = 0)
    (hmiddleLast : overlapDet v q p middle last = 0) :
    overlapDet v q p first last = 0 := by
  have hidentity :
      v middle * overlapDet v q p first last =
        v first * overlapDet v q p middle last +
          v last * overlapDet v q p first middle := by
    simp [overlapDet]
    ring
  rw [hfirstMiddle, hmiddleLast, mul_zero, mul_zero, add_zero] at hidentity
  exact (mul_eq_zero.mp hidentity).resolve_left hmiddle

/-- The zero-color relation on a nonzero speed family is a setoid. -/
def overlapZeroSetoid
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0) (q p : ℕ) :
    Setoid (Fin 13) where
  r first second := overlapDet v q p first second = 0
  iseqv := ⟨
    overlapDet_zero_refl v q p,
    fun hzero => overlapDet_zero_symm v q p _ _ hzero,
    fun hfirstMiddle hmiddleLast =>
      overlapDet_zero_trans v q p _ _ _ (hv _) hfirstMiddle hmiddleLast⟩

/-- A path of zero-color edges has zero-color endpoints.  Thus every connected
component of the zero-color graph is a clique. -/
theorem overlapDet_zero_of_reflTransGen
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0) (q p : ℕ)
    (first last : Fin 13)
    (hpath : Relation.ReflTransGen
      (fun left right => overlapDet v q p left right = 0) first last) :
    overlapDet v q p first last = 0 := by
  exact Relation.ReflTransGen.head_induction_on hpath
    (overlapDet_zero_refl v q p last)
    (fun hfirstMiddle _hpath hinduction =>
      overlapDet_zero_trans v q p _ _ last
        (hv _) hfirstMiddle hinduction)

/-- Pairwise zero color glues a nonempty stalk to one rational parameter. -/
theorem exists_common_witnessSlope_of_pairwise_zero
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0) (q p : ℕ)
    (stalk : Finset (Fin 13)) (hstalk : stalk.Nonempty)
    (hzero : Set.Pairwise (stalk : Set (Fin 13))
      fun first second => overlapDet v q p first second = 0) :
    ∃ slope : ℚ, ∀ index ∈ stalk,
      witnessSlope v q p index = slope := by
  obtain ⟨root, hroot⟩ := hstalk
  refine ⟨witnessSlope v q p root, ?_⟩
  intro index hindex
  by_cases hindexRoot : index = root
  · subst index
    rfl
  · apply (overlapDet_eq_zero_iff_witnessSlope_eq
      v q p index root (hv index) (hv root)).mp
    exact hzero hindex hroot hindexRoot

/-- Connected zero color suffices: transitivity first upgrades every root path
to a root edge, after which all vertices share the root slope. -/
theorem exists_common_witnessSlope_of_zero_connected
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0) (q p : ℕ)
    (stalk : Finset (Fin 13)) (root : Fin 13) (_hroot : root ∈ stalk)
    (hconnected : ∀ index ∈ stalk,
      Relation.ReflTransGen
        (fun left right => overlapDet v q p left right = 0) root index) :
    ∃ slope : ℚ, ∀ index ∈ stalk,
      witnessSlope v q p index = slope := by
  refine ⟨witnessSlope v q p root, ?_⟩
  intro index hindex
  have hzero := overlapDet_zero_of_reflTransGen
    v hv q p root index (hconnected index hindex)
  exact (overlapDet_eq_zero_iff_witnessSlope_eq
    v q p root index (hv root) (hv index)).mp hzero |>.symm

/-- A common rational slope has one primitive integer numerator/denominator.
The denominator divides every speed in the stalk and every witness is the
same numerator times the reduced speed.  This is the exact local parameter
needed by the aligned-resonance count. -/
theorem exists_common_primitiveWitnessParameter_of_zero_connected
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0) (q p : ℕ)
    (stalk : Finset (Fin 13)) (root : Fin 13) (hroot : root ∈ stalk)
    (hconnected : ∀ index ∈ stalk,
      Relation.ReflTransGen
        (fun left right => overlapDet v q p left right = 0) root index) :
    ∃ denominator : ℕ, ∃ numerator : ℤ,
      0 < denominator ∧ numerator.natAbs.Coprime denominator ∧
      ∀ index ∈ stalk,
        (denominator : ℤ) ∣ v index ∧
        failWitness v q p index =
          (v index / (denominator : ℤ)) * numerator := by
  obtain ⟨slope, hslope⟩ :=
    exists_common_witnessSlope_of_zero_connected
      v hv q p stalk root hroot hconnected
  refine ⟨slope.den, slope.num, slope.den_pos, slope.reduced, ?_⟩
  intro index hindex
  have hslopeDivInt :
      slope = Rat.divInt (failWitness v q p index) (v index) := by
    rw [Rat.divInt_eq_div]
    exact (hslope index hindex).symm
  obtain ⟨scale, hwitness, hspeed⟩ :=
    Rat.num_den_mk (hv index) hslopeDivInt
  have hdenominator : (slope.den : ℤ) ≠ 0 := by
    exact_mod_cast slope.den_nz
  have hdivides : (slope.den : ℤ) ∣ v index := by
    refine ⟨scale, ?_⟩
    simpa [mul_comm] using hspeed
  refine ⟨hdivides, ?_⟩
  have hquotient : v index / (slope.den : ℤ) = scale :=
    Int.ediv_eq_of_eq_mul_left hdenominator hspeed
  rw [hquotient]
  exact hwitness

/-- One bad runner on a primitive common-slope stalk forces exact resonance
once its reduced speed reaches the sharp integrality window. -/
theorem primitiveWitnessParameter_resonates_of_bad
    (v : Fin 13 → ℤ) (q p : ℕ) (index : Fin 13)
    (denominator : ℕ) (numerator : ℤ)
    (hvelocity : v index ≠ 0)
    (_hdenominator : 0 < denominator)
    (hdivides : (denominator : ℤ) ∣ v index)
    (hwitness : failWitness v q p index =
      (v index / (denominator : ℤ)) * numerator)
    (hbad : 14 *
      |v index * (p : ℤ) - failWitness v q p index * (q : ℤ)| <
        (q : ℤ))
    (hwindow : (q : ℤ) ≤
      14 * |v index / (denominator : ℤ)|) :
    (denominator : ℤ) * (p : ℤ) = numerator * (q : ℤ) := by
  let reducedSpeed : ℤ := v index / (denominator : ℤ)
  have hfactor : v index = (denominator : ℤ) * reducedSpeed := by
    exact (Int.mul_ediv_cancel' hdivides).symm
  have herror :
      v index * (p : ℤ) - failWitness v q p index * (q : ℤ) =
        reducedSpeed *
          ((denominator : ℤ) * p - numerator * q) := by
    rw [hfactor, hwitness]
    dsimp [reducedSpeed]
    ring
  rw [herror, abs_mul] at hbad
  have hbad' :
      14 * |reducedSpeed| *
          |(denominator : ℤ) * p - numerator * q| < (q : ℤ) := by
    simpa [mul_assoc] using hbad
  have hscalePos : (0 : ℤ) < 14 * |reducedSpeed| := by
    have hreducedNe : reducedSpeed ≠ 0 := by
      intro hzero
      apply hvelocity
      rw [hfactor, hzero, mul_zero]
    exact mul_pos (by omega) (abs_pos.mpr hreducedNe)
  have hscaled :
      14 * |reducedSpeed| *
          |(denominator : ℤ) * p - numerator * q| <
        14 * |reducedSpeed| * 1 :=
    lt_of_lt_of_le hbad' (by simpa [reducedSpeed] using hwindow)
  have habs :
      |(denominator : ℤ) * p - numerator * q| < 1 :=
    (mul_lt_mul_iff_right₀ hscalePos).mp hscaled
  exact sub_eq_zero.mp (Int.abs_lt_one_iff.mp habs)

/-- Without a window assumption, badness has one exact alternative: either
the primitive parameter resonates, or the reduced speed lies below `q / 14`.
The threshold is sharp because a nonzero integer cross-product can have
absolute value one. -/
theorem primitiveWitnessParameter_resonates_or_reducedSpeed_small_of_bad
    (v : Fin 13 → ℤ) (q p : ℕ) (index : Fin 13)
    (denominator : ℕ) (numerator : ℤ)
    (hvelocity : v index ≠ 0)
    (hdivides : (denominator : ℤ) ∣ v index)
    (hwitness : failWitness v q p index =
      (v index / (denominator : ℤ)) * numerator)
    (hbad : 14 *
      |v index * (p : ℤ) - failWitness v q p index * (q : ℤ)| <
        (q : ℤ)) :
    (denominator : ℤ) * (p : ℤ) = numerator * (q : ℤ) ∨
      14 * |v index / (denominator : ℤ)| < (q : ℤ) := by
  let reducedSpeed : ℤ := v index / (denominator : ℤ)
  have hfactor : v index = (denominator : ℤ) * reducedSpeed := by
    exact (Int.mul_ediv_cancel' hdivides).symm
  have hreducedNe : reducedSpeed ≠ 0 := by
    intro hzero
    apply hvelocity
    rw [hfactor, hzero, mul_zero]
  have herror :
      v index * (p : ℤ) - failWitness v q p index * (q : ℤ) =
        reducedSpeed *
          ((denominator : ℤ) * p - numerator * q) := by
    rw [hfactor, hwitness]
    dsimp [reducedSpeed]
    ring
  rw [herror, abs_mul] at hbad
  by_cases hresonance :
      (denominator : ℤ) * (p : ℤ) = numerator * (q : ℤ)
  · exact Or.inl hresonance
  · right
    have herrorNe :
        (denominator : ℤ) * (p : ℤ) - numerator * (q : ℤ) ≠ 0 :=
      sub_ne_zero.mpr hresonance
    have hone : (1 : ℤ) ≤
        |(denominator : ℤ) * p - numerator * q| :=
      Int.one_le_abs herrorNe
    have hscaleNonneg : (0 : ℤ) ≤ 14 * |reducedSpeed| := by positivity
    have hlower : 14 * |reducedSpeed| ≤
        14 * |reducedSpeed| *
          |(denominator : ℤ) * p - numerator * q| := by
      nlinarith
    dsimp [reducedSpeed] at hbad ⊢
    exact lt_of_le_of_lt hlower (by simpa [mul_assoc] using hbad)

/-- An exact primitive resonance carries the reduced multiplier modulus into
the primitive denominator. -/
theorem reduced_modulus_dvd_primitive_denominator_of_resonance
    (p q denominator : ℕ) (numerator : ℤ) (hq : 0 < q)
    (hresonance : (denominator : ℤ) * (p : ℤ) =
      numerator * (q : ℤ)) :
    q / Nat.gcd p q ∣ denominator := by
  have hqDvdInt : (q : ℤ) ∣ (p : ℤ) * (denominator : ℤ) := by
    refine ⟨numerator, ?_⟩
    calc
      (p : ℤ) * denominator = denominator * p := by ring
      _ = numerator * q := hresonance
      _ = q * numerator := by ring
  have hqDvdNat : q ∣ p * denominator := by
    exact_mod_cast hqDvdInt
  exact (LRC14Grand.dvd_mul_iff_reduced_modulus_dvd
    p q denominator hq).mp hqDvdNat

/-- A primitive witness parameter in exact resonance is the multiplier
fraction in lowest terms, including its numerator. -/
theorem primitiveWitnessParameter_eq_reducedMultiplier_of_resonance
    (p q denominator : ℕ) (numerator : ℤ) (hq : 0 < q)
    (hprimitive : numerator.natAbs.Coprime denominator)
    (hresonance : (denominator : ℤ) * (p : ℤ) =
      numerator * (q : ℤ)) :
    denominator = q / Nat.gcd p q ∧
      numerator = (p / Nat.gcd p q : ℕ) := by
  let common := Nat.gcd p q
  have hcommonPos : 0 < common := by
    simpa [common] using Nat.gcd_pos_of_pos_right p hq
  have hcommonDvdP : common ∣ p := by
    simpa [common] using Nat.gcd_dvd_left p q
  have hcommonDvdQ : common ∣ q := by
    simpa [common] using Nat.gcd_dvd_right p q
  have hqReducedPos : 0 < q / common :=
    Nat.div_pos (Nat.le_of_dvd hq hcommonDvdQ) hcommonPos
  have hreducedCoprime : (p / common).Coprime (q / common) := by
    simpa [common] using Nat.coprime_div_gcd_div_gcd hcommonPos
  have hcommonInt : (common : ℤ) ≠ 0 := by
    exact_mod_cast (ne_of_gt hcommonPos)
  have hpFactor : p = common * (p / common) := by
    exact (Nat.mul_div_cancel' hcommonDvdP).symm
  have hqFactor : q = common * (q / common) := by
    exact (Nat.mul_div_cancel' hcommonDvdQ).symm
  have hpCastFactor : (p : ℤ) = (common : ℤ) * (p / common : ℕ) := by
    exact_mod_cast hpFactor
  have hqCastFactor : (q : ℤ) = (common : ℤ) * (q / common : ℕ) := by
    exact_mod_cast hqFactor
  have hreducedResonance :
      (denominator : ℤ) * (p / common : ℕ) =
        numerator * (q / common : ℕ) := by
    have hwithCommon :
        ((denominator : ℤ) * (p / common : ℕ)) * (common : ℤ) =
          (numerator * (q / common : ℕ)) * (common : ℤ) := by
      calc
        ((denominator : ℤ) * (p / common : ℕ)) * common =
            (denominator : ℤ) * p := by rw [hpCastFactor]; ring
        _ = numerator * q := hresonance
        _ = (numerator * (q / common : ℕ)) * common := by
          rw [hqCastFactor]
          ring
    exact mul_right_cancel₀ hcommonInt hwithCommon
  have hqReducedDvdDenominator : q / common ∣ denominator := by
    simpa [hreducedCoprime.gcd_eq_one] using
      reduced_modulus_dvd_primitive_denominator_of_resonance
        (p / common) (q / common) denominator numerator hqReducedPos
          hreducedResonance
  have hdenominatorDvdInt :
      (denominator : ℤ) ∣ numerator * (q / common : ℕ) := by
    exact ⟨(p / common : ℕ), hreducedResonance.symm⟩
  have hdenominatorDvdNat :
      denominator ∣ numerator.natAbs * (q / common) := by
    have hnatAbs := (Int.natAbs_dvd_natAbs).mpr hdenominatorDvdInt
    simpa [Int.natAbs_mul] using hnatAbs
  have hdenominatorDvdQReduced : denominator ∣ q / common :=
    hprimitive.symm.dvd_of_dvd_mul_left hdenominatorDvdNat
  have hdenominatorEq : denominator = q / common :=
    Nat.dvd_antisymm hdenominatorDvdQReduced hqReducedDvdDenominator
  refine ⟨by simpa [common] using hdenominatorEq, ?_⟩
  have hqReducedInt : ((q / common : ℕ) : ℤ) ≠ 0 := by
    exact_mod_cast (ne_of_gt hqReducedPos)
  have hnumerator : numerator = (p / common : ℕ) := by
    apply Eq.symm
    apply mul_right_cancel₀ hqReducedInt
    calc
      ((p / common : ℕ) : ℤ) * (q / common : ℕ) =
          (q / common : ℕ) * (p / common : ℕ) := by ring
      _ = numerator * (q / common : ℕ) := by
        simpa [hdenominatorEq] using hreducedResonance
  simpa [common] using hnumerator

/-- At a coprime multiplier, both primitive parameter coordinates agree
exactly with the multiplier numerator and denominator. -/
theorem primitiveWitnessParameter_eq_of_coprime_resonance
    (p q denominator : ℕ) (numerator : ℤ) (hq : 0 < q)
    (hpq : p.Coprime q)
    (hprimitive : numerator.natAbs.Coprime denominator)
    (hresonance : (denominator : ℤ) * (p : ℤ) =
      numerator * (q : ℤ)) :
    denominator = q ∧ numerator = (p : ℤ) := by
  simpa [hpq.gcd_eq_one] using
    primitiveWitnessParameter_eq_reducedMultiplier_of_resonance
      p q denominator numerator hq hprimitive hresonance

/-- At a coprime multiplier, an exact primitive resonance forces the whole
modulus into the primitive denominator. -/
theorem modulus_dvd_primitive_denominator_of_coprime_resonance
    (p q denominator : ℕ) (numerator : ℤ) (hq : 0 < q)
    (hpq : p.Coprime q)
    (hresonance : (denominator : ℤ) * (p : ℤ) =
      numerator * (q : ℤ)) :
    q ∣ denominator := by
  simpa [hpq.gcd_eq_one] using
    reduced_modulus_dvd_primitive_denominator_of_resonance
      p q denominator numerator hq hresonance

/-- A connected zero-color bad stalk in the sharp primitive window carries
the reduced multiplier modulus in every stalk speed. -/
theorem reduced_modulus_dvd_stalk_speeds_of_zero_connected_bad
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0)
    (q p : ℕ) (hq : 0 < q)
    (stalk : Finset (Fin 13)) (root pivot : Fin 13)
    (hroot : root ∈ stalk) (hpivot : pivot ∈ stalk)
    (hconnected : ∀ index ∈ stalk,
      Relation.ReflTransGen
        (fun left right => overlapDet v q p left right = 0) root index)
    (hpivotBad : ¬ inBand v q p pivot)
    (hwindow : ∀ denominator : ℕ,
      (denominator : ℤ) ∣ v pivot →
      (q : ℤ) ≤ 14 * |v pivot / (denominator : ℤ)|) :
    ∀ index ∈ stalk,
      ((q / Nat.gcd p q : ℕ) : ℤ) ∣ v index := by
  obtain ⟨denominator, numerator, hdenominator, hprimitive,
      hparameter⟩ :=
    exists_common_primitiveWitnessParameter_of_zero_connected
      v hv q p stalk root hroot hconnected
  have hpivotParameter := hparameter pivot hpivot
  have hresonance := primitiveWitnessParameter_resonates_of_bad
    v q p pivot denominator numerator (hv pivot) hdenominator
      hpivotParameter.1 hpivotParameter.2
      (bad_at_witness v q p pivot hq hpivotBad)
      (hwindow denominator hpivotParameter.1)
  have hparameterEq :=
    primitiveWitnessParameter_eq_reducedMultiplier_of_resonance
      p q denominator numerator hq hprimitive hresonance
  have hreducedDenominatorInt :
      ((q / Nat.gcd p q : ℕ) : ℤ) = (denominator : ℤ) := by
    exact_mod_cast hparameterEq.1.symm
  intro index hindex
  rw [hreducedDenominatorInt]
  exact (hparameter index hindex).1

/-- A connected zero-color bad stalk at a coprime multiplier cannot remain a
merely local alignment once the primitive window applies: the modulus divides
every speed in the stalk. -/
theorem modulus_dvd_stalk_speeds_of_zero_connected_bad_coprime
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0)
    (q p : ℕ) (hq : 0 < q) (hpq : p.Coprime q)
    (stalk : Finset (Fin 13)) (root pivot : Fin 13)
    (hroot : root ∈ stalk) (hpivot : pivot ∈ stalk)
    (hconnected : ∀ index ∈ stalk,
      Relation.ReflTransGen
        (fun left right => overlapDet v q p left right = 0) root index)
    (hpivotBad : ¬ inBand v q p pivot)
    (hwindow : ∀ denominator : ℕ,
      (denominator : ℤ) ∣ v pivot →
      (q : ℤ) ≤ 14 * |v pivot / (denominator : ℤ)|) :
    ∀ index ∈ stalk, (q : ℤ) ∣ v index := by
  simpa [hpq.gcd_eq_one] using
    reduced_modulus_dvd_stalk_speeds_of_zero_connected_bad
      v hv q p hq stalk root pivot hroot hpivot hconnected hpivotBad hwindow

#print axioms overlapDet_eq_zero_iff_witnessSlope_eq
#print axioms overlapDet_zero_trans
#print axioms overlapZeroSetoid
#print axioms overlapDet_zero_of_reflTransGen
#print axioms exists_common_witnessSlope_of_pairwise_zero
#print axioms exists_common_witnessSlope_of_zero_connected
#print axioms exists_common_primitiveWitnessParameter_of_zero_connected
#print axioms primitiveWitnessParameter_resonates_of_bad
#print axioms primitiveWitnessParameter_resonates_or_reducedSpeed_small_of_bad
#print axioms reduced_modulus_dvd_primitive_denominator_of_resonance
#print axioms primitiveWitnessParameter_eq_reducedMultiplier_of_resonance
#print axioms primitiveWitnessParameter_eq_of_coprime_resonance
#print axioms modulus_dvd_primitive_denominator_of_coprime_resonance
#print axioms reduced_modulus_dvd_stalk_speeds_of_zero_connected_bad
#print axioms modulus_dvd_stalk_speeds_of_zero_connected_bad_coprime

end LRC14Concrete
end LonelyRunner
