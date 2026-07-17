import TournamentH7.LRCZeroColorGluing

/-!
# The exact zero-color stalk fork

A connected zero-color bad stalk has one primitive witness parameter.  Applying
the sharp integrality fork to every vertex, not only to one chosen pivot, gives
the useful global alternative: either that parameter is exactly the reduced
multiplier `p/q`, so the reduced modulus divides every stalk speed, or every
speed after division by the common parameter denominator lies below `q/14`.

The sharp window is not automatic from the usual short-modulus hypothesis.
For `q = 99`, `p = 1`, denominator `98`, numerator `1`, and speeds `98s` with
`1 ≤ s ≤ 7`, the witnesses `s` form an all-bad zero-color stalk while
`14s < 99`; the parameter is nonresonant even though `q ≤ 7 * 98`.

The faithful carrier is the connected zero-color component together with its
primitive rational label.  Orienting runners as a tournament forgets the label
and the shared resonance equation; there is no meaningful tournament switch in
this exact dichotomy.
-/

namespace LonelyRunner
namespace LRCZeroColorStalkFork

open LRC14Concrete
open scoped Classical

/-- Seven distinct nonzero speeds cannot all have distinct positive reduced
magnitudes below `7`.  This is the exact pigeonhole that turns the sharp
stalk fork into an unconditional resonance statement for `q ≤ 98`. -/
theorem not_all_reduced_small_of_seven_distinct
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0)
    (hdistinct : ∀ first second, first ≠ second →
      |v first| ≠ |v second|)
    (q denominator : ℕ) (hq : q ≤ 98)
    (stalk : Finset (Fin 13)) (hcard : stalk.card = 7)
    (hdiv : ∀ index ∈ stalk, (denominator : ℤ) ∣ v index)
    (hsmall : ∀ index ∈ stalk,
      14 * |v index / (denominator : ℤ)| < (q : ℤ)) : False := by
  let reducedMagnitude : Fin 13 → ℕ := fun index =>
    (v index / (denominator : ℤ)).natAbs
  have hspeedMagnitude : ∀ index ∈ stalk,
      (v index).natAbs = denominator * reducedMagnitude index := by
    intro index hindex
    have hfactor :
        (denominator : ℤ) * (v index / (denominator : ℤ)) = v index :=
      Int.mul_ediv_cancel' (hdiv index hindex)
    calc
      (v index).natAbs =
          ((denominator : ℤ) *
            (v index / (denominator : ℤ))).natAbs := by rw [hfactor]
      _ = denominator * reducedMagnitude index := by
        rw [Int.natAbs_mul, Int.natAbs_natCast]
  have hpositive : ∀ index ∈ stalk, 0 < reducedMagnitude index := by
    intro index hindex
    rw [Nat.pos_iff_ne_zero]
    intro hzero
    have hspeedZero : (v index).natAbs = 0 := by
      rw [hspeedMagnitude index hindex, hzero, Nat.mul_zero]
    exact hv index (Int.natAbs_eq_zero.mp hspeedZero)
  have hltSeven : ∀ index ∈ stalk, reducedMagnitude index < 7 := by
    intro index hindex
    have hsmallIndex := hsmall index hindex
    have habsCast :
        |v index / (denominator : ℤ)| =
          (reducedMagnitude index : ℤ) := by
      simp [reducedMagnitude]
    rw [habsCast] at hsmallIndex
    omega
  have hinjective : Set.InjOn reducedMagnitude stalk := by
    intro first hfirst second hsecond hequal
    by_contra hne
    apply hdistinct first second hne
    have hnatAbs : (v first).natAbs = (v second).natAbs := by
      rw [hspeedMagnitude first hfirst, hspeedMagnitude second hsecond,
        hequal]
    have hcast : ((v first).natAbs : ℤ) = (v second).natAbs := by
      exact_mod_cast hnatAbs
    simpa only [Int.natCast_natAbs] using hcast
  have hsubset : stalk.image reducedMagnitude ⊆ Finset.Icc 1 6 := by
    intro magnitude hmagnitude
    rw [Finset.mem_image] at hmagnitude
    obtain ⟨index, hindex, rfl⟩ := hmagnitude
    rw [Finset.mem_Icc]
    have hpos := hpositive index hindex
    have hlt := hltSeven index hindex
    omega
  have hcardImage : (stalk.image reducedMagnitude).card = 7 := by
    rw [Finset.card_image_iff.mpr hinjective, hcard]
  have hcardBound := Finset.card_le_card hsubset
  rw [hcardImage, Nat.card_Icc] at hcardBound
  omega

/-- A connected zero-color bad stalk is either the exact reduced multiplier
fiber or is uniformly below the sharp reduced-speed window. -/
theorem exists_primitiveParameter_resonance_or_all_reduced_small
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0)
    (q p : ℕ) (hq : 0 < q)
    (stalk : Finset (Fin 13)) (root : Fin 13) (hroot : root ∈ stalk)
    (hconnected : ∀ index ∈ stalk,
      Relation.ReflTransGen
        (fun left right => overlapDet v q p left right = 0) root index)
    (hbad : ∀ index ∈ stalk, ¬ inBand v q p index) :
    ∃ denominator : ℕ, ∃ numerator : ℤ,
      0 < denominator ∧ numerator.natAbs.Coprime denominator ∧
      (∀ index ∈ stalk,
        (denominator : ℤ) ∣ v index ∧
        failWitness v q p index =
          (v index / (denominator : ℤ)) * numerator) ∧
      (((denominator = q / Nat.gcd p q ∧
          numerator = ((p / Nat.gcd p q : ℕ) : ℤ)) ∧
          ∀ index ∈ stalk,
            ((q / Nat.gcd p q : ℕ) : ℤ) ∣ v index) ∨
        ∀ index ∈ stalk,
          14 * |v index / (denominator : ℤ)| < (q : ℤ)) := by
  obtain ⟨denominator, numerator, hdenominator, hprimitive, hparameter⟩ :=
    exists_common_primitiveWitnessParameter_of_zero_connected
      v hv q p stalk root hroot hconnected
  refine ⟨denominator, numerator, hdenominator, hprimitive, hparameter, ?_⟩
  by_cases hresonance :
      (denominator : ℤ) * (p : ℤ) = numerator * (q : ℤ)
  · left
    have hnormalized :=
      primitiveWitnessParameter_eq_reducedMultiplier_of_resonance
        p q denominator numerator hq hprimitive hresonance
    refine ⟨⟨hnormalized.1, hnormalized.2⟩, ?_⟩
    intro index hindex
    have hreducedDenominator :
        ((q / Nat.gcd p q : ℕ) : ℤ) = (denominator : ℤ) := by
      exact_mod_cast hnormalized.1.symm
    rw [hreducedDenominator]
    exact (hparameter index hindex).1
  · right
    intro index hindex
    have hfork :=
      primitiveWitnessParameter_resonates_or_reducedSpeed_small_of_bad
        v q p index denominator numerator (hv index)
          (hparameter index hindex).1 (hparameter index hindex).2
          (bad_at_witness v q p index hq (hbad index hindex))
    exact hfork.resolve_left hresonance

/-- Coprime multipliers give the cleanest form: either the common primitive
parameter is exactly `(p,q)` and `q` divides every stalk speed, or every
reduced stalk speed lies below `q/14`. -/
theorem exists_primitiveParameter_coprime_resonance_or_all_reduced_small
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0)
    (q p : ℕ) (hq : 0 < q) (hpq : p.Coprime q)
    (stalk : Finset (Fin 13)) (root : Fin 13) (hroot : root ∈ stalk)
    (hconnected : ∀ index ∈ stalk,
      Relation.ReflTransGen
        (fun left right => overlapDet v q p left right = 0) root index)
    (hbad : ∀ index ∈ stalk, ¬ inBand v q p index) :
    ∃ denominator : ℕ, ∃ numerator : ℤ,
      0 < denominator ∧ numerator.natAbs.Coprime denominator ∧
      (∀ index ∈ stalk,
        (denominator : ℤ) ∣ v index ∧
        failWitness v q p index =
          (v index / (denominator : ℤ)) * numerator) ∧
      (((denominator = q ∧ numerator = (p : ℤ)) ∧
          ∀ index ∈ stalk, (q : ℤ) ∣ v index) ∨
        ∀ index ∈ stalk,
          14 * |v index / (denominator : ℤ)| < (q : ℤ)) := by
  simpa [hpq.gcd_eq_one] using
    exists_primitiveParameter_resonance_or_all_reduced_small
      v hv q p hq stalk root hroot hconnected hbad

/-- On a seven-runner stalk with distinct absolute speeds, the nonresonant
side of the sharp fork is impossible for every modulus `q ≤ 98`.  Hence the
shared primitive parameter is exactly the reduced multiplier, and the reduced
modulus divides every stalk speed. -/
theorem exists_primitiveParameter_resonance_of_seven_stalk_q_le_98
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0)
    (hdistinct : ∀ first second, first ≠ second →
      |v first| ≠ |v second|)
    (q p : ℕ) (hq : 0 < q) (hq98 : q ≤ 98)
    (stalk : Finset (Fin 13)) (hcard : stalk.card = 7)
    (root : Fin 13) (hroot : root ∈ stalk)
    (hconnected : ∀ index ∈ stalk,
      Relation.ReflTransGen
        (fun left right => overlapDet v q p left right = 0) root index)
    (hbad : ∀ index ∈ stalk, ¬ inBand v q p index) :
    ∃ denominator : ℕ, ∃ numerator : ℤ,
      0 < denominator ∧ numerator.natAbs.Coprime denominator ∧
      (∀ index ∈ stalk,
        (denominator : ℤ) ∣ v index ∧
        failWitness v q p index =
          (v index / (denominator : ℤ)) * numerator) ∧
      (denominator = q / Nat.gcd p q ∧
        numerator = ((p / Nat.gcd p q : ℕ) : ℤ)) ∧
      ∀ index ∈ stalk,
        ((q / Nat.gcd p q : ℕ) : ℤ) ∣ v index := by
  obtain ⟨denominator, numerator, hdenominator, hprimitive,
      hparameter, hresonance | hsmall⟩ :=
    exists_primitiveParameter_resonance_or_all_reduced_small
      v hv q p hq stalk root hroot hconnected hbad
  · exact ⟨denominator, numerator, hdenominator, hprimitive,
      hparameter, hresonance.1, hresonance.2⟩
  · exact False.elim
      (not_all_reduced_small_of_seven_distinct
        v hv hdistinct q denominator hq98 stalk hcard
          (fun index hindex => (hparameter index hindex).1) hsmall)

/-- Coprime multipliers in the same finite window give the strongest divisor
form: the common primitive parameter is `(p,q)` and `q` divides all seven
stalk speeds. -/
theorem exists_primitiveParameter_coprime_resonance_of_seven_stalk_q_le_98
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0)
    (hdistinct : ∀ first second, first ≠ second →
      |v first| ≠ |v second|)
    (q p : ℕ) (hq : 0 < q) (hq98 : q ≤ 98) (hpq : p.Coprime q)
    (stalk : Finset (Fin 13)) (hcard : stalk.card = 7)
    (root : Fin 13) (hroot : root ∈ stalk)
    (hconnected : ∀ index ∈ stalk,
      Relation.ReflTransGen
        (fun left right => overlapDet v q p left right = 0) root index)
    (hbad : ∀ index ∈ stalk, ¬ inBand v q p index) :
    ∃ denominator : ℕ, ∃ numerator : ℤ,
      0 < denominator ∧ numerator.natAbs.Coprime denominator ∧
      (∀ index ∈ stalk,
        (denominator : ℤ) ∣ v index ∧
        failWitness v q p index =
          (v index / (denominator : ℤ)) * numerator) ∧
      (denominator = q ∧ numerator = (p : ℤ)) ∧
      ∀ index ∈ stalk, (q : ℤ) ∣ v index := by
  simpa [hpq.gcd_eq_one] using
    exists_primitiveParameter_resonance_of_seven_stalk_q_le_98
      v hv hdistinct q p hq hq98 stalk hcard root hroot hconnected hbad

#print axioms not_all_reduced_small_of_seven_distinct
#print axioms exists_primitiveParameter_resonance_or_all_reduced_small
#print axioms exists_primitiveParameter_coprime_resonance_or_all_reduced_small
#print axioms exists_primitiveParameter_resonance_of_seven_stalk_q_le_98
#print axioms exists_primitiveParameter_coprime_resonance_of_seven_stalk_q_le_98

end LRCZeroColorStalkFork
end LonelyRunner
