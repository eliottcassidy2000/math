import Mathlib.Algebra.CharP.Basic
import Mathlib.Data.Nat.Prime.Infinite
import Mathlib.NumberTheory.NumberField.Basic
import Mathlib.RingTheory.Ideal.GoingUp
import Mathlib.RingTheory.Ideal.Int
import Mathlib.RingTheory.Ideal.Norm.AbsNorm
import Mathlib.RingTheory.Ideal.Quotient.HasFiniteQuotients
import Mathlib.RingTheory.Localization.AtPrime.Basic
import Mathlib.RingTheory.Localization.FractionRing
import Mathlib.RingTheory.Spectrum.Maximal.Basic

/-!
# Good finite-place reduction for the GMC(2) argument

This module isolates two algebraic facts needed by the Frobenius
lowest-balanced-face proof.

* A finite family of nonzero algebraic integers admits an arbitrarily large
  rational prime `p` and a maximal ideal above `p` at which every member of the
  family is a local unit.  The residue field is finite of characteristic `p`,
  and every chosen element has nonzero residue.
* A common nonzero factor must be cancelled in the domain before passage to a
  quotient.  A nonzero reduced normalized sum therefore certifies that the raw
  sum was nonzero.

The first theorem is a genuine good-prime existence result, not an assumed
interface.  It is stated for algebraic integers.  The second existence theorem
extends it to arbitrary number-field elements by explicitly choosing finitely
many integral numerators and denominators and avoiding all of them at one
finite place.
-/

open Finset
open scoped NumberField

namespace GMC2GoodReduction

/-- An algebraic integer whose absolute norm is smaller than the rational
prime below `P` cannot lie in `P`.

The proof uses the prime-power absolute norm of a maximal ideal and the fact
that `P` lies over `(p)`. -/
theorem not_mem_prime_over_of_norm_lt
    {K : Type*} [Field K] [NumberField K]
    (p : ℕ) (hp : p.Prime) (P : Ideal (𝓞 K))
    (hPmax : P.IsMaximal)
    (hPover : P.comap (algebraMap ℤ (𝓞 K)) = Ideal.span {(p : ℤ)})
    {x : 𝓞 K} (hx : x ≠ 0)
    (hnorm : (Algebra.norm ℤ x).natAbs < p) :
    x ∉ P := by
  letI : P.IsMaximal := hPmax
  intro hxP
  obtain ⟨q, n, hn, hqP, hq, hPnorm⟩ :=
    Ideal.exists_prime_and_absNorm_eq_pow P
  have hqUnder : (q : ℤ) ∈ P.comap (algebraMap ℤ (𝓞 K)) := by
    rw [Ideal.mem_comap]
    have hcast : algebraMap ℤ (𝓞 K) (q : ℤ) = (q : 𝓞 K) := by simp
    rw [hcast]
    exact hqP
  rw [hPover, Ideal.mem_span_singleton] at hqUnder
  have hpq : p ∣ q := Int.natCast_dvd_natCast.mp hqUnder
  have hpqEq : p = q := (Nat.prime_dvd_prime_iff_eq hp hq).mp hpq
  have hpAbsNorm : p ∣ P.absNorm := by
    rw [hPnorm, hpqEq]
    exact dvd_pow_self q (Nat.ne_of_gt hn)
  have hpNormInt : (p : ℤ) ∣ Algebra.norm ℤ x :=
    (Int.natCast_dvd_natCast.mpr hpAbsNorm).trans
      (Ideal.absNorm_dvd_norm_of_mem hxP)
  have hpNorm : p ∣ (Algebra.norm ℤ x).natAbs :=
    Int.natCast_dvd.mp hpNormInt
  have hnormPos : 0 < (Algebra.norm ℤ x).natAbs := by
    rw [Int.natAbs_pos]
    exact (Algebra.norm_ne_zero_iff.mpr hx)
  exact (not_le_of_gt hnorm) (Nat.le_of_dvd hnormPos hpNorm)

/-- **Finite-family good reduction over a number field.**

Given finitely many nonzero algebraic integers and a lower bound, there is a
larger rational prime `p` and a maximal ideal `P` above `(p)` such that:

* the residue field `(𝓞 K) ⧸ P` is finite of characteristic `p`;
* every selected algebraic integer has nonzero residue;
* every selected algebraic integer is a unit in the local ring `(𝓞 K)_P`.

The prime is chosen larger than every absolute field norm in the family.
Consequently no prime above it can contain a selected element. -/
theorem exists_good_reduction_for_integral_family
    {K ι : Type*} [Field K] [NumberField K]
    (bound : ℕ) (S : Finset ι) (x : ι → 𝓞 K)
    (hx : ∀ i ∈ S, x i ≠ 0) :
    ∃ p : ℕ, ∃ P : MaximalSpectrum (𝓞 K),
      p.Prime ∧ bound < p ∧
      P.asIdeal.comap (algebraMap ℤ (𝓞 K)) = Ideal.span {(p : ℤ)} ∧
      Finite ((𝓞 K) ⧸ P.asIdeal) ∧ CharP ((𝓞 K) ⧸ P.asIdeal) p ∧
      (∀ i ∈ S, Ideal.Quotient.mk P.asIdeal (x i) ≠ 0) ∧
      ∀ i ∈ S,
        IsUnit (algebraMap (𝓞 K) (Localization.AtPrime P.asIdeal) (x i)) := by
  let normBound : ℕ := S.sup fun i => (Algebra.norm ℤ (x i)).natAbs
  obtain ⟨p, hpLower, hp⟩ :=
    Nat.exists_infinite_primes (max bound normBound + 1)
  have hbound : bound < p := by
    exact lt_of_lt_of_le (Nat.lt_succ_iff.mpr (Nat.le_max_left _ _)) hpLower
  have hnormBound : normBound < p := by
    exact lt_of_lt_of_le (Nat.lt_succ_iff.mpr (Nat.le_max_right _ _)) hpLower
  letI : Fact p.Prime := ⟨hp⟩
  have hker : RingHom.ker (algebraMap ℤ (𝓞 K)) = ⊥ :=
    (RingHom.injective_iff_ker_eq_bot _).mp (RingHom.injective_int _)
  obtain ⟨P0, hPmax, hPover⟩ :=
    Ideal.exists_ideal_over_maximal_of_isIntegral
      (R := ℤ) (S := 𝓞 K) (Ideal.span {(p : ℤ)}) (hker ▸ bot_le)
  let P : MaximalSpectrum (𝓞 K) := ⟨P0, hPmax⟩
  have hPover' :
      P.asIdeal.comap (algebraMap ℤ (𝓞 K)) = Ideal.span {(p : ℤ)} := by
    simpa only [P] using hPover
  have hPne : P.asIdeal ≠ ⊥ := NeZero.ne P.asIdeal
  have hfinite : Finite ((𝓞 K) ⧸ P.asIdeal) :=
    Ring.HasFiniteQuotients.finiteQuotient hPne
  have hpMem : (p : 𝓞 K) ∈ P.asIdeal := by
    have : (p : ℤ) ∈ P.asIdeal.comap (algebraMap ℤ (𝓞 K)) := by
      rw [hPover']
      exact Ideal.mem_span_singleton_self (p : ℤ)
    rw [Ideal.mem_comap] at this
    have hcast : algebraMap ℤ (𝓞 K) (p : ℤ) = (p : 𝓞 K) := by simp
    rwa [hcast] at this
  have hpZero : (p : (𝓞 K) ⧸ P.asIdeal) = 0 := by
    change Ideal.Quotient.mk P.asIdeal (p : 𝓞 K) = 0
    exact Ideal.Quotient.eq_zero_iff_mem.mpr hpMem
  have hchar : CharP ((𝓞 K) ⧸ P.asIdeal) p :=
    (CharP.charP_iff_prime_eq_zero hp).2 hpZero
  have havoid : ∀ i ∈ S, x i ∉ P.asIdeal := by
    intro i hiS
    apply not_mem_prime_over_of_norm_lt p hp P.asIdeal P.isMaximal hPover' (hx i hiS)
    exact (Finset.le_sup (s := S) (f := fun j => (Algebra.norm ℤ (x j)).natAbs)
      hiS).trans_lt hnormBound
  refine ⟨p, P, hp, hbound, hPover', hfinite, hchar, ?_, ?_⟩
  · intro i hiS
    exact (Ideal.Quotient.eq_zero_iff_mem.not).2 (havoid i hiS)
  · intro i hiS
    exact IsLocalization.map_units (Localization.AtPrime P.asIdeal)
      (⟨x i, havoid i hiS⟩ : P.asIdeal.primeCompl)

/-- Reduction of a fraction of algebraic integers at a bundled maximal ideal.
The field structure is installed only inside this definition. -/
noncomputable def fractionResidue
    {K : Type*} [Field K] [NumberField K]
    (P : MaximalSpectrum (𝓞 K)) (num den : 𝓞 K) :
    (𝓞 K) ⧸ P.asIdeal := by
  letI := Ideal.Quotient.field P.asIdeal
  exact Ideal.Quotient.mk P.asIdeal num / Ideal.Quotient.mk P.asIdeal den

theorem fractionResidue_ne_zero
    {K : Type*} [Field K] [NumberField K]
    (P : MaximalSpectrum (𝓞 K)) {num den : 𝓞 K}
    (hnum : Ideal.Quotient.mk P.asIdeal num ≠ 0)
    (hden : Ideal.Quotient.mk P.asIdeal den ≠ 0) :
    fractionResidue P num den ≠ 0 := by
  letI := Ideal.Quotient.field P.asIdeal
  exact div_ne_zero hnum hden

/-- Fraction reduction is independent of the chosen integral presentation as
long as both denominators survive at `P`.  This is the precise local
specialization statement replacing the impossible mixed-characteristic field
homomorphism `K → (𝓞 K) ⧸ P`. -/
theorem fractionResidue_eq_of_eq
    {K : Type*} [Field K] [NumberField K]
    (P : MaximalSpectrum (𝓞 K))
    {num₁ den₁ num₂ den₂ : 𝓞 K}
    (hden₁ : Ideal.Quotient.mk P.asIdeal den₁ ≠ 0)
    (hden₂ : Ideal.Quotient.mk P.asIdeal den₂ ≠ 0)
    (hfrac :
      algebraMap (𝓞 K) K num₁ / algebraMap (𝓞 K) K den₁ =
        algebraMap (𝓞 K) K num₂ / algebraMap (𝓞 K) K den₂) :
    fractionResidue P num₁ den₁ = fractionResidue P num₂ den₂ := by
  letI := Ideal.Quotient.field P.asIdeal
  change Ideal.Quotient.mk P.asIdeal num₁ / Ideal.Quotient.mk P.asIdeal den₁ =
    Ideal.Quotient.mk P.asIdeal num₂ / Ideal.Quotient.mk P.asIdeal den₂
  apply (div_eq_div_iff hden₁ hden₂).2
  have hden₁O : den₁ ≠ 0 := by
    intro hzero
    apply hden₁
    simp [hzero]
  have hden₂O : den₂ ≠ 0 := by
    intro hzero
    apply hden₂
    simp [hzero]
  have hden₁K : algebraMap (𝓞 K) K den₁ ≠ 0 :=
    IsFractionRing.to_map_ne_zero_of_mem_nonZeroDivisors
      (mem_nonZeroDivisors_iff_ne_zero.mpr hden₁O)
  have hden₂K : algebraMap (𝓞 K) K den₂ ≠ 0 :=
    IsFractionRing.to_map_ne_zero_of_mem_nonZeroDivisors
      (mem_nonZeroDivisors_iff_ne_zero.mpr hden₂O)
  have hcrossK :
      algebraMap (𝓞 K) K (num₁ * den₂) =
        algebraMap (𝓞 K) K (num₂ * den₁) := by
    simpa only [map_mul] using (div_eq_div_iff hden₁K hden₂K).1 hfrac
  have hcrossO : num₁ * den₂ = num₂ * den₁ :=
    IsFractionRing.injective (𝓞 K) K hcrossK
  simpa only [map_mul] using
    congrArg (Ideal.Quotient.mk P.asIdeal) hcrossO

/-- Good reduction for a finite family of arbitrary nonzero number-field
elements.

The theorem chooses integral numerator/denominator presentations, then chooses
one good prime avoiding every numerator and denominator.  The resulting chosen
fraction residue is nonzero for each original field element.  Both integral
pieces are also exhibited as units in the corresponding local ring.  This does
not assert a mixed-characteristic field homomorphism from all of `K`. -/
theorem exists_good_reduction_for_number_field_family
    {K ι : Type*} [Field K] [NumberField K]
    (bound : ℕ) (S : Finset ι) (c : ι → K)
    (hc : ∀ i ∈ S, c i ≠ 0) :
    ∃ p : ℕ, ∃ P : MaximalSpectrum (𝓞 K),
      ∃ num den : ι → 𝓞 K,
      p.Prime ∧ bound < p ∧
      P.asIdeal.comap (algebraMap ℤ (𝓞 K)) = Ideal.span {(p : ℤ)} ∧
      Finite ((𝓞 K) ⧸ P.asIdeal) ∧ CharP ((𝓞 K) ⧸ P.asIdeal) p ∧
      (∀ i ∈ S,
        algebraMap (𝓞 K) K (num i) / algebraMap (𝓞 K) K (den i) = c i) ∧
      (∀ i ∈ S, fractionResidue P (num i) (den i) ≠ 0) ∧
      ∀ i ∈ S,
        IsUnit (algebraMap (𝓞 K) (Localization.AtPrime P.asIdeal) (num i)) ∧
        IsUnit (algebraMap (𝓞 K) (Localization.AtPrime P.asIdeal) (den i)) := by
  classical
  choose num den hden hfrac using fun i =>
    IsFractionRing.div_surjective (𝓞 K) (c i)
  have hnumNe : ∀ i ∈ S, num i ≠ 0 := by
    intro i hiS hzero
    apply hc i hiS
    rw [← hfrac i, hzero, map_zero, zero_div]
  have hdenNe : ∀ i ∈ S, den i ≠ 0 := by
    intro i hiS
    exact mem_nonZeroDivisors_iff_ne_zero.mp (hden i)
  let T : Finset (ι ⊕ ι) := S.image Sum.inl ∪ S.image Sum.inr
  let z : ι ⊕ ι → 𝓞 K := Sum.elim num den
  have hz : ∀ j ∈ T, z j ≠ 0 := by
    intro j hj
    rcases Finset.mem_union.1 hj with hj | hj
    · rcases Finset.mem_image.1 hj with ⟨i, hiS, rfl⟩
      exact hnumNe i hiS
    · rcases Finset.mem_image.1 hj with ⟨i, hiS, rfl⟩
      exact hdenNe i hiS
  obtain ⟨p, P, hp, hbound, hPover, hfinite, hchar, hresidue, hlocal⟩ :=
    exists_good_reduction_for_integral_family bound T z hz
  have hnumResidue : ∀ i ∈ S, Ideal.Quotient.mk P.asIdeal (num i) ≠ 0 := by
    intro i hiS
    exact hresidue (Sum.inl i) (Finset.mem_union_left _ <| Finset.mem_image_of_mem _ hiS)
  have hdenResidue : ∀ i ∈ S, Ideal.Quotient.mk P.asIdeal (den i) ≠ 0 := by
    intro i hiS
    exact hresidue (Sum.inr i) (Finset.mem_union_right _ <| Finset.mem_image_of_mem _ hiS)
  have hnumLocal : ∀ i ∈ S,
      IsUnit (algebraMap (𝓞 K) (Localization.AtPrime P.asIdeal) (num i)) := by
    intro i hiS
    exact hlocal (Sum.inl i) (Finset.mem_union_left _ <| Finset.mem_image_of_mem _ hiS)
  have hdenLocal : ∀ i ∈ S,
      IsUnit (algebraMap (𝓞 K) (Localization.AtPrime P.asIdeal) (den i)) := by
    intro i hiS
    exact hlocal (Sum.inr i) (Finset.mem_union_right _ <| Finset.mem_image_of_mem _ hiS)
  refine ⟨p, P, num, den, hp, hbound, hPover, hfinite, hchar, ?_, ?_, ?_⟩
  · intro i hiS
    exact hfrac i
  · intro i hiS
    exact fractionResidue_ne_zero P (hnumResidue i hiS) (hdenResidue i hiS)
  · intro i hiS
    exact ⟨hnumLocal i hiS, hdenLocal i hiS⟩

/-- If every raw summand has a common factor `d`, a vanishing raw sum can be
cancelled in a domain before reduction.  No inverse of the residue of `d` is
used or needed. -/
theorem normalized_sum_eq_zero_of_raw_sum_eq_zero
    {A ι : Type*} [CommRing A] [IsDomain A]
    (S : Finset ι) (d : A) (raw normalized : ι → A)
    (hd : d ≠ 0) (hterm : ∀ i ∈ S, raw i = d * normalized i)
    (hraw : ∑ i ∈ S, raw i = 0) :
    ∑ i ∈ S, normalized i = 0 := by
  have hfactor : (∑ i ∈ S, raw i) = d * ∑ i ∈ S, normalized i := by
    calc
      (∑ i ∈ S, raw i) = ∑ i ∈ S, d * normalized i := by
        apply Finset.sum_congr rfl
        intro i hiS
        exact hterm i hiS
      _ = d * ∑ i ∈ S, normalized i := by
        rw [Finset.mul_sum]
  rw [hfactor] at hraw
  exact (mul_eq_zero.mp hraw).resolve_left hd

/-- Consumer form of cancellation: if the normalized sum has nonzero image in
any quotient, then the raw sum was nonzero.  This is the safe formal pattern
for cancelling `(p*A₀)!` before reducing modulo the good prime. -/
theorem raw_sum_ne_zero_of_normalized_reduction_ne_zero
    {A ι : Type*} [CommRing A] [IsDomain A]
    (I : Ideal A) (S : Finset ι) (d : A)
    (raw normalized : ι → A)
    (hd : d ≠ 0) (hterm : ∀ i ∈ S, raw i = d * normalized i)
    (hreduced : Ideal.Quotient.mk I (∑ i ∈ S, normalized i) ≠ 0) :
    ∑ i ∈ S, raw i ≠ 0 := by
  intro hraw
  have hnormalized := normalized_sum_eq_zero_of_raw_sum_eq_zero
    S d raw normalized hd hterm hraw
  apply hreduced
  rw [hnormalized, map_zero]

end GMC2GoodReduction

#print axioms GMC2GoodReduction.exists_good_reduction_for_integral_family
#print axioms GMC2GoodReduction.fractionResidue_eq_of_eq
#print axioms GMC2GoodReduction.exists_good_reduction_for_number_field_family
#print axioms GMC2GoodReduction.raw_sum_ne_zero_of_normalized_reduction_ne_zero
