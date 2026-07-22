import TournamentH7.GMC2Henselian

/-!
# Reciprocal Hensel lifting for the small roots of `X^M - t R(X)`

After the ramified substitution `t = s^M`, `X = s Z`, the small-root
equation is

`Z^M - R(s Z) = 0` over `F⟦s⟧`.

This polynomial need not be monic: its degree can be larger than `M`, while
its reduction modulo `s` has degree `M`.  Its constant coefficient is a unit,
however.  Reversing the polynomial and scaling by the inverse of that
constant coefficient produces a monic polynomial.  A simple nonzero root of
the reversed residue polynomial therefore lifts by ordinary monic Hensel;
the lift remains a unit, and taking its reciprocal returns a root of the
original nonmonic polynomial.

The two theorems below isolate that reusable mechanism.  The second theorem
specializes the local-ring bookkeeping to `F⟦s⟧`: it reduces a small-root
construction to two explicit constant-coefficient checks for the monicized
reverse.  Those checks are the finite polynomial algebra left in the
specialization to `Z^M - R(s Z)`.
-/

open Polynomial

namespace GMC2ReciprocalSmallRoots

section LocalRing

variable {A : Type*} [CommRing A] [HenselianLocalRing A]

/-- Scale the reverse of `p` by the inverse of its unit constant coefficient.
The resulting polynomial is monic. -/
noncomputable def monicReverse (p : A[X]) (h0 : IsUnit (p.coeff 0)) : A[X] :=
  Polynomial.C (↑(h0.unit⁻¹) : A) * p.reverse

/-- `monicReverse` is genuinely monic; this is the step that removes the
degree-drop obstruction from the public Henselian API. -/
theorem monic_monicReverse (p : A[X]) (h0 : IsUnit (p.coeff 0)) :
    (monicReverse p h0).Monic := by
  rw [Polynomial.Monic.def, monicReverse,
    Polynomial.leadingCoeff_C_mul_of_isUnit (h0.unit⁻¹).isUnit,
    Polynomial.reverse_leadingCoeff,
    Polynomial.trailingCoeff_eq_coeff_zero h0.ne_zero]
  exact h0.val_inv_mul

/-- **Reciprocal simple-root Hensel lift.**

Suppose the monicized reverse of `p` has a simple approximate root `b₀`
modulo the maximal ideal, and `b₀` is a unit.  Monic Hensel produces a root
`b`.  Since `b ≡ b₀` modulo the maximal ideal, `b` is still a unit; its
reciprocal `z` is then a root of the original polynomial `p`.

Returning both `b` and `z`, together with `b*z=1`, avoids hiding a choice of
inverse in the statement and is convenient for residue-field calculations. -/
theorem exists_unit_reciprocal_root
    (p : A[X]) (h0 : IsUnit (p.coeff 0))
    (b₀ : A) (hb₀ : IsUnit b₀)
    (heval : (monicReverse p h0).eval b₀ ∈ IsLocalRing.maximalIdeal A)
    (hderiv : IsUnit ((monicReverse p h0).derivative.eval b₀)) :
    ∃ b z : A,
      IsUnit b ∧ b * z = 1 ∧
      (monicReverse p h0).IsRoot b ∧ p.IsRoot z ∧
      b - b₀ ∈ IsLocalRing.maximalIdeal A := by
  obtain ⟨b, hbroot, hbclose⟩ :=
    HenselianLocalRing.is_henselian
      (monicReverse p h0) (monic_monicReverse p h0) b₀ heval hderiv
  have hb : IsUnit b := by
    rw [← IsLocalRing.notMem_maximalIdeal]
    intro hbmem
    have hb₀mem : b₀ ∈ IsLocalRing.maximalIdeal A := by
      have hsub := (IsLocalRing.maximalIdeal A).sub_mem hbmem hbclose
      simpa using hsub
    exact (IsLocalRing.notMem_maximalIdeal.mpr hb₀) hb₀mem
  have hreverse : p.reverse.IsRoot b := by
    rw [Polynomial.IsRoot.def] at hbroot ⊢
    have hmul : (↑(h0.unit⁻¹) : A) * p.reverse.eval b = 0 := by
      simpa [monicReverse] using hbroot
    exact ((h0.unit⁻¹).isUnit.mul_right_eq_zero).mp hmul
  letI : Invertible b := hb.invertible
  refine ⟨b, ⅟ b, hb, by simp, hbroot, ?_, hbclose⟩
  rw [Polynomial.IsRoot.def]
  apply (Polynomial.eval₂_reverse_eq_zero_iff (RingHom.id A) (⅟ b) p).mp
  simpa using hreverse.eq_zero

end LocalRing

section PowerSeries

variable {F : Type*} [Field F]

/-- Over `F⟦s⟧`, reciprocal Hensel reduces root construction to two
constant-coefficient checks.  If the monicized reverse vanishes at
`b₀ = C(a⁻¹)` modulo `s` and its derivative is nonzero there modulo `s`,
then the original polynomial has a root with constant coefficient `a`.

For `p(Z)=Z^M-R(sZ)`, the two hypotheses reduce respectively to
`a^M=R(0)` and `(M:F)≠0`; this finite specialization is deliberately kept
separate from the local-ring mechanism. -/
theorem exists_root_with_constantCoeff_of_monicReverse
    (p : (PowerSeries F)[X]) (h0 : IsUnit (p.coeff 0))
    (a : F) (ha : a ≠ 0)
    (heval : PowerSeries.constantCoeff
        ((monicReverse p h0).eval (PowerSeries.C a⁻¹)) = 0)
    (hderiv : PowerSeries.constantCoeff
        ((monicReverse p h0).derivative.eval (PowerSeries.C a⁻¹)) ≠ 0) :
    ∃ z : PowerSeries F, p.IsRoot z ∧ PowerSeries.constantCoeff z = a := by
  have hb₀ : IsUnit (PowerSeries.C a⁻¹) := by
    rw [PowerSeries.isUnit_iff_constantCoeff, PowerSeries.constantCoeff_C]
    exact isUnit_iff_ne_zero.mpr (inv_ne_zero ha)
  have hevalIdeal :
      (monicReverse p h0).eval (PowerSeries.C a⁻¹) ∈
        IsLocalRing.maximalIdeal (PowerSeries F) := by
    rw [PowerSeries.maximalIdeal_eq_span_X, Ideal.mem_span_singleton,
      PowerSeries.X_dvd_iff]
    exact heval
  have hderivUnit :
      IsUnit ((monicReverse p h0).derivative.eval (PowerSeries.C a⁻¹)) := by
    rw [PowerSeries.isUnit_iff_constantCoeff]
    exact isUnit_iff_ne_zero.mpr hderiv
  obtain ⟨b, z, _hb, hbz, _hbroot, hpz, hbclose⟩ :=
    exists_unit_reciprocal_root p h0 (PowerSeries.C a⁻¹) hb₀
      hevalIdeal hderivUnit
  have hcbsub : PowerSeries.constantCoeff (b - PowerSeries.C a⁻¹) = 0 := by
    rw [← PowerSeries.X_dvd_iff, ← Ideal.mem_span_singleton,
      ← PowerSeries.maximalIdeal_eq_span_X]
    exact hbclose
  have hcb : PowerSeries.constantCoeff b = a⁻¹ := by
    rw [map_sub, PowerSeries.constantCoeff_C, sub_eq_zero] at hcbsub
    exact hcbsub
  have hczProduct : a⁻¹ * PowerSeries.constantCoeff z = 1 := by
    have h := congrArg PowerSeries.constantCoeff hbz
    simpa [map_mul, hcb] using h
  have hcz : PowerSeries.constantCoeff z = a := by
    have h := congrArg (fun x : F ↦ a * x) hczProduct
    simpa [mul_assoc, ha] using h
  exact ⟨z, hpz, hcz⟩

end PowerSeries

end GMC2ReciprocalSmallRoots

#print axioms GMC2ReciprocalSmallRoots.monic_monicReverse
#print axioms GMC2ReciprocalSmallRoots.exists_unit_reciprocal_root
#print axioms GMC2ReciprocalSmallRoots.exists_root_with_constantCoeff_of_monicReverse
