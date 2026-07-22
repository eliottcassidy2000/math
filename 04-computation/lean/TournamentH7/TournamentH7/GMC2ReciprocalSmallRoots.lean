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

/-- The polynomial in the rescaled root `Z` obtained from
`X^M - t R(X)` after `t = s^M` and `X = s Z`, with the common factor
`s^M` removed.  Thus this is literally `Z^M - R(s Z)`. -/
noncomputable def ramifiedPolynomial (M : ℕ) (R : F[X]) :
    (PowerSeries F)[X] :=
  X ^ M - (R.map PowerSeries.C).comp (C PowerSeries.X * X)

/-- Reduction modulo the power-series parameter forgets every positive
`Z`-coefficient contributed by `R(sZ)`. -/
theorem map_constantCoeff_ramifiedPolynomial (M : ℕ) (R : F[X]) :
    (ramifiedPolynomial M R).map PowerSeries.constantCoeff =
      X ^ M - C (R.coeff 0) := by
  simp [ramifiedPolynomial, Polynomial.map_comp, Polynomial.coeff_zero_eq_eval_zero]

theorem constantCoeff_coeff_zero_ramifiedPolynomial
    (M : ℕ) (R : F[X]) (hM : M ≠ 0) :
    PowerSeries.constantCoeff ((ramifiedPolynomial M R).coeff 0) =
      -R.coeff 0 := by
  have h := congrArg (fun q : F[X] ↦ q.coeff 0)
    (map_constantCoeff_ramifiedPolynomial M R)
  simpa [hM, Ne.symm hM] using h

theorem isUnit_coeff_zero_ramifiedPolynomial
    (M : ℕ) (R : F[X]) (hM : M ≠ 0) (hr₀ : R.coeff 0 ≠ 0) :
    IsUnit ((ramifiedPolynomial M R).coeff 0) := by
  rw [PowerSeries.isUnit_iff_constantCoeff,
    constantCoeff_coeff_zero_ramifiedPolynomial M R hM]
  exact isUnit_iff_ne_zero.mpr (neg_ne_zero.mpr hr₀)

theorem natDegree_ramifiedPolynomial
    (M : ℕ) (R : F[X]) (hdeg : M < R.natDegree) :
    (ramifiedPolynomial M R).natDegree = R.natDegree := by
  have hmap : (R.map PowerSeries.C).natDegree = R.natDegree :=
    Polynomial.natDegree_map_eq_of_injective PowerSeries.C_injective R
  have hX : (PowerSeries.X : PowerSeries F) ≠ 0 := by simp
  have hinner : (C (PowerSeries.X : PowerSeries F) * X : (PowerSeries F)[X]).natDegree = 1 :=
    Polynomial.natDegree_C_mul_X _ hX
  have hcomp :
      ((R.map PowerSeries.C).comp
        (C (PowerSeries.X : PowerSeries F) * X)).natDegree = R.natDegree := by
    rw [Polynomial.natDegree_comp, hmap, hinner, mul_one]
  unfold ramifiedPolynomial
  rw [← hcomp]
  apply Polynomial.natDegree_sub_eq_right_of_natDegree_lt
  rw [Polynomial.natDegree_X_pow, hcomp]
  exact hdeg

/-- Before monic scaling, reversing the genuine degree-`d` ramified
polynomial and then reducing modulo `s` exposes the factor
`W^(d-M) (W^M-r₀)`.  This identity is where keeping the pre-reduction
degree matters. -/
theorem map_constantCoeff_reverse_ramifiedPolynomial
    (M : ℕ) (R : F[X]) (hdeg : M < R.natDegree) :
    (ramifiedPolynomial M R).reverse.map PowerSeries.constantCoeff =
      X ^ (R.natDegree - M) - C (R.coeff 0) * X ^ R.natDegree := by
  rw [Polynomial.reverse, natDegree_ramifiedPolynomial M R hdeg,
    ← Polynomial.reflect_map, map_constantCoeff_ramifiedPolynomial,
    Polynomial.reflect_sub, Polynomial.reflect_monomial, Polynomial.reflect_C,
    Polynomial.revAt_le (Nat.le_of_lt hdeg)]

/-- The nonzero residue root `a⁻¹` of the reversed polynomial. -/
theorem eval_reverseResidue_at_inv
    (M d : ℕ) (r₀ a : F) (hMd : M ≤ d) (ha : a ≠ 0)
    (hpow : a ^ M = r₀) :
    (X ^ (d - M) - C r₀ * X ^ d : F[X]).eval a⁻¹ = 0 := by
  have hunit : r₀ * (a⁻¹) ^ M = 1 := by
    rw [← hpow, inv_pow]
    exact mul_inv_cancel₀ (pow_ne_zero M ha)
  have hpowd : (a⁻¹) ^ d = (a⁻¹) ^ (d - M) * (a⁻¹) ^ M := by
    calc
      (a⁻¹) ^ d = (a⁻¹) ^ ((d - M) + M) := by
        exact congrArg (fun n : ℕ ↦ (a⁻¹) ^ n) (by omega)
      _ = (a⁻¹) ^ (d - M) * (a⁻¹) ^ M := by rw [pow_add]
  rw [eval_sub, eval_pow, eval_X, eval_mul, eval_C, eval_pow, eval_X, hpowd]
  calc
    (a⁻¹) ^ (d - M) - r₀ * ((a⁻¹) ^ (d - M) * (a⁻¹) ^ M) =
        (a⁻¹) ^ (d - M) * (1 - r₀ * (a⁻¹) ^ M) := by ring
    _ = 0 := by rw [hunit]; simp

/-- The derivative at the same residue root is exactly
`-M a⁻⁽ᵈ⁻ᴹ⁻¹⁾`.  This makes separability depend only on the expected
characteristic condition and `a ≠ 0`. -/
theorem derivative_eval_reverseResidue_at_inv
    (M d : ℕ) (r₀ a : F) (hMd : M < d) (ha : a ≠ 0)
    (hpow : a ^ M = r₀) :
    (X ^ (d - M) - C r₀ * X ^ d : F[X]).derivative.eval a⁻¹ =
      -(M : F) * (a⁻¹) ^ (d - M - 1) := by
  have hk : 0 < d - M := Nat.sub_pos_of_lt hMd
  have hsum : d - 1 = (d - M - 1) + M := by omega
  have hunit : r₀ * (a⁻¹) ^ M = 1 := by
    rw [← hpow, inv_pow]
    exact mul_inv_cancel₀ (pow_ne_zero M ha)
  simp only [Polynomial.derivative_sub, Polynomial.derivative_pow,
    Polynomial.derivative_X, mul_one, Polynomial.derivative_mul,
    Polynomial.derivative_C, zero_mul, zero_add, eval_sub, eval_mul,
    eval_pow, eval_X, eval_C]
  rw [hsum, pow_add]
  calc
    (((d - M : ℕ) : F)) * (a⁻¹) ^ (d - M - 1) -
        r₀ * ((d : F) * ((a⁻¹) ^ (d - M - 1) * (a⁻¹) ^ M)) =
      ((((d - M : ℕ) : F) - (d : F) * (r₀ * (a⁻¹) ^ M))) *
        (a⁻¹) ^ (d - M - 1) := by ring
    _ = (((d - M : ℕ) : F) - (d : F)) * (a⁻¹) ^ (d - M - 1) := by
      rw [hunit, mul_one]
    _ = -(M : F) * (a⁻¹) ^ (d - M - 1) := by
      have hd : d - M + M = d := Nat.sub_add_cancel (Nat.le_of_lt hMd)
      have hcast := congrArg (fun n : ℕ ↦ (n : F)) hd
      push_cast at hcast
      rw [← hcast]
      ring

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

/-- **Concrete small-root specialization.**

If `R` has genuine degree larger than `M`, nonzero constant coefficient,
`a^M = R(0)`, and the characteristic does not divide `M`, then the
ramified polynomial `Z^M-R(sZ)` has a power-series root reducing to `a`.
No degree-preserving Hensel factorization is used: the proof applies simple
root Hensel to the monicized reverse. -/
theorem exists_root_ramifiedPolynomial
    (M : ℕ) (R : F[X]) (hdeg : M < R.natDegree)
    (hr₀ : R.coeff 0 ≠ 0) (a : F) (ha : a ≠ 0)
    (hpow : a ^ M = R.coeff 0) (hMchar : (M : F) ≠ 0) :
    ∃ z : PowerSeries F,
      (ramifiedPolynomial M R).IsRoot z ∧
        PowerSeries.constantCoeff z = a := by
  have hM : M ≠ 0 := by
    intro h
    subst M
    simp at hMchar
  let h₀ : IsUnit ((ramifiedPolynomial M R).coeff 0) :=
    isUnit_coeff_zero_ramifiedPolynomial M R hM hr₀
  have hreverseEval :
      PowerSeries.constantCoeff
          ((ramifiedPolynomial M R).reverse.eval (PowerSeries.C a⁻¹)) = 0 := by
    calc
      PowerSeries.constantCoeff
          ((ramifiedPolynomial M R).reverse.eval (PowerSeries.C a⁻¹)) =
          ((ramifiedPolynomial M R).reverse.map PowerSeries.constantCoeff).eval a⁻¹ := by
            symm
            simpa using
              (Polynomial.eval_map_apply
                (p := (ramifiedPolynomial M R).reverse)
                PowerSeries.constantCoeff (PowerSeries.C a⁻¹))
      _ = 0 := by
        rw [map_constantCoeff_reverse_ramifiedPolynomial M R hdeg]
        exact eval_reverseResidue_at_inv M R.natDegree (R.coeff 0) a
          (Nat.le_of_lt hdeg) ha hpow
  have hreverseDeriv :
      PowerSeries.constantCoeff
          ((ramifiedPolynomial M R).reverse.derivative.eval
            (PowerSeries.C a⁻¹)) ≠ 0 := by
    have heq :
        PowerSeries.constantCoeff
            ((ramifiedPolynomial M R).reverse.derivative.eval
              (PowerSeries.C a⁻¹)) =
          -(M : F) * (a⁻¹) ^ (R.natDegree - M - 1) := by
      calc
        PowerSeries.constantCoeff
            ((ramifiedPolynomial M R).reverse.derivative.eval
              (PowerSeries.C a⁻¹)) =
            ((ramifiedPolynomial M R).reverse.derivative.map
              PowerSeries.constantCoeff).eval a⁻¹ := by
                symm
                simpa using
                  (Polynomial.eval_map_apply
                    (p := (ramifiedPolynomial M R).reverse.derivative)
                    PowerSeries.constantCoeff (PowerSeries.C a⁻¹))
        _ = (((ramifiedPolynomial M R).reverse.map
              PowerSeries.constantCoeff).derivative).eval a⁻¹ := by
                rw [Polynomial.derivative_map]
        _ = -(M : F) * (a⁻¹) ^ (R.natDegree - M - 1) := by
          rw [map_constantCoeff_reverse_ramifiedPolynomial M R hdeg]
          exact derivative_eval_reverseResidue_at_inv M R.natDegree
            (R.coeff 0) a hdeg ha hpow
    rw [heq]
    exact mul_ne_zero (neg_ne_zero.mpr hMchar) <|
      pow_ne_zero _ (inv_ne_zero ha)
  have heval : PowerSeries.constantCoeff
      ((monicReverse (ramifiedPolynomial M R) h₀).eval
        (PowerSeries.C a⁻¹)) = 0 := by
    simp only [monicReverse, eval_mul, eval_C, map_mul, hreverseEval, mul_zero]
  have hscale :
      PowerSeries.constantCoeff (↑(h₀.unit⁻¹) : PowerSeries F) ≠ 0 := by
    exact ((h₀.unit⁻¹).isUnit.map PowerSeries.constantCoeff).ne_zero
  have hderiv : PowerSeries.constantCoeff
      ((monicReverse (ramifiedPolynomial M R) h₀).derivative.eval
        (PowerSeries.C a⁻¹)) ≠ 0 := by
    simp only [monicReverse, Polynomial.derivative_mul, Polynomial.derivative_C,
      zero_mul, zero_add, eval_mul, eval_C, map_mul]
    exact mul_ne_zero hscale hreverseDeriv
  exact exists_root_with_constantCoeff_of_monicReverse
    (ramifiedPolynomial M R) h₀ a ha heval hderiv

end PowerSeries

end GMC2ReciprocalSmallRoots

#print axioms GMC2ReciprocalSmallRoots.monic_monicReverse
#print axioms GMC2ReciprocalSmallRoots.exists_unit_reciprocal_root
#print axioms GMC2ReciprocalSmallRoots.exists_root_with_constantCoeff_of_monicReverse
#print axioms GMC2ReciprocalSmallRoots.exists_root_ramifiedPolynomial
