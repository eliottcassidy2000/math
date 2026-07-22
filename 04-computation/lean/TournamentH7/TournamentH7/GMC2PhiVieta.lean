import Mathlib

/-!
# Vieta valuation-0 input for the THM-2067 wrapper

`GMC2Thm2067Wrapper.thm2067_contradiction` consumes two number-theoretic inputs about the roots of
`Φ(X) = Xᴹ − t·R(X)` over `F(t)`:

* `hS` (THM-1550): the small-root product equals `c·t` (`t`-adic valuation 1) — the deep analytic gap;
* `hΩ` (**Vieta**): the full-root product is a nonzero *constant* `d ∈ F` (`t`-adic valuation 0).

This module supplies the number-theoretic content of the second input.  By Vieta's formula the product
of the roots of `Φ` (over a splitting field) is `(−1)^{deg} · Φ.coeff 0 / Φ.leadingCoeff`, so the point
is that this ratio is *constant in `t`*: `Φ.coeff 0 = −t·r₀` and `Φ.leadingCoeff = −t·lc(R)`, and the
`t` cancels, leaving `r₀/lc(R) ∈ F`.  Hence the root product has `t`-adic valuation `0`.
-/

open Polynomial

namespace GMC2PhiVieta

variable {F : Type*} [Field F]

/-- `Φ = Xᴹ − t·R` over `F(t)` (`t = RatFunc.X`, `R` embedded as constants). -/
noncomputable def Phi (R : F[X]) (M : ℕ) : (RatFunc F)[X] :=
  X ^ M - C (RatFunc.X) * R.map (algebraMap F (RatFunc F))

theorem natDegree_t_Rmap (R : F[X]) :
    (C (RatFunc.X) * R.map (algebraMap F (RatFunc F))).natDegree = R.natDegree := by
  rw [Polynomial.natDegree_C_mul (by simpa using RatFunc.X_ne_zero),
    Polynomial.natDegree_map_eq_of_injective (algebraMap F (RatFunc F)).injective]

/-- With `M < deg R`, the `−t·R` term dominates, so `natDegree Φ = deg R`. -/
theorem natDegree_Phi (R : F[X]) (M : ℕ) (hMd : M < R.natDegree) :
    (Phi R M).natDegree = R.natDegree := by
  unfold Phi
  rw [sub_eq_add_neg, Polynomial.natDegree_add_eq_right_of_natDegree_lt, natDegree_neg,
    natDegree_t_Rmap R]
  rw [natDegree_neg, natDegree_t_Rmap R, Polynomial.natDegree_X_pow]; exact hMd

/-- `Φ.coeff 0 = −t·r₀`. -/
theorem coeff_zero_Phi (R : F[X]) (M : ℕ) (hM : 1 ≤ M) :
    (Phi R M).coeff 0 = - RatFunc.X * (algebraMap F (RatFunc F)) (R.coeff 0) := by
  unfold Phi
  rw [Polynomial.coeff_sub, Polynomial.coeff_X_pow, if_neg (by omega),
    Polynomial.coeff_C_mul, Polynomial.coeff_map]; ring

/-- `Φ.leadingCoeff = −t·lc(R)`. -/
theorem leadingCoeff_Phi (R : F[X]) (M : ℕ) (hMd : M < R.natDegree) :
    (Phi R M).leadingCoeff = - RatFunc.X * (algebraMap F (RatFunc F)) R.leadingCoeff := by
  simp only [Polynomial.leadingCoeff]
  rw [natDegree_Phi R M hMd]
  unfold Phi
  rw [Polynomial.coeff_sub, Polynomial.coeff_X_pow, if_neg (by omega),
    Polynomial.coeff_C_mul, Polynomial.coeff_map]; ring

/-- **The Vieta core (valuation-0 fact).**  `Φ.coeff 0 / Φ.leadingCoeff = r₀/lc(R)` is the image of a
CONSTANT of `F` (the `t` cancels), so the product of the roots of `Φ` (= `±` this ratio) has `t`-adic
valuation `0` — the `hΩ` input of `GMC2Thm2067Wrapper.thm2067_contradiction`. -/
theorem coeff_ratio_Phi_eq_const (R : F[X]) (M : ℕ) (hM : 1 ≤ M) (hR : R ≠ 0)
    (hMd : M < R.natDegree) :
    (Phi R M).coeff 0 / (Phi R M).leadingCoeff
      = (algebraMap F (RatFunc F)) (R.coeff 0 / R.leadingCoeff) := by
  rw [coeff_zero_Phi R M hM, leadingCoeff_Phi R M hMd, map_div₀]
  have hX : (RatFunc.X : RatFunc F) ≠ 0 := RatFunc.X_ne_zero
  have hlc : (algebraMap F (RatFunc F)) R.leadingCoeff ≠ 0 := by
    simp only [Ne, (algebraMap F (RatFunc F)).injective.eq_iff' (map_zero _)]
    exact Polynomial.leadingCoeff_ne_zero.mpr hR
  field_simp

end GMC2PhiVieta

#print axioms GMC2PhiVieta.coeff_ratio_Phi_eq_const
