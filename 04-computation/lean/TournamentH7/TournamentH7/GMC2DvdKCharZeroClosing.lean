import Mathlib

/-!
# A power series with vanishing formal derivative is constant (characteristic zero)

Over a characteristic-zero field, `PowerSeries.derivativeFun f = 0` forces `f = C (constantCoeff f)`.
Mathlib has the forward direction (`PowerSeries.derivativeFun_C`, `derivativeFun_one`) but **not** this
converse; `coeff_eq_zero_of_derivativeFun_eq_zero` / `eq_C_of_derivativeFun_eq_zero` supply it.

## Role in GMC(2)/DvdK

This is the **char-0 back half** of the exp/log-free route to the last DvdK identity
`h(0,t) = exp(-∑ D_m tᵐ/m)`.  Viewing `Φ = Xᴹ − t·R` in `(F⟦t⟧)⟦X⟧`, the Weierstrass factorization
`Φ = P·h` (mac-mini) and the `[X⁰]`-Laurent log-derivative identity give
`d_t(h.coeff₀) / h.coeff₀ = -∑_{m≥1} D_m t^{m-1}` in `F⟦t⟧`.  Under the DvdK vanishing hypothesis
`D_m = 0 ∀ m ≥ 1` (boxeph's `generatingFunction_eq_one` supplies `F(t)=1`), the right side is `0`, so
`d_t(h.coeff₀) = 0`.  `eq_one_of_derivativeFun_eq_zero` then turns this — with `h.coeff₀(0) = 1` — into
`h.coeff₀ = 1` **without any exp/log or Puiseux**, purely by "zero derivative ⟹ constant in char 0".
`factorCoeff0_eq_of_unit_eq_one` performs the final trivial substitution
`P.coeff 0 = [X⁰]Φ = -(t·r₀)`, giving `Π = (-1)ᴹ P.coeff 0 = c·t`.

The module is **frame-agnostic**: it consumes the log-derivative output (`derivativeFun g = 0`) as a
hypothesis, so it composes with the `[X⁰]`-Laurent frame the moment that lands, and it touches no
Weierstrass internals.  Kernel-pure (`#print axioms` below).
-/

open PowerSeries

namespace GMC2DvdKCharZeroClosing

variable {F : Type*} [Field F] [CharZero F]

/-- **Core.** Over a characteristic-zero field, if the formal derivative `derivativeFun f` vanishes
then every positive-index coefficient of `f` is zero.  (`coeff_derivativeFun` gives
`coeff n (derivativeFun f) = coeff (n+1) f · (n+1)`; char zero makes `(n+1)` a nonzero factor.) -/
theorem coeff_eq_zero_of_derivativeFun_eq_zero {f : PowerSeries F}
    (hf : derivativeFun f = 0) {n : ℕ} (hn : 1 ≤ n) : coeff n f = 0 := by
  obtain ⟨k, rfl⟩ : ∃ k, n = k + 1 := ⟨n - 1, by omega⟩
  have hk := coeff_derivativeFun f k
  rw [hf] at hk
  simp only [map_zero] at hk
  -- hk : 0 = coeff (k+1) f * (↑k + 1)
  have hne : (k + 1 : F) ≠ 0 := by exact_mod_cast Nat.succ_ne_zero k
  exact (mul_eq_zero.mp hk.symm).resolve_right hne

/-- Over a characteristic-zero field, a power series with vanishing formal derivative equals the
constant series at its constant term.  (Converse of `PowerSeries.derivativeFun_C`; not in Mathlib.) -/
theorem eq_C_of_derivativeFun_eq_zero {f : PowerSeries F} (hf : derivativeFun f = 0) :
    f = C (constantCoeff f) := by
  ext n
  rcases Nat.eq_zero_or_pos n with hn | hn
  · subst hn
    rw [coeff_zero_eq_constantCoeff, constantCoeff_C]
  · rw [coeff_eq_zero_of_derivativeFun_eq_zero hf hn, coeff_C, if_neg (by omega)]

/-- **Assembly wrapper (exp/log-free DvdK closing).**  The constant-in-`X` coefficient
`g = h.coeff₀ : F⟦t⟧` of the Weierstrass unit: if its `t`-derivative vanishes (the log-derivative
identity under `D_m = 0`) and its value at `t = 0` is `1`, then `g = 1` — no exp/log, no Puiseux. -/
theorem eq_one_of_derivativeFun_eq_zero {g : PowerSeries F}
    (hg0 : constantCoeff g = 1) (hderiv : derivativeFun g = 0) : g = 1 := by
  rw [eq_C_of_derivativeFun_eq_zero hderiv, hg0, map_one]

omit [CharZero F] in
/-- **Final closing step both DvdK routes consume.**  With the Weierstrass unit's constant-in-`X`
coefficient collapsed to `1`, the small-root factor's constant-in-`X` coefficient `P₀` equals the
`[X⁰]` of `Φ = Xᴹ − t·R`.  Stated frame-agnostically: from `P₀ · g = rhs` and `g = 1`, get
`P₀ = rhs`.  (With `rhs = -(t·r₀)` this is `P.coeff 0 = -(t·r₀)`, hence `Π = (-1)ᴹ P₀ = c·t`.) -/
theorem factorCoeff0_eq_of_unit_eq_one {P₀ g rhs : PowerSeries F}
    (hg : g = 1) (hfact : P₀ * g = rhs) : P₀ = rhs := by
  rwa [hg, mul_one] at hfact

end GMC2DvdKCharZeroClosing

#print axioms GMC2DvdKCharZeroClosing.coeff_eq_zero_of_derivativeFun_eq_zero
#print axioms GMC2DvdKCharZeroClosing.eq_C_of_derivativeFun_eq_zero
#print axioms GMC2DvdKCharZeroClosing.eq_one_of_derivativeFun_eq_zero
#print axioms GMC2DvdKCharZeroClosing.factorCoeff0_eq_of_unit_eq_one
