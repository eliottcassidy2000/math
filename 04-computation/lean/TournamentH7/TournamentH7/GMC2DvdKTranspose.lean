import Mathlib

/-!
# The transpose hom `φ : (PowerSeries F)⟦X⟧ → PowerSeries (LaurentSeries F)` (the shared GMC(2) glue)

The sole remaining glue of the whole GMC(2)/hderiv formalization is mapping the Weierstrass factors
`P, h` of `Φ = xᴹ − t·R` (living in `F⟦t⟧⟦x⟧`) into death-star's frame `(F⸨x⸩)⟦t⟧`.  This module
builds that ring hom from scratch (no Mathlib nested-power-series curry iso exists):

* **`tau`** — the swap of nested power-series variables `F⟦t⟧⟦x⟧ → F⟦x⟧⟦t⟧`, `(τf)`'s `(t^a,x^b)`
  coefficient is `f`'s `(x^b,t^a)` coefficient.  It is a ring hom; multiplicativity is the double-sum
  reorder `Finset.sum_comm` on the two convolution antidiagonals (`tau_mul`).
* **`phi := map (HahnSeries.ofPowerSeries) ∘ tauHom`** — post-composing the coefficient inclusion
  `F⟦x⟧ ↪ F⸨x⸩` lands in the frame.  `phi X = C(single 1 1)` = the frame's `x` (`phi_X`).

Because `phi` sends the `x`-power-series `h` into the `x`-support-`≥0` part of the frame with the right
`x⁰`-coefficient, it supplies exactly the `Φ = P·h` factorization and `xCoeff0(h)=unitCoeff0` that
`GMC2DvdKHderiv.hderiv_of_frame`'s `(a)` needs.
-/

open PowerSeries

namespace GMC2DvdKTranspose

variable {F : Type*} [Field F]

/-- The swap of nested power-series variables `τ : F⟦t⟧⟦X⟧ → F⟦X⟧⟦t⟧`:
`(τ f)`'s `(t^a, X^b)` coefficient is `f`'s `(X^b, t^a)` coefficient. -/
noncomputable def tau (f : PowerSeries (PowerSeries F)) : PowerSeries (PowerSeries F) :=
  PowerSeries.mk fun k => PowerSeries.mk fun n =>
    PowerSeries.coeff (R := F) k (PowerSeries.coeff (R := PowerSeries F) n f)

@[simp] theorem coeff_coeff_tau (f : PowerSeries (PowerSeries F)) (a b : ℕ) :
    PowerSeries.coeff (R := F) b (PowerSeries.coeff (R := PowerSeries F) a (tau f))
      = PowerSeries.coeff (R := F) a (PowerSeries.coeff (R := PowerSeries F) b f) := by
  simp [tau, PowerSeries.coeff_mk]

theorem tau_add (f g : PowerSeries (PowerSeries F)) : tau (f + g) = tau f + tau g := by
  refine PowerSeries.ext fun a => PowerSeries.ext fun b => ?_
  simp [coeff_coeff_tau]

theorem tau_one : tau (1 : PowerSeries (PowerSeries F)) = 1 := by
  refine PowerSeries.ext fun a => PowerSeries.ext fun b => ?_
  rw [coeff_coeff_tau]
  by_cases ha : a = 0 <;> by_cases hb : b = 0 <;>
    simp [PowerSeries.coeff_one, PowerSeries.coeff_zero_eq_constantCoeff, ha, hb]

theorem tau_mul (f g : PowerSeries (PowerSeries F)) : tau (f * g) = tau f * tau g := by
  refine PowerSeries.ext fun a => PowerSeries.ext fun b => ?_
  rw [coeff_coeff_tau]
  simp only [PowerSeries.coeff_mul, map_sum, coeff_coeff_tau]
  exact Finset.sum_comm


theorem tau_zero : tau (0 : PowerSeries (PowerSeries F)) = 0 := by
  refine PowerSeries.ext fun a => PowerSeries.ext fun b => ?_
  simp [coeff_coeff_tau]

/-- The swap bundled as a ring hom. -/
noncomputable def tauHom : PowerSeries (PowerSeries F) →+* PowerSeries (PowerSeries F) where
  toFun := tau
  map_one' := tau_one
  map_mul' := tau_mul
  map_zero' := tau_zero
  map_add' := tau_add

/-- The transpose `φ : (PowerSeries F)⟦X⟧ → PowerSeries (LaurentSeries F)`. -/
noncomputable def phi : PowerSeries (PowerSeries F) →+* PowerSeries (LaurentSeries F) :=
  (PowerSeries.map (HahnSeries.ofPowerSeries ℤ F)).comp tauHom

theorem tau_X : tau (PowerSeries.X : PowerSeries (PowerSeries F))
    = PowerSeries.C (PowerSeries.X : PowerSeries F) := by
  refine PowerSeries.ext fun a => PowerSeries.ext fun b => ?_
  rw [coeff_coeff_tau]
  by_cases hb : b = 1 <;> by_cases ha : a = 0 <;>
    simp [PowerSeries.coeff_X, PowerSeries.coeff_C, PowerSeries.coeff_one, ha, hb]

theorem phi_X : phi (PowerSeries.X : PowerSeries (PowerSeries F))
    = PowerSeries.C (HahnSeries.single (1 : ℤ) (1 : F)) := by
  rw [phi, RingHom.comp_apply, show tauHom (PowerSeries.X : PowerSeries (PowerSeries F))
    = PowerSeries.C (PowerSeries.X : PowerSeries F) from tau_X, PowerSeries.map_C,
    HahnSeries.ofPowerSeries_X]

end GMC2DvdKTranspose

#print axioms GMC2DvdKTranspose.phi
#print axioms GMC2DvdKTranspose.phi_X
