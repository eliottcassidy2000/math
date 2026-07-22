import Mathlib
import TournamentH7.GMC2DvdKFrame

/-!
# DvdK `hderiv`, leg (c): the coefficient-extraction `xCoeff0(logDeriv Φ)`

Builds on `GMC2DvdKFrame` (death-star's `(LaurentSeries F)⟦t⟧` frame).  This module supplies the
**extraction leg** of the `hderiv` identity: computing `logDeriv Φ` (`= −R/Φ`) in closed form and
reading off its `[x⁰]`.

The key algebra (this section, general over any `CommRing`): the geometric inverse
`(1 − w·X)⁻¹ = ∑ wⁿ Xⁿ`, which is **not in Mathlib** in this form.  Then, in the frame,
`Φ = C(xᴹ)·(1 − X·C(w))` with `w = Rl·x⁻ᴹ`, so `logDeriv Φ = −mk (n ↦ wⁿ⁺¹)`, whose `xCoeff0`
is `−mk (n ↦ (Rlⁿ⁺¹).coeff (M·(n+1)))` — the DvdK moments `D_{m}` (Check A).
-/

open PowerSeries

namespace GMC2DvdKFrameExtraction

/-! ## The geometric inverse `(1 − w·X)⁻¹ = ∑ wⁿ Xⁿ` (general, reusable — Mathlib gap) -/

section Geometric

variable {R : Type*} [CommRing R]

/-- `(1 − C w · X) · (∑ wⁿ Xⁿ) = 1`: the geometric series is the inverse of `1 − w·X`. -/
theorem oneSubCX_mul_mkGeom (w : R) :
    (1 - C w * X) * mk (fun n => w ^ n) = 1 := by
  rw [sub_mul, one_mul, mul_assoc]
  ext n
  rw [map_sub, coeff_C_mul, coeff_mk]
  cases n with
  | zero => simp [coeff_zero_X_mul]
  | succ k =>
    rw [coeff_succ_X_mul, coeff_mk, coeff_one, if_neg (Nat.succ_ne_zero k), pow_succ]
    ring

/-- `1 − C w · X` as a unit of `R⟦X⟧`, with inverse `∑ wⁿ Xⁿ`. -/
noncomputable def unitOneSubCX (w : R) : (PowerSeries R)ˣ :=
  ⟨1 - C w * X, mk (fun n => w ^ n), oneSubCX_mul_mkGeom w,
    by rw [mul_comm]; exact oneSubCX_mul_mkGeom w⟩

theorem isUnit_oneSubCX (w : R) : IsUnit (1 - C w * X) :=
  ⟨unitOneSubCX w, rfl⟩

/-- `Ring.inverse (1 − C w · X) = ∑ wⁿ Xⁿ`. -/
theorem inverse_oneSubCX (w : R) :
    Ring.inverse (1 - C w * X) = mk (fun n => w ^ n) := by
  rw [show (1 - C w * X) = (unitOneSubCX w : PowerSeries R) from rfl, Ring.inverse_unit]
  rfl

end Geometric

end GMC2DvdKFrameExtraction

#print axioms GMC2DvdKFrameExtraction.oneSubCX_mul_mkGeom
#print axioms GMC2DvdKFrameExtraction.inverse_oneSubCX
