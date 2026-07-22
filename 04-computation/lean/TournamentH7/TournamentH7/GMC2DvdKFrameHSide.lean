import Mathlib
import TournamentH7.GMC2DvdKFrame

/-!
# The `hderiv` h-side (a): `[x⁰](h_t/h) = g_t/g` on the disk

The last-but-one open input to `hderiv`, discharged via one structural idea (reflection
`hderiv-disk-annulus-split-hside-and-transpose-kps-S128c153`):

**`xCoeff0` is a ring homomorphism on the disk subring `F⟦x⟧⟦t⟧`** (the constant-`x` term of a product of
*power* series is the product of constant terms — false for Laurent series, verified: `1 + t(x+x⁻¹)` breaks
it).  The Weierstrass unit `h` is a genuine power series in `x`, so it lives on the disk, and the h-side is
just **`logDeriv` commuting with a ring homomorphism**.

* `logDeriv_map` (general, reusable — Mathlib gap): `map ψ (logDeriv u) = logDeriv (map ψ u)` for any ring
  hom `ψ` and unit `u`, from `derivativeFun_map` + ring homs preserving unit inverses.
* `xCoeff0_map_ofPowerSeries`: `[x⁰] ∘ ofPowerSeries = constantCoeff`, i.e. on the disk `xCoeff0` is the
  ring hom `map constantCoeff`.
* `xCoeff0_logDeriv_map_ofPowerSeries` = death-star's `hderiv_of_frame` hypothesis `ha`, for
  `hfr = map ofPowerSeries H` (the form the transpose delivers).
-/

open PowerSeries

namespace GMC2DvdKFrameHSide

/-! ## General: `logDeriv` commutes with ring homomorphisms -/

section General

variable {A B : Type*} [CommRing A] [CommRing B] (ψ : A →+* B)

/-- `PowerSeries.map ψ` commutes with the formal derivative `derivativeFun`. -/
theorem derivativeFun_map (f : PowerSeries A) :
    PowerSeries.map ψ (derivativeFun f) = derivativeFun (PowerSeries.map ψ f) := by
  ext n
  rw [coeff_map, coeff_derivativeFun, coeff_derivativeFun, coeff_map, map_mul, map_add,
    map_natCast, map_one]

/-- A ring homomorphism preserves the `Ring.inverse` of a unit. -/
theorem map_ringInverse_unit {u : PowerSeries A} (hu : IsUnit u) :
    PowerSeries.map ψ (Ring.inverse u) = Ring.inverse (PowerSeries.map ψ u) := by
  have hu' : IsUnit (PowerSeries.map ψ u) := hu.map (PowerSeries.map ψ)
  have h1 : PowerSeries.map ψ u * PowerSeries.map ψ (Ring.inverse u) = 1 := by
    rw [← map_mul, Ring.mul_inverse_cancel u hu, map_one]
  calc PowerSeries.map ψ (Ring.inverse u)
      = (Ring.inverse (PowerSeries.map ψ u) * PowerSeries.map ψ u)
          * PowerSeries.map ψ (Ring.inverse u) := by rw [Ring.inverse_mul_cancel _ hu', one_mul]
    _ = Ring.inverse (PowerSeries.map ψ u) * 1 := by rw [mul_assoc, h1]
    _ = Ring.inverse (PowerSeries.map ψ u) := mul_one _

/-- **`logDeriv` commutes with a ring homomorphism** (on units). -/
theorem logDeriv_map {u : PowerSeries A} (hu : IsUnit u) :
    PowerSeries.map ψ (GMC2DvdKFrame.logDeriv u) = GMC2DvdKFrame.logDeriv (PowerSeries.map ψ u) := by
  rw [GMC2DvdKFrame.logDeriv, GMC2DvdKFrame.logDeriv, map_mul, derivativeFun_map,
    map_ringInverse_unit ψ hu]

end General

/-! ## Disk specialization: `xCoeff0` on `F⟦x⟧⟦t⟧` is the ring hom `map constantCoeff`

On the **disk subring** `F⟦x⟧⟦t⟧ ↪ F⸨x⸩⟦t⟧` (image of `map (ofPowerSeries ℤ F)`), `xCoeff0` is the ring
homomorphism `map constantCoeff` — because `[x⁰]` is multiplicative on *power* series (false on Laurent
series, e.g. `1 + t(x+x⁻¹)`).  The Weierstrass unit lands here, so the h-side is `logDeriv_map` twice. -/

section Disk

variable {F : Type*} [Field F]

/-- **`[x⁰]` on the disk subring is `constantCoeff`.** -/
theorem xCoeff0_map_ofPowerSeries (H : PowerSeries (PowerSeries F)) :
    GMC2DvdKFrame.xCoeff0 (PowerSeries.map (HahnSeries.ofPowerSeries ℤ F) H)
      = PowerSeries.map (PowerSeries.constantCoeff (R := F)) H := by
  ext k
  rw [GMC2DvdKFrame.coeff_xCoeff0, coeff_map, coeff_map]
  simpa using HahnSeries.ofPowerSeries_apply_coeff (PowerSeries.coeff k H) 0

/-- **The h-side (a), on the disk.**  For a unit `H : F⟦x⟧⟦t⟧`, with `hfr = map ofPowerSeries H` (the form
the transpose delivers), `xCoeff0(logDeriv hfr) = g_t · g⁻¹` with `g = xCoeff0 hfr`.  Exactly
`GMC2DvdKHderiv.hderiv_of_frame`'s hypothesis `ha`. -/
theorem xCoeff0_logDeriv_map_ofPowerSeries {H : PowerSeries (PowerSeries F)} (hH : IsUnit H) :
    GMC2DvdKFrame.xCoeff0 (GMC2DvdKFrame.logDeriv (PowerSeries.map (HahnSeries.ofPowerSeries ℤ F) H))
      = derivativeFun (GMC2DvdKFrame.xCoeff0 (PowerSeries.map (HahnSeries.ofPowerSeries ℤ F) H))
        * Ring.inverse (GMC2DvdKFrame.xCoeff0 (PowerSeries.map (HahnSeries.ofPowerSeries ℤ F) H)) := by
  simp only [xCoeff0_map_ofPowerSeries]
  rw [← logDeriv_map (HahnSeries.ofPowerSeries ℤ F) hH, xCoeff0_map_ofPowerSeries,
    logDeriv_map (PowerSeries.constantCoeff (R := F)) hH, GMC2DvdKFrame.logDeriv]

end Disk

end GMC2DvdKFrameHSide

#print axioms GMC2DvdKFrameHSide.logDeriv_map
#print axioms GMC2DvdKFrameHSide.xCoeff0_map_ofPowerSeries
#print axioms GMC2DvdKFrameHSide.xCoeff0_logDeriv_map_ofPowerSeries
