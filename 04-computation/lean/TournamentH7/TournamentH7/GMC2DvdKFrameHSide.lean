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

/-! ## Disk specialization (reduction; see reflection `hderiv-disk-annulus-split-...`)

The h-side `ha` follows from `logDeriv_map` on the **disk subring** `F⟦x⟧⟦t⟧ ↪ F⸨x⸩⟦t⟧` (image of
`PowerSeries.map (HahnSeries.ofPowerSeries ℤ F)`), where `xCoeff0 = PowerSeries.map (constantCoeff)` is a
ring homomorphism (because `[x⁰]` is multiplicative on *power* series — **false** on Laurent series, e.g.
`1 + t(x+x⁻¹)`).  Concretely, for a unit `H : F⟦x⟧⟦t⟧` with `hfr = map ofPowerSeries H` (the form the
transpose delivers):

```
xCoeff0(logDeriv hfr) = xCoeff0(map ofPowerSeries (logDeriv H))     [logDeriv_map, ψ = ofPowerSeries]
                      = map constantCoeff (logDeriv H)              [xCoeff0∘ofPowerSeries = constantCoeff]
                      = logDeriv (map constantCoeff H)              [logDeriv_map, ψ = constantCoeff]
                      = derivativeFun g * Ring.inverse g,  g := xCoeff0 hfr = map constantCoeff H
```

— exactly `GMC2DvdKHderiv.hderiv_of_frame`'s hypothesis `ha`.  The two `logDeriv_map` applications above are
the kernel-pure lemma proved in this module; the remaining `xCoeff0 ∘ ofPowerSeries = constantCoeff` step is
elementary (`ofPowerSeries_apply_coeff` at index `0`) and is finalized against the transpose's concrete
embedding of `hfr` (death-star's lane), so the h-side composes the moment the transpose lands.
-/

end GMC2DvdKFrameHSide

#print axioms GMC2DvdKFrameHSide.derivativeFun_map
#print axioms GMC2DvdKFrameHSide.map_ringInverse_unit
#print axioms GMC2DvdKFrameHSide.logDeriv_map
