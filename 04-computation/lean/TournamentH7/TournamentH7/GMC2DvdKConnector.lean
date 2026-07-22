import TournamentH7.GMC2DvdKTranspose
import TournamentH7.GMC2DvdKWeierstrass
import TournamentH7.GMC2DvdKFrame

/-!
# The connector `φ(Φ_weier) = PhiFrame` — the transpose meets the frame

The transpose `φ` (GMC2DvdKTranspose) carries the Weierstrass world `F⟦t⟧⟦x⟧` into death-star's frame
`(F⸨x⸩)⟦t⟧`.  This module proves it carries `Φ_weier = xᴹ − t·R` to `PhiFrame`, so that `φ(Φ = P·h)`
gives the frame factorization `PhiFrame = φ(P)·φ(h)` feeding `GMC2DvdKHderiv.hderiv_of_frame`.
-/

open PowerSeries GMC2DvdKTranspose GMC2DvdKFrame

namespace GMC2DvdKConnector

variable {F : Type*} [Field F]

theorem phi_R_map (R : Polynomial F) :
    phi ((R.map (algebraMap F (PowerSeries F)) : PowerSeries (PowerSeries F)))
      = PowerSeries.C (Polynomial.aeval (HahnSeries.single (1 : ℤ) (1 : F)) R) := by
  have key : (phi.comp (Polynomial.coeToPowerSeries.ringHom)).comp
        (Polynomial.mapRingHom (algebraMap F (PowerSeries F)))
      = (PowerSeries.C (R := LaurentSeries F)).comp
          (Polynomial.aeval (HahnSeries.single (1 : ℤ) (1 : F))).toRingHom := by
    apply Polynomial.ringHom_ext
    · intro a
      simp only [RingHom.comp_apply, Polynomial.coe_mapRingHom, Polynomial.map_C,
        Polynomial.coeToPowerSeries.ringHom_apply, Polynomial.coe_C, AlgHom.toRingHom_eq_coe,
        RingHom.coe_coe, Polynomial.aeval_C]
      rw [show (algebraMap F (PowerSeries F)) a = PowerSeries.C a from rfl, phi_C_C,
        LaurentSeries.algebraMap_apply]
    · simp only [RingHom.comp_apply, Polynomial.coe_mapRingHom, Polynomial.map_X,
        Polynomial.coeToPowerSeries.ringHom_apply, Polynomial.coe_X, AlgHom.toRingHom_eq_coe,
        RingHom.coe_coe, Polynomial.aeval_X, phi_X]
  have h := DFunLike.congr_fun key R
  simpa only [RingHom.comp_apply, Polynomial.coe_mapRingHom,
    Polynomial.coeToPowerSeries.ringHom_apply, AlgHom.toRingHom_eq_coe, RingHom.coe_coe] using h

/-- **The connector: `φ(Φ_weier) = PhiFrame`.** -/
theorem phi_Phi (R : Polynomial F) (M : ℕ) :
    phi (GMC2DvdKWeierstrass.Phi R M)
      = PhiFrame (Polynomial.aeval (HahnSeries.single (1 : ℤ) (1 : F)) R) M := by
  rw [GMC2DvdKWeierstrass.Phi, PhiFrame, map_sub, map_pow, phi_X, map_mul, phi_C_X, phi_R_map,
    ← map_pow, HahnSeries.single_pow, one_pow]

end GMC2DvdKConnector

#print axioms GMC2DvdKConnector.phi_Phi
