import Mathlib

/-!
# `F⟦X⟧` is a Henselian local ring — foundation for the THM-1550 unramified-Hensel lift

The DvdK residual of the GMC(2) formalization (THM-2067) reduces to THM-1550: the small-root
product `Π(t) = c·t` of `Φ(X) = X^M - t·R(X)`.  death-star-S106's route substitutes `X = sZ`,
`t = sᴹ`, so `Φ(sZ) = sᴹ·ψ(Z)` with `ψ(Z) = Zᴹ - R(sZ)` over `F⟦s⟧`, and lifts the small factor
by Hensel over the complete local ring `F⟦s⟧`.  This module supplies the missing foundation
(boxeph-S234 obstacle (i)): `F⟦X⟧` is a Henselian local ring.  It is local with maximal ideal
`(X)` and `(X)`-adically complete, hence Henselian at its maximal ideal.
-/

open PowerSeries

namespace GMC2Henselian

variable {F : Type*} [Field F]

/-- The Henselian-at-`(X)` structure on `F⟦X⟧` is free from adic completeness. -/
example : HenselianRing (PowerSeries F) (Ideal.span {(X : PowerSeries F)}) := inferInstance

/-- **`F⟦X⟧` is a Henselian local ring.**  (boxeph-S234 obstacle (i) for the THM-1550 Hensel lift.) -/
instance powerSeries_henselianLocalRing : HenselianLocalRing (PowerSeries F) where
  is_henselian := by
    intro f hf a₀ hfa hunit
    have hHR : HenselianRing (PowerSeries F) (IsLocalRing.maximalIdeal (PowerSeries F)) := by
      rw [maximalIdeal_eq_span_X]; infer_instance
    exact hHR.is_henselian f hf a₀ hfa (hunit.map _)

end GMC2Henselian

#print axioms GMC2Henselian.powerSeries_henselianLocalRing
