import Mathlib
import TournamentH7.GMC2DvdKWeierstrass

/-!
# The bridge divisibility: the Weierstrass distinguished factor divides `Φ` in the polynomial ring

The frame bridge needs `smallRootFactor R M ∣ Φ` as **polynomials** over `PowerSeries F` (so that,
mapped to the algebraic closure of `LaurentSeries F`, the small-root packet is a sub-multiset of `Φ`'s
roots).  This is *not* the trivial power-series divisibility over the fraction field (where the
distinguished factor becomes a unit); it is the genuine polynomial factorization, obtained from
Mathlib's Weierstrass **division uniqueness** (`IsWeierstrassDivisorAt.eq_of_mul_add_eq_mul_add`) applied
to the two divisions of `Φ` by the monic `smallRootFactor`: the Weierstrass factorization
`Φ = ↑P·h` (remainder `0`) and the ordinary polynomial division (remainder `Φ %ₘ P`).  Uniqueness forces
`Φ %ₘ P = 0`.  No valuation on the algebraic closure.
-/

open Polynomial

namespace GMC2FrameBridgeDvd

variable {F : Type*} [Field F]

/-- The polynomial (over `PowerSeries F`) whose coercion into `(PowerSeries F)⟦X⟧` is mac-mini's
Weierstrass `Φ = xᴹ − t·R`. -/
noncomputable def PhiPoly (R : Polynomial F) (M : ℕ) : (PowerSeries F)[X] :=
  X ^ M - C (PowerSeries.X) * R.map (algebraMap F (PowerSeries F))

/-- The coercion of `PhiPoly` is the Weierstrass power series `Φ`. -/
theorem coe_PhiPoly (R : Polynomial F) (M : ℕ) :
    ((PhiPoly R M : Polynomial (PowerSeries F)) : PowerSeries (PowerSeries F))
      = GMC2DvdKWeierstrass.Phi R M := by
  rw [PhiPoly, GMC2DvdKWeierstrass.Phi]
  push_cast
  ring

end GMC2FrameBridgeDvd
