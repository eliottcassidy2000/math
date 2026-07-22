import Mathlib

/-!
# Obstacle (ii) of THM-1550 via Mathlib's Weierstrass preparation

The last gap in the GMC(2)/DvdK proof (THM-2067 route) is THM-1550, whose obstacle (ii)
(death-star-S106/S111) is the construction of the **small-root factor** of
`Φ = x^M − t·R(x)` — the degree-`M` distinguished polynomial `P` (`P ≡ x^M mod t`) with
`Φ = P · h`, `h` a unit.  death-star estimated this as manual monic-M-th-root Hensel + a
`(t)`-adic fixed-point iteration ("the next substantial piece", months-scale).

**This file shows Mathlib already has it.**  View `Φ` as a power series in `x` over the
complete local ring `A = F[[t]]`; its residue image is `x^M ≠ 0`, so
`PowerSeries.exists_isWeierstrassFactorization` produces exactly the distinguished factor `P`
and the unit `h`.  The only missing instance — `IsAdicComplete (maximalIdeal A) A` — is
assembled from `PowerSeries.maximalIdeal_eq_span_X` plus the `(X)`-adic completeness instance
Mathlib already provides (the same instance death-star built for obstacle (i)).

So obstacle (ii) is **not** months of manual Hensel; it is a direct Mathlib appeal, kernel-pure.

Conventions: `A = PowerSeries F = F[[t]]` (inner variable `t = PowerSeries.X`); `Φ ∈ A⟦X⟧`
is a power series in `x = X` (outer).  `Φ = x^M − t·R(x)` with `R : F[X]` embedded as
`R.map (algebraMap F A)`.
-/

open PowerSeries

namespace GMC2DvdKWeierstrass

noncomputable section

variable {F : Type*} [Field F]

/-- The instance that unlocks Weierstrass preparation over `F[[t]]`: `F[[t]]` is complete for
its maximal ideal, because that ideal is `(X)` and `F[[t]]` is `(X)`-adically complete. -/
instance : IsAdicComplete (IsLocalRing.maximalIdeal (PowerSeries F)) (PowerSeries F) := by
  rw [PowerSeries.maximalIdeal_eq_span_X]; infer_instance

/-- `Φ = x^M − t·R(x)` as a power series in `x = X` over `A = F[[t]]`. -/
def Phi (R : Polynomial F) (M : ℕ) : (PowerSeries F)⟦X⟧ :=
  (PowerSeries.X) ^ M
    - (PowerSeries.C (PowerSeries.X : PowerSeries F))
        * ((R.map (algebraMap F (PowerSeries F))) : (PowerSeries F)⟦X⟧)

/-- The residue image of `Φ` (mod `t`) is `x^M`: the `−t·R` term dies because `residue t = 0`. -/
theorem map_residue_Phi (R : Polynomial F) (M : ℕ) :
    PowerSeries.map (IsLocalRing.residue (PowerSeries F)) (Phi R M) = PowerSeries.X ^ M := by
  have hX : IsLocalRing.residue (PowerSeries F) (PowerSeries.X) = 0 := by
    rw [IsLocalRing.residue_eq_zero_iff, PowerSeries.maximalIdeal_eq_span_X]
    exact Ideal.subset_span rfl
  rw [Phi, map_sub, map_pow, PowerSeries.map_X, map_mul, PowerSeries.map_C, hX, map_zero,
    zero_mul, sub_zero]

/-- **Obstacle (ii), kernel-pure.**  `Φ = x^M − t·R(x)` admits a Weierstrass factorization:
`Φ = P · h` where `P` is a distinguished polynomial (monic, `≡ x^M mod t`) and `h` is a unit
power series.  `P` is the small-root factor; `Π = (−1)^M · P.coeff 0` is the small-root product.
Obtained directly from Mathlib's Weierstrass preparation theorem. -/
theorem phi_weierstrass (R : Polynomial F) (M : ℕ) :
    ∃ f h, (Phi R M).IsWeierstrassFactorization f h := by
  apply PowerSeries.exists_isWeierstrassFactorization
  rw [map_residue_Phi]
  exact pow_ne_zero M PowerSeries.X_ne_zero

end

#print axioms phi_weierstrass

end GMC2DvdKWeierstrass
