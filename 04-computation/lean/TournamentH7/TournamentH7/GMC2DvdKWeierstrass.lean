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

/-- The residue image of `Φ` is nonzero (`= x^M`), the hypothesis for Weierstrass preparation. -/
theorem phi_residue_ne_zero (R : Polynomial F) (M : ℕ) :
    PowerSeries.map (IsLocalRing.residue (PowerSeries F)) (Phi R M) ≠ 0 := by
  rw [map_residue_Phi]; exact pow_ne_zero M PowerSeries.X_ne_zero

/-- **Obstacle (ii), kernel-pure.**  `Φ = x^M − t·R(x)` admits a Weierstrass factorization:
`Φ = P · h` where `P` is a distinguished polynomial (monic, `≡ x^M mod t`) and `h` is a unit
power series.  `P` is the small-root factor; `Π = (−1)^M · P.coeff 0` is the small-root product.
Obtained directly from Mathlib's Weierstrass preparation theorem. -/
theorem phi_weierstrass (R : Polynomial F) (M : ℕ) :
    ∃ f h, (Phi R M).IsWeierstrassFactorization f h :=
  PowerSeries.exists_isWeierstrassFactorization (phi_residue_ne_zero R M)

/-- The **small-root factor** `P` of `Φ`: the distinguished polynomial in its Weierstrass
factorization.  Monic, `P ≡ x^M mod t`, and `Φ = P · (unit)`. -/
def smallRootFactor (R : Polynomial F) (M : ℕ) : Polynomial (PowerSeries F) :=
  PowerSeries.weierstrassDistinguished (Phi R M) (phi_residue_ne_zero R M)

/-- `P` has degree **exactly `M`** — there are exactly `M` small roots, so
`Π = (−1)^M · P.coeff 0` has `t`-valuation `1`. -/
theorem smallRootFactor_natDegree (R : Polynomial F) (M : ℕ) :
    (smallRootFactor R M).natDegree = M := by
  have H := PowerSeries.isWeierstrassFactorization_weierstrassDistinguished_weierstrassUnit
    (phi_residue_ne_zero R M)
  rw [smallRootFactor, H.natDegree_eq_toNat_order_map, map_residue_Phi, PowerSeries.order_X_pow,
    ENat.toNat_coe]

/-- `P` is monic. -/
theorem smallRootFactor_monic (R : Polynomial F) (M : ℕ) : (smallRootFactor R M).Monic :=
  (PowerSeries.isDistinguishedAt_weierstrassDistinguished (phi_residue_ne_zero R M)).monic

/-- `Φ = P · h` with `h` a unit: the small-root factor divides `Φ`. -/
theorem phi_eq_smallRootFactor_mul (R : Polynomial F) (M : ℕ) :
    ∃ h, IsUnit h ∧ Phi R M = (smallRootFactor R M : (PowerSeries F)⟦X⟧) * h := by
  have H := PowerSeries.isWeierstrassFactorization_weierstrassDistinguished_weierstrassUnit
    (phi_residue_ne_zero R M)
  exact ⟨_, H.isUnit, H.eq_mul⟩

/-- The `x`-constant term of `Φ` is `−t·r₀` (`r₀ = R.coeff 0`, the lowest coefficient). -/
theorem constantCoeff_Phi (R : Polynomial F) (M : ℕ) (hM : 1 ≤ M) :
    PowerSeries.constantCoeff (R := PowerSeries F) (Phi R M)
      = - PowerSeries.X * (algebraMap F (PowerSeries F)) (R.coeff 0) := by
  have hcR : PowerSeries.constantCoeff (R := PowerSeries F)
      ((R.map (algebraMap F (PowerSeries F))) : (PowerSeries F)⟦X⟧)
      = (algebraMap F (PowerSeries F)) (R.coeff 0) := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff, Polynomial.coeff_coe, Polynomial.coeff_map]
  rw [Phi, map_sub, map_pow, map_mul, PowerSeries.constantCoeff_X,
    zero_pow (by omega : M ≠ 0), PowerSeries.constantCoeff_C, hcR]
  ring

/-- **The multiplicative THM-1550 crux reduces to a single scalar identity.**
From `Φ = P·h` (Weierstrass), taking the `x`-constant term gives
`P.coeff 0 · h(0) = −t·r₀`, where `h(0) := constantCoeff h` is a unit of `F[[t]]` (`= 1 mod t`).
Hence the small-root product `Π = (−1)^M·P.coeff 0` satisfies `Π · h(0) = c·t` with
`c = (−1)^{M+1} r₀`, so **`Π = c·t ⟺ h(0,t) = 1`**.  This isolates the sole remaining analytic
input as exactly `h(0,t) = 1` under `D_m = 0` (equivalently `h(0,t) = exp(−∑ D_m tᵐ/m)`). -/
theorem coeff_zero_smallRootFactor_mul_unit (R : Polynomial F) (M : ℕ) (hM : 1 ≤ M) :
    (smallRootFactor R M).coeff 0
      * PowerSeries.constantCoeff (R := PowerSeries F)
          (PowerSeries.weierstrassUnit (Phi R M) (phi_residue_ne_zero R M))
      = - PowerSeries.X * (algebraMap F (PowerSeries F)) (R.coeff 0) := by
  have H := PowerSeries.isWeierstrassFactorization_weierstrassDistinguished_weierstrassUnit
    (phi_residue_ne_zero R M)
  have key := congrArg (PowerSeries.constantCoeff (R := PowerSeries F)) H.eq_mul
  rw [map_mul, ← PowerSeries.coeff_zero_eq_constantCoeff, Polynomial.coeff_coe,
    PowerSeries.coeff_zero_eq_constantCoeff, constantCoeff_Phi R M hM] at key
  exact key.symm

/-- Reducing mod `t` (the ring hom `constantCoeff : F[[t]] → F`) sends `Φ` to `x^M`. -/
theorem map_constantCoeff_Phi (R : Polynomial F) (M : ℕ) (hM : 1 ≤ M) :
    PowerSeries.map (PowerSeries.constantCoeff (R := F)) (Phi R M) = PowerSeries.X ^ M := by
  rw [Phi, map_sub, map_pow, PowerSeries.map_X, map_mul, PowerSeries.map_C,
    PowerSeries.constantCoeff_X, map_zero, zero_mul, sub_zero]

/-- The distinguished small-root factor reduces mod `t` to `x^M`. -/
theorem smallRootFactor_map_constantCoeff (R : Polynomial F) (M : ℕ) :
    (smallRootFactor R M).map (PowerSeries.constantCoeff (R := F)) = (Polynomial.X : Polynomial F) ^ M := by
  have hd := PowerSeries.isDistinguishedAt_weierstrassDistinguished (phi_residue_ne_zero R M)
  have hdeg : (smallRootFactor R M).natDegree = M := smallRootFactor_natDegree R M
  ext k
  rw [Polynomial.coeff_map, Polynomial.coeff_X_pow]
  by_cases hk : k = M
  · rw [if_pos hk, hk]
    have h1 : (smallRootFactor R M).coeff M = 1 := by
      have hmon := (smallRootFactor_monic R M).coeff_natDegree
      rwa [hdeg] at hmon
    rw [h1, map_one]
  · rw [if_neg hk]
    rcases lt_or_gt_of_ne hk with hlt | hgt
    · have hmem : (smallRootFactor R M).coeff k ∈ IsLocalRing.maximalIdeal (PowerSeries F) :=
        hd.toIsWeaklyEisensteinAt.mem (lt_of_lt_of_eq hlt hdeg.symm)
      rw [PowerSeries.maximalIdeal_eq_span_X, Ideal.mem_span_singleton] at hmem
      obtain ⟨c, hc⟩ := hmem
      rw [hc, map_mul, PowerSeries.constantCoeff_X, zero_mul]
    · rw [Polynomial.coeff_eq_zero_of_natDegree_lt (lt_of_eq_of_lt hdeg hgt), map_zero]

/-- **`hconst` for kps's char-0 closing: `h(0,0) = 1`.**  The Weierstrass unit's `x`-constant
term `h(0,t)` has `t`-constant term `1`, because reducing `Φ = P·h` mod `t` gives `x^M = x^M·(h mod t)`,
so `h ≡ 1 mod t`. -/
theorem constantCoeff_constantCoeff_weierstrassUnit (R : Polynomial F) (M : ℕ) (hM : 1 ≤ M) :
    PowerSeries.constantCoeff (R := F)
        (PowerSeries.constantCoeff (R := PowerSeries F)
          (PowerSeries.weierstrassUnit (Phi R M) (phi_residue_ne_zero R M))) = 1 := by
  have H := PowerSeries.isWeierstrassFactorization_weierstrassDistinguished_weierstrassUnit
    (phi_residue_ne_zero R M)
  -- `map(constantCoeff) ↑P = x^M`
  have hPm : PowerSeries.map (PowerSeries.constantCoeff (R := F))
      (↑(smallRootFactor R M) : (PowerSeries F)⟦X⟧) = (PowerSeries.X : F⟦X⟧) ^ M := by
    rw [← Polynomial.polynomial_map_coe, smallRootFactor_map_constantCoeff, Polynomial.coe_pow,
      Polynomial.coe_X]
  -- map(constantCoeff) h = 1
  have hh : PowerSeries.map (PowerSeries.constantCoeff (R := F))
      (PowerSeries.weierstrassUnit (Phi R M) (phi_residue_ne_zero R M)) = 1 := by
    have key := congrArg (PowerSeries.map (PowerSeries.constantCoeff (R := F))) H.eq_mul
    rw [map_mul, map_constantCoeff_Phi R M hM,
      show (↑(PowerSeries.weierstrassDistinguished (Phi R M) (phi_residue_ne_zero R M))
        : (PowerSeries F)⟦X⟧) = ↑(smallRootFactor R M) from rfl, hPm] at key
    exact (mul_left_cancel₀ (pow_ne_zero M PowerSeries.X_ne_zero)
      (by rw [mul_one]; exact key)).symm
  -- take `coeff 0` (in x) of `map(constantCoeff) h = 1`
  have e := congrArg (PowerSeries.coeff (R := F) 0) hh
  simpa only [PowerSeries.coeff_map, PowerSeries.coeff_zero_eq_constantCoeff, map_one] using e

end

#print axioms phi_weierstrass
#print axioms smallRootFactor_natDegree
#print axioms coeff_zero_smallRootFactor_mul_unit
#print axioms constantCoeff_constantCoeff_weierstrassUnit

end GMC2DvdKWeierstrass
