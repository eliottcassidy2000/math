import Mathlib

/-!
# The unified `(LaurentSeries F)⟦t⟧` frame for the DvdK log-derivative identity (`hderiv`)

The last open lemma of the whole GMC(2) formalization is `hderiv`: `d_t(h(0,t)) = 0` under the DvdK
vanishing `D_m = 0`, where `h` is the Weierstrass unit of `Φ = xᴹ − t·R` over `F⟦t⟧`.  mac-mini
flagged the obstacle as *"taking `[x⁰]` consistently across two completions"* — `h(0,t)` is `x`-adic
(`h ∈ F⟦t⟧⟦x⟧`), while the generating function `F(t) = [x⁰](xᴹ/Φ)` is `t`-adic.

**This module supplies the resolution: work in the single ring `PowerSeries (LaurentSeries F)`
(`= (F⸨x⸩)⟦t⟧`).**  Because `LaurentSeries F = F⸨x⸩` is a **field** (`IsFractionRing F⟦x⟧ F⸨x⸩`), a
`t`-power-series is a unit as soon as its constant-in-`t` coefficient is nonzero.  So `Φ` (const-`t`
coeff `xᴹ ≠ 0`), the distinguished factor `P` (const-`t` coeff `xᴹ`), and the unit `h` (const-`t`
coeff `1`) are **all units in one ring** — no fraction field, honest division, and `h ∈ F⟦t⟧⟦x⟧`
embeds here (its `x`-support is bounded below).  The two completions collapse to one.

The `[x⁰]` extraction is then the `F⟦t⟧`-additive map `xCoeff0 : (F⸨x⸩)⟦t⟧ →+ F⟦t⟧` reading the `x⁰`
Hahn-coefficient of each `t`-coefficient.  In this frame the `hderiv` identity is
`xCoeff0 (−R/Φ) = xCoeff0 (P_t/P) + xCoeff0 (h_t/h)`, with `xCoeff0 (P_t/P) = 0` (P monic of `x`-degree
`M`, `P_t` of `x`-degree `< M`) and `xCoeff0 (h_t/h) = d_t(h(0,t))/h(0,t)`.
-/

open PowerSeries

namespace GMC2DvdKFrame

variable {F : Type*} [Field F]

/-- **A `t`-power-series over the field `F⸨x⸩` is a unit iff its constant-in-`t` term is nonzero.**
This is the whole point of the frame: in a field, "nonzero constant term" is all it takes. -/
theorem isUnit_iff_constantCoeff_ne_zero (φ : PowerSeries (LaurentSeries F)) :
    IsUnit φ ↔ PowerSeries.constantCoeff (R := LaurentSeries F) φ ≠ 0 := by
  rw [PowerSeries.isUnit_iff_constantCoeff, isUnit_iff_ne_zero]

/-- `Φ = xᴹ − t·R` in the frame (`t = PowerSeries.X`, inner `x` the Laurent-series variable
`single 1 1`, `R` a Laurent series embedded constant-in-`t`). -/
noncomputable def PhiFrame (Rl : LaurentSeries F) (M : ℕ) : PowerSeries (LaurentSeries F) :=
  PowerSeries.C ((HahnSeries.single (1 : ℤ) (1 : F)) ^ M) - PowerSeries.X * PowerSeries.C Rl

theorem constantCoeff_PhiFrame (Rl : LaurentSeries F) (M : ℕ) :
    PowerSeries.constantCoeff (R := LaurentSeries F) (PhiFrame Rl M)
      = (HahnSeries.single (1 : ℤ) (1 : F)) ^ M := by
  simp [PhiFrame]

/-- **`Φ` is a unit in `(LaurentSeries F)⟦t⟧`.**  (`xᴹ = single 1 1 ^ M ≠ 0` in the field `F⸨x⸩`.) -/
theorem isUnit_PhiFrame (Rl : LaurentSeries F) (M : ℕ) : IsUnit (PhiFrame Rl M) := by
  rw [isUnit_iff_constantCoeff_ne_zero, constantCoeff_PhiFrame]
  exact pow_ne_zero M (by simp)

/-- **The `[x⁰]` functional** `(F⸨x⸩)⟦t⟧ →+ F⟦t⟧`: read the `x⁰` Hahn-coefficient of each
`t`-coefficient.  It is additive (`x⁰`-coefficient is `F`-linear on `F⸨x⸩`), the operator that turns
the frame identity into a `t`-series identity. -/
noncomputable def xCoeff0 : PowerSeries (LaurentSeries F) →+ PowerSeries F where
  toFun φ := PowerSeries.mk fun k => (PowerSeries.coeff (R := LaurentSeries F) k φ).coeff 0
  map_zero' := by ext k; simp
  map_add' a b := by
    ext k
    simp only [PowerSeries.coeff_mk, map_add, HahnSeries.coeff_add]

@[simp] theorem coeff_xCoeff0 (φ : PowerSeries (LaurentSeries F)) (k : ℕ) :
    PowerSeries.coeff (R := F) k (xCoeff0 φ)
      = (PowerSeries.coeff (R := LaurentSeries F) k φ).coeff 0 := by
  simp [xCoeff0, PowerSeries.coeff_mk]

/-! ## The logarithmic derivative (general, reusable)

`logDeriv φ = φ'/φ` for a formal power series.  The one structural fact behind the `hderiv` identity
is that on a *product of units* it is additive: `logDeriv (Φ) = logDeriv P + logDeriv h` becomes
`−R/Φ = P_t/P + h_t/h` when `Φ = P·h`.  This is a clean group-hom-like statement (units → additive),
absent from Mathlib, valid over any commutative ring. -/

section LogDeriv

variable {R : Type*} [CommRing R]

/-- The logarithmic derivative `φ'/φ` of a formal power series (via `Ring.inverse`; the honest `φ'·φ⁻¹`
whenever `φ` is a unit). -/
noncomputable def logDeriv (φ : R⟦X⟧) : R⟦X⟧ := PowerSeries.derivativeFun φ * Ring.inverse φ

/-- **The logarithmic derivative is additive on products of units** (the algebraic core of `hderiv`).
`Φ = P·h ⟹ logDeriv Φ = logDeriv P + logDeriv h`, i.e. `−R/Φ = P_t/P + h_t/h`. -/
theorem logDeriv_mul {φ ψ : R⟦X⟧} (hφ : IsUnit φ) (hψ : IsUnit ψ) :
    logDeriv (φ * ψ) = logDeriv φ + logDeriv ψ := by
  have hφi : φ * Ring.inverse φ = 1 := Ring.mul_inverse_cancel φ hφ
  have hψi : ψ * Ring.inverse ψ = 1 := Ring.mul_inverse_cancel ψ hψ
  simp only [logDeriv, PowerSeries.derivativeFun_mul, Ring.mul_inverse_rev]
  linear_combination (PowerSeries.derivativeFun φ * Ring.inverse φ) * hψi
    + (PowerSeries.derivativeFun ψ * Ring.inverse ψ) * hφi

/-- Log-derivative of a unit times the unit recovers the derivative: `(logDeriv φ)·φ = φ'`. -/
theorem logDeriv_mul_self {φ : R⟦X⟧} (hφ : IsUnit φ) :
    logDeriv φ * φ = PowerSeries.derivativeFun φ := by
  rw [logDeriv, mul_assoc, Ring.inverse_mul_cancel φ hφ, mul_one]

end LogDeriv

/-- **The `hderiv` identity assembles for free in the frame.**  For `Φ = φ·ψ` with `φ, ψ` units in
`(LaurentSeries F)⟦t⟧`, the `[x⁰]` of the log-derivative splits additively:
`[x⁰](logDeriv Φ) = [x⁰](logDeriv φ) + [x⁰](logDeriv ψ)` — immediate from `logDeriv_mul` (Stage 2) plus
additivity of `xCoeff0` (Stage 1).  Instantiated at the Weierstrass `Φ = P·h`, the left side is
`[x⁰](−R/Φ) = −∑_{m≥1} D_m t^{m-1}` (the generating function) and the right side is
`[x⁰](P_t/P) + [x⁰](h_t/h) = 0 + d_t(h(0,t))/h(0,t)` — so under `D_m = 0`, `d_t(h(0,t)) = 0 = hderiv`.
The three remaining inputs are all frame-local: `[x⁰](logDeriv P) = 0` (P monic of x-degree M),
`[x⁰](−R/Φ) = −∑ D_m t^{m-1}` (geometric series), and the transpose embedding
`F⟦t⟧⟦x⟧ ↪ (F⸨x⸩)⟦t⟧` carrying the Weierstrass `P, h`. -/
theorem xCoeff0_logDeriv_mul {φ ψ : PowerSeries (LaurentSeries F)}
    (hφ : IsUnit φ) (hψ : IsUnit ψ) :
    xCoeff0 (logDeriv (φ * ψ)) = xCoeff0 (logDeriv φ) + xCoeff0 (logDeriv ψ) := by
  rw [logDeriv_mul hφ hψ, map_add]

end GMC2DvdKFrame

#print axioms GMC2DvdKFrame.isUnit_PhiFrame
#print axioms GMC2DvdKFrame.xCoeff0
#print axioms GMC2DvdKFrame.logDeriv_mul
#print axioms GMC2DvdKFrame.xCoeff0_logDeriv_mul
