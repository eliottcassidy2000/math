import Mathlib

open PowerSeries

namespace ScratchFrame2

variable {F : Type*} [Field F]

-- LaurentSeries F = F⸨X⸩ is a field (fraction field of F⟦X⟧).
example : Field (LaurentSeries F) := inferInstance

-- The Laurent-series variable x = (single 1 1) is nonzero, so x^M is a unit (field).
example (M : ℕ) : IsUnit ((HahnSeries.single (1 : ℤ) (1 : F)) ^ M) := by
  apply IsUnit.pow
  rw [isUnit_iff_ne_zero]
  simp [HahnSeries.single_ne_zero]

/-- `Φ = x^M − t·R` in the UNIFIED frame `(LaurentSeries F)⟦t⟧` (`t = PowerSeries.X`, inner `x` the
Laurent-series variable, `R` embedded constant-in-`t`).  Because `LaurentSeries F` is a FIELD, the
constant-in-`t` coefficient `x^M ≠ 0` is a unit, so **Φ is a unit** — no fraction field needed, and
`h ∈ F⟦t⟧⟦x⟧` embeds here too (its x-support is bounded below), dissolving the two-completion bridge. -/
noncomputable def PhiFrame (Rl : LaurentSeries F) (M : ℕ) : PowerSeries (LaurentSeries F) :=
  PowerSeries.C ((HahnSeries.single (1 : ℤ) (1 : F)) ^ M) - PowerSeries.X * PowerSeries.C Rl

theorem constantCoeff_PhiFrame (Rl : LaurentSeries F) (M : ℕ) :
    PowerSeries.constantCoeff (R := LaurentSeries F) (PhiFrame Rl M)
      = (HahnSeries.single (1 : ℤ) (1 : F)) ^ M := by
  simp [PhiFrame]

/-- **`Φ` is a unit in `(LaurentSeries F)⟦t⟧`.** -/
theorem isUnit_PhiFrame (Rl : LaurentSeries F) (M : ℕ) : IsUnit (PhiFrame Rl M) := by
  rw [PowerSeries.isUnit_iff_constantCoeff, constantCoeff_PhiFrame]
  apply IsUnit.pow
  rw [isUnit_iff_ne_zero]
  simp [HahnSeries.single_ne_zero]

end ScratchFrame2
