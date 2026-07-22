import Mathlib

open PowerSeries LaurentPolynomial

namespace ScratchFrame

variable {F : Type*} [Field F]

/-- `x^M = T M` is a unit in `LaurentPolynomial F`. -/
theorem isUnit_TM (M : ℕ) : IsUnit (LaurentPolynomial.T (M : ℤ) : LaurentPolynomial F) :=
  ⟨⟨LaurentPolynomial.T (M : ℤ), LaurentPolynomial.T (-(M : ℤ)),
      by rw [← LaurentPolynomial.T_add]; simp, by rw [← LaurentPolynomial.T_add]; simp⟩, rfl⟩

/-- `Φ = x^M − t·R` in the unified frame `(LaurentPolynomial F)⟦t⟧` (`t = PowerSeries.X`,
`x^M = T M`, `R` a Laurent polynomial embedded constant-in-`t`). -/
noncomputable def PhiFrame (Rl : LaurentPolynomial F) (M : ℕ) : PowerSeries (LaurentPolynomial F) :=
  PowerSeries.C (LaurentPolynomial.T (M : ℤ) : LaurentPolynomial F)
    - PowerSeries.X * PowerSeries.C Rl

/-- Its constant-in-`t` coefficient is `x^M = T M`. -/
theorem constantCoeff_PhiFrame (Rl : LaurentPolynomial F) (M : ℕ) :
    PowerSeries.constantCoeff (R := LaurentPolynomial F) (PhiFrame Rl M)
      = LaurentPolynomial.T (M : ℤ) := by
  simp [PhiFrame]

/-- **The load-bearing fact: `Φ` is a UNIT in `(LaurentPolynomial F)⟦t⟧`.**  So `R/Φ`, `P_t/P`,
`h_t/h` are honest (no fraction field), and the whole log-derivative identity lives in one ring. -/
theorem isUnit_PhiFrame (Rl : LaurentPolynomial F) (M : ℕ) : IsUnit (PhiFrame Rl M) := by
  rw [PowerSeries.isUnit_iff_constantCoeff, constantCoeff_PhiFrame]
  exact isUnit_TM M

end ScratchFrame
