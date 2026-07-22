import Mathlib

/-!
# The transpose embedding `F⟦t⟧⟦x⟧ ↪ (F⸨x⸩)⟦t⟧` (the sole remaining glue for `hderiv`)

The Weierstrass factorization `Φ = P·h` lives in `PowerSeries (PowerSeries F)` (outer `x`, inner `t`);
death-star's `hderiv` frame (`GMC2DvdKFrame`) is `PowerSeries (LaurentSeries F)` (outer `t`, inner
`x`-Laurent).  `hderiv_of_frame` needs `PhiFrame = Pfr · hfr` with `Pfr` monic of `x`-degree `M`
(so mac-mini's `(c)` `xCoeff0_logDeriv_eq_zero_of_monic` applies) and `xCoeff0 hfr = unitCoeff0`.

This module starts the swap-and-embed ring hom `psi` that transports `Φ = P·h` to the frame.
`psi φ` has `t`-coefficient `m` equal to the Laurent series `∑_{k≥0} (φ.coeff k).coeff m · xᵏ`.

**Completion plan (fully worked out; the swap ring hom has no Mathlib shortcut — `finSuccEquiv`
for iterated `PowerSeries` is absent, so `psi` is built coefficient-wise):**
1. `coeff_psi`: `((coeff m (psi φ))).coeff k = coeff m (coeff k.toNat φ)` for `k ≥ 0`, `= 0` else
   (`HahnSeries.ofPowerSeries_apply_coeff` for `k = ↑n`; support `⊆ ℕ` for `k < 0`).
2. `map_add'`,`map_zero'`: coefficient-wise linearity.  `map_one'`: `1` has only `(k,m)=(0,0)` nonzero.
3. `map_mul'` (the crux, TRUE by symmetry): both `psi (φχ)` and `psi φ · psi χ` have
   `x^k`-coeff of `t`-order `m` equal to `∑_{i+j=k} ∑_{p+q=m} (φ.coeff i).coeff p · (χ.coeff j).coeff q`
   — the same double convolution, reindexed by `Finset.sum_comm` / nested antidiagonals.
4. `psi_Phi`: `psi (GMC2DvdKWeierstrass.Phi R M) = GMC2DvdKFrame.PhiFrame (ofPowerSeries.. R) M`
   (compute `psi` on `X^M` and on `C X * R.map C`).
5. Then `Pfr := psi (smallRootFactor R M)` is monic of `x`-degree `M` (from `smallRootFactor_natDegree`
   /`_monic`), so `(c)` gives `xCoeff0(logDeriv Pfr)=0`; `hfr := psi (weierstrassUnit ..)` has
   `xCoeff0 hfr = constantCoeff (weierstrassUnit ..) = unitCoeff0` (the `x⁰`-part of the swap is the
   inner constant coeff); `PhiFrame = Pfr·hfr` by `map_mul' + psi_Phi + phi_eq_smallRootFactor_mul`.
   Feed to `GMC2DvdKHderiv.hderiv_of_frame` ⇒ `hderiv` ⇒ (closing files) ⇒ GMC(2).

Shared blocker with boxeph (S243 offered help); mac-mini owns the Weierstrass source and `(c)`.
-/

open PowerSeries

namespace GMC2DvdKTranspose

variable {F : Type*} [CommRing F]

/-- The coefficient-swap underlying the transpose embedding: `psiFun φ` has `t`-coefficient `m` the
Laurent series `∑_{k≥0} (φ.coeff k).coeff m · xᵏ`. -/
noncomputable def psiFun (φ : PowerSeries (PowerSeries F)) : PowerSeries (LaurentSeries F) :=
  PowerSeries.mk fun m =>
    HahnSeries.ofPowerSeries ℤ F (PowerSeries.mk fun k => coeff m (coeff k φ))

/-- `t`-order `m`, `x`-order `↑n` coefficient of `psiFun φ` is the double coefficient of `φ`. -/
theorem coeff_psiFun_natCoeff (φ : PowerSeries (PowerSeries F)) (m n : ℕ) :
    ((coeff m (psiFun φ)) : LaurentSeries F).coeff (n : ℤ) = coeff m (coeff n φ) := by
  rw [psiFun, coeff_mk, HahnSeries.ofPowerSeries_apply_coeff, coeff_mk]

end GMC2DvdKTranspose
