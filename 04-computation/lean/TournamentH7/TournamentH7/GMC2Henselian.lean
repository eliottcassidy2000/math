import Mathlib

/-!
# `F⟦X⟧` is Henselian, with `M`-th roots of units — foundations for the THM-1550 small-root lift

The DvdK residual of the GMC(2) formalization (codex THM-2067) reduces to THM-1550: the small-root
product `Π(t) = c·t` of `Φ(X) = X^M − t·R(X)`.  The death-star-S106 route substitutes `X = sZ`,
`t = sᴹ`, so `Φ(sZ) = sᴹ·ψ(Z)` with `ψ(Z) = Zᴹ − R(sZ)` over `F⟦s⟧`, and lifts the small factor by
Hensel over the complete local ring `F⟦s⟧`.

This module supplies two foundations that boxeph-S234 named as blockers:

* **Obstacle (i):** `F⟦X⟧` is a Henselian local ring (`powerSeries_henselianLocalRing`).  Free from
  adic completeness: `F⟦X⟧` is local with maximal ideal `(X)` and `(X)`-adically complete.

* **A refinement of obstacle (ii):** the naive small-factor extraction of `ψ` fails because `ψ` is
  **non-monic** (`ψ mod s = Zᴹ − r₀` drops the degree `d → M`), and Mathlib's Henselian API only
  lifts *simple roots of monic polynomials* (`HenselianLocalRing.TFAE` has no factorization item),
  not degree-dropping factorizations.  The fix is to build the `M` small roots individually as
  `Z = a·Y` with `Yᴹ =` a principal unit — and `M`-th roots of units *do* come from **monic** Hensel,
  since `Zᴹ − C u` is monic.  `exists_pow_eq_of_constantCoeff_pow` is that step, the key sub-lemma of
  the refined route (the remaining pieces are the fixed-point iteration `Y ↦ (R(s·a·Y)/r₀)^{1/M}` and
  THM-1550's Wiener–Hopf product identity).
-/

open Polynomial

namespace GMC2Henselian

variable {F : Type*} [Field F]

/-- The Henselian-at-`(X)` structure on `F⟦X⟧` is free from adic completeness. -/
example : HenselianRing (PowerSeries F) (Ideal.span {(PowerSeries.X : PowerSeries F)}) :=
  inferInstance

/-- **`F⟦X⟧` is a Henselian local ring.**  (boxeph-S234 obstacle (i) for the THM-1550 Hensel lift.)
It is local with maximal ideal `(X)` and `(X)`-adically complete, hence Henselian at its maximal
ideal; the derivative unit condition transfers along `Ideal.Quotient.mk`. -/
instance powerSeries_henselianLocalRing : HenselianLocalRing (PowerSeries F) where
  is_henselian := by
    intro f hf a₀ hfa hunit
    have hHR : HenselianRing (PowerSeries F) (IsLocalRing.maximalIdeal (PowerSeries F)) := by
      rw [PowerSeries.maximalIdeal_eq_span_X]; infer_instance
    exact hHR.is_henselian f hf a₀ hfa (hunit.map _)

/-- **`M`-th roots of a unit via monic Hensel.**  If `a₀^M` is the constant term of `u`, then `u` has
an `M`-th root in `F⟦X⟧` with constant term `a₀`.  The monic polynomial `Xᴹ − C u` has `a₀` as a
simple root modulo `(X)` (its derivative `M·a₀^{M-1}` is a unit since `(M : F) ≠ 0`, `a₀ ≠ 0`), so the
Henselian simple-root property lifts it.  This is the degree-drop-free ingredient of the refined
THM-1550 small-root construction. -/
theorem exists_pow_eq_of_constantCoeff_pow (M : ℕ) (hM : (M : F) ≠ 0)
    (u : PowerSeries F) (a₀ : F) (ha₀ : a₀ ≠ 0)
    (h : a₀ ^ M = PowerSeries.constantCoeff u) :
    ∃ Y : PowerSeries F, Y ^ M = u ∧ PowerSeries.constantCoeff Y = a₀ := by
  have hM0 : M ≠ 0 := by rintro rfl; simp at hM
  have hmon : (X ^ M - C u : (PowerSeries F)[X]).Monic := monic_X_pow_sub_C u hM0
  have heval : ∀ z : PowerSeries F, (X ^ M - C u : (PowerSeries F)[X]).eval z = z ^ M - u := by
    intro z; simp
  have hmem : (X ^ M - C u : (PowerSeries F)[X]).eval (PowerSeries.C a₀)
      ∈ IsLocalRing.maximalIdeal (PowerSeries F) := by
    rw [PowerSeries.maximalIdeal_eq_span_X, Ideal.mem_span_singleton, PowerSeries.X_dvd_iff,
      heval, map_sub, map_pow, PowerSeries.constantCoeff_C, h, sub_self]
  have hunitD : IsUnit ((X ^ M - C u : (PowerSeries F)[X]).derivative.eval (PowerSeries.C a₀)) := by
    have hd : (X ^ M - C u : (PowerSeries F)[X]).derivative.eval (PowerSeries.C a₀)
        = (M : PowerSeries F) * (PowerSeries.C a₀) ^ (M - 1) := by
      simp [derivative_X_pow]
    rw [hd, PowerSeries.isUnit_iff_constantCoeff, map_mul, map_pow, PowerSeries.constantCoeff_C,
      map_natCast]
    exact isUnit_iff_ne_zero.mpr (mul_ne_zero hM (pow_ne_zero _ ha₀))
  obtain ⟨Y, hYroot, hYclose⟩ :=
    HenselianLocalRing.is_henselian _ hmon (PowerSeries.C a₀) hmem hunitD
  refine ⟨Y, sub_eq_zero.mp (by rw [← heval Y]; exact hYroot), ?_⟩
  have h0 : PowerSeries.constantCoeff (Y - PowerSeries.C a₀) = 0 := by
    rw [← PowerSeries.X_dvd_iff, ← Ideal.mem_span_singleton, ← PowerSeries.maximalIdeal_eq_span_X]
    exact hYclose
  rw [map_sub, PowerSeries.constantCoeff_C, sub_eq_zero] at h0
  exact h0

end GMC2Henselian

#print axioms GMC2Henselian.powerSeries_henselianLocalRing
#print axioms GMC2Henselian.exists_pow_eq_of_constantCoeff_pow
