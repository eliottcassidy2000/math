import Mathlib

/-!
# Irreducibility of the DvdK parameter polynomial `Y^M - t·R(Y)` over `F(t)`

This is **death-star-S106 target (2)** for the GMC(2)/DvdK Galois route (THM-2067):
the polynomial `Φ = Y^M - t·R(Y)` is irreducible over the rational function field `F(t)`,
whenever `R(0) ≠ 0` and `M ≥ 1`.  It is the input that makes the Galois action on the roots
of `Φ` transitive (a single orbit), which the orbit-product lemma then uses.

The proof is elementary and self-contained (imports only Mathlib), kernel-pure, and stated for
an arbitrary field `F` — so it is directly Mathlib-PR-shaped:

* `Φ` is **linear in the parameter `t`** (the inner variable `X`), so under the bivariate swap
  `X ↔ Y` it becomes `C(X^M) - C(R)·Y`, a degree-one polynomial in the outer variable whose two
  coefficients are `X^M` and `-R`.  These are **coprime** in `F[X]` because `R(0) ≠ 0` forces
  `X ∤ R`, so the swapped polynomial is primitive, hence irreducible (degree one + primitive).
* Irreducibility transfers back along the swap (a ring automorphism), and then to `F(t) = F(X)`
  by Gauss's lemma (`IsPrimitive.irreducible_iff_irreducible_map_fraction_map`).

Conventions (matching `Mathlib.Algebra.Polynomial.Bivariate`): `F[X][Y] = Polynomial (Polynomial F)`,
the inner variable `X` is the DvdK parameter `t`, and the outer variable `Y` is the polynomial
variable of `Φ`.  `R : F[X]` is used as `R.map C : F[X][Y]` — i.e. `R` in the outer variable `Y`.
-/

open Polynomial
open scoped Polynomial.Bivariate

namespace GMC2DvdKParameterIrreducible

variable {R : Type*} [CommRing R] [IsDomain R]

/-- A degree-one primitive polynomial over an integral domain is irreducible. -/
theorem irreducible_of_natDegree_eq_one_isPrimitive
    {p : R[X]} (hdeg : p.natDegree = 1) (hprim : p.IsPrimitive) : Irreducible p := by
  have hp0 : p ≠ 0 := by
    intro h; rw [h, natDegree_zero] at hdeg; exact one_ne_zero hdeg.symm
  refine ⟨?_, ?_⟩
  · intro hu
    have := natDegree_eq_zero_of_isUnit hu
    rw [hdeg] at this; exact one_ne_zero this
  · intro a b hab
    have ha : a ≠ 0 := fun h => hp0 (by rw [hab, h, zero_mul])
    have hb : b ≠ 0 := fun h => hp0 (by rw [hab, h, mul_zero])
    have hdegs : a.natDegree + b.natDegree = 1 := by
      rw [← natDegree_mul ha hb, ← hab, hdeg]
    rcases Nat.add_eq_one_iff.mp hdegs with ⟨ha0, _⟩ | ⟨_, hb0⟩
    · obtain ⟨c, hc⟩ := natDegree_eq_zero.mp ha0
      left
      have hCdvd : C c ∣ p := ⟨b, by rw [hab, ← hc]⟩
      exact (hc ▸ (hprim c hCdvd).map C)
    · obtain ⟨c, hc⟩ := natDegree_eq_zero.mp hb0
      right
      have hCdvd : C c ∣ p := ⟨a, by rw [hab, ← hc]; ring⟩
      exact (hc ▸ (hprim c hCdvd).map C)

section Field
variable {F : Type*} [Field F]

/-- **One-variable Duistermaat–van der Kallen, irreducibility half** (death-star-S106 target 2).
For `Rp : F[X]` with `Rp(0) ≠ 0` and `M ≥ 1`, the polynomial `Y^M - t·Rp(Y)` is irreducible over
the rational function field `F(t)`.  Here `t` is the inner variable `X`, `Y` is the outer
polynomial variable, and the coefficients live in `F[X] = F[t]`. -/
theorem irreducible_map_ratfunc
    (Rp : F[X]) (M : ℕ) (hM : 1 ≤ M) (hR0 : Rp.eval 0 ≠ 0) :
    Irreducible
      (((Y : F[X][Y]) ^ M - C (Polynomial.X) * Rp.map C).map
        (algebraMap F[X] (RatFunc F))) := by
  set Φ : F[X][Y] := (Y : F[X][Y]) ^ M - C (Polynomial.X) * Rp.map C with hΦ
  have hRp0 : Rp ≠ 0 := fun h => hR0 (by rw [h, eval_zero])
  -- swap `X ↔ Y`: `Φ` becomes the degree-one `g = C(X^M) - C(Rp)·Y`
  have hswap : Bivariate.swap Φ = C (Polynomial.X ^ M) - C Rp * Y := by
    rw [hΦ, map_sub, map_pow, map_mul, Bivariate.swap_X, Bivariate.swap_Y]
    have hsw : Bivariate.swap (Rp.map C) = C Rp := by
      have h1 : Bivariate.swap (C Rp) = Rp.map C := Bivariate.swap_C Rp
      have h2 := congrArg Bivariate.swap h1
      rw [Bivariate.swap_swap_apply] at h2
      exact h2.symm
    rw [hsw, map_pow]; ring
  set g : F[X][Y] := C (Polynomial.X ^ M) - C Rp * Y with hg
  have hg1 : g.coeff 1 = -Rp := by
    rw [hg, coeff_sub, coeff_C, coeff_C_mul, coeff_X, if_neg (by norm_num), if_pos rfl,
      mul_one, zero_sub]
  have hg0 : g.coeff 0 = Polynomial.X ^ M := by
    rw [hg, coeff_sub, coeff_C, coeff_C_mul, coeff_X, if_pos rfl, if_neg (by norm_num),
      mul_zero, sub_zero]
  have hgdeg : g.natDegree = 1 := by
    have hupper : g.natDegree ≤ 1 := by
      rw [hg]
      refine (natDegree_sub_le _ _).trans ?_
      simp only [natDegree_C, Nat.zero_le, sup_le_iff, true_and]
      calc (C Rp * Y).natDegree ≤ (C Rp).natDegree + (Y : F[X][Y]).natDegree := natDegree_mul_le
        _ ≤ 0 + 1 := by rw [natDegree_C, natDegree_X]
        _ = 1 := by ring
    have hlower : 1 ≤ g.natDegree :=
      le_natDegree_of_ne_zero (by rw [hg1]; exact neg_ne_zero.mpr hRp0)
    omega
  -- `X ∤ Rp` (from `Rp(0) ≠ 0`) makes the coefficients `X^M`, `-Rp` coprime, so `g` is primitive
  have hXnotdvd : ¬ (Polynomial.X ∣ Rp) := by
    rw [Polynomial.X_dvd_iff, Polynomial.coeff_zero_eq_eval_zero]; exact hR0
  have hcop : IsCoprime (Polynomial.X ^ M : F[X]) Rp :=
    ((Polynomial.irreducible_X).coprime_iff_not_dvd.mpr hXnotdvd).pow_left
  have hgprim : g.IsPrimitive := by
    rw [isPrimitive_iff_isUnit_of_C_dvd]
    intro c hc
    rw [Polynomial.C_dvd_iff_dvd_coeff] at hc
    have hd0 : c ∣ Polynomial.X ^ M := by have := hc 0; rwa [hg0] at this
    have hd1 : c ∣ Rp := by have := hc 1; rw [hg1] at this; exact (dvd_neg).mp this
    exact hcop.isUnit_of_dvd' hd0 hd1
  have hgirr : Irreducible g := irreducible_of_natDegree_eq_one_isPrimitive hgdeg hgprim
  -- transfer irreducibility back through the swap, then to `F(t)` by Gauss's lemma
  have hΦirr : Irreducible Φ :=
    (MulEquiv.irreducible_iff Bivariate.swap).mp (by rw [hswap]; exact hgirr)
  have hΦnd : Φ.natDegree ≠ 0 := by
    have hcM : Φ.coeff M = 1 - Polynomial.X * C (Rp.coeff M) := by
      rw [hΦ, coeff_sub, coeff_C_mul, coeff_map, coeff_X_pow, if_pos rfl]
    have hne : Φ.coeff M ≠ 0 := by
      rw [hcM]; intro h; have := congrArg (Polynomial.eval (0 : F)) h; simp at this
    have := le_natDegree_of_ne_zero hne; omega
  exact ((hΦirr.isPrimitive hΦnd).irreducible_iff_irreducible_map_fraction_map).mp hΦirr

end Field

#print axioms irreducible_map_ratfunc

end GMC2DvdKParameterIrreducible
