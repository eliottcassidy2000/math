import Mathlib
import TournamentH7.GMC2DvdKFrame

/-!
# The x-degree toolkit for the frame (the (c) degree lemma `xCoeff0(logDeriv P) = 0`)

`xDegLE φ d` says every `t`-coefficient of a frame element `φ ∈ (F⸨x⸩)⟦t⟧` has `x`-support in
`(-∞, d]`.  The payoff is `xCoeff0_eq_zero_of_xDegLE_neg`: a strictly-negative x-degree kills the `x⁰`
coefficient.  Closure under `+`, `×` (x-degrees add, via `HahnSeries.support_mul_subset`), and `∂_t`
(preserved) reduces the degree lemma `(c)` to the structural x-degrees of `P` and `1/P`:
`∂_t P` has x-degree `≤ M-1`, `1/P` has x-degree `≤ -M`, so `logDeriv P = ∂_t P / P` has x-degree
`≤ -1 < 0`, hence `xCoeff0(logDeriv P) = 0`.  Reusable, kernel-pure, no residues/Puiseux.
-/

open PowerSeries GMC2DvdKFrame

namespace GMC2DvdKXDeg

variable {F : Type*} [Field F]

/-- `xDegLE φ d` : every `t`-coefficient of `φ` has `x`-support in `(-∞, d]` (x-degree ≤ d). -/
def xDegLE (φ : PowerSeries (LaurentSeries F)) (d : ℤ) : Prop :=
  ∀ k : ℕ, (PowerSeries.coeff (R := LaurentSeries F) k φ).support ⊆ Set.Iic d

theorem xCoeff0_eq_zero_of_xDegLE_neg {φ : PowerSeries (LaurentSeries F)} {d : ℤ}
    (h : xDegLE φ d) (hd : d < 0) : xCoeff0 φ = 0 := by
  ext k
  rw [coeff_xCoeff0, map_zero]
  by_contra hne
  exact absurd (h k hne) (by simp only [Set.mem_Iic, not_le]; omega)

theorem xDegLE_add {φ ψ : PowerSeries (LaurentSeries F)} {d : ℤ}
    (hφ : xDegLE φ d) (hψ : xDegLE ψ d) : xDegLE (φ + ψ) d := by
  intro k
  rw [map_add]
  exact (HahnSeries.support_add_subset _ _).trans (Set.union_subset (hφ k) (hψ k))

/-- Laurent-series helper: `x`-degrees add under multiplication. -/
theorem support_mul_Iic {a b : LaurentSeries F} {d e : ℤ}
    (ha : a.support ⊆ Set.Iic d) (hb : b.support ⊆ Set.Iic e) :
    (a * b).support ⊆ Set.Iic (d + e) := by
  intro n hn
  obtain ⟨m, hm, p, hp, rfl⟩ := Set.mem_add.mp (HahnSeries.support_mul_subset hn)
  exact Set.mem_Iic.mpr (add_le_add (Set.mem_Iic.mp (ha hm)) (Set.mem_Iic.mp (hb hp)))

theorem xDegLE_mul {φ ψ : PowerSeries (LaurentSeries F)} {d e : ℤ}
    (hφ : xDegLE φ d) (hψ : xDegLE ψ e) : xDegLE (φ * ψ) (d + e) := by
  intro k n hn
  rw [HahnSeries.mem_support, PowerSeries.coeff_mul, HahnSeries.coeff_sum] at hn
  obtain ⟨p, _, hp⟩ := Finset.exists_ne_zero_of_sum_ne_zero hn
  exact support_mul_Iic (hφ p.1) (hψ p.2) hp

/-- The nat-cast constant `(↑m : F⸨x⸩)` has `x`-support `⊆ {0}`. -/
theorem support_natCast_Iic (m : ℕ) : ((m : LaurentSeries F)).support ⊆ Set.Iic 0 := by
  intro n hn; rw [Set.mem_Iic]; by_contra hle; push_neg at hle
  refine hn ?_
  rw [← map_natCast (HahnSeries.C (R := F) (Γ := ℤ)) m, HahnSeries.C_apply,
    HahnSeries.coeff_single_of_ne (by omega : n ≠ 0)]

/-- `∂_t` preserves x-degree (its `t`-coefficient is `φ.coeff (k+1)` times a constant). -/
theorem xDegLE_derivativeFun {φ : PowerSeries (LaurentSeries F)} {d : ℤ}
    (hφ : xDegLE φ d) : xDegLE (PowerSeries.derivativeFun φ) d := by
  intro k
  rw [PowerSeries.coeff_derivativeFun]
  have := support_mul_Iic (hφ (k + 1)) (support_natCast_Iic (F := F) (k + 1))
  rwa [add_zero, Nat.cast_add, Nat.cast_one] at this

/-- The constant `1` has x-degree `≤ 0`. -/
theorem xDegLE_one : xDegLE (1 : PowerSeries (LaurentSeries F)) 0 := by
  intro k n hn
  rw [Set.mem_Iic]; by_contra hle; push_neg at hle
  refine hn ?_
  rcases Nat.eq_zero_or_pos k with hk | hk
  · subst hk
    rw [show (PowerSeries.coeff (R := LaurentSeries F) 0 1) = 1 by simp]
    rw [← HahnSeries.single_zero_one, HahnSeries.coeff_single_of_ne (by omega : n ≠ 0)]
  · rw [PowerSeries.coeff_one, if_neg (by omega)]; simp

/-- The monomial `C(x^m)` (constant in `t`) has x-degree `≤ m`. -/
theorem xDegLE_C_xpow (m : ℤ) :
    xDegLE (PowerSeries.C (HahnSeries.single m (1 : F))) m := by
  intro k n hn
  rw [Set.mem_Iic]; by_contra hle; push_neg at hle
  refine hn ?_
  rcases Nat.eq_zero_or_pos k with hk | hk
  · subst hk
    rw [PowerSeries.coeff_zero_eq_constantCoeff, PowerSeries.constantCoeff_C,
      HahnSeries.coeff_single_of_ne (by omega : n ≠ m)]
  · rw [PowerSeries.coeff_C, if_neg (by omega)]; simp

end GMC2DvdKXDeg

#print axioms GMC2DvdKXDeg.xCoeff0_eq_zero_of_xDegLE_neg
#print axioms GMC2DvdKXDeg.xDegLE_mul
#print axioms GMC2DvdKXDeg.xDegLE_derivativeFun
