import Mathlib

/-!
# Irreducibility of `Φ = X^M − t·R(X)` over `F(t)` — the transitivity input to THM-2067

For `R ∈ F[X]` with `R(0) ≠ 0`, the polynomial `Φ = X^M − t·R(X)` is irreducible over `F(t)`, so its
Galois group acts transitively on its roots (`Polynomial.Gal.galAction_isPretransitive`) — the input
consumed by the orbit-product core (`GMC2OrbitProduct`).

Strategy (Gauss): `Φ` is *linear* in `t` over `F[X]`, with coprime coefficients `X^M` and `−R`
(coprime because `X ∤ R`, i.e. `R(0) ≠ 0`), hence irreducible in `F[X][t]`.  This module proves that
first step (`phi_t_irreducible`); the transport across the `F[X][t] ≅ F[t][X]` swap and Gauss to
`F(t)[X]` is the remaining step.
-/

open Polynomial

namespace GMC2PhiIrreducible

variable {F : Type*} [Field F]

/-- `X^M` and `R` are relatively prime in `F[X]` when `R(0) ≠ 0` (so `X ∤ R`). -/
theorem isRelPrime_X_pow_R (R : F[X]) (M : ℕ) (hR0 : R.coeff 0 ≠ 0) :
    IsRelPrime ((X : F[X]) ^ M) R := by
  have hXR : IsCoprime (X : F[X]) R :=
    (prime_X.coprime_iff_not_dvd).mpr (by rw [Polynomial.X_dvd_iff]; exact hR0)
  exact (hXR.pow_left).isRelPrime

/-- **`Φ` is irreducible as a linear polynomial in `t` over `F[X]`.**  Here the outer variable `X` of
`F[X][t]` plays the role of `t`; `Φ = C(Xᴹ) − C(R)·t`. -/
theorem phi_t_irreducible (R : F[X]) (M : ℕ) (hR0 : R.coeff 0 ≠ 0) :
    Irreducible
      (Polynomial.C ((X : F[X]) ^ M) - Polynomial.C R * X : Polynomial (Polynomial F)) := by
  have hRne : R ≠ 0 := fun h => hR0 (by rw [h]; simp)
  have hrw : (Polynomial.C ((X : F[X]) ^ M) - Polynomial.C R * X : Polynomial (Polynomial F))
      = Polynomial.C (-R) * X + Polynomial.C ((X : F[X]) ^ M) := by
    simp only [map_neg]; ring
  rw [hrw]
  apply Polynomial.irreducible_of_degree_eq_one_of_isRelPrime_coeff
  · -- degree of the linear form `C(-R)·X + C(Xᴹ)` is 1
    rw [Polynomial.degree_linear (neg_ne_zero.mpr hRne)]
  · -- coeff 0 = Xᴹ, coeff 1 = −R, and these are relatively prime
    have hc0 : (Polynomial.C (-R) * X + Polynomial.C ((X : F[X]) ^ M) :
        Polynomial (Polynomial F)).coeff 0 = (X : F[X]) ^ M := by
      rw [Polynomial.coeff_add, Polynomial.coeff_C_zero, Polynomial.coeff_C_mul,
        Polynomial.coeff_X_zero, mul_zero, zero_add]
    have hc1 : (Polynomial.C (-R) * X + Polynomial.C ((X : F[X]) ^ M) :
        Polynomial (Polynomial F)).coeff 1 = -R := by
      rw [Polynomial.coeff_add, Polynomial.coeff_C_mul, Polynomial.coeff_X_one, mul_one,
        Polynomial.coeff_C, if_neg (by norm_num : ¬ ((1 : ℕ) = 0)), add_zero]
    rw [hc0, hc1]
    exact (isRelPrime_X_pow_R R M hR0).neg_right

end GMC2PhiIrreducible

#print axioms GMC2PhiIrreducible.phi_t_irreducible
