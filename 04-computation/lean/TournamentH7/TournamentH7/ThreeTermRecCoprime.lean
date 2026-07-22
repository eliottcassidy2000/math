import Mathlib

/-!
# Coprimality of consecutive polynomials of a three-term recurrence

A *monic three-term recurrence* over a commutative ring `R` is determined by two coefficient
sequences `a b : ℕ → R` and produces the monic polynomials
`p 0 = 1`, `p 1 = X - C (a 0)`, `p (n+2) = (X - C (a (n+1))) * p (n+1) - C (b (n+1)) * p n`.

This is the shape of every family of orthogonal polynomials (Favard's theorem: a sequence of monic
polynomials is orthogonal for some moment functional iff it satisfies such a recurrence with
`b n ≠ 0`) and of Sturm sequences. The main result `ThreeTermRec.isCoprime_succ` is that when every
off-diagonal coefficient `b n` is a unit, consecutive members `p n` and `p (n+1)` are coprime; over
a field this is the classical fact that consecutive orthogonal polynomials share no root, uniformly
for Hermite, Legendre, Laguerre, Chebyshev, Gegenbauer, ….

## Main results
* `ThreeTermRec.isCoprime_succ` — consecutive members are coprime (units off-diagonal).
* `ThreeTermRec.isCoprime_succ_of_ne_zero` — the field specialization (`b n ≠ 0`).
* `ThreeTermRec.noCommonRoot` — consecutive members have no common root.
-/

open Polynomial

/-- Coefficient data of a monic three-term recurrence over a commutative ring `R`. -/
structure ThreeTermRec (R : Type*) [CommRing R] where
  /-- Diagonal (recentring) coefficients. -/
  a : ℕ → R
  /-- Off-diagonal coefficients; units (nonzero, over a field) in the main results. -/
  b : ℕ → R

namespace ThreeTermRec

variable {R : Type*} [CommRing R]

/-- The monic polynomial sequence of the recurrence:
`p 0 = 1`, `p 1 = X - C (a 0)`, `p (n+2) = (X - C (a (n+1))) * p (n+1) - C (b (n+1)) * p n`. -/
noncomputable def p (T : ThreeTermRec R) : ℕ → R[X]
  | 0       => 1
  | 1       => X - C (T.a 0)
  | (n + 2) => (X - C (T.a (n + 1))) * p T (n + 1) - C (T.b (n + 1)) * p T n

variable (T : ThreeTermRec R)

@[simp] theorem p_zero : T.p 0 = 1 := rfl

theorem p_one : T.p 1 = X - C (T.a 0) := rfl

theorem p_add_two (n : ℕ) :
    T.p (n + 2) = (X - C (T.a (n + 1))) * T.p (n + 1) - C (T.b (n + 1)) * T.p n := rfl

/-- Each member `p n` is monic of degree `n`. -/
theorem monic_and_natDegree [Nontrivial R] (n : ℕ) :
    (T.p n).Monic ∧ (T.p n).natDegree = n := by
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    match n with
    | 0 => exact ⟨monic_one, by simp⟩
    | 1 => exact ⟨monic_X_sub_C _, by rw [p_one, natDegree_X_sub_C]⟩
    | (k + 2) =>
        obtain ⟨hm1, hd1⟩ := ih (k + 1) (by omega)
        obtain ⟨_, hdk⟩ := ih k (by omega)
        have hA : ((X - C (T.a (k + 1))) * T.p (k + 1)).Monic := (monic_X_sub_C _).mul hm1
        have hAd : ((X - C (T.a (k + 1))) * T.p (k + 1)).natDegree = k + 2 := by
          rw [(monic_X_sub_C (T.a (k + 1))).natDegree_mul hm1, natDegree_X_sub_C, hd1]
          omega
        have hBk : (C (T.b (k + 1)) * T.p k).natDegree ≤ k :=
          le_trans (natDegree_C_mul_le _ _) hdk.le
        have hBlt : (C (T.b (k + 1)) * T.p k).degree
            < ((X - C (T.a (k + 1))) * T.p (k + 1)).degree := by
          rw [degree_eq_natDegree hA.ne_zero, hAd]
          calc (C (T.b (k + 1)) * T.p k).degree
              ≤ ((C (T.b (k + 1)) * T.p k).natDegree : WithBot ℕ) := degree_le_natDegree
            _ ≤ (k : WithBot ℕ) := by exact_mod_cast hBk
            _ < ((k + 2 : ℕ) : WithBot ℕ) := by exact_mod_cast (by omega : k < k + 2)
        have hpe : T.p (k + 2)
            = (X - C (T.a (k + 1))) * T.p (k + 1) + -(C (T.b (k + 1)) * T.p k) := by
          rw [p_add_two]; ring
        refine ⟨?_, ?_⟩
        · rw [hpe]; exact hA.add_of_left (by rwa [degree_neg])
        · rw [hpe, natDegree_eq_of_degree_eq
              (degree_add_eq_left_of_degree_lt (by rwa [degree_neg])), hAd]

/-- Each member `p n` is monic. -/
theorem monic [Nontrivial R] (n : ℕ) : (T.p n).Monic := (T.monic_and_natDegree n).1

/-- Each member `p n` has degree `n`. -/
theorem natDegree_p [Nontrivial R] (n : ℕ) : (T.p n).natDegree = n := (T.monic_and_natDegree n).2

/-- **Consecutive members of a three-term recurrence are coprime** whenever every off-diagonal
coefficient `b n` is a unit. Holds over any commutative ring.

The proof is a one-line Bézout update: if `u·pₙ + v·pₙ₊₁ = 1`, then using
`pₙ = b⁻¹((X - a)·pₙ₊₁ - pₙ₊₂)` one gets an explicit Bézout identity for `(pₙ₊₁, pₙ₊₂)`. -/
theorem isCoprime_succ (hb : ∀ n, IsUnit (T.b (n + 1))) :
    ∀ n, IsCoprime (T.p n) (T.p (n + 1)) := by
  intro n
  induction n with
  | zero => exact isCoprime_one_left
  | succ k ih =>
      obtain ⟨u, v, huv⟩ := ih
      have hCb : IsUnit (C (T.b (k + 1))) := isUnit_C.mpr (hb k)
      set binv := Ring.inverse (C (T.b (k + 1))) with hbdef
      have hbinv : binv * C (T.b (k + 1)) = 1 := Ring.inverse_mul_cancel _ hCb
      have hrec : T.p (k + 2)
          = (X - C (T.a (k + 1))) * T.p (k + 1) - C (T.b (k + 1)) * T.p k := rfl
      refine ⟨u * binv * (X - C (T.a (k + 1))) + v, -(u * binv), ?_⟩
      linear_combination (-(u * binv)) * hrec + huv + (u * T.p k) * hbinv

/-- Field specialization: over a field, `b (n+1) ≠ 0` suffices. -/
theorem isCoprime_succ_of_ne_zero {K : Type*} [Field K] (T : ThreeTermRec K)
    (hb : ∀ n, T.b (n + 1) ≠ 0) (n : ℕ) : IsCoprime (T.p n) (T.p (n + 1)) :=
  T.isCoprime_succ (fun n => isUnit_iff_ne_zero.mpr (hb n)) n

/-- **Consecutive members share no root** (off-diagonal coefficients units). -/
theorem noCommonRoot [Nontrivial R] (hb : ∀ n, IsUnit (T.b (n + 1))) (n : ℕ) (x : R) :
    ¬ ((T.p n).eval x = 0 ∧ (T.p (n + 1)).eval x = 0) := by
  rintro ⟨h1, h2⟩
  have hco := (T.isCoprime_succ hb n).map (evalRingHom x)
  simp only [coe_evalRingHom, h1, h2] at hco
  exact not_isUnit_zero (isCoprime_zero_left.mp hco)

/-! ### The probabilists' Hermite family as an instance -/

/-- The probabilists' Hermite recurrence `Heₙ₊₁ = X·Heₙ − n·Heₙ₋₁`: here `a ≡ 0` and `b k = k`
(only `b (n+1) = n+1 ≠ 0` is used). -/
noncomputable def hermite : ThreeTermRec ℝ where
  a := fun _ => 0
  b := fun n => (n : ℝ)

theorem hermite_b_isUnit (n : ℕ) : IsUnit (hermite.b (n + 1)) :=
  isUnit_iff_ne_zero.mpr (by simp [hermite]; positivity)

/-- Consecutive probabilists' Hermite polynomials are coprime. -/
theorem hermite_isCoprime_succ (n : ℕ) :
    IsCoprime (hermite.p n) (hermite.p (n + 1)) :=
  hermite.isCoprime_succ hermite_b_isUnit n

end ThreeTermRec

/-
All public results depend only on `[propext, Classical.choice, Quot.sound]`
(no `sorry`, no `native_decide`). Verified via `#print axioms ThreeTermRec.isCoprime_succ` etc.
-/
