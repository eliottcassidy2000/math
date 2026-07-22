import Mathlib.RingTheory.Polynomial.Hermite.Basic
import Mathlib.Tactic
import TournamentH7.ThreeTermRecurrence

/-!
# Mathlib's Hermite polynomials as a three-term recurrence, and their no-common-root property

Mathlib defines the probabilists' Hermite polynomials `Polynomial.hermite : ℕ → ℤ[X]` and proves
only the *derivative-form* recurrence `hermite_succ : hermite (n+1) = X * hermite n - (hermite n)'`.
This file adds:

* `derivative_hermite_succ` — the ladder relation `(hermite (n+1))' = (n+1) • hermite n`
  (Mathlib has no such lemma);
* `hermite_recurrence` — the classical **three-term recurrence**
  `hermite (n+2) = X * hermite (n+1) - (n+1) * hermite n`;
* `hermite_no_common_root` — **no two consecutive Hermite polynomials share a root**,
  obtained by exhibiting `Polynomial.hermite` as an instance of the abstract
  `ThreeTermRecurrence.ThreeTerm` family and invoking its general no-common-root theorem.

Kernel-pure (`propext, Classical.choice, Quot.sound`).
-/

set_option autoImplicit false

namespace Polynomial

open Polynomial

/-- The ladder relation for probabilists' Hermite polynomials: `He_{n+1}' = (n+1) He_n`.
Proved from Mathlib's `hermite_succ` alone; Mathlib does not otherwise record it. -/
theorem derivative_hermite_succ (n : ℕ) :
    derivative (hermite (n + 1)) = ((n : ℤ) + 1) • hermite n := by
  induction n with
  | zero => simp [hermite_zero]
  | succ k ih =>
    have hk : hermite (k + 1) = X * hermite k - derivative (hermite k) := hermite_succ k
    rw [hermite_succ (k + 1), derivative_sub, derivative_mul, derivative_X, ih]
    simp only [map_smul, mul_smul_comm, one_mul]
    rw [hk]
    module

/-- The classical **three-term recurrence** for probabilists' Hermite polynomials:
`He_{n+2} = X · He_{n+1} - (n+1) · He_n`. -/
theorem hermite_recurrence (n : ℕ) :
    hermite (n + 2) = X * hermite (n + 1) - C ((n : ℤ) + 1) * hermite n := by
  rw [hermite_succ (n + 1), derivative_hermite_succ,
    show ((n : ℤ) + 1) • hermite n = C ((n : ℤ) + 1) * hermite n from smul_eq_C_mul _]

end Polynomial

namespace ThreeTermRecurrence

open Polynomial

/-- Mathlib's Hermite family packaged as an abstract three-term recurrence over `ℤ`:
recentering `a ≡ 0`, off-diagonal `b n = n` (so `b (n+1) = n+1 ≠ 0`). -/
def hermiteZ : ThreeTerm ℤ where
  a := fun _ => 0
  b := fun m => (m : ℤ)
  hb := fun n => by positivity

/-- The abstract polynomial sequence of `hermiteZ` is exactly Mathlib's `Polynomial.hermite`. -/
theorem poly_eq_hermite (n : ℕ) : hermiteZ.poly n = hermite n := by
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    match n with
    | 0 => simp [ThreeTerm.poly]
    | 1 => simp [ThreeTerm.poly, hermiteZ]
    | (k + 2) =>
      rw [ThreeTerm.poly_succ_succ, hermite_recurrence, ih (k + 1) (by omega), ih k (by omega)]
      simp [hermiteZ]

/-- **No two consecutive Hermite polynomials share a root.** A direct application of the
general `ThreeTerm.no_common_root_poly` to the Hermite instance `hermiteZ`. -/
theorem hermite_no_common_root (n : ℕ) (x : ℤ)
    (h1 : (hermite n).IsRoot x) (h2 : (hermite (n + 1)).IsRoot x) : False :=
  hermiteZ.no_common_root_poly n x
    (by rw [poly_eq_hermite]; exact h1) (by rw [poly_eq_hermite]; exact h2)

end ThreeTermRecurrence
