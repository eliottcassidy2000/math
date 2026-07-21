import Mathlib

/-!
# No common root for monic three-term recurrences

A monic three-term recurrence over a commutative ring `R`,
`p 0 = 1`,  `p 1 x = x - a 0`,  `p (n+2) x = (x - a (n+1)) * p (n+1) x - b (n+1) * p n x`,
with nonvanishing off-diagonal `b n ≠ 0`, produces a polynomial sequence in which
**no two consecutive members share a root** (over an integral domain), and hence no
point is a root of the whole family.

This is the abstract core of the classical fact that orthogonal polynomials
(Favard families — Legendre, Hermite, Chebyshev, …) have no common roots; it is proved
from the recurrence alone, needing only `b n ≠ 0` and no zero divisors.

Provenance: extracted and generalized (`ℝ → any integral domain`) from the `GMC2Hermite`
development. Kernel-pure (`propext, Classical.choice, Quot.sound`).
-/

set_option autoImplicit false

namespace ThreeTermRecurrence

/-- Data of a monic three-term recurrence with nonvanishing off-diagonal coefficients. -/
structure ThreeTerm (R : Type*) [Zero R] where
  /-- Diagonal (recentering) coefficients. -/
  a : ℕ → R
  /-- Off-diagonal coefficients. -/
  b : ℕ → R
  /-- The off-diagonal coefficients never vanish. -/
  hb : ∀ n, b n ≠ 0

variable {R : Type*} [CommRing R]

/-- The monic polynomial sequence attached to a three-term recurrence, evaluated at `x`:
`p 0 = 1`, `p 1 = x - a 0`, `p (n+2) = (x - a (n+1)) * p (n+1) - b (n+1) * p n`. -/
def ThreeTerm.p (T : ThreeTerm R) : ℕ → R → R
  | 0,       _ => 1
  | 1,       x => x - T.a 0
  | (n + 2), x => (x - T.a (n + 1)) * T.p (n + 1) x - T.b (n + 1) * T.p n x

@[simp] theorem ThreeTerm.p_zero (T : ThreeTerm R) (x : R) : T.p 0 x = 1 := rfl

theorem ThreeTerm.p_succ_succ (T : ThreeTerm R) (n : ℕ) (x : R) :
    T.p (n + 2) x = (x - T.a (n + 1)) * T.p (n + 1) x - T.b (n + 1) * T.p n x := rfl

/-- **Consecutive members of a monic three-term recurrence with `b ≠ 0` have no common
root**, over any integral domain. -/
theorem ThreeTerm.no_common_root [IsDomain R] (T : ThreeTerm R) :
    ∀ (n : ℕ) (x : R), T.p n x = 0 → T.p (n + 1) x = 0 → False := by
  intro n
  induction n with
  | zero =>
      intro x h _
      rw [T.p_zero] at h
      exact one_ne_zero h
  | succ k ih =>
      intro x h1 h2
      have hrec := T.p_succ_succ k x
      rw [h2, h1] at hrec
      have hbpk : T.b (k + 1) * T.p k x = 0 := by linear_combination hrec
      have hk : T.p k x = 0 := (mul_eq_zero.mp hbpk).resolve_left (T.hb (k + 1))
      exact ih x hk h1

/-- No point is a root of every member of a three-term family (over an integral domain). -/
theorem ThreeTerm.exists_nonvanishing [IsDomain R] (T : ThreeTerm R) (x : R) :
    ∃ n : ℕ, T.p n x ≠ 0 := by
  by_cases h : T.p 0 x = 0
  · exact ⟨1, fun h1 => T.no_common_root 0 x h h1⟩
  · exact ⟨0, h⟩

/-- The probabilists' Hermite family as an instance over `ℝ`: `a ≡ 0`, `b n = n + 1 ≠ 0`.
(Any characteristic-zero field works; `ℝ` chosen for a concrete witness.) -/
noncomputable def hermiteReal : ThreeTerm ℝ where
  a := fun _ => 0
  b := fun n => (n : ℝ) + 1
  hb := fun n => by positivity

end ThreeTermRecurrence
