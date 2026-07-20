/-
# GMC(2), the {-1,0,1} stratum: no point kills every Hermite polynomial

kind-pasteur-2026-07-20-S128c120.  Companion to THM-1585 and
`02-court/active/CASE-gamma-bridge-domination-step.md`.

## Why this lemma

klein-S351's "Gamma Bridge" proposed to derive NC2 (hence GMC(2)) from the toral
statement by arguing that `E_r[psi_m] = sum_k c_k * k!` is *dominated by its top term*.
THM-1585 refutes that: measured in exact arithmetic, the ratio of the second term to the
top term GROWS linearly in `m` (to 45x), and the top term's share of the sum falls to
0.04%.  The mass sits at an INTERIOR index -- a saddle -- so no domination argument can
work, and the bridge as published is broken.

The repair is not a better estimate; it is a change of object.  On the `{-1,0,1}`
stratum at `M = 1`, Lagrange-Buermann gives the exact closed form

  `psi_m = (1/m) * sum_k  m!/(k!^2 (m-2k)!) * W^k * B^{m-2k}`,   `B = b`, `W = r*a*c`,

and for CONSTANT `a, b, c` the Gaussian radial average `E_r[r^k] = k!` cancels one `k!`:

  `m * E_r[psi_m] = sum_k m!/(k!(m-2k)!) * w^k * b^{m-2k} = s^m * He_m(b/s)`,
  `w = a*c`,  `s = sqrt(-2w)`,

where `He` is the probabilists' Hermite polynomial.  So `E_r[psi_m] = 0` iff `b/s` is a
root of `He_m`, and `P` lies in the nullcone iff `b/s` is a root of `He_m` for EVERY `m`.

**That is impossible, and this file proves it.**  Consecutive Hermite polynomials share
no root, so no point can be a root of all of them.  Hence for `a*c != 0` -- i.e. whenever
both extreme charges are present, i.e. whenever `P` is two-sided -- some `E_r[psi_m]` is
nonzero and `P` is not in the nullcone.  That is the one-sided conjecture on this
stratum, with no asymptotics, no domination, and no saddle-point analysis.

Everything below is `sorry`-free and uses no `native_decide`.

## Scope, stated honestly

This file formalises the step that replaces the false one: *no common root*.  It does
NOT formalise the Lagrange-Buermann closed form or the identity
`sum_k m!/(k!(m-2k)!) (-1/2)^k b^{m-2k} = He_m(b)`; both are verified exactly in
`04-computation/hermite_closure_gmc2_kps_S128c120.py` but are not proved here.  So this
is one link of the chain, machine-checked; the chain itself is not yet fully formal.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic

namespace GMC2Hermite

/-- Probabilists' Hermite polynomials, as functions `ℝ → ℝ`, defined by the three-term
recurrence `He 0 = 1`, `He 1 = X`, `He (n+2) = X * He (n+1) - (n+1) * He n`. -/
def He : ℕ → ℝ → ℝ
  | 0,       _ => 1
  | 1,       x => x
  | (n + 2), x => x * He (n + 1) x - ((n : ℝ) + 1) * He n x

@[simp] theorem He_zero (x : ℝ) : He 0 x = 1 := rfl

@[simp] theorem He_one (x : ℝ) : He 1 x = x := rfl

theorem He_succ_succ (n : ℕ) (x : ℝ) :
    He (n + 2) x = x * He (n + 1) x - ((n : ℝ) + 1) * He n x := rfl

theorem He_two (x : ℝ) : He 2 x = x ^ 2 - 1 := by
  rw [He_succ_succ]; simp; ring

theorem He_three (x : ℝ) : He 3 x = x ^ 3 - 3 * x := by
  rw [He_succ_succ, He_two]; simp; ring

/-- **Consecutive Hermite polynomials have no common root.**

The recurrence `He (n+2) = X * He (n+1) - (n+1) * He n` turns a common root of
`He (n+1)` and `He (n+2)` into a root of `He n` (using `n + 1 ≠ 0`, i.e. characteristic
zero), so the descent reaches `He 0 = 1`, which has no root. -/
theorem no_common_root : ∀ (n : ℕ) (x : ℝ), He n x = 0 → He (n + 1) x = 0 → False := by
  intro n
  induction n with
  | zero =>
      intro x h _
      rw [He_zero] at h
      exact one_ne_zero h
  | succ k ih =>
      intro x h1 h2
      -- `h1 : He (k+1) x = 0`, `h2 : He (k+2) x = 0`
      have hrec : He (k + 2) x = x * He (k + 1) x - ((k : ℝ) + 1) * He k x :=
        He_succ_succ k x
      rw [h2, h1] at hrec
      -- now `hrec : 0 = x * 0 - (k+1) * He k x`
      have hk1 : ((k : ℝ) + 1) * He k x = 0 := by linarith
      have hne : ((k : ℝ) + 1) ≠ 0 := by positivity
      have hk : He k x = 0 := by
        rcases mul_eq_zero.mp hk1 with h | h
        · exact absurd h hne
        · exact h
      exact ih x hk h1

/-- **No real number is a root of every Hermite polynomial of positive degree.**

This is the statement that closes the sign-mixed branch of the `{-1,0,1}` stratum: since
`m * E_r[psi_m] = s^m * He m (b/s)`, membership of the nullcone would require `b/s` to be
a common root of all `He m`, which cannot happen. -/
theorem exists_nonvanishing (x : ℝ) : ∃ n : ℕ, 1 ≤ n ∧ He n x ≠ 0 := by
  by_cases h : He 1 x = 0
  · exact ⟨2, by norm_num, fun h2 => no_common_root 1 x h h2⟩
  · exact ⟨1, le_rfl, h⟩

/-- Sharper form: among any two consecutive degrees, at least one is nonvanishing. -/
theorem nonvanishing_of_consecutive (n : ℕ) (x : ℝ) :
    He n x ≠ 0 ∨ He (n + 1) x ≠ 0 := by
  by_cases h : He n x = 0
  · exact Or.inr (fun h2 => no_common_root n x h h2)
  · exact Or.inl h

/-!
## The same closure covers mac-mini-S140's degree-1 layer

THM-1600 (mac-mini-S140) computes the Laplace/GMC(1) layer at degree 1 and finds it is an
IDENTITY: `L((a*v + b)^m) = m! * a^m * e_m(b/a)`, where `e_m` is the TRUNCATED EXPONENTIAL
`e_m(z) = sum_{k <= m} z^k / k!`.  That is the same shape as the degree-2 answer above --
the Gamma average of an `m`-th power is a classical polynomial sequence evaluated at one
fixed point -- so the nullcone condition again asks for a COMMON ROOT of the whole
sequence, and again there is none.

The proof is even shorter than the Hermite one: consecutive truncated exponentials differ
by a monomial, `e_{n+1}(z) - e_n(z) = z^{n+1}/(n+1)!`, so a common root forces `z = 0`,
where `e_n(0) = 1 != 0`.

Two different classical families, one argument shape.  Neither needs an estimate.
-/

/-- Truncated exponential `trExp n z = sum_{k = 0}^{n} z^k / k!`. -/
noncomputable def trExp : ℕ → ℝ → ℝ
  | 0,       _ => 1
  | (n + 1), z => trExp n z + z ^ (n + 1) / (Nat.factorial (n + 1) : ℝ)

@[simp] theorem trExp_zero_arg : ∀ n : ℕ, trExp n 0 = 1
  | 0 => rfl
  | (n + 1) => by simp [trExp, trExp_zero_arg n]

theorem trExp_succ_sub (n : ℕ) (z : ℝ) :
    trExp (n + 1) z - trExp n z = z ^ (n + 1) / (Nat.factorial (n + 1) : ℝ) := by
  simp [trExp]

/-- **Consecutive truncated exponentials have no common root** — the degree-1 analogue of
`no_common_root`, closing mac-mini-S140/THM-1600's layer by the same algebraic route. -/
theorem trExp_no_common_root (n : ℕ) (z : ℝ) :
    trExp n z = 0 → trExp (n + 1) z = 0 → False := by
  intro h0 h1
  have hd : z ^ (n + 1) / (Nat.factorial (n + 1) : ℝ) = 0 := by
    rw [← trExp_succ_sub n z, h0, h1]; ring
  have hfac : (Nat.factorial (n + 1) : ℝ) ≠ 0 := by
    exact_mod_cast Nat.factorial_ne_zero (n + 1)
  have hz : z ^ (n + 1) = 0 := by
    rcases div_eq_zero_iff.mp hd with h | h
    · exact h
    · exact absurd h hfac
  have : z = 0 := pow_eq_zero_iff (Nat.succ_ne_zero n) |>.mp hz
  rw [this] at h0
  rw [trExp_zero_arg] at h0
  exact one_ne_zero h0

/-- No real number is a root of every truncated exponential. -/
theorem trExp_exists_nonvanishing (z : ℝ) : ∃ n : ℕ, trExp n z ≠ 0 := by
  by_cases h : trExp 0 z = 0
  · exact ⟨1, fun h1 => trExp_no_common_root 0 z h h1⟩
  · exact ⟨0, h⟩

end GMC2Hermite
