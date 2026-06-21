/-
  TournamentH7.LRCFactorialAtom -- finite factorial atom identities for LRC14.

  This module isolates the algebra behind HYP-2720..HYP-2728.  It does not
  formalize analytic LRC estimates; it proves the seven-coordinate binomial
  identities that make the `Q_0` origin atom an alternating factorial boundary.
-/

namespace LonelyRunner
namespace FactorialAtom

/-- LRC14 has six inner sectors, so the missed-count atom has coordinates 0..6. -/
def atomCount : Nat := 7

/-- A signed missed-count atom vector `q_t`, `t=0..6`. -/
abbrev Atom := Fin atomCount -> Int

/-- Tiny local binomial coefficient, avoiding a Mathlib dependency for this
finite audit module. -/
def choose : Nat -> Nat -> Nat
  | _, 0 => 1
  | 0, _ + 1 => 0
  | n + 1, k + 1 => choose n k + choose n (k + 1)

/-- Natural binomial coefficients, cast to `Int`. -/
def chooseInt (n k : Nat) : Int := (choose n k : Int)

/-- Sum over the seven atom coordinates. -/
def sum7 (f : Nat -> Int) : Int :=
  f 0 + f 1 + f 2 + f 3 + f 4 + f 5 + f 6

/-- The alternating sign `(-1)^k` as an integer. -/
def altSign (k : Nat) : Int := if k % 2 = 0 then 1 else -1

/-- The finite-difference atom packet created by a unit factorial moment `W_j`. -/
def basis (j : Fin atomCount) : Atom :=
  fun t =>
    if t.val <= j.val then altSign (j.val - t.val) * chooseInt j.val t.val else 0

/-- Factorial moment readout `W_j(q)=sum_{t>=j} binom(t,j) q_t`. -/
def moment (q : Atom) (j : Fin atomCount) : Int :=
  sum7 fun t =>
    if h : t < atomCount then
      if j.val <= t then chooseInt t j.val * q ⟨t, h⟩ else 0
    else
      0

/-- The origin atom coordinate. -/
def q0 (q : Atom) : Int := q ⟨0, by decide⟩

/-- Coefficient of atom coordinate `q_t` in the alternating boundary of moments. -/
def originCoeff (t : Fin atomCount) : Int :=
  sum7 fun j =>
    if j <= t.val then altSign j * chooseInt t.val j else 0

/-- The Bonferroni4 high-tail readout `U4=q0+q5+5q6`. -/
def U4 (q : Atom) : Int :=
  q ⟨0, by decide⟩ + q ⟨5, by decide⟩ + 5 * q ⟨6, by decide⟩

/-- Low leakage readout `W1+W2` used by HYP-2722/HYP-2726. -/
def low12 (q : Atom) : Int :=
  moment q ⟨1, by decide⟩ + moment q ⟨2, by decide⟩

/-- The finite-difference packets are exactly dual to the factorial moments. -/
theorem basis_moment_delta :
    ∀ i j : Fin atomCount, moment (basis j) i = if i = j then 1 else 0 := by
  native_decide

/-- The `Q_0` coordinate is the alternating boundary of factorial moments. -/
theorem originCoeff_delta :
    ∀ t : Fin atomCount, originCoeff t = if t.val = 0 then 1 else 0 := by
  native_decide

/-- On a unit packet, `Q_0(B_j)=(-1)^j`. -/
theorem basis_q0_sign :
    ∀ j : Fin atomCount, q0 (basis j) = altSign j.val := by
  native_decide

/-- `U4` sees packets through degree four and is blind to `B_5,B_6`. -/
theorem U4_basis :
    ∀ j : Fin atomCount, U4 (basis j) = if j.val <= 4 then altSign j.val else 0 := by
  native_decide

/-- The low `W1+W2` leakage sees exactly the `B_1` and `B_2` packet axes. -/
theorem low12_basis :
    ∀ j : Fin atomCount, low12 (basis j) =
      if j.val = 1 ∨ j.val = 2 then 1 else 0 := by
  native_decide

/-! ### Axiom audit

These are closed finite identities over seven coordinates.  In a working Lean
environment, the axiom output should be empty or contain only Lean foundations.
-/

#print axioms basis_moment_delta
#print axioms originCoeff_delta
#print axioms basis_q0_sign
#print axioms U4_basis
#print axioms low12_basis

end FactorialAtom
end LonelyRunner
