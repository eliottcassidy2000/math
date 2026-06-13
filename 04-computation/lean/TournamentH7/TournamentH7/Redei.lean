/-
  TournamentH7.Redei — Rédei's classical theorems on H(T)

  Reference:
    · Rédei, L. (1934), *Ein kombinatorischer Satz*, Acta Litt. Sci. Szeged
      7, 39–43.

  Rédei proved two classical results about the number of Hamiltonian paths
  in a tournament:

    (R1) **Existence.**  Every (finite) tournament has at least one
         Hamiltonian path. Equivalently `H(T) ≥ 1`.

    (R2) **Parity.**  The number of Hamiltonian paths in a tournament is
         odd. Equivalently `Odd (H T)`.

  Both are classical 1934 results, and are recorded as axioms here. A
  fully formal Lean proof of (R1) (induction by inserting a vertex into
  an existing path of T\{v}) is straightforward but lengthy; (R2) is the
  more subtle of the two and uses a *swap involution* on pairs of paths
  differing by an adjacent swap that creates / destroys an inversion.
-/

import TournamentH7.SCC
import Mathlib.Algebra.Ring.Parity

namespace Tournament

variable {n : ℕ}

/-! ### (R1) Rédei: every tournament has a Hamiltonian path

    Statement as axiom; for n = 0 the count of Hamiltonian paths is 1
    (the empty path on the empty vertex set) under standard conventions. -/

/-- **Axiom (R1 — Rédei 1934).** Every (nonempty) tournament has at least
    one directed Hamiltonian path:  `H(T) ≥ 1` for `n ≥ 1`. -/
axiom redei_existence {n : ℕ} (hn : 1 ≤ n) (T : Tournament n) : 1 ≤ H T

/-! ### (R2) Rédei: the count is odd -/

/-- **Axiom (R2 — Rédei 1934, parity).** The number of Hamiltonian paths
    in any tournament is odd. -/
axiom redei_parity {n : ℕ} (hn : 1 ≤ n) (T : Tournament n) : Odd (H T)

/-! ### Consequences -/

/-- H(T) ≠ 0 for nonempty tournaments. -/
theorem H_pos {n : ℕ} (hn : 1 ≤ n) (T : Tournament n) : H T ≠ 0 := by
  intro h
  have h1 : 1 ≤ H T := redei_existence hn T
  omega

/-- H(T) ≠ 2 for any tournament: by parity, H is odd. -/
theorem H_ne_two {n : ℕ} (hn : 1 ≤ n) (T : Tournament n) : H T ≠ 2 := by
  intro h
  have hodd : Odd (H T) := redei_parity hn T
  rw [h] at hodd
  -- Odd 2 is False
  exact (by decide : ¬ Odd 2) hodd

/-- H(T) ≠ 4, 6, 8, … — any even value. -/
theorem H_ne_even {n : ℕ} (hn : 1 ≤ n) (T : Tournament n) (m : ℕ) (hm : Even m) :
    H T ≠ m := by
  intro h
  have hodd : Odd (H T) := redei_parity hn T
  rw [h] at hodd
  -- Odd m and Even m gives False
  exact (Nat.not_odd_iff_even.mpr hm) hodd

end Tournament
