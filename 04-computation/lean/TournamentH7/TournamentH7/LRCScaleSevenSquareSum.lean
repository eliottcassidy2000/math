/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: codex-c7-lean (LRC multi-agent project, 2026-07-17)
-/
import Mathlib

/-!
# The scale-seven square-sum obstruction

This file isolates the terminal algebraic contradiction in the scale-`c = 7`,
Hamming-six branch.  Let `z₀, ..., z₅` be nonzero squares in `ZMod 7`, let

`S = z₀ + ... + z₅`,

and suppose each opposite pair has the same prescribed sum

`2 S = zᵢ + zᵢ₊₃`.

Adding the equations for `i = 0, 1, 2` gives `6 S = S`, hence `S = 0` in
characteristic seven.  Thus every opposite pair sums to zero.  But two nonzero
squares modulo seven cannot sum to zero: this would make `-1` a square, while
the nonzero squares are exactly `1, 2, 4`.

The first theorem below keeps the linear part independent of the square
assumption, so it can be reused by later reductions.  The only finite check is
the seven-element residue lemma `nonzero_squares_do_not_sum_to_zero`, proved by
kernel `decide`.  No `sorry`; no `native_decide`.
-/

namespace LonelyRunner
namespace ScaleSevenSquareSum

/-- The index opposite `i` on the cyclic six-set. -/
def antipode (i : Fin 6) : Fin 6 :=
  ⟨((i : ℕ) + 3) % 6, by omega⟩

/-- Two nonzero squares in `ZMod 7` cannot add to zero. -/
theorem nonzero_squares_do_not_sum_to_zero :
    ∀ a b : ZMod 7,
      a ≠ 0 → IsSquare a → b ≠ 0 → IsSquare b → a + b ≠ 0 := by
  decide

/-- The three opposite-pair equations already force the total sum to vanish.
No quadratic-residue hypothesis is needed for this linear step. -/
theorem sum_eq_zero_of_opposite_pair_equations
    (z : Fin 6 → ZMod 7)
    (hpair : ∀ i : Fin 6,
      2 * (∑ j : Fin 6, z j) = z i + z (antipode i)) :
    (∑ j : Fin 6, z j) = 0 := by
  let S : ZMod 7 := ∑ j : Fin 6, z j
  have h0 := hpair (0 : Fin 6)
  have h1 := hpair (1 : Fin 6)
  have h2 := hpair (2 : Fin 6)
  change 2 * S = z 0 + z 3 at h0
  change 2 * S = z 1 + z 4 at h1
  change 2 * S = z 2 + z 5 at h2
  have hS : S = z 0 + z 1 + z 2 + z 3 + z 4 + z 5 := by
    simp [S, Fin.sum_univ_succ]
    ring
  have hfive : (5 : ZMod 7) * S = 0 := by
    linear_combination h0 + h1 + h2 - hS
  have hunit : (3 : ZMod 7) * 5 = 1 := by decide
  calc
    S = ((3 : ZMod 7) * 5) * S := by rw [hunit]; ring
    _ = 3 * (5 * S) := by ring
    _ = 0 := by rw [hfive]; ring

/-- There is no six-tuple of nonzero squares in `ZMod 7` satisfying all six
opposite-pair equations `2 S = zᵢ + zᵢ₊₃`. -/
theorem no_nonzero_square_solution
    (z : Fin 6 → ZMod 7)
    (hnz : ∀ i : Fin 6, z i ≠ 0)
    (hsq : ∀ i : Fin 6, IsSquare (z i))
    (hpair : ∀ i : Fin 6,
      2 * (∑ j : Fin 6, z j) = z i + z (antipode i)) :
    False := by
  have hS := sum_eq_zero_of_opposite_pair_equations z hpair
  have h03 := hpair (0 : Fin 6)
  change 2 * (∑ j : Fin 6, z j) = z 0 + z 3 at h03
  rw [hS] at h03
  norm_num at h03
  exact nonzero_squares_do_not_sum_to_zero (z 0) (z 3)
    (hnz 0) (hsq 0) (hnz 3) (hsq 3) h03.symm

end ScaleSevenSquareSum
end LonelyRunner
