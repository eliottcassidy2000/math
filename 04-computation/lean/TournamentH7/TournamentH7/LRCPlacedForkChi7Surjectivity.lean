/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Arithmetic core for THM-1260's placed-fork / chi7 guardrail

The paper theorem constructs a local marked blocker path and sharp two-wall
toothpick fork for every compact rung, binary digit, and seven-label speed
colour word.  This module checks the four decisive residue rows, the exact
rung and seam identities, the half-gap binary address split, and the
continued exact blocker clock.
-/

namespace LRC14
namespace PlacedForkChi7Surjectivity

def chi7Residue (n : ℕ) : ℤ :=
  if n % 7 = 0 then 0
  else if n % 7 = 1 ∨ n % 7 = 2 ∨ n % 7 = 4 then 1
  else -1

/-- The sharp normalized seam equation `H = (7r-1)D+1` realizes all four
ordered target/provider quadratic-colour pairs already modulo seven. -/
theorem sharp_pair_colour_table :
    (chi7Residue 4, chi7Residue (1 + 7 - 4)) = (1, 1) ∧
    (chi7Residue 2, chi7Residue (1 + 7 - 2)) = (1, -1) ∧
    (chi7Residue 6, chi7Residue (1 + 7 - 6)) = (-1, 1) ∧
    (chi7Residue 3, chi7Residue (1 + 7 - 3)) = (-1, -1) := by
  norm_num [chi7Residue]

/-- The provider is one gcd-sheet above the resonant toothpick multiplier. -/
theorem sharp_rung_detuning (r D sheet : ℤ) :
    sheet * ((7 * r - 1) * D + 1) -
      (7 * r - 1) * (sheet * D) = sheet := by
  ring

/-- Left and right wall overlaps in the symmetric placed fork are equal to
the sharp gcd/lcm quantum. -/
theorem sharp_left_wall_quantum
    {j h r sheet : ℚ} (hj : j ≠ 0) (hh : h ≠ 0)
    (hdef : h = (7 * r - 1) * j + sheet) :
    (1 / 2 - (7 * r - 1) / (14 * h)) -
        (1 / 2 - 1 / (14 * j)) = sheet / (14 * j * h) := by
  field_simp [hj, hh]
  linarith

theorem sharp_right_wall_quantum
    {j h r sheet : ℚ} (hj : j ≠ 0) (hh : h ≠ 0)
    (hdef : h = (7 * r - 1) * j + sheet) :
    (1 / 2 + 1 / (14 * j)) -
        (1 / 2 + (7 * r - 1) / (14 * h)) =
      sheet / (14 * j * h) := by
  field_simp [hj, hh]
  linarith

/-- For odd carrier `2C+1` and even target `2J`, the two nearest centered
target numerators are exactly the binary lower/upper addresses. -/
theorem binary_half_gap_address_split (C J : ℤ) :
    2 * (C + J) - ((2 * C + 1) + 2 * J) = -1 ∧
    2 * (C + J + 1) - ((2 * C + 1) + 2 * J) = 1 := by
  constructor <;> ring

/-- A multiple of the target centered clock is an exact next blocker at the
same marked phase. -/
theorem continued_blocker_clock (m P Q : ℤ) :
    (m * Q) * P = (m * P) * Q := by
  ring

/-- The worst compact rung follows from the same margin used in the paper
construction; smaller rungs only decrease the provider. -/
theorem terminal_compact_margin :
    (7 * 335 - 1) * 100004 + 2 < 2345 * (100004 - 13) := by
  norm_num

#print axioms sharp_pair_colour_table
#print axioms sharp_rung_detuning
#print axioms sharp_left_wall_quantum
#print axioms sharp_right_wall_quantum
#print axioms binary_half_gap_address_split
#print axioms continued_blocker_clock
#print axioms terminal_compact_margin

end PlacedForkChi7Surjectivity
end LRC14
