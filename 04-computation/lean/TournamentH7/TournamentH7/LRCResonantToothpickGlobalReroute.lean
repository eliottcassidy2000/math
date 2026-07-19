/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Arithmetic core for THM-1243 resonant-toothpick global reroute

The paper theorem reduces the thirteen circle distances to two parity ledgers.
This module checks the master-clock congruence, every ledger lower bound, the
strict `3/28` depth, and the separation from the locally erased carrier gap.
-/

namespace LRC14
namespace ResonantToothpickGlobalReroute

/-- On the master clock `q=2a+7`, an even numerator `p=2s` turns the eight
consecutive speeds `a+r` into the odd-multiplier word `(2r-7)s`. -/
theorem consecutive_packet_modEq
    (a q p s r : ℤ) (hq : q = 2 * a + 7) (hp : p = 2 * s) :
    (a + r) * p ≡ (2 * r - 7) * s [ZMOD q] := by
  rw [Int.modEq_iff_dvd]
  refine ⟨-s, ?_⟩
  rw [hq, hp]
  ring

/-- The resonant blocker `14a` is the terminal multiplier `-98s`. -/
theorem terminal_blocker_modEq
    (a q p s : ℤ) (hq : q = 2 * a + 7) (hp : p = 2 * s) :
    14 * a * p ≡ -98 * s [ZMOD q] := by
  rw [Int.modEq_iff_dvd]
  refine ⟨-14 * s, ?_⟩
  rw [hq, hp]
  ring

/-- The second THM-1239 blocker `7a+4` is the terminal multiplier `-41s`. -/
theorem second_blocker_modEq
    (a q p s : ℤ) (hq : q = 2 * a + 7) (hp : p = 2 * s) :
    (7 * a + 4) * p ≡ -41 * s [ZMOD q] := by
  rw [Int.modEq_iff_dvd]
  refine ⟨-7 * s, ?_⟩
  rw [hq, hp]
  ring

/-- Every distinct least-residue numerator in the even branch is at least
the central value `s=3h+2`. -/
theorem even_ledger_lower
    {h x : ℤ} (hh : 14 ≤ h)
    (hx : x = 6 * h + 4 ∨ x = 12 * h + 8 ∨ x = 10 * h - 3 ∨
      x = 4 * h - 7 ∨ x = 7 * h - 5 ∨ x = 13 * h - 1 ∨
      x = 9 * h + 6 ∨ x = 3 * h + 2 ∨ x = 14 * h - 97 ∨
      x = 11 * h + 46) :
    3 * h + 2 ≤ x := by
  rcases hx with h₁ | h₂ | h₃ | h₄ | h₅ | h₆ | h₇ | h₈ | h₉ | h₁₀
  all_goals subst x <;> omega

/-- Every distinct least-residue numerator in the odd branch is at least
the central value `s=3h+4`. -/
theorem odd_ledger_lower
    {h x : ℤ} (hh : 13 ≤ h)
    (hx : x = 6 * h + 8 ∨ x = 12 * h + 16 ∨ x = 10 * h - 1 ∨
      x = 4 * h - 9 ∨ x = 7 * h - 5 ∨ x = 13 * h + 3 ∨
      x = 9 * h + 12 ∨ x = 3 * h + 4 ∨ x = 14 * h - 139 ∨
      x = 11 * h + 72 ∨ x = 17 * h - 49) :
    3 * h + 4 ≤ x := by
  rcases hx with h₁ | h₂ | h₃ | h₄ | h₅ | h₆ | h₇ | h₈ | h₉ | h₁₀ | h₁₁
  all_goals subst x <;> omega

/-- Even rows have exact numerator excess `29` above depth `3/28`. -/
theorem even_depth_excess (h : ℤ) :
    28 * (3 * h + 2) - 3 * (28 * h + 9) = 29 := by
  ring

/-- Odd rows have exact numerator excess `43` above depth `3/28`. -/
theorem odd_depth_excess (h : ℤ) :
    28 * (3 * h + 4) - 3 * (28 * h + 23) = 43 := by
  ring

/-- The even parity certificate is strictly deeper than `3/28`, and hence
than the LRC(14) threshold. -/
theorem even_depth_certificate {h : ℤ} (hh : 14 ≤ h) :
    3 * (28 * h + 9) < 28 * (3 * h + 2) ∧
      28 * h + 9 < 14 * (3 * h + 2) := by
  constructor <;> omega

/-- The odd parity certificate is strictly deeper than `3/28`, and hence
than the LRC(14) threshold. -/
theorem odd_depth_certificate {h : ℤ} (hh : 13 ≤ h) :
    3 * (28 * h + 23) < 28 * (3 * h + 4) ∧
      28 * h + 23 < 14 * (3 * h + 4) := by
  constructor <;> omega

/-- Both parity numerators lie strictly to the right of phase `1/5`. -/
theorem parity_phase_right_of_one_fifth {m e : ℤ}
    (hm : 0 ≤ m) (he : 0 ≤ e) :
    14 * m + 9 < 5 * (3 * m + 4 + e) := by
  omega

/-- The THM-1239 selected carrier gap lies strictly to the left of `1/6`
once `m≥5`. -/
theorem selected_gap_left_of_one_sixth {m : ℤ} (hm : 5 ≤ m) :
    6 * (14 * m + 13) < 14 * (7 * m + 1) := by
  omega

#print axioms consecutive_packet_modEq
#print axioms terminal_blocker_modEq
#print axioms second_blocker_modEq
#print axioms even_ledger_lower
#print axioms odd_ledger_lower
#print axioms even_depth_excess
#print axioms odd_depth_excess
#print axioms even_depth_certificate
#print axioms odd_depth_certificate
#print axioms parity_phase_right_of_one_fifth
#print axioms selected_gap_left_of_one_sixth

end ResonantToothpickGlobalReroute
end LRC14
