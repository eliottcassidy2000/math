/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Arithmetic consumers for THM-1242 seven-wall beat sunflower law

This module checks the common-zero cardinality threshold, its exact
seven/eight-wall criticality, and the concrete six-petal full clock on
`Z/15Z`.  The speed-to-mask reduction is supplied by THM-1217.
-/

namespace LRC14
namespace SevenWallBeatSunflower

/-- Six masks of size at most `A`, with the five units of overlap supplied by
a common point, have union size at most `6A-5`. -/
theorem six_common_zero_union_cap
    {A U s₁ s₂ s₃ s₄ s₅ s₆ : ℕ}
    (hcredit : U + 5 ≤ s₁ + s₂ + s₃ + s₄ + s₅ + s₆)
    (h₁ : s₁ ≤ A) (h₂ : s₂ ≤ A) (h₃ : s₃ ≤ A)
    (h₄ : s₄ ≤ A) (h₅ : s₅ ≤ A) (h₆ : s₆ ≤ A) :
    U ≤ 6 * A - 5 := by
  omega

/-- If `h=ceil(Q/14)`, encoded by the lower shell inequality, the seven-wall
threshold `12h-10` never exceeds `Q` for `Q≥2`. -/
theorem seven_wall_threshold_all_Q
    {Q h : ℕ} (hQ : 2 ≤ Q) (hh : 1 ≤ h)
    (hlower : 14 * h ≤ Q + 13) :
    12 * h ≤ Q + 10 := by
  omega

/-- At eight walls the analogous common-zero threshold misses the first
modulus of every noninitial ceiling shell by exactly one. -/
theorem eight_wall_first_shell_failure
    {h Q : ℕ} (hh : 2 ≤ h) (hQ : Q = 14 * h - 13) :
    14 * h - 12 = Q + 1 := by
  omega

/-- The common-zero cap is one below the escape threshold. -/
theorem escape_threshold_from_union_cap
    {A U : ℕ} (hA : 1 ≤ A) (hU : U ≤ 6 * A - 5) :
    U < 6 * A - 4 := by
  omega

def q15Pair : Finset (Fin 15) := {0, 1, 14}
def q15Five : Finset (Fin 15) := {0, 3, 6, 9, 12}
def q15Three : Finset (Fin 15) := {0, 5, 10}
def q15Two : Finset (Fin 15) := {0, 7, 8}
def q15Four : Finset (Fin 15) := {0, 4, 11}
def q15Seven : Finset (Fin 15) := {0, 2, 13}

def q15Mask : Fin 6 → Finset (Fin 15)
  | 0 => q15Pair
  | 1 => q15Five
  | 2 => q15Three
  | 3 => q15Two
  | 4 => q15Four
  | 5 => q15Seven

/-- The six masks cover the complete master clock. -/
theorem q15_full_clock :
    q15Pair ∪ q15Five ∪ q15Three ∪ q15Two ∪ q15Four ∪ q15Seven =
      Finset.univ := by
  decide

/-- Distinct petals meet exactly in their common zero. -/
theorem q15_sunflower_intersections :
    ∀ i j : Fin 6, i ≠ j → q15Mask i ∩ q15Mask j = {0} := by
  decide

/-- Exact mask-size and saturated Hunter-tree ledger. -/
theorem q15_cardinality_ledger :
    q15Pair.card = 3 ∧ q15Five.card = 5 ∧ q15Three.card = 3 ∧
      q15Two.card = 3 ∧ q15Four.card = 3 ∧ q15Seven.card = 3 ∧
      q15Pair.card + q15Five.card + q15Three.card + q15Two.card +
        q15Four.card + q15Seven.card - 5 = 15 := by
  decide

/-- The nonzero petal sizes partition the fourteen nonzero residues. -/
theorem q15_nonzero_petal_sizes :
    q15Pair.card - 1 = 2 ∧ q15Five.card - 1 = 4 ∧
      q15Three.card - 1 = 2 ∧ q15Two.card - 1 = 2 ∧
      q15Four.card - 1 = 2 ∧ q15Seven.card - 1 = 2 ∧
      (q15Pair.card - 1) + (q15Five.card - 1) +
        (q15Three.card - 1) + (q15Two.card - 1) +
        (q15Four.card - 1) + (q15Seven.card - 1) = 14 := by
  decide

#print axioms six_common_zero_union_cap
#print axioms seven_wall_threshold_all_Q
#print axioms eight_wall_first_shell_failure
#print axioms escape_threshold_from_union_cap
#print axioms q15_full_clock
#print axioms q15_sunflower_intersections
#print axioms q15_cardinality_ledger
#print axioms q15_nonzero_petal_sizes

end SevenWallBeatSunflower
end LRC14
