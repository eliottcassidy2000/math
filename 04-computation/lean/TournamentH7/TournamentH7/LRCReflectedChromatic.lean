/-
  TournamentH7.LRCReflectedChromatic -- finite graph core of the reflected
  signed-pair chromatic theorem candidate.

  The analytic signed-pair floors and the 3,003-body graph census remain
  exact-prose/Python providers.  This module kernel-checks the finite
  combinatorics after that census:

  * six vertices colored by five levels repeat on K6;
  * each of the two exceptional good graphs has chromatic number exactly four;
  * spread four does NOT force all five offsets to occur, witnessed by a
    proper four-color word on each exceptional graph.

  The last item is the formal guardrail against the rejected Delta=1 D=4
  shortcut.  No `native_decide` and no `sorry`.
-/
import Mathlib

namespace LonelyRunner.LRC14.ReflectedChromatic

/-- The three missing edges on the first exceptional body form the path
`2--3--4--5` in slot coordinates. -/
def pathBad (i j : Fin 6) : Prop :=
  (i = 2 ∧ j = 3) ∨ (i = 3 ∧ j = 2) ∨
  (i = 3 ∧ j = 4) ∨ (i = 4 ∧ j = 3) ∨
  (i = 4 ∧ j = 5) ∨ (i = 5 ∧ j = 4)

/-- The two missing edges on the second exceptional body form a matching. -/
def matchingBad (i j : Fin 6) : Prop :=
  (i = 2 ∧ j = 4) ∨ (i = 4 ∧ j = 2) ∨
  (i = 3 ∧ j = 5) ∨ (i = 5 ∧ j = 3)

def pathGood (i j : Fin 6) : Prop := i ≠ j ∧ ¬ pathBad i j
def matchingGood (i j : Fin 6) : Prop := i ≠ j ∧ ¬ matchingBad i j

def ProperOn {k : ℕ} (edge : Fin 6 → Fin 6 → Prop)
    (color : Fin 6 → Fin k) : Prop :=
  ∀ i j, edge i j → color i ≠ color j

/-- The K6 pigeonhole used on 3,001 bodies in the D=4 sector. -/
theorem six_into_five_repeat (color : Fin 6 → Fin 5) :
    ∃ i j, i ≠ j ∧ color i = color j := by
  obtain ⟨i, j, hne, heq⟩ :=
    Fintype.exists_ne_map_eq_of_card_lt color (by norm_num)
  exact ⟨i, j, hne, heq⟩

/-- The explicit K4s used by the exact graph census. -/
def pathClique : Fin 4 → Fin 6 := ![0, 1, 2, 4]
def matchingClique : Fin 4 → Fin 6 := ![0, 1, 2, 3]

theorem path_clique_good :
    ∀ i j, i ≠ j → pathGood (pathClique i) (pathClique j) := by
  decide

theorem matching_clique_good :
    ∀ i j, i ≠ j → matchingGood (matchingClique i) (matchingClique j) := by
  decide

/-- Neither exceptional graph admits a proper three-coloring: restrict to
the displayed K4 and apply four-into-three pigeonhole. -/
theorem path_not_three_colorable
    (color : Fin 6 → Fin 3) : ¬ ProperOn pathGood color := by
  intro hproper
  obtain ⟨i, j, hne, heq⟩ :=
    Fintype.exists_ne_map_eq_of_card_lt (fun k => color (pathClique k)) (by norm_num)
  exact (hproper _ _ (path_clique_good i j hne)) heq

theorem matching_not_three_colorable
    (color : Fin 6 → Fin 3) : ¬ ProperOn matchingGood color := by
  intro hproper
  obtain ⟨i, j, hne, heq⟩ :=
    Fintype.exists_ne_map_eq_of_card_lt (fun k => color (matchingClique k)) (by norm_num)
  exact (hproper _ _ (matching_clique_good i j hne)) heq

/-- Explicit proper four-colorings prove the matching upper bound. -/
def pathFourColor : Fin 6 → Fin 4 := ![0, 1, 2, 2, 3, 3]
def matchingFourColor : Fin 6 → Fin 4 := ![0, 1, 2, 3, 2, 3]

theorem path_four_colorable : ProperOn pathGood pathFourColor := by decide
theorem matching_four_colorable : ProperOn matchingGood matchingFourColor := by decide

/-- Proper D=4 words using only four of the five offsets.  These are the
minimal witnesses to the failed implication `spread = 4 -> all offsets occur`. -/
def pathSpreadFourWitness : Fin 6 → Fin 5 := ![0, 1, 2, 2, 4, 4]
def matchingSpreadFourWitness : Fin 6 → Fin 5 := ![0, 1, 2, 4, 2, 4]

theorem path_spread_four_witness_proper :
    ProperOn pathGood pathSpreadFourWitness := by decide

theorem matching_spread_four_witness_proper :
    ProperOn matchingGood matchingSpreadFourWitness := by decide

theorem path_spread_four_witness_values :
    pathSpreadFourWitness 0 = 0 ∧ pathSpreadFourWitness 1 = 1 ∧
    pathSpreadFourWitness 2 = 2 ∧ pathSpreadFourWitness 3 = 2 ∧
    pathSpreadFourWitness 4 = 4 ∧ pathSpreadFourWitness 5 = 4 := by
  decide

theorem matching_spread_four_witness_values :
    matchingSpreadFourWitness 0 = 0 ∧ matchingSpreadFourWitness 1 = 1 ∧
    matchingSpreadFourWitness 2 = 2 ∧ matchingSpreadFourWitness 3 = 4 ∧
    matchingSpreadFourWitness 4 = 2 ∧ matchingSpreadFourWitness 5 = 4 := by
  decide

#print axioms six_into_five_repeat
#print axioms path_not_three_colorable
#print axioms matching_not_three_colorable
#print axioms path_four_colorable
#print axioms matching_four_colorable
#print axioms path_spread_four_witness_proper
#print axioms matching_spread_four_witness_proper

end LonelyRunner.LRC14.ReflectedChromatic
