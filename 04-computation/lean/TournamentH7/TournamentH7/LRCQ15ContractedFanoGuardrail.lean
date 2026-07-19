/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Finite core for THM-1247 q=15 contracted-Fano guardrail

This module kernel-checks the six q=15 masks, two explicit invariant Fano
planes and their common contraction, and the blocker-complete lonely witness.
The hand proofs classify the two invariant planes and establish CRT freedom of
the speed `chi_7` colouring.
-/

namespace LRC14
namespace Q15ContractedFanoGuardrail

def normNum (q v p : ℕ) : ℕ :=
  min ((v * p) % q) (q - (v * p) % q)

def dangerMask15 (v : ℕ) : Finset (Fin 15) :=
  Finset.univ.filter fun p => 14 * normNum 15 v p.val < 15

def sixMaskClasses : Finset (Finset (Fin 15)) :=
  (Finset.univ.erase (0 : Fin 15)).image fun s => dangerMask15 s.val

theorem exactly_six_mask_classes : sixMaskClasses.card = 6 := by
  decide

def petalWord : List (Fin 6) :=
  [0, 1, 2, 3, 4, 2, 5, 5, 2, 4, 3, 2, 1, 0]

theorem petal_word_palindromic : petalWord.reverse = petalWord := by
  decide

def fano0 : Finset (Finset (Fin 7)) :=
  {{0, 1, 2}, {0, 3, 4}, {0, 5, 6}, {1, 3, 5},
   {1, 4, 6}, {2, 3, 6}, {2, 4, 5}}

def fano1 : Finset (Finset (Fin 7)) :=
  {{0, 1, 5}, {0, 2, 6}, {0, 3, 4}, {1, 2, 3},
   {1, 4, 6}, {2, 4, 5}, {3, 5, 6}}

def timesTwo : Fin 7 → Fin 7
  | 0 => 1
  | 1 => 3
  | 2 => 5
  | 3 => 6
  | 4 => 4
  | 5 => 2
  | 6 => 0

def transportPlane (P : Finset (Finset (Fin 7))) : Finset (Finset (Fin 7)) :=
  P.image fun line => line.image timesTwo

theorem fano_planes_pair_unique :
    (∀ i j : Fin 7, i ≠ j →
      (fano0.filter fun line => i ∈ line ∧ j ∈ line).card = 1) ∧
    (∀ i j : Fin 7, i ≠ j →
      (fano1.filter fun line => i ∈ line ∧ j ∈ line).card = 1) := by
  decide

theorem fano_planes_invariant :
    transportPlane fano0 = fano0 ∧ transportPlane fano1 = fano1 := by
  decide

def contractPoint : Fin 7 → Fin 6
  | 0 => 0
  | 1 => 1
  | 2 => 2
  | 3 => 3
  | 4 => 4
  | 5 => 2
  | 6 => 5

def contractedPlane (P : Finset (Finset (Fin 7))) :
    Multiset (Multiset (Fin 6)) :=
  P.1.map fun line => line.1.map contractPoint

theorem invariant_planes_same_contraction :
    contractedPlane fano0 = contractedPlane fano1 := by
  decide

def negativeLine : Finset (Fin 7) := {2, 4, 5}

theorem negative_line_degenerates :
    negativeLine ∈ fano0 ∧ negativeLine ∈ fano1 ∧
      negativeLine.1.map contractPoint = ({2, 2, 4} : Multiset (Fin 6)) := by
  decide

def witnessSpeed : Fin 7 → ℕ
  | 0 => 1
  | 1 => 2
  | 2 => 7
  | 3 => 8
  | 4 => 18
  | 5 => 19
  | 6 => 20

theorem witness_q15_full_clock :
    ∀ p : Fin 15, ∃ i : Fin 7, p ∈ dangerMask15 (witnessSpeed i) := by
  decide

theorem defining_pair_mask_equal : dangerMask15 7 = dangerMask15 8 := by
  decide

theorem witness_lonely_at_twelve :
    ∀ i : Fin 7, 12 < 14 * normNum 12 (witnessSpeed i) 1 := by
  decide

def fastSpeed : Fin 6 → ℕ
  | 0 => 2
  | 1 => 7
  | 2 => 8
  | 3 => 18
  | 4 => 19
  | 5 => 20

def spokeQ : Fin 6 → ℕ
  | 0 => 3
  | 1 => 8
  | 2 => 9
  | 3 => 19
  | 4 => 20
  | 5 => 21

def spokeP : Fin 6 → ℕ
  | 0 => 2
  | 1 => 4
  | 2 => 5
  | 3 => 10
  | 4 => 10
  | 5 => 11

def spokeBlocker : Fin 6 → Fin 6
  | 0 => 3
  | 1 => 0
  | 2 => 3
  | 3 => 0
  | 4 => 0
  | 5 => 0

theorem every_carrier_spoke_deep_and_blocked :
    ∀ i : Fin 6,
      4 * normNum (spokeQ i) 1 (spokeP i) > spokeQ i ∧
      4 * normNum (spokeQ i) (fastSpeed i) (spokeP i) > spokeQ i ∧
      spokeBlocker i ≠ i ∧
      14 * normNum (spokeQ i) (fastSpeed (spokeBlocker i)) (spokeP i) <
        spokeQ i := by
  decide

structure SumBeatRow where
  left : Fin 6
  right : Fin 6
  numerator : ℕ
  blocker : Fin 6

def sumBeatRow : Fin 15 → SumBeatRow
  | 0 => ⟨0, 1, 1, 3⟩
  | 1 => ⟨0, 2, 1, 5⟩
  | 2 => ⟨0, 3, 2, 5⟩
  | 3 => ⟨0, 4, 3, 1⟩
  | 4 => ⟨0, 5, 3, 1⟩
  | 5 => ⟨1, 2, 3, 5⟩
  | 6 => ⟨1, 3, 3, 2⟩
  | 7 => ⟨1, 4, 13, 0⟩
  | 8 => ⟨1, 5, 3, 3⟩
  | 9 => ⟨2, 3, 11, 1⟩
  | 10 => ⟨2, 4, 3, 3⟩
  | 11 => ⟨2, 5, 3, 4⟩
  | 12 => ⟨3, 4, 5, 1⟩
  | 13 => ⟨3, 5, 4, 4⟩
  | 14 => ⟨4, 5, 5, 2⟩

theorem every_fast_sum_beat_positive_and_blocked :
    ∀ row : Fin 15,
      let R := sumBeatRow row
      let x := fastSpeed R.left
      let y := fastSpeed R.right
      let q := x + y
      let D := normNum q x R.numerator
      14 * R.numerator > q ∧ 14 * R.numerator < 13 * q ∧
      D = normNum q y R.numerator ∧ q < 14 * D ∧
      R.blocker ≠ R.left ∧ R.blocker ≠ R.right ∧
      14 * normNum q (fastSpeed R.blocker) R.numerator < q := by
  decide

#print axioms exactly_six_mask_classes
#print axioms petal_word_palindromic
#print axioms fano_planes_pair_unique
#print axioms fano_planes_invariant
#print axioms invariant_planes_same_contraction
#print axioms negative_line_degenerates
#print axioms witness_q15_full_clock
#print axioms defining_pair_mask_equal
#print axioms witness_lonely_at_twelve
#print axioms every_carrier_spoke_deep_and_blocked
#print axioms every_fast_sum_beat_positive_and_blocked

end Q15ContractedFanoGuardrail
end LRC14

