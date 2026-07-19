/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import TournamentH7.LRCBinaryPhaseWordLanding
import TournamentH7.LRCCenteredSurvivorProtrusion
import TournamentH7.LRCClosestReturnLeafPaidStar

/-!
# Centered-tail return-or-terminal arithmetic (THM-1274)

This module checks the abstract arithmetic and finite-label consumers behind
the protrusion-side continuation of a slowest-rooted blocker two-cycle.

The paper topology remains explicit: identification of the protruding safe
component, placement of the reverse marked tooth outside that component,
alignment of phase and chronological order, extraction of the actual
five-owner subword, and interpretation of private-tooth counts.  No covering
or interval-extraction assertion is hidden in a definition or axiom here.
-/

namespace LRC14
namespace CenteredTailReturnOrTerminal

/-! ## The protruding side fixes phase order and the binary descent digit -/

/-- If the safe component protrudes through the right endpoint of the carrier
gap, a reverse phase which lies in the carrier gap but outside the component
must lie to the left of the centered slow phase. -/
theorem right_protrusion_forces_reverse_phase_left
    {gLeft gRight sLeft sRight tSlow tReverse : ℚ}
    (hComponentRightOutside : gRight < sRight)
    (hSlowInside : sLeft < tSlow ∧ tSlow < sRight)
    (hReverseInGap : gLeft < tReverse ∧ tReverse < gRight)
    (hReverseOutsideComponent : tReverse < sLeft ∨ sRight < tReverse) :
    tReverse < tSlow := by
  rcases hReverseOutsideComponent with hleft | hright
  · linarith [hSlowInside.1]
  · linarith [hReverseInGap.2, hComponentRightOutside]

/-- The reflected statement: a left-protruding safe component forces the
reverse marked phase to the right of the centered slow phase. -/
theorem left_protrusion_forces_reverse_phase_right
    {gLeft gRight sLeft sRight tSlow tReverse : ℚ}
    (hComponentLeftOutside : sLeft < gLeft)
    (hSlowInside : sLeft < tSlow ∧ tSlow < sRight)
    (hReverseInGap : gLeft < tReverse ∧ tReverse < gRight)
    (hReverseOutsideComponent : tReverse < sLeft ∨ sRight < tReverse) :
    tSlow < tReverse := by
  rcases hReverseOutsideComponent with hleft | hright
  · linarith [hComponentLeftOutside, hReverseInGap.1]
  · linarith [hSlowInside.2]

/-- On a right protrusion the aligned descent difference is negative, so a
binary relative digit must be one.  In the application
`diff = Q * (tReverse - tSlow)`. -/
theorem right_protrusion_forces_descent_digit_one
    {gLeft gRight sLeft sRight tSlow tReverse Q theta : ℚ}
    {delta : ℤ}
    (hComponentRightOutside : gRight < sRight)
    (hSlowInside : sLeft < tSlow ∧ tSlow < sRight)
    (hReverseInGap : gLeft < tReverse ∧ tReverse < gRight)
    (hReverseOutsideComponent : tReverse < sLeft ∨ sRight < tReverse)
    (hQ : 0 < Q) (hTheta : 0 < theta)
    (hDifference : Q * (tReverse - tSlow) = theta - (delta : ℚ))
    (hBinary : delta = 0 ∨ delta = 1) :
    delta = 1 := by
  have hOrder := right_protrusion_forces_reverse_phase_left
    hComponentRightOutside hSlowInside hReverseInGap
    hReverseOutsideComponent
  rcases hBinary with rfl | rfl
  · norm_num at hDifference
    nlinarith
  · rfl

/-- On a left protrusion the aligned descent difference is positive, so a
binary relative digit must be zero. -/
theorem left_protrusion_forces_descent_digit_zero
    {gLeft gRight sLeft sRight tSlow tReverse Q theta : ℚ}
    {delta : ℤ}
    (hComponentLeftOutside : sLeft < gLeft)
    (hSlowInside : sLeft < tSlow ∧ tSlow < sRight)
    (hReverseInGap : gLeft < tReverse ∧ tReverse < gRight)
    (hReverseOutsideComponent : tReverse < sLeft ∨ sRight < tReverse)
    (hQ : 0 < Q) (hTheta : theta < 1)
    (hDifference : Q * (tReverse - tSlow) = theta - (delta : ℚ))
    (hBinary : delta = 0 ∨ delta = 1) :
    delta = 0 := by
  have hOrder := left_protrusion_forces_reverse_phase_right
    hComponentLeftOutside hSlowInside hReverseInGap
    hReverseOutsideComponent
  rcases hBinary with rfl | rfl
  · rfl
  · norm_num at hDifference
    nlinarith

/-! ## Five-owner closest returns -/

/-- A globally closest repeated pair in a word over five labels has at most
five intervening edges.  The word is represented by its index type; the
topological extraction of this owner map remains a paper input. -/
theorem closest_repeat_over_five_labels
    {N : ℕ} (owner : Fin N → Fin 5) (p q : Fin N)
    (hpq : p.val < q.val) (hRepeat : owner p = owner q)
    (hClosest : ∀ i j : Fin N, i.val < j.val → owner i = owner j →
      q.val - p.val ≤ j.val - i.val) :
    q.val - p.val ≤ 5 := by
  have hSelf := hClosest p q hpq hRepeat
  by_cases hSmall : N ≤ 5
  · have hp := p.isLt
    have hq := q.isLt
    omega
  · have hSix : 6 ≤ N := by omega
    let embed : Fin 6 → Fin N := fun i =>
      ⟨i.val, lt_of_lt_of_le i.isLt hSix⟩
    obtain ⟨i, j, hij, hOwner⟩ :=
      Fintype.exists_ne_map_eq_of_card_lt
        (fun i : Fin 6 => owner (embed i)) (by norm_num)
    have hval : i.val ≠ j.val := by
      intro h
      apply hij
      exact Fin.ext h
    rcases lt_or_gt_of_ne hval with hlt | hgt
    · have hBound := hClosest (embed i) (embed j) (by simpa [embed] using hlt) hOwner
      dsimp [embed] at hBound
      have hi := i.isLt
      have hj := j.isLt
      have hGap : j.val - i.val ≤ 5 := by omega
      exact hBound.trans hGap
    · have hBound := hClosest (embed j) (embed i) (by simpa [embed] using hgt) hOwner.symm
      dsimp [embed] at hBound
      have hi := i.isLt
      have hj := j.isLt
      have hGap : i.val - j.val ≤ 5 := by omega
      exact hBound.trans hGap

/-- The return-polygon inequality over at most five edges gives the strict
factor `3/2`, improving the six-owner `6/5` leaf when the slowest owner has
been excluded from the protrusion-side word. -/
theorem five_owner_return_forces_three_halves
    {a b ratioSum R : ℚ} {m : ℕ}
    (ha : 0 < a) (hb : 0 < b)
    (hmLower : 2 ≤ m) (hmUpper : m ≤ 5) (hR : 1 ≤ R)
    (hReturn : 7 * R - 1 < ratioSum)
    (hMinimum : ratioSum ≤ (m - 1 : ℕ) * (a / b)) :
    3 * b < 2 * a := by
  have hCount : m - 1 ≤ 4 := by omega
  have hCountQ : ((m - 1 : ℕ) : ℚ) ≤ 4 := by exact_mod_cast hCount
  have hRatioNonnegative : 0 ≤ a / b := by positivity
  have hCountBound : ((m - 1 : ℕ) : ℚ) * (a / b) ≤ 4 * (a / b) := by
    exact mul_le_mul_of_nonneg_right hCountQ hRatioNonnegative
  have hSix : (6 : ℚ) ≤ 7 * R - 1 := by linarith
  have hRatio : (6 : ℚ) < 4 * (a / b) :=
    lt_of_le_of_lt hSix (lt_of_lt_of_le hReturn (hMinimum.trans hCountBound))
  have hCancel : (a / b) * b = a := by field_simp
  nlinarith

/-- A repetition-free chronological word over the five non-slowest labels
has at most five tooth occurrences. -/
theorem nodup_five_owner_word_length (word : List (Fin 5))
    (hNodup : word.Nodup) :
    word.length ≤ 5 := by
  simpa using hNodup.length_le_card

/-! ## Six-owner terminal path and private-count rigidity -/

/-- A repetition-free path containing six occurrences already contains every
one of the six labels, so any additional occurrence repeats a path label. -/
theorem six_label_path_plus_extra_repeats
    (word : List (Fin 6)) (hNodup : word.Nodup)
    (hLength : word.length = 6) (extra : Fin 6) :
    extra ∈ word := by
  by_contra hMissing
  have hExtended : (extra :: word).Nodup := by
    simp [hMissing, hNodup]
  have hBound := hExtended.length_le_card
  norm_num [hLength] at hBound

/-- Seven occurrences over six labels cannot be repetition-free.  Selecting
a closest repeated pair is the remaining finite combinatorial operation. -/
theorem seven_occurrences_over_six_labels_have_repeat
    (word : List (Fin 6)) (hLength : 7 ≤ word.length) :
    ¬word.Nodup := by
  intro hNodup
  have hBound := hNodup.length_le_card
  norm_num at hBound
  omega

/-- Six positive private-tooth counts with total at most six are all exactly
one. -/
theorem six_private_counts_are_one
    (count : Fin 6 → ℕ)
    (hPositive : ∀ i, 1 ≤ count i)
    (hTotal : (∑ i, count i) ≤ 6) :
    ∀ i, count i = 1 := by
  intro i
  have hLower : (∑ j : Fin 6, (1 : ℕ)) ≤ ∑ j, count j := by
    exact Finset.sum_le_sum fun j _ => hPositive j
  have hOnes : (∑ _j : Fin 6, (1 : ℕ)) = 6 := by norm_num
  have hSumEquality : (∑ j : Fin 6, (1 : ℕ)) = ∑ j, count j := by
    omega
  have hPoint : (1 : ℕ) = count i :=
    (Finset.sum_eq_sum_iff_of_le (fun j _ => hPositive j)).mp
      hSumEquality i (Finset.mem_univ i)
  exact hPoint.symm

/-- Abstract ceiling-count consumer.  The geometric count provider only has
to show that a speed above `7c` would create a second private tooth.  Count
rigidity then forces every one of the six speeds below that threshold. -/
theorem private_count_rigidity_forces_speed_cap
    (c : ℕ) (speed count : Fin 6 → ℕ)
    (hPositive : ∀ i, 1 ≤ count i)
    (hTotal : (∑ i, count i) ≤ 6)
    (hLargeHasTwo : ∀ i, 7 * c < speed i → 2 ≤ count i) :
    ∀ i, speed i ≤ 7 * c := by
  have hUnit := six_private_counts_are_one count hPositive hTotal
  intro i
  by_contra hLarge
  have hTwo := hLargeHasTwo i (by omega)
  have hOne := hUnit i
  omega

#print axioms right_protrusion_forces_reverse_phase_left
#print axioms left_protrusion_forces_reverse_phase_right
#print axioms right_protrusion_forces_descent_digit_one
#print axioms left_protrusion_forces_descent_digit_zero
#print axioms closest_repeat_over_five_labels
#print axioms five_owner_return_forces_three_halves
#print axioms nodup_five_owner_word_length
#print axioms six_label_path_plus_extra_repeats
#print axioms seven_occurrences_over_six_labels_have_repeat
#print axioms six_private_counts_are_one
#print axioms private_count_rigidity_forces_speed_cap

end CenteredTailReturnOrTerminal
end LRC14
