import TournamentH7.RatIntervals

/-!
Exact rational overlap for two normalized circular intervals.  The proof
handles all start-order and wrap/no-wrap cases and identifies the quadratic
clip sum with the two-tile overlap tent.  Comb specialization and cyclic
reindexing deliberately live in a downstream module.
-/

namespace LonelyRunner
namespace LRCB5WrappedToothOverlap

open RatIntervals

set_option maxHeartbeats 300000

abbrev Interval := ℚ × ℚ

/-- One half-open circular interval with normalized start and length. -/
def circleInterval (start length : ℚ) : Region :=
  if start + length ≤ 1 then
    [(start, start + length)]
  else
    [(start, 1), (0, start + length - 1)]

/-- Rational version of the two-tile overlap tent. -/
def pairOverlapQ (length₁ length₂ shift : ℚ) : ℚ :=
  max 0 (min length₁ (length₂ - shift)) +
    max 0 (min length₁ (1 - shift + length₂) - (1 - shift))

theorem pairOverlapQ_swap_pos
    (length₁ length₂ shift : ℚ)
    (hlength₁0 : 0 ≤ length₁) (hlength₁1 : length₁ ≤ 1)
    (hlength₂0 : 0 ≤ length₂) (hlength₂1 : length₂ ≤ 1)
    (hshift0 : 0 < shift) (hshift1 : shift < 1) :
    pairOverlapQ length₁ length₂ shift =
      pairOverlapQ length₂ length₁ (1 - shift) := by
  unfold pairOverlapQ
  simp only [max_def, min_def]
  split_ifs <;> linarith

theorem pairOverlapQ_swap_zero
    (length₁ length₂ : ℚ)
    (hlength₁0 : 0 ≤ length₁) (hlength₁1 : length₁ ≤ 1)
    (hlength₂0 : 0 ≤ length₂) (hlength₂1 : length₂ ≤ 1) :
    pairOverlapQ length₁ length₂ 0 =
      pairOverlapQ length₂ length₁ 0 := by
  unfold pairOverlapQ
  simp only [sub_zero]
  rw [max_eq_right (le_min hlength₁0 hlength₂0),
    max_eq_left (by linarith [min_le_left length₁ (1 + length₂)]),
    max_eq_right (le_min hlength₂0 hlength₁0),
    max_eq_left (by linarith [min_le_left length₂ (1 + length₁)]),
    min_comm]

theorem fract_sub_eq_of_le
    {first second : ℚ}
    (hfirst0 : 0 ≤ first) (hfirst1 : first < 1)
    (hsecond0 : 0 ≤ second) (hsecond1 : second < 1)
    (hsecondFirst : second ≤ first) :
    Int.fract (first - second) = first - second := by
  apply Int.fract_eq_self.mpr
  constructor <;> linarith

theorem fract_sub_eq_add_one_of_lt
    {first second : ℚ}
    (hfirst0 : 0 ≤ first) (hfirst1 : first < 1)
    (hsecond0 : 0 ≤ second) (hsecond1 : second < 1)
    (hfirstSecond : first < second) :
    Int.fract (first - second) = first - second + 1 := by
  unfold Int.fract
  have hfloor : ⌊first - second⌋ = -1 := by
    rw [Int.floor_eq_iff]
    push_cast
    constructor <;> linarith
  rw [hfloor]
  push_cast
  ring

theorem overlap_le_nowrap_nowrap
    (start₁ start₂ length₁ length₂ : ℚ)
    (hstart₁0 : 0 ≤ start₁) (hstart₁1 : start₁ < 1)
    (hstart₂0 : 0 ≤ start₂) (hstart₂1 : start₂ < 1)
    (hlength₁0 : 0 ≤ length₁) (hlength₁1 : length₁ ≤ 1)
    (hlength₂0 : 0 ≤ length₂) (hlength₂1 : length₂ ≤ 1)
    (horder : start₂ ≤ start₁)
    (hwrap₁ : start₁ + length₁ ≤ 1)
    (hwrap₂ : start₂ + length₂ ≤ 1) :
    max 0 (min (start₁ + length₁) (start₂ + length₂) - max start₁ start₂) =
      pairOverlapQ length₁ length₂ (start₁ - start₂) := by
  unfold pairOverlapQ
  simp only [max_def, min_def]
  split_ifs <;> linarith

theorem overlap_le_wrap_nowrap
    (start₁ start₂ length₁ length₂ : ℚ)
    (hstart₁0 : 0 ≤ start₁) (hstart₁1 : start₁ < 1)
    (hstart₂0 : 0 ≤ start₂) (hstart₂1 : start₂ < 1)
    (hlength₁0 : 0 ≤ length₁) (hlength₁1 : length₁ ≤ 1)
    (hlength₂0 : 0 ≤ length₂) (hlength₂1 : length₂ ≤ 1)
    (horder : start₂ ≤ start₁)
    (hwrap₁ : ¬start₁ + length₁ ≤ 1)
    (hwrap₂ : start₂ + length₂ ≤ 1) :
    max 0 (min 1 (start₂ + length₂) - max start₁ start₂) +
        max 0 (min (start₁ + length₁ - 1) (start₂ + length₂) - max 0 start₂) =
      pairOverlapQ length₁ length₂ (start₁ - start₂) := by
  unfold pairOverlapQ
  simp only [max_def, min_def]
  split_ifs <;> linarith

theorem overlap_le_nowrap_wrap
    (start₁ start₂ length₁ length₂ : ℚ)
    (hstart₁0 : 0 ≤ start₁) (hstart₁1 : start₁ < 1)
    (hstart₂0 : 0 ≤ start₂) (hstart₂1 : start₂ < 1)
    (hlength₁0 : 0 ≤ length₁) (hlength₁1 : length₁ ≤ 1)
    (hlength₂0 : 0 ≤ length₂) (hlength₂1 : length₂ ≤ 1)
    (horder : start₂ ≤ start₁)
    (hwrap₁ : start₁ + length₁ ≤ 1)
    (hwrap₂ : ¬start₂ + length₂ ≤ 1) :
    max 0 (min (start₁ + length₁) 1 - max start₁ start₂) +
        max 0 (min (start₁ + length₁) (start₂ + length₂ - 1) - max start₁ 0) =
      pairOverlapQ length₁ length₂ (start₁ - start₂) := by
  unfold pairOverlapQ
  have hsum₂ : 1 < start₂ + length₂ := lt_of_not_ge hwrap₂
  have hoverflow₂le : start₂ + length₂ - 1 ≤ start₂ := by linarith
  have hcrossZero :
      min (start₁ + length₁) (start₂ + length₂ - 1) - start₁ ≤ 0 := by
    linarith [min_le_right (start₁ + length₁) (start₂ + length₂ - 1)]
  have hcover : length₁ ≤ length₂ - (start₁ - start₂) := by
    linarith
  have hsecondZero :
      min length₁ (1 - (start₁ - start₂) + length₂) -
          (1 - (start₁ - start₂)) ≤ 0 := by
    linarith [min_le_left length₁ (1 - (start₁ - start₂) + length₂)]
  rw [min_eq_left hwrap₁, max_eq_left horder,
    show start₁ + length₁ - start₁ = length₁ by ring,
    max_eq_right hlength₁0, max_eq_left hstart₁0,
    max_eq_left hcrossZero, add_zero,
    min_eq_left hcover, max_eq_right hlength₁0,
    max_eq_left hsecondZero, add_zero]

theorem overlap_le_wrap_wrap
    (start₁ start₂ length₁ length₂ : ℚ)
    (hstart₁0 : 0 ≤ start₁) (hstart₁1 : start₁ < 1)
    (hstart₂0 : 0 ≤ start₂) (hstart₂1 : start₂ < 1)
    (hlength₁0 : 0 ≤ length₁) (hlength₁1 : length₁ ≤ 1)
    (hlength₂0 : 0 ≤ length₂) (hlength₂1 : length₂ ≤ 1)
    (horder : start₂ ≤ start₁)
    (hwrap₁ : ¬start₁ + length₁ ≤ 1)
    (hwrap₂ : ¬start₂ + length₂ ≤ 1) :
    max 0 (min 1 1 - max start₁ start₂) +
        max 0 (min 1 (start₂ + length₂ - 1) - max start₁ 0) +
        max 0 (min (start₁ + length₁ - 1) 1 - max 0 start₂) +
        max 0 (min (start₁ + length₁ - 1) (start₂ + length₂ - 1) - max 0 0) =
      pairOverlapQ length₁ length₂ (start₁ - start₂) := by
  unfold pairOverlapQ
  have hsum₁ : 1 < start₁ + length₁ := lt_of_not_ge hwrap₁
  have hsum₂ : 1 < start₂ + length₂ := lt_of_not_ge hwrap₂
  let overflow₁ := start₁ + length₁ - 1
  let overflow₂ := start₂ + length₂ - 1
  have hoverflow₁0 : 0 < overflow₁ := by dsimp [overflow₁]; linarith
  have hoverflow₂0 : 0 < overflow₂ := by dsimp [overflow₂]; linarith
  have hoverflow₁le : overflow₁ ≤ start₁ := by dsimp [overflow₁]; linarith
  have hoverflow₂le : overflow₂ ≤ start₂ := by dsimp [overflow₂]; linarith
  have hoverflow₁one : overflow₁ ≤ 1 := hoverflow₁le.trans hstart₁1.le
  have hoverflow₂one : overflow₂ ≤ 1 := hoverflow₂le.trans hstart₂1.le
  have hcrossZero : overflow₂ - start₁ ≤ 0 := by linarith
  have hminimumNonneg : 0 ≤ min overflow₁ overflow₂ :=
    le_min hoverflow₁0.le hoverflow₂0.le
  have hfirstNonneg : 0 ≤ (1 - start₁) + min overflow₁ overflow₂ := by
    linarith
  have hsecondLarge :
      length₁ ≤ 1 - (start₁ - start₂) + length₂ := by
    linarith
  have hfirstMinimum :
      min length₁ (length₂ - (start₁ - start₂)) =
        (1 - start₁) + min overflow₁ overflow₂ := by
    by_cases hoverflow : overflow₁ ≤ overflow₂
    · rw [min_eq_left hoverflow]
      have hlength : length₁ ≤ length₂ - (start₁ - start₂) := by
        dsimp [overflow₁, overflow₂] at hoverflow
        linarith
      rw [min_eq_left hlength]
      dsimp [overflow₁]
      ring
    · have hoverflow' : overflow₂ ≤ overflow₁ := le_of_not_ge hoverflow
      rw [min_eq_right hoverflow']
      have hlength : length₂ - (start₁ - start₂) ≤ length₁ := by
        dsimp [overflow₁, overflow₂] at hoverflow'
        linarith
      rw [min_eq_right hlength]
      dsimp [overflow₂]
      ring
  rw [min_self, max_eq_left horder,
    max_eq_right (by linarith : (0 : ℚ) ≤ 1 - start₁),
    min_eq_right hoverflow₂one,
    max_eq_left hstart₁0,
    max_eq_left hcrossZero,
    min_eq_left hoverflow₁one,
    max_eq_right hstart₂0,
    max_self, sub_zero,
    show min (start₁ + length₁ - 1) (start₂ + length₂ - 1) =
        min overflow₁ overflow₂ by rfl,
    max_eq_right hminimumNonneg,
    min_eq_left hsecondLarge, hfirstMinimum,
    max_eq_right hfirstNonneg]
  dsimp [overflow₁]
  ring

theorem overlap_lt_nowrap_nowrap
    (start₁ start₂ length₁ length₂ : ℚ)
    (hstart₁0 : 0 ≤ start₁) (hstart₁1 : start₁ < 1)
    (hstart₂0 : 0 ≤ start₂) (hstart₂1 : start₂ < 1)
    (hlength₁0 : 0 ≤ length₁) (hlength₁1 : length₁ ≤ 1)
    (hlength₂0 : 0 ≤ length₂) (hlength₂1 : length₂ ≤ 1)
    (horder : start₁ < start₂)
    (hwrap₁ : start₁ + length₁ ≤ 1)
    (hwrap₂ : start₂ + length₂ ≤ 1) :
    max 0 (min (start₁ + length₁) (start₂ + length₂) - max start₁ start₂) =
      pairOverlapQ length₁ length₂ (start₁ - start₂ + 1) := by
  unfold pairOverlapQ
  simp only [max_def, min_def]
  split_ifs <;> linarith

theorem overlap_lt_wrap_nowrap
    (start₁ start₂ length₁ length₂ : ℚ)
    (hstart₁0 : 0 ≤ start₁) (hstart₁1 : start₁ < 1)
    (hstart₂0 : 0 ≤ start₂) (hstart₂1 : start₂ < 1)
    (hlength₁0 : 0 ≤ length₁) (hlength₁1 : length₁ ≤ 1)
    (hlength₂0 : 0 ≤ length₂) (hlength₂1 : length₂ ≤ 1)
    (horder : start₁ < start₂)
    (hwrap₁ : ¬start₁ + length₁ ≤ 1)
    (hwrap₂ : start₂ + length₂ ≤ 1) :
    max 0 (min 1 (start₂ + length₂) - max start₁ start₂) +
        max 0 (min (start₁ + length₁ - 1) (start₂ + length₂) - max 0 start₂) =
      pairOverlapQ length₁ length₂ (start₁ - start₂ + 1) := by
  unfold pairOverlapQ
  simp only [max_def, min_def]
  split_ifs <;> linarith

theorem overlap_lt_nowrap_wrap
    (start₁ start₂ length₁ length₂ : ℚ)
    (hstart₁0 : 0 ≤ start₁) (hstart₁1 : start₁ < 1)
    (hstart₂0 : 0 ≤ start₂) (hstart₂1 : start₂ < 1)
    (hlength₁0 : 0 ≤ length₁) (hlength₁1 : length₁ ≤ 1)
    (hlength₂0 : 0 ≤ length₂) (hlength₂1 : length₂ ≤ 1)
    (horder : start₁ < start₂)
    (hwrap₁ : start₁ + length₁ ≤ 1)
    (hwrap₂ : ¬start₂ + length₂ ≤ 1) :
    max 0 (min (start₁ + length₁) 1 - max start₁ start₂) +
        max 0 (min (start₁ + length₁) (start₂ + length₂ - 1) - max start₁ 0) =
      pairOverlapQ length₁ length₂ (start₁ - start₂ + 1) := by
  have hswapped := overlap_le_wrap_nowrap
    start₂ start₁ length₂ length₁
    hstart₂0 hstart₂1 hstart₁0 hstart₁1
    hlength₂0 hlength₂1 hlength₁0 hlength₁1 horder.le
    hwrap₂ hwrap₁
  have hshift0 : 0 < start₁ - start₂ + 1 := by linarith
  have hshift1 : start₁ - start₂ + 1 < 1 := by linarith
  calc
    max 0 (min (start₁ + length₁) 1 - max start₁ start₂) +
        max 0 (min (start₁ + length₁) (start₂ + length₂ - 1) - max start₁ 0) =
        pairOverlapQ length₂ length₁ (start₂ - start₁) := by
      simpa [min_comm, max_comm, add_comm] using hswapped
    _ = pairOverlapQ length₁ length₂ (start₁ - start₂ + 1) := by
      rw [show start₂ - start₁ = 1 - (start₁ - start₂ + 1) by ring]
      exact (pairOverlapQ_swap_pos length₁ length₂
        (start₁ - start₂ + 1)
        hlength₁0 hlength₁1 hlength₂0 hlength₂1 hshift0 hshift1).symm

theorem overlap_lt_wrap_wrap
    (start₁ start₂ length₁ length₂ : ℚ)
    (hstart₁0 : 0 ≤ start₁) (hstart₁1 : start₁ < 1)
    (hstart₂0 : 0 ≤ start₂) (hstart₂1 : start₂ < 1)
    (hlength₁0 : 0 ≤ length₁) (hlength₁1 : length₁ ≤ 1)
    (hlength₂0 : 0 ≤ length₂) (hlength₂1 : length₂ ≤ 1)
    (horder : start₁ < start₂)
    (hwrap₁ : ¬start₁ + length₁ ≤ 1)
    (hwrap₂ : ¬start₂ + length₂ ≤ 1) :
    max 0 (min 1 1 - max start₁ start₂) +
        max 0 (min 1 (start₂ + length₂ - 1) - max start₁ 0) +
        max 0 (min (start₁ + length₁ - 1) 1 - max 0 start₂) +
        max 0 (min (start₁ + length₁ - 1) (start₂ + length₂ - 1) - max 0 0) =
      pairOverlapQ length₁ length₂ (start₁ - start₂ + 1) := by
  have hswapped := overlap_le_wrap_wrap
    start₂ start₁ length₂ length₁
    hstart₂0 hstart₂1 hstart₁0 hstart₁1
    hlength₂0 hlength₂1 hlength₁0 hlength₁1 horder.le
    hwrap₂ hwrap₁
  have hshift0 : 0 < start₁ - start₂ + 1 := by linarith
  have hshift1 : start₁ - start₂ + 1 < 1 := by linarith
  calc
    max 0 (min 1 1 - max start₁ start₂) +
        max 0 (min 1 (start₂ + length₂ - 1) - max start₁ 0) +
        max 0 (min (start₁ + length₁ - 1) 1 - max 0 start₂) +
        max 0 (min (start₁ + length₁ - 1) (start₂ + length₂ - 1) - max 0 0) =
        pairOverlapQ length₂ length₁ (start₂ - start₁) := by
      simpa [min_comm, max_comm, add_comm, add_left_comm, add_assoc] using hswapped
    _ = pairOverlapQ length₁ length₂ (start₁ - start₂ + 1) := by
      rw [show start₂ - start₁ = 1 - (start₁ - start₂ + 1) by ring]
      exact (pairOverlapQ_swap_pos length₁ length₂
        (start₁ - start₂ + 1)
        hlength₁0 hlength₁1 hlength₂0 hlength₂1 hshift0 hshift1).symm

/-- Exact clip-sum formula for two normalized circular intervals. -/
theorem length_inter_circleInterval_eq_pairOverlapQ
    (start₁ start₂ length₁ length₂ : ℚ)
    (hstart₁0 : 0 ≤ start₁) (hstart₁1 : start₁ < 1)
    (hstart₂0 : 0 ≤ start₂) (hstart₂1 : start₂ < 1)
    (hlength₁0 : 0 ≤ length₁) (hlength₁1 : length₁ ≤ 1)
    (hlength₂0 : 0 ≤ length₂) (hlength₂1 : length₂ ≤ 1) :
    RatIntervals.length
        (inter (circleInterval start₁ length₁)
          (circleInterval start₂ length₂)) =
      pairOverlapQ length₁ length₂ (Int.fract (start₁ - start₂)) := by
  by_cases horder : start₂ ≤ start₁
  · rw [fract_sub_eq_of_le hstart₁0 hstart₁1 hstart₂0 hstart₂1 horder]
    unfold circleInterval inter RatIntervals.length clip
    split_ifs <;> simp only [List.flatMap_cons, List.flatMap_nil,
      List.map_cons, List.map_nil, List.append_nil, List.map_append,
      List.sum_cons, List.sum_nil, add_zero]
    all_goals
      first
      | simpa [add_assoc] using
          overlap_le_nowrap_nowrap start₁ start₂ length₁ length₂
            hstart₁0 hstart₁1 hstart₂0 hstart₂1
            hlength₁0 hlength₁1 hlength₂0 hlength₂1 horder
            (by assumption) (by assumption)
      | simpa [add_assoc] using
          overlap_le_wrap_nowrap start₁ start₂ length₁ length₂
            hstart₁0 hstart₁1 hstart₂0 hstart₂1
            hlength₁0 hlength₁1 hlength₂0 hlength₂1 horder
            (by assumption) (by assumption)
      | simpa [add_assoc] using
          overlap_le_nowrap_wrap start₁ start₂ length₁ length₂
            hstart₁0 hstart₁1 hstart₂0 hstart₂1
            hlength₁0 hlength₁1 hlength₂0 hlength₂1 horder
            (by assumption) (by assumption)
      | simpa [add_assoc] using
          overlap_le_wrap_wrap start₁ start₂ length₁ length₂
            hstart₁0 hstart₁1 hstart₂0 hstart₂1
            hlength₁0 hlength₁1 hlength₂0 hlength₂1 horder
            (by assumption) (by assumption)
  · have horder' : start₁ < start₂ := lt_of_not_ge horder
    rw [fract_sub_eq_add_one_of_lt hstart₁0 hstart₁1 hstart₂0 hstart₂1 horder']
    unfold circleInterval inter RatIntervals.length clip
    split_ifs <;> simp only [List.flatMap_cons, List.flatMap_nil,
      List.map_cons, List.map_nil, List.append_nil, List.map_append,
      List.sum_cons, List.sum_nil, add_zero]
    all_goals
      first
      | simpa [add_assoc] using
          overlap_lt_nowrap_nowrap start₁ start₂ length₁ length₂
            hstart₁0 hstart₁1 hstart₂0 hstart₂1
            hlength₁0 hlength₁1 hlength₂0 hlength₂1 horder'
            (by assumption) (by assumption)
      | simpa [add_assoc] using
          overlap_lt_wrap_nowrap start₁ start₂ length₁ length₂
            hstart₁0 hstart₁1 hstart₂0 hstart₂1
            hlength₁0 hlength₁1 hlength₂0 hlength₂1 horder'
            (by assumption) (by assumption)
      | simpa [add_assoc] using
          overlap_lt_nowrap_wrap start₁ start₂ length₁ length₂
            hstart₁0 hstart₁1 hstart₂0 hstart₂1
            hlength₁0 hlength₁1 hlength₂0 hlength₂1 horder'
            (by assumption) (by assumption)
      | simpa [add_assoc] using
          overlap_lt_wrap_wrap start₁ start₂ length₁ length₂
            hstart₁0 hstart₁1 hstart₂0 hstart₂1
            hlength₁0 hlength₁1 hlength₂0 hlength₂1 horder'
            (by assumption) (by assumption)

#print axioms length_inter_circleInterval_eq_pairOverlapQ

end LRCB5WrappedToothOverlap
end LonelyRunner
