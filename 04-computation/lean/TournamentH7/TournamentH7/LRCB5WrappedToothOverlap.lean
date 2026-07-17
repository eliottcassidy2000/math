import TournamentH7.LRCRationalOpenComb

namespace LonelyRunner
namespace LRCB5WrappedToothOverlap

open RatIntervals LRCRationalOpenComb

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
  simp only [max_def, min_def]
  split_ifs <;> linarith

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
  rw [max_eq_right horder, min_self,
    max_eq_right (by linarith : (0 : ℚ) ≤ 1 - start₁),
    min_eq_right hoverflow₂one,
    max_eq_right hstart₁0,
    max_eq_left hcrossZero,
    min_eq_left hoverflow₁one,
    max_eq_right hstart₂0,
    max_eq_right (sub_nonneg.mpr (le_min hoverflow₁0.le hoverflow₂0.le)),
    max_self, min_eq_left hsecondLarge, hfirstMinimum]
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
  unfold pairOverlapQ
  simp only [max_def, min_def]
  split_ifs <;> linarith

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
  unfold pairOverlapQ
  simp only [max_def, min_def]
  split_ifs <;> linarith

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

def rawToothStart (w k : ℕ) : ℚ :=
  ((k : ℚ) - 1 / 14) / w

def toothLengthQ (w : ℕ) : ℚ :=
  1 / (7 * w)

def circularTooth (w k : ℕ) : Region :=
  circleInterval (Int.fract (rawToothStart w k)) (toothLengthQ w)

def circularCombRegion (w : ℕ) : Region :=
  (List.range w).flatMap (circularTooth w)

theorem fract_rawToothStart_zero
    (w : ℕ) (hw : 0 < w) :
    Int.fract (rawToothStart w 0) = 1 - 1 / (14 * w : ℚ) := by
  have hwQ : (0 : ℚ) < w := by exact_mod_cast hw
  unfold rawToothStart Int.fract
  have hfloor : ⌊(((0 : ℚ) - 1 / 14) / w)⌋ = -1 := by
    rw [Int.floor_eq_iff]
    push_cast
    constructor
    · rw [le_div_iff₀ hwQ]
      linarith [show (0 : ℚ) < 1 / 14 by norm_num]
    · positivity
  rw [hfloor]
  push_cast
  field_simp
  ring

theorem fract_rawToothStart_internal
    (w k : ℕ) (hw : 0 < w) (hk0 : 0 < k) (hkw : k < w) :
    Int.fract (rawToothStart w k) = rawToothStart w k := by
  apply Int.fract_eq_self.mpr
  have hwQ : (0 : ℚ) < w := by exact_mod_cast hw
  have hk0Q : (1 : ℚ) ≤ k := by exact_mod_cast hk0
  have hkwQ : (k : ℚ) < w := by exact_mod_cast hkw
  unfold rawToothStart
  constructor
  · exact div_nonneg (by linarith [show (0 : ℚ) < 1 / 14 by norm_num]) hwQ.le
  · rw [div_lt_one hwQ]
    linarith

theorem circularTooth_zero_eq_boundary
    (w : ℕ) (hw : 0 < w) :
    circularTooth w 0 =
      [ratOpenCombInterval w w, ratOpenCombInterval w 0] := by
  have hwQ : (0 : ℚ) < w := by exact_mod_cast hw
  rw [circularTooth, fract_rawToothStart_zero w hw]
  unfold circleInterval toothLengthQ ratOpenCombInterval
  have hsmall : (0 : ℚ) < 1 / (14 * w : ℚ) := by positivity
  have hsmallOne : 1 / (14 * w : ℚ) < 1 := by
    rw [div_lt_one (by positivity : (0 : ℚ) < 14 * w)]
    have hwOne : (1 : ℚ) ≤ w := by exact_mod_cast hw
    linarith
  have hwrap : ¬1 - 1 / (14 * w : ℚ) + 1 / (7 * w : ℚ) ≤ 1 := by
    field_simp
    positivity
  rw [if_neg hwrap]
  congr 1
  · apply Prod.ext <;> simp only [Prod.fst, Prod.snd]
    · rw [max_eq_right]
      · field_simp
        ring
      · rw [div_nonneg_iff]
        exact Or.inl ⟨by norm_num, hwQ.le⟩
    · rw [min_eq_left]
      field_simp
      linarith
  · congr 1
    apply Prod.ext <;> simp only [Prod.fst, Prod.snd]
    · rw [max_eq_left]
      · rfl
      · positivity
    · rw [min_eq_right]
      · field_simp
        ring
      · exact hsmallOne.le

theorem circularTooth_internal_eq_interval
    (w k : ℕ) (hw : 0 < w) (hk0 : 0 < k) (hkw : k < w) :
    circularTooth w k = [ratOpenCombInterval w k] := by
  have hwQ : (0 : ℚ) < w := by exact_mod_cast hw
  have hk0Q : (1 : ℚ) ≤ k := by exact_mod_cast hk0
  have hkwQ : (k : ℚ) < w := by exact_mod_cast hkw
  rw [circularTooth, fract_rawToothStart_internal w k hw hk0 hkw]
  unfold circleInterval toothLengthQ rawToothStart ratOpenCombInterval
  have hnowrap :
      ((k : ℚ) - 1 / 14) / w + 1 / (7 * w : ℚ) ≤ 1 := by
    rw [div_add_div_same, div_le_one hwQ]
    linarith
  rw [if_pos hnowrap]
  congr 1
  apply Prod.ext <;> simp only [Prod.fst, Prod.snd]
  · rw [max_eq_right]
    exact div_nonneg (by linarith) hwQ.le
  · rw [min_eq_right]
    · field_simp
      ring
    · rw [div_le_one hwQ]
      linarith

theorem ratOpenCombTail_succ_eq_append
    (w start count : ℕ) :
    ratOpenCombTail w start (count + 1) =
      ratOpenCombTail w start count ++
        [ratOpenCombInterval w (start + count)] := by
  induction count generalizing start with
  | zero => simp [ratOpenCombTail]
  | succ count inductionHypothesis =>
      simp only [Nat.succ_eq_add_one]
      rw [show count + 1 + 1 = (count + 1) + 1 by rfl,
        ratOpenCombTail, ratOpenCombTail,
        inductionHypothesis (start + 1)]
      simp only [List.cons_append, List.cons.injEq, true_and]
      congr 2
      omega

theorem ratOpenCombRegion_eq_tail_append_boundary (w : ℕ) :
    ratOpenCombRegion w =
      ratOpenCombTail w 0 w ++ [ratOpenCombInterval w w] := by
  unfold ratOpenCombRegion
  simpa using ratOpenCombTail_succ_eq_append w 0 w

theorem circularInternalFlatMap_eq_tail
    (w start count : ℕ) (hw : 0 < w) (hstart : 0 < start)
    (hbound : start + count ≤ w) :
    (List.range' start count).flatMap (circularTooth w) =
      ratOpenCombTail w start count := by
  induction count generalizing start with
  | zero => simp [ratOpenCombTail]
  | succ count inductionHypothesis =>
      rw [show count + 1 = count + 1 by rfl, List.range'_succ]
      simp only [List.flatMap_cons]
      rw [circularTooth_internal_eq_interval w start hw hstart (by omega),
        ratOpenCombTail]
      congr 1
      exact inductionHypothesis (start + 1) (by omega) (by omega)

theorem circularCombRegion_eq_boundary_cons_tail
    (w : ℕ) (hw : 0 < w) :
    circularCombRegion w =
      ratOpenCombInterval w w :: ratOpenCombTail w 0 w := by
  obtain ⟨count, rfl⟩ := Nat.exists_eq_succ_of_ne_zero hw.ne'
  unfold circularCombRegion
  rw [List.range_eq_range', List.range'_succ]
  simp only [List.flatMap_cons]
  rw [circularTooth_zero_eq_boundary (count + 1) (by omega),
    circularInternalFlatMap_eq_tail (count + 1) 1 count
      (by omega) (by omega) (by omega),
    ratOpenCombTail]
  rfl

theorem circularCombRegion_perm_ratOpenCombRegion
    (w : ℕ) (hw : 0 < w) :
    (circularCombRegion w).Perm (ratOpenCombRegion w) := by
  rw [circularCombRegion_eq_boundary_cons_tail w hw,
    ratOpenCombRegion_eq_tail_append_boundary]
  simpa using
    (List.perm_append_comm (l₁ := [ratOpenCombInterval w w])
      (l₂ := ratOpenCombTail w 0 w))

theorem length_eq_of_perm {left right : Region} (hperm : left.Perm right) :
    RatIntervals.length left = RatIntervals.length right := by
  unfold RatIntervals.length
  exact (hperm.map fun interval => max 0 (interval.2 - interval.1)).sum_eq

theorem length_inter_eq_of_perm_left
    {first first' second : Region} (hperm : first.Perm first') :
    RatIntervals.length (inter first second) =
      RatIntervals.length (inter first' second) := by
  unfold inter
  apply length_eq_of_perm
  exact hperm.flatMap fun _ _ => List.Perm.refl _

theorem length_inter_eq_of_perm_right
    {first second second' : Region} (hperm : second.Perm second') :
    RatIntervals.length (inter first second) =
      RatIntervals.length (inter first second') := by
  rw [RatIntervals.length_inter_comm first second,
    RatIntervals.length_inter_comm first second']
  exact length_inter_eq_of_perm_left hperm

theorem length_inter_ratOpenCombRegion_eq_circularCombRegion
    (first second : ℕ) (hfirst : 0 < first) (hsecond : 0 < second) :
    RatIntervals.length
        (inter (ratOpenCombRegion first) (ratOpenCombRegion second)) =
      RatIntervals.length
        (inter (circularCombRegion first) (circularCombRegion second)) := by
  calc
    RatIntervals.length
        (inter (ratOpenCombRegion first) (ratOpenCombRegion second)) =
        RatIntervals.length
          (inter (circularCombRegion first) (ratOpenCombRegion second)) :=
      length_inter_eq_of_perm_left
        (circularCombRegion_perm_ratOpenCombRegion first hfirst).symm
    _ = RatIntervals.length
          (inter (circularCombRegion first) (circularCombRegion second)) :=
      length_inter_eq_of_perm_right
        (circularCombRegion_perm_ratOpenCombRegion second hsecond).symm

#print axioms length_inter_circleInterval_eq_pairOverlapQ

end LRCB5WrappedToothOverlap
end LonelyRunner
