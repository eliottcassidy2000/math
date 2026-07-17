/-
  TournamentH7.LRCCanonicalCircleAtlas

  Canonical adjacent-interval coalescing for the rational circle atlas.
-/

import TournamentH7.LRCLocalDensityBlockGluing

namespace LRC14.LocalDensityBlockGluing

open LonelyRunner.RatIntervals

/-- A sorted region has no artificial internal seam when every earlier piece
ends strictly before every later piece starts. -/
def RegionSeamFree (region : Region) : Prop :=
  region.Pairwise fun left right => left.2 < right.1

/-- Prepend one interval to an already coalesced tail, merging the unique
adjacent head when it exists. -/
def prependCoalesced (interval : ℚ × ℚ) : Region → Region
  | [] => [interval]
  | next :: remaining =>
      if interval.2 = next.1 then
        (interval.1, next.2) :: remaining
      else
        interval :: next :: remaining

/-- Merge every chain of adjacent intervals in a sorted region. -/
def coalesceAdjacent : Region → Region
  | [] => []
  | interval :: remaining =>
      prependCoalesced interval (coalesceAdjacent remaining)

theorem mem_cons_region_iff {x : ℚ} {interval : ℚ × ℚ}
    {region : Region} :
    mem x (interval :: region) ↔
      (interval.1 ≤ x ∧ x < interval.2) ∨ mem x region := by
  unfold mem
  constructor
  · rintro ⟨source, hsource, hx⟩
    rcases List.mem_cons.mp hsource with rfl | hsource
    · exact Or.inl hx
    · exact Or.inr ⟨source, hsource, hx⟩
  · rintro (hx | ⟨source, hsource, hx⟩)
    · exact ⟨interval, List.mem_cons_self .., hx⟩
    · exact ⟨source, List.mem_cons_of_mem _ hsource, hx⟩

theorem mem_merge_adjacent_iff {x : ℚ}
    {left right : ℚ × ℚ} {remaining : Region}
    (hleftLive : left.1 < left.2) (hrightLive : right.1 < right.2)
    (hadjacent : left.2 = right.1) :
    mem x ((left.1, right.2) :: remaining) ↔
      mem x (left :: right :: remaining) := by
  rw [mem_cons_region_iff, mem_cons_region_iff, mem_cons_region_iff]
  constructor
  · rintro (hmerged | hremaining)
    · by_cases hleft : x < left.2
      · exact Or.inl ⟨hmerged.1, hleft⟩
      · exact Or.inr (Or.inl
          ⟨by rw [← hadjacent]; exact le_of_not_gt hleft, hmerged.2⟩)
    · exact Or.inr (Or.inr hremaining)
  · rintro (hleft | hright | hremaining)
    · exact Or.inl ⟨hleft.1, lt_trans hleft.2 (by
          rw [hadjacent]
          exact hrightLive)⟩
    · exact Or.inl ⟨by
          dsimp
          linarith [hleftLive, hright.1], hright.2⟩
    · exact Or.inr hremaining

theorem coalesceAdjacent_cons_head
    (interval : ℚ × ℚ) (region : Region) :
    ∃ head remaining,
      coalesceAdjacent (interval :: region) = head :: remaining ∧
        head.1 = interval.1 := by
  induction region generalizing interval with
  | nil =>
      exact ⟨interval, [], rfl, rfl⟩
  | cons next region inductionHypothesis =>
      obtain ⟨head, remaining, hcoalesce, hstart⟩ :=
        inductionHypothesis next
      rw [coalesceAdjacent, hcoalesce]
      unfold prependCoalesced
      by_cases hadjacent : interval.2 = head.1
      · exact ⟨(interval.1, head.2), remaining, by simp [hadjacent], rfl⟩
      · exact ⟨interval, head :: remaining, by simp [hadjacent], rfl⟩

/-- Coalescing a `Norm` region preserves its exact half-open carrier and
produces a normalized presentation with no adjacent internal seams. -/
theorem coalesceAdjacent_spec : ∀ {region : Region}, Norm region →
    Norm (coalesceAdjacent region) ∧
      RegionSeamFree (coalesceAdjacent region) ∧
      ∀ x, mem x (coalesceAdjacent region) ↔ mem x region := by
  intro region
  induction region with
  | nil =>
      intro _
      simp [coalesceAdjacent, RegionSeamFree, mem,
        LonelyRunner.RatIntervals.Norm]
  | cons interval region inductionHypothesis =>
      intro hnorm
      cases region with
      | nil =>
          exact ⟨by simpa [coalesceAdjacent, prependCoalesced] using hnorm,
            by simp [coalesceAdjacent, prependCoalesced, RegionSeamFree],
            by intro x; simp [coalesceAdjacent, prependCoalesced]⟩
      | cons next remaining =>
          have htailNorm : Norm (next :: remaining) := norm_tail hnorm
          obtain ⟨hcoalescedNorm, hcoalescedSeamFree, hcarrier⟩ :=
            inductionHypothesis htailNorm
          obtain ⟨head, tail, hcoalesce, hheadStart⟩ :=
            coalesceAdjacent_cons_head next remaining
          rw [hcoalesce] at hcoalescedNorm hcoalescedSeamFree hcarrier
          rw [coalesceAdjacent, hcoalesce]
          unfold prependCoalesced
          have hintervalLive : interval.1 < interval.2 := norm_head_lt hnorm
          have hheadLive : head.1 < head.2 :=
            norm_head_lt hcoalescedNorm
          have hintervalOrder : interval.2 ≤ head.1 := by
            rw [hheadStart]
            exact hnorm.2.1
          by_cases hadjacent : interval.2 = head.1
          · simp only [hadjacent, if_true]
            have htailNorm : Norm tail := norm_tail hcoalescedNorm
            have hmergedNorm : Norm ((interval.1, head.2) :: tail) := by
              apply norm_cons_region (by linarith) htailNorm
              intro later hlater
              exact norm_head_le hcoalescedNorm later hlater
            have hmergedSeamFree :
                RegionSeamFree ((interval.1, head.2) :: tail) := by
              unfold RegionSeamFree at hcoalescedSeamFree ⊢
              rw [List.pairwise_cons]
              exact ⟨(List.pairwise_cons.mp hcoalescedSeamFree).1,
                (List.pairwise_cons.mp hcoalescedSeamFree).2⟩
            refine ⟨hmergedNorm, hmergedSeamFree, ?_⟩
            intro x
            calc
              mem x ((interval.1, head.2) :: tail) ↔
                  mem x (interval :: head :: tail) :=
                mem_merge_adjacent_iff hintervalLive hheadLive hadjacent
              _ ↔
                  (interval.1 ≤ x ∧ x < interval.2) ∨
                    mem x (head :: tail) := mem_cons_region_iff
              _ ↔
                  (interval.1 ≤ x ∧ x < interval.2) ∨
                    mem x (next :: remaining) :=
                or_congr Iff.rfl (hcarrier x)
              _ ↔ mem x (interval :: next :: remaining) :=
                mem_cons_region_iff.symm
          · simp only [hadjacent, if_false]
            have hstrict : interval.2 < head.1 := lt_of_le_of_ne
              hintervalOrder hadjacent
            have hnewNorm : Norm (interval :: head :: tail) := by
              apply norm_cons_region hintervalLive hcoalescedNorm
              intro later hlater
              rcases List.mem_cons.mp hlater with rfl | hlater
              · exact hstrict.le
              · have hheadOrder := norm_head_le hcoalescedNorm later hlater
                linarith
            have hnewSeamFree :
                RegionSeamFree (interval :: head :: tail) := by
              unfold RegionSeamFree at hcoalescedSeamFree ⊢
              rw [List.pairwise_cons]
              constructor
              · intro later hlater
                rcases List.mem_cons.mp hlater with rfl | hlater
                · exact hstrict
                · have hheadOrder := norm_head_le hcoalescedNorm later hlater
                  linarith
              · exact hcoalescedSeamFree
            refine ⟨hnewNorm, hnewSeamFree, ?_⟩
            intro x
            calc
              mem x (interval :: head :: tail) ↔
                  (interval.1 ≤ x ∧ x < interval.2) ∨
                    mem x (head :: tail) := mem_cons_region_iff
              _ ↔
                  (interval.1 ≤ x ∧ x < interval.2) ∨
                    mem x (next :: remaining) :=
                or_congr Iff.rfl (hcarrier x)
              _ ↔ mem x (interval :: next :: remaining) :=
                mem_cons_region_iff.symm

theorem norm_coalesceAdjacent {region : Region} (hnorm : Norm region) :
    Norm (coalesceAdjacent region) :=
  (coalesceAdjacent_spec hnorm).1

theorem regionSeamFree_coalesceAdjacent {region : Region}
    (hnorm : Norm region) :
    RegionSeamFree (coalesceAdjacent region) :=
  (coalesceAdjacent_spec hnorm).2.1

theorem mem_coalesceAdjacent_iff {region : Region} (hnorm : Norm region)
    (x : ℚ) :
    mem x (coalesceAdjacent region) ↔ mem x region :=
  (coalesceAdjacent_spec hnorm).2.2 x

theorem regionInUnit_coalesceAdjacent {region : Region}
    (hunit : RegionInUnit region) :
    RegionInUnit (coalesceAdjacent region) := by
  induction region with
  | nil => simp [coalesceAdjacent, RegionInUnit]
  | cons interval region inductionHypothesis =>
      have hintervalUnit := hunit interval (List.mem_cons_self ..)
      have htailUnit : RegionInUnit region := by
        intro piece hpiece
        exact hunit piece (List.mem_cons_of_mem _ hpiece)
      have hcoalescedTailUnit := inductionHypothesis htailUnit
      cases region with
      | nil =>
          simpa [coalesceAdjacent, prependCoalesced, RegionInUnit] using
            hintervalUnit
      | cons next remaining =>
          obtain ⟨head, tail, hcoalesce, _⟩ :=
            coalesceAdjacent_cons_head next remaining
          rw [hcoalesce] at hcoalescedTailUnit
          rw [coalesceAdjacent, hcoalesce]
          unfold prependCoalesced
          by_cases hadjacent : interval.2 = head.1
          · simp only [hadjacent, if_true]
            intro piece hpiece
            rcases List.mem_cons.mp hpiece with rfl | hpiece
            · exact ⟨hintervalUnit.1,
                (hcoalescedTailUnit head (List.mem_cons_self ..)).2⟩
            · exact hcoalescedTailUnit piece
                (List.mem_cons_of_mem _ hpiece)
          · simp only [hadjacent, if_false]
            intro piece hpiece
            rcases List.mem_cons.mp hpiece with rfl | hpiece
            · exact hintervalUnit
            · exact hcoalescedTailUnit piece hpiece

/-- Symmetric formulation of seam-freeness, useful when two pieces are
obtained existentially rather than in list order. -/
def RegionHasNoAdjacentSeams (region : Region) : Prop :=
  ∀ left ∈ region, ∀ right ∈ region, left ≠ right →
    left.2 ≠ right.1

theorem regionHasNoAdjacentSeams_of_norm_seamFree :
    ∀ {region : Region}, Norm region → RegionSeamFree region →
      RegionHasNoAdjacentSeams region := by
  intro region
  induction region with
  | nil => simp [RegionHasNoAdjacentSeams]
  | cons head tail inductionHypothesis =>
      intro hnorm hseamFree left hleft right hright hne
      have hheadLive := norm_head_lt hnorm
      have hpairwise := List.pairwise_cons.mp hseamFree
      rcases List.mem_cons.mp hleft with hleftHead | hleftTail
      · subst left
        rcases List.mem_cons.mp hright with hrightHead | hrightTail
        · subst right
          exact (hne rfl).elim
        · have hstrict := hpairwise.1 right hrightTail
          linarith
      · rcases List.mem_cons.mp hright with hrightHead | hrightTail
        · subst right
          have hstrict := hpairwise.1 left hleftTail
          have hleftLive :=
            regionLive_of_norm (norm_tail hnorm) left hleftTail
          linarith
        · exact inductionHypothesis (norm_tail hnorm) hpairwise.2
            left hleftTail right hrightTail hne

theorem regionHasNoAdjacentSeams_coalesceAdjacent {region : Region}
    (hnorm : Norm region) :
    RegionHasNoAdjacentSeams (coalesceAdjacent region) :=
  regionHasNoAdjacentSeams_of_norm_seamFree
    (norm_coalesceAdjacent hnorm)
    (regionSeamFree_coalesceAdjacent hnorm)

def rotationSplitIndicator (shift : ℚ) (interval : ℚ × ℚ) : ℕ :=
  wrapSplitIndicator (interval.1 + shift, interval.2 + shift)

theorem rotationSplitCount_eq_indicatorSum (shift : ℚ) (region : Region) :
    rotationSplitCount shift region =
      (region.map (rotationSplitIndicator shift)).sum := by
  unfold rotationSplitCount wrapSplitCount translate
  rw [List.map_map]
  rfl

theorem rotationSplitIndicator_eq_zero_or_one
    (shift : ℚ) (interval : ℚ × ℚ) :
    rotationSplitIndicator shift interval = 0 ∨
      rotationSplitIndicator shift interval = 1 := by
  unfold rotationSplitIndicator wrapSplitIndicator
  split <;> simp_all

theorem exists_integer_of_rotationSplitIndicator_eq_one
    (shift : ℚ) (interval : ℚ × ℚ)
    (hsplit : rotationSplitIndicator shift interval = 1) :
    ∃ cutInteger : ℤ,
      interval.1 + shift < (cutInteger : ℚ) ∧
        (cutInteger : ℚ) < interval.2 + shift := by
  unfold rotationSplitIndicator wrapSplitIndicator at hsplit
  by_cases hend :
      interval.2 + shift - (⌊interval.1 + shift⌋ : ℚ) ≤ 1
  · simp [hend] at hsplit
  · refine ⟨⌊interval.1 + shift⌋ + 1, ?_, ?_⟩
    · push_cast
      exact Int.lt_floor_add_one (interval.1 + shift)
    · push_cast
      have hoverhang :
          1 < interval.2 + shift -
            (⌊interval.1 + shift⌋ : ℚ) := lt_of_not_ge hend
      linarith

theorem rotationSplitIndicator_eq_zero_of_order
    (shift : ℚ) {left right : ℚ × ℚ}
    (hleftUnit : 0 ≤ left.1 ∧ left.2 ≤ 1)
    (hrightUnit : 0 ≤ right.1 ∧ right.2 ≤ 1)
    (horder : left.2 ≤ right.1)
    (hleftSplit : rotationSplitIndicator shift left = 1) :
    rotationSplitIndicator shift right = 0 := by
  rcases rotationSplitIndicator_eq_zero_or_one shift right with hzero | hsplit
  · exact hzero
  · obtain ⟨leftInteger, hleftLower, hleftUpper⟩ :=
      exists_integer_of_rotationSplitIndicator_eq_one shift left hleftSplit
    obtain ⟨rightInteger, hrightLower, hrightUpper⟩ :=
      exists_integer_of_rotationSplitIndicator_eq_one shift right hsplit
    have hintegerOrderRat :
        (leftInteger : ℚ) < (rightInteger : ℚ) := by
      linarith
    have hintegerOrder : leftInteger < rightInteger := by
      exact_mod_cast hintegerOrderRat
    have hintegerGapRat :
        (rightInteger : ℚ) - (leftInteger : ℚ) < 1 := by
      linarith [hleftUnit.1, hrightUnit.2]
    have hintegerGap : rightInteger - leftInteger < 1 := by
      exact_mod_cast hintegerGapRat
    omega

theorem rotationSplitCount_le_one
    (shift : ℚ) {region : Region}
    (hnorm : Norm region) (hunit : RegionInUnit region) :
    rotationSplitCount shift region ≤ 1 := by
  induction region with
  | nil => simp [rotationSplitCount_eq_indicatorSum]
  | cons interval region inductionHypothesis =>
      have htailNorm := norm_tail hnorm
      have htailUnit : RegionInUnit region := by
        intro piece hpiece
        exact hunit piece (List.mem_cons_of_mem _ hpiece)
      rcases rotationSplitIndicator_eq_zero_or_one shift interval with
        hzero | hsplit
      · rw [rotationSplitCount_eq_indicatorSum]
        simp only [List.map_cons, List.sum_cons, hzero, zero_add]
        simpa [rotationSplitCount_eq_indicatorSum] using
          inductionHypothesis htailNorm htailUnit
      · have htailZero :
            ∀ later ∈ region,
              rotationSplitIndicator shift later = 0 := by
          intro later hlater
          apply rotationSplitIndicator_eq_zero_of_order shift
          · exact hunit interval (List.mem_cons_self ..)
          · exact hunit later (List.mem_cons_of_mem _ hlater)
          · exact norm_head_le hnorm later hlater
          · exact hsplit
        have hsumZero :
            (region.map (rotationSplitIndicator shift)).sum = 0 := by
          rw [List.sum_eq_zero_iff]
          intro value hvalue
          rw [List.mem_map] at hvalue
          obtain ⟨later, hlater, rfl⟩ := hvalue
          exact htailZero later hlater
        rw [rotationSplitCount_eq_indicatorSum]
        simp [hsplit, hsumZero]

theorem rotationSplitIndicator_eq_zero_of_count_zero
    (shift : ℚ) {region : Region}
    (hzero : rotationSplitCount shift region = 0)
    {interval : ℚ × ℚ} (hinterval : interval ∈ region) :
    rotationSplitIndicator shift interval = 0 := by
  have hsumZero :
      (region.map (rotationSplitIndicator shift)).sum = 0 := by
    simpa [rotationSplitCount_eq_indicatorSum] using hzero
  apply List.sum_eq_zero_iff.mp hsumZero
  exact List.mem_map.mpr ⟨interval, hinterval, rfl⟩

theorem exists_rotationSplitIndicator_eq_one_of_count_eq_one
    (shift : ℚ) {region : Region}
    (hone : rotationSplitCount shift region = 1) :
    ∃ interval ∈ region, rotationSplitIndicator shift interval = 1 := by
  have hsumPos :
      0 < (region.map (rotationSplitIndicator shift)).sum := by
    rw [← rotationSplitCount_eq_indicatorSum, hone]
    omega
  obtain ⟨value, hvalue, hpositive⟩ :=
    List.sum_pos_iff_exists_pos_nat.mp hsumPos
  rw [List.mem_map] at hvalue
  obtain ⟨interval, hinterval, rfl⟩ := hvalue
  refine ⟨interval, hinterval, ?_⟩
  rcases rotationSplitIndicator_eq_zero_or_one shift interval with
    hzero | hone
  · omega
  · exact hone

theorem wrapOne_rotation_eq_pair_of_split
    (shift : ℚ) (interval : ℚ × ℚ)
    (hsplit : rotationSplitIndicator shift interval = 1) :
    wrapOne (interval.1 + shift, interval.2 + shift) =
      [((interval.1 + shift) - ⌊interval.1 + shift⌋, 1),
        (0, (interval.2 + shift) - ⌊interval.1 + shift⌋ - 1)] := by
  have hoverhang :
      ¬ (interval.2 + shift -
        (⌊interval.1 + shift⌋ : ℚ) ≤ 1) := by
    intro hend
    simp [rotationSplitIndicator, wrapSplitIndicator, hend] at hsplit
  simp [wrapOne, hoverhang]

theorem boundaryPiecesJoined_translateCirc_eq_true_of_split
    (shift : ℚ) {region : Region}
    (hnorm : Norm region) (hunit : RegionInUnit region)
    {interval : ℚ × ℚ} (hinterval : interval ∈ region)
    (hsplit : rotationSplitIndicator shift interval = 1) :
    boundaryPiecesJoined (translateCirc shift region) = true := by
  have hcountPositive : 0 < rotationSplitCount shift region := by
    by_contra hnot
    have hzero : rotationSplitCount shift region = 0 := by omega
    have := rotationSplitIndicator_eq_zero_of_count_zero
      shift hzero hinterval
    omega
  have hcountOne : rotationSplitCount shift region = 1 := by
    have hle := rotationSplitCount_le_one shift hnorm hunit
    omega
  have hregionNonempty : 1 ≤ intervalComponentCount region := by
    unfold intervalComponentCount
    have := List.length_pos_of_mem hinterval
    omega
  have htwo :
      2 ≤ intervalComponentCount (translateCirc shift region) := by
    rw [intervalComponentCount_translateCirc, hcountOne]
    omega
  let endPiece : ℚ × ℚ :=
    ((interval.1 + shift) - ⌊interval.1 + shift⌋, 1)
  let startPiece : ℚ × ℚ :=
    (0, (interval.2 + shift) - ⌊interval.1 + shift⌋ - 1)
  have hendPiece : endPiece ∈ translateCirc shift region := by
    unfold translateCirc wrap translate
    rw [List.mem_flatMap]
    refine ⟨(interval.1 + shift, interval.2 + shift), ?_, ?_⟩
    · exact List.mem_map.mpr ⟨interval, hinterval, rfl⟩
    · rw [wrapOne_rotation_eq_pair_of_split shift interval hsplit]
      exact List.mem_cons_self
  have hstartPiece : startPiece ∈ translateCirc shift region := by
    unfold translateCirc wrap translate
    rw [List.mem_flatMap]
    refine ⟨(interval.1 + shift, interval.2 + shift), ?_, ?_⟩
    · exact List.mem_map.mpr ⟨interval, hinterval, rfl⟩
    · rw [wrapOne_rotation_eq_pair_of_split shift interval hsplit]
      exact List.mem_cons_of_mem _ (List.mem_singleton_self _)
  have hends :
      (translateCirc shift region).any
          (fun piece => piece.2 == 1) = true :=
    List.any_eq_true.mpr ⟨endPiece, hendPiece, by simp [endPiece]⟩
  have hstarts :
      (translateCirc shift region).any
          (fun piece => piece.1 == 0) = true :=
    List.any_eq_true.mpr ⟨startPiece, hstartPiece, by simp [startPiece]⟩
  simp [boundaryPiecesJoined, htwo, hstarts, hends]

theorem exists_integer_eq_left_of_no_split_start_zero
    (shift : ℚ) (interval piece : ℚ × ℚ)
    (hnoSplit : rotationSplitIndicator shift interval = 0)
    (hpiece : piece ∈ wrapOne (interval.1 + shift, interval.2 + shift))
    (hstart : piece.1 = 0) :
    ∃ integer : ℤ, interval.1 + shift = (integer : ℚ) := by
  have hend :
      interval.2 + shift -
        (⌊interval.1 + shift⌋ : ℚ) ≤ 1 := by
    by_contra hnot
    have : rotationSplitIndicator shift interval = 1 := by
      simp [rotationSplitIndicator, wrapSplitIndicator, hnot]
    omega
  unfold wrapOne at hpiece
  simp only [hend, if_true] at hpiece
  rcases List.mem_singleton.mp hpiece with rfl
  refine ⟨⌊interval.1 + shift⌋, ?_⟩
  change interval.1 + shift -
    (⌊interval.1 + shift⌋ : ℚ) = 0 at hstart
  linarith

theorem exists_integer_eq_right_of_no_split_end_one
    (shift : ℚ) (interval piece : ℚ × ℚ)
    (hnoSplit : rotationSplitIndicator shift interval = 0)
    (hpiece : piece ∈ wrapOne (interval.1 + shift, interval.2 + shift))
    (hendPoint : piece.2 = 1) :
    ∃ integer : ℤ, interval.2 + shift = (integer : ℚ) := by
  have hend :
      interval.2 + shift -
        (⌊interval.1 + shift⌋ : ℚ) ≤ 1 := by
    by_contra hnot
    have : rotationSplitIndicator shift interval = 1 := by
      simp [rotationSplitIndicator, wrapSplitIndicator, hnot]
    omega
  unfold wrapOne at hpiece
  simp only [hend, if_true] at hpiece
  rcases List.mem_singleton.mp hpiece with rfl
  refine ⟨⌊interval.1 + shift⌋ + 1, ?_⟩
  dsimp at hendPoint
  push_cast
  linarith

theorem rotationSplitCount_ne_zero_of_boundary_join
    (shift : ℚ) {region : Region}
    (hnorm : Norm region) (hunit : RegionInUnit region)
    (hnoSeams : RegionHasNoAdjacentSeams region)
    (hsourceBoundary : boundaryPiecesJoined region = false)
    (htargetBoundary :
      boundaryPiecesJoined (translateCirc shift region) = true) :
    rotationSplitCount shift region ≠ 0 := by
  intro hsplitZero
  have htargetData :
      2 ≤ intervalComponentCount (translateCirc shift region) ∧
      (translateCirc shift region).any
          (fun piece => piece.1 == 0) = true ∧
      (translateCirc shift region).any
          (fun piece => piece.2 == 1) = true := by
    have htargetData := htargetBoundary
    simp only [boundaryPiecesJoined, Bool.and_eq_true,
      decide_eq_true_eq] at htargetData
    exact ⟨htargetData.1.1, htargetData.1.2, htargetData.2⟩
  obtain ⟨startPiece, hstartPiece, hstartEq⟩ :=
    List.any_eq_true.mp htargetData.2.1
  obtain ⟨endPiece, hendPiece, hendEq⟩ :=
    List.any_eq_true.mp htargetData.2.2
  have hstartCoordinate : startPiece.1 = 0 := by
    simpa using hstartEq
  have hendCoordinate : endPiece.2 = 1 := by
    simpa using hendEq
  unfold translateCirc wrap translate at hstartPiece hendPiece
  rw [List.mem_flatMap] at hstartPiece hendPiece
  obtain ⟨startTranslated, hstartTranslated, hstartWrap⟩ := hstartPiece
  obtain ⟨endTranslated, hendTranslated, hendWrap⟩ := hendPiece
  rw [List.mem_map] at hstartTranslated hendTranslated
  obtain ⟨startSource, hstartSource, hstartTranslatedEq⟩ := hstartTranslated
  obtain ⟨endSource, hendSource, hendTranslatedEq⟩ := hendTranslated
  subst startTranslated
  subst endTranslated
  have hstartNoSplit :=
    rotationSplitIndicator_eq_zero_of_count_zero
      shift hsplitZero hstartSource
  have hendNoSplit :=
    rotationSplitIndicator_eq_zero_of_count_zero
      shift hsplitZero hendSource
  obtain ⟨startInteger, hstartInteger⟩ :=
    exists_integer_eq_left_of_no_split_start_zero
      shift startSource startPiece hstartNoSplit hstartWrap hstartCoordinate
  obtain ⟨endInteger, hendInteger⟩ :=
    exists_integer_eq_right_of_no_split_end_one
      shift endSource endPiece hendNoSplit hendWrap hendCoordinate
  have hstartUnit := hunit startSource hstartSource
  have hendUnit := hunit endSource hendSource
  have hstartLive := regionLive_of_norm hnorm startSource hstartSource
  have hendLive := regionLive_of_norm hnorm endSource hendSource
  have hintegerGapLowerRat :
      (-1 : ℚ) < (endInteger : ℚ) - (startInteger : ℚ) := by
    linarith
  have hintegerGapUpperRat :
      (endInteger : ℚ) - (startInteger : ℚ) ≤ 1 := by
    linarith
  have hintegerGapLower : -1 < endInteger - startInteger := by
    exact_mod_cast hintegerGapLowerRat
  have hintegerGapUpper : endInteger - startInteger ≤ 1 := by
    exact_mod_cast hintegerGapUpperRat
  have hintegerCases :
      endInteger - startInteger = 0 ∨
        endInteger - startInteger = 1 := by
    omega
  rcases hintegerCases with hsameInteger | hnextInteger
  · have hsameIntegerRat :
        (endInteger : ℚ) - (startInteger : ℚ) = 0 := by
      exact_mod_cast hsameInteger
    have hadjacent : endSource.2 = startSource.1 := by
      linarith
    have hdistinct : endSource ≠ startSource := by
      intro hequal
      subst endSource
      linarith
    exact hnoSeams endSource hendSource startSource hstartSource
      hdistinct hadjacent
  · have hnextIntegerRat :
        (endInteger : ℚ) - (startInteger : ℚ) = 1 := by
      exact_mod_cast hnextInteger
    have hstartAtZero : startSource.1 = 0 := by
      linarith
    have hendAtOne : endSource.2 = 1 := by
      linarith
    have hsourceTwo : 2 ≤ intervalComponentCount region := by
      have hcount := intervalComponentCount_translateCirc shift region
      rw [hsplitZero, Nat.add_zero] at hcount
      linarith
    have hsourceStarts :
        region.any (fun piece => piece.1 == 0) = true :=
      List.any_eq_true.mpr
        ⟨startSource, hstartSource, by simp [hstartAtZero]⟩
    have hsourceEnds :
        region.any (fun piece => piece.2 == 1) = true :=
      List.any_eq_true.mpr
        ⟨endSource, hendSource, by simp [hendAtOne]⟩
    have : boundaryPiecesJoined region = true := by
      simp [boundaryPiecesJoined, hsourceTwo, hsourceStarts, hsourceEnds]
    rw [hsourceBoundary] at this
    contradiction

/-- A normalized unit-chart presentation with neither internal adjacent seams
nor an old boundary join satisfies the exact moving-cut balance automatically. -/
theorem boundaryFaithfulRotation_of_noAdjacentSeams
    (shift : ℚ) {region : Region}
    (hnorm : Norm region) (hunit : RegionInUnit region)
    (hnoSeams : RegionHasNoAdjacentSeams region)
    (hsourceBoundary : boundaryPiecesJoined region = false) :
    BoundaryFaithfulRotation shift region := by
  unfold BoundaryFaithfulRotation
  have hle := rotationSplitCount_le_one shift hnorm hunit
  rcases Nat.le_one_iff_eq_zero_or_eq_one.mp hle with hzero | hone
  · have htargetBoundary :
        boundaryPiecesJoined (translateCirc shift region) = false := by
      by_contra hnot
      have htrue :
          boundaryPiecesJoined (translateCirc shift region) = true := by
        exact Bool.eq_true_of_not_eq_false hnot
      exact rotationSplitCount_ne_zero_of_boundary_join shift hnorm hunit
        hnoSeams hsourceBoundary htrue hzero
    simp [hzero, boundaryCorrection, htargetBoundary]
  · obtain ⟨interval, hinterval, hsplit⟩ :=
      exists_rotationSplitIndicator_eq_one_of_count_eq_one shift hone
    have htargetBoundary :=
      boundaryPiecesJoined_translateCirc_eq_true_of_split
        shift hnorm hunit hinterval hsplit
    simp [hone, boundaryCorrection, htargetBoundary]

theorem regionSeamFree_deleteAnchoredTooth
    (region : Region) (tooth : AnchoredCircularTooth)
    (hunit : RegionInUnit region) (hseamFree : RegionSeamFree region) :
    RegionSeamFree (deleteAnchoredTooth region tooth) := by
  unfold deleteAnchoredTooth diff1F RegionSeamFree
  rw [List.pairwise_flatMap]
  constructor
  · intro interval hinterval
    have hintervalUnit := hunit interval hinterval
    unfold cutF cut
    have hleftDead :
        ¬ (interval.1 < min interval.2 0) := by
      linarith [min_le_right interval.2 0]
    rw [List.filter_cons_of_neg (by simpa using hleftDead)]
    by_cases hrightLive : max interval.1 tooth.width < interval.2
    · rw [List.filter_cons_of_pos (by simpa using hrightLive),
          List.filter_nil]
      simp
    · rw [List.filter_cons_of_neg (by simpa using hrightLive),
          List.filter_nil]
      simp
  · apply hseamFree.imp
    intro left right hstrict leftPiece hleftPiece rightPiece hrightPiece
    have hleftBounds := cutF_piece_bounds hleftPiece
    have hrightBounds := cutF_piece_bounds hrightPiece
    linarith

def RegionStartsPositive (region : Region) : Prop :=
  ∀ interval ∈ region, 0 < interval.1

theorem regionStartsPositive_deleteAnchoredTooth
    (region : Region) (tooth : AnchoredCircularTooth)
    (hunit : RegionInUnit region) :
    RegionStartsPositive (deleteAnchoredTooth region tooth) := by
  intro interval hinterval
  exact deleteAnchoredTooth_piece_start_pos region tooth hunit hinterval

theorem exists_integer_eq_source_right_of_piece_end_lt_one
    {source piece : ℚ × ℚ}
    (hpiece : piece ∈ wrapOne source) (hend : piece.2 < 1) :
    ∃ integer : ℤ, source.2 = piece.2 + integer := by
  unfold wrapOne at hpiece
  by_cases hoverhang : source.2 - (⌊source.1⌋ : ℚ) ≤ 1
  · simp only [hoverhang, if_true] at hpiece
    rcases List.mem_singleton.mp hpiece with rfl
    exact ⟨⌊source.1⌋, by dsimp; ring⟩
  · simp only [hoverhang, if_false] at hpiece
    rcases List.mem_cons.mp hpiece with rfl | hpiece
    · simp at hend
    · rcases List.mem_singleton.mp hpiece with rfl
      refine ⟨⌊source.1⌋ + 1, ?_⟩
      push_cast
      ring

theorem exists_integer_eq_source_left_of_piece_start_pos
    {source piece : ℚ × ℚ}
    (hpiece : piece ∈ wrapOne source) (hstart : 0 < piece.1) :
    ∃ integer : ℤ, source.1 = piece.1 + integer := by
  unfold wrapOne at hpiece
  by_cases hoverhang : source.2 - (⌊source.1⌋ : ℚ) ≤ 1
  · simp only [hoverhang, if_true] at hpiece
    rcases List.mem_singleton.mp hpiece with rfl
    exact ⟨⌊source.1⌋, (Int.fract_add_floor source.1).symm⟩
  · simp only [hoverhang, if_false] at hpiece
    rcases List.mem_cons.mp hpiece with rfl | hpiece
    · exact ⟨⌊source.1⌋, (Int.fract_add_floor source.1).symm⟩
    · rcases List.mem_singleton.mp hpiece with rfl
      simp at hstart

theorem regionHasNoAdjacentSeams_translateCirc
    (shift : ℚ) {region : Region}
    (hnorm : Norm region) (hunit : RegionInUnit region)
    (hnoSeams : RegionHasNoAdjacentSeams region)
    (hpositive : RegionStartsPositive region) :
    RegionHasNoAdjacentSeams (translateCirc shift region) := by
  intro leftPiece hleftPiece rightPiece hrightPiece hpChanged hadjacent
  have htranslatedUnit := regionInUnit_translateCirc shift region hunit
  have htranslatedLive :=
    regionLive_translateCirc shift region (regionLive_of_norm hnorm)
  have hleftPieceUnit := htranslatedUnit leftPiece hleftPiece
  have hrightPieceUnit := htranslatedUnit rightPiece hrightPiece
  have hleftPieceLive := htranslatedLive leftPiece hleftPiece
  have hrightPieceLive := htranslatedLive rightPiece hrightPiece
  have hsharedPositive : 0 < leftPiece.2 := by
    linarith
  have hsharedLtOne : rightPiece.1 < 1 := by
    linarith
  unfold translateCirc wrap translate at hleftPiece hrightPiece
  rw [List.mem_flatMap] at hleftPiece hrightPiece
  obtain ⟨leftTranslated, hleftTranslated, hleftWrap⟩ := hleftPiece
  obtain ⟨rightTranslated, hrightTranslated, hrightWrap⟩ := hrightPiece
  rw [List.mem_map] at hleftTranslated hrightTranslated
  obtain ⟨leftSource, hleftSource, hleftTranslatedEq⟩ := hleftTranslated
  obtain ⟨rightSource, hrightSource, hrightTranslatedEq⟩ := hrightTranslated
  subst leftTranslated
  subst rightTranslated
  obtain ⟨leftInteger, hleftInteger⟩ :=
    exists_integer_eq_source_right_of_piece_end_lt_one
      hleftWrap (by linarith)
  obtain ⟨rightInteger, hrightInteger⟩ :=
    exists_integer_eq_source_left_of_piece_start_pos
      hrightWrap (by linarith)
  have hleftUnit := hunit leftSource hleftSource
  have hrightUnit := hunit rightSource hrightSource
  have hleftLive := regionLive_of_norm hnorm leftSource hleftSource
  have hrightLive := regionLive_of_norm hnorm rightSource hrightSource
  have hintegerGapLowerRat :
      (-1 : ℚ) ≤ (rightInteger : ℚ) - (leftInteger : ℚ) := by
    linarith
  have hintegerGapUpperRat :
      (rightInteger : ℚ) - (leftInteger : ℚ) < 1 := by
    linarith
  have hintegerGapLower : -1 ≤ rightInteger - leftInteger := by
    exact_mod_cast hintegerGapLowerRat
  have hintegerGapUpper : rightInteger - leftInteger < 1 := by
    exact_mod_cast hintegerGapUpperRat
  have hintegerCases :
      rightInteger - leftInteger = -1 ∨
        rightInteger - leftInteger = 0 := by
    omega
  rcases hintegerCases with hboundaryInteger | hsameInteger
  · have hboundaryIntegerRat :
        (rightInteger : ℚ) - (leftInteger : ℚ) = -1 := by
      exact_mod_cast hboundaryInteger
    have hrightAtZero : rightSource.1 = 0 := by
      linarith
    exact (ne_of_gt (hpositive rightSource hrightSource)) hrightAtZero
  · have hsameIntegerRat :
        (rightInteger : ℚ) - (leftInteger : ℚ) = 0 := by
      exact_mod_cast hsameInteger
    have hsourceAdjacent : leftSource.2 = rightSource.1 := by
      linarith
    have hsourceDistinct : leftSource ≠ rightSource := by
      intro hequal
      subst rightSource
      linarith
    exact hnoSeams leftSource hleftSource rightSource hrightSource
      hsourceDistinct hsourceAdjacent

theorem regionHasNoAdjacentSeams_of_perm {left right : Region}
    (hperm : left.Perm right) (hnoSeams : RegionHasNoAdjacentSeams right) :
    RegionHasNoAdjacentSeams left := by
  intro first hfirst second hsecond hne
  exact hnoSeams first (hperm.mem_iff.mp hfirst)
    second (hperm.mem_iff.mp hsecond) hne

theorem regionSeamFree_of_norm_noAdjacentSeams :
    ∀ {region : Region}, Norm region → RegionHasNoAdjacentSeams region →
      RegionSeamFree region := by
  intro region
  induction region with
  | nil => simp [RegionSeamFree]
  | cons head tail inductionHypothesis =>
      intro hnorm hnoSeams
      unfold RegionSeamFree
      rw [List.pairwise_cons]
      constructor
      · intro later hlater
        have horder := norm_head_le hnorm later hlater
        have hdistinct : head ≠ later := by
          intro hequal
          subst later
          linarith [norm_head_lt hnorm]
        have hnotAdjacent := hnoSeams head (List.mem_cons_self ..)
          later (List.mem_cons_of_mem _ hlater) hdistinct
        exact lt_of_le_of_ne horder hnotAdjacent
      · apply inductionHypothesis (norm_tail hnorm)
        intro left hleft right hright hne
        exact hnoSeams left (List.mem_cons_of_mem _ hleft)
          right (List.mem_cons_of_mem _ hright) hne

theorem regionSeamFree_sortedTranslateCirc
    (shift : ℚ) {region : Region}
    (hnorm : Norm region) (hunit : RegionInUnit region)
    (hnoSeams : RegionHasNoAdjacentSeams region)
    (hpositive : RegionStartsPositive region) :
    RegionSeamFree (sortedTranslateCirc shift region) := by
  have htranslatedSeparated :=
    regionSeparated_translateCirc_of_norm shift region hnorm hunit
  have hsortedNorm := norm_sortedTranslateCirc_of_separated shift region
    (regionLive_of_norm hnorm) htranslatedSeparated
  apply regionSeamFree_of_norm_noAdjacentSeams hsortedNorm
  apply regionHasNoAdjacentSeams_of_perm (sortedTranslateCirc_perm shift region)
  exact regionHasNoAdjacentSeams_translateCirc shift hnorm hunit
    hnoSeams hpositive

/-- The existing survivor recursion is already canonical after the first
anchored deletion: survivors and charts have no artificial internal seams,
and every positive survivor starts strictly to the right of the chart cut. -/
theorem rationalCircle_seamFree_invariants
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth) :
    ∀ toothCount,
      RegionSeamFree (rationalCircleChart shift tooth toothCount) ∧
      RegionSeamFree (rationalCircleSurvivor shift tooth toothCount) ∧
      (1 ≤ toothCount →
        RegionStartsPositive
          (rationalCircleSurvivor shift tooth toothCount)) := by
  intro toothCount
  induction toothCount with
  | zero =>
      simp [rationalCircleChart, rationalCircleRechart,
        rationalCircleSurvivor, RegionSeamFree]
  | succ toothCount inductionHypothesis =>
      have hsurvivorSeamFree :
          RegionSeamFree
            (rationalCircleSurvivor shift tooth (toothCount + 1)) := by
        rw [rationalCircleSurvivor_succ]
        exact regionSeamFree_deleteAnchoredTooth _ _
          (regionInUnit_rationalCircleChart shift tooth toothCount)
          inductionHypothesis.1
      have hsurvivorPositive :
          RegionStartsPositive
            (rationalCircleSurvivor shift tooth (toothCount + 1)) := by
        rw [rationalCircleSurvivor_succ]
        exact regionStartsPositive_deleteAnchoredTooth _ _
          (regionInUnit_rationalCircleChart shift tooth toothCount)
      have hsurvivorNoSeams :
          RegionHasNoAdjacentSeams
            (rationalCircleSurvivor shift tooth (toothCount + 1)) :=
        regionHasNoAdjacentSeams_of_norm_seamFree
          (norm_rationalCircleSurvivor shift tooth (toothCount + 1))
          hsurvivorSeamFree
      have hchartSeamFree :
          RegionSeamFree
            (rationalCircleChart shift tooth (toothCount + 1)) := by
        change RegionSeamFree
          (sortedTranslateCirc (shift (toothCount + 1))
            (rationalCircleSurvivor shift tooth (toothCount + 1)))
        exact regionSeamFree_sortedTranslateCirc _
          (norm_rationalCircleSurvivor shift tooth (toothCount + 1))
          (regionInUnit_rationalCircleSurvivor shift tooth (toothCount + 1))
          hsurvivorNoSeams hsurvivorPositive
      exact ⟨hchartSeamFree, hsurvivorSeamFree,
        fun _ => hsurvivorPositive⟩

theorem regionSeamFree_rationalCircleSurvivor
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    RegionSeamFree (rationalCircleSurvivor shift tooth toothCount) :=
  (rationalCircle_seamFree_invariants shift tooth toothCount).2.1

theorem regionHasNoAdjacentSeams_rationalCircleSurvivor
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    RegionHasNoAdjacentSeams
      (rationalCircleSurvivor shift tooth toothCount) :=
  regionHasNoAdjacentSeams_of_norm_seamFree
    (norm_rationalCircleSurvivor shift tooth toothCount)
    (regionSeamFree_rationalCircleSurvivor shift tooth toothCount)

theorem boundaryFaithfulRotation_rationalCircleSurvivor
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) (hpositive : 1 ≤ toothCount) :
    BoundaryFaithfulRotation (shift toothCount)
      (rationalCircleSurvivor shift tooth toothCount) :=
  boundaryFaithfulRotation_of_noAdjacentSeams _
    (norm_rationalCircleSurvivor shift tooth toothCount)
    (regionInUnit_rationalCircleSurvivor shift tooth toothCount)
    (regionHasNoAdjacentSeams_rationalCircleSurvivor
      shift tooth toothCount)
    (boundaryPiecesJoined_rationalCircleSurvivor_eq_false
      shift tooth toothCount hpositive)

theorem positiveRotationTopologyCertificate_rationalCircleSurvivor
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) (hpositive : 1 ≤ toothCount) :
    PositiveRotationTopologyCertificate (shift toothCount)
      (rationalCircleSurvivor shift tooth toothCount) :=
  (positiveRotationTopologyCertificate_rationalCircleSurvivor_iff
    shift tooth toothCount).mpr
      (boundaryFaithfulRotation_rationalCircleSurvivor
        shift tooth toothCount hpositive)

/-- Hypothesis-free rational circle atlas: canonical seam preservation closes
the last external topology certificate. -/
def canonicalRationalCircularToothAtlas
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth) :
    CircularToothAtlas :=
  rationalCircularToothAtlasOfBoundaryFaithfulRotations shift tooth
    fun toothCount hpositive =>
      boundaryFaithfulRotation_rationalCircleSurvivor
        shift tooth toothCount hpositive

/-- The concrete rational survivor component cap now has no topology
hypothesis. -/
theorem canonical_rational_circle_component_count_le_tooth_count
    (shift : ℕ → ℚ) (tooth : ℕ → AnchoredCircularTooth)
    (toothCount : ℕ) :
    1 ≤ toothCount →
      circularComponentCount
          (rationalCircleSurvivor shift tooth toothCount) ≤ toothCount :=
  circular_component_count_le_tooth_count
    (canonicalRationalCircularToothAtlas shift tooth) toothCount

/-! ## Axiom audit -/

#print axioms coalesceAdjacent_spec
#print axioms boundaryFaithfulRotation_of_noAdjacentSeams
#print axioms rationalCircle_seamFree_invariants
#print axioms boundaryFaithfulRotation_rationalCircleSurvivor
#print axioms canonical_rational_circle_component_count_le_tooth_count

end LRC14.LocalDensityBlockGluing
