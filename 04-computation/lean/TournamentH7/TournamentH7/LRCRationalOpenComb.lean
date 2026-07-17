import TournamentH7.LRCOpenDangerComb
import TournamentH7.LRCStrictInterMerge

/-!
Normalized rational strict-open danger combs and their exact two-pointer pair
intersection.  The rational endpoints expose exact grid arithmetic while the
cast lemmas connect the lists to the real comb carrier.
-/

namespace LonelyRunner
namespace LRCRationalOpenComb

open RatIntervals LRCStrictInterMerge LRCOpenDangerComb MeasureTheory

noncomputable section

abbrev Interval := ℚ × ℚ

/-- Rational version of one clipped strict tooth. -/
def ratOpenCombInterval (w k : ℕ) : Interval :=
  (max 0 ((((k : ℕ) : ℚ) - (1 : ℚ) / 14) / (w : ℚ)),
    min 1 ((((k : ℕ) : ℚ) + (1 : ℚ) / 14) / (w : ℚ)))

/-- Consecutive rational teeth, starting at `start`. -/
def ratOpenCombTail (w : ℕ) : ℕ → ℕ → Region
  | _, 0 => []
  | start, count + 1 =>
      ratOpenCombInterval w start :: ratOpenCombTail w (start + 1) count

/-- The full `w+1`-tooth strict rational comb on `(0,1)`. -/
def ratOpenCombRegion (w : ℕ) : Region :=
  ratOpenCombTail w 0 (w + 1)

theorem ratOpenCombInterval_live
    (w : ℕ) (hw : 0 < w) (k : ℕ) (hk : k ≤ w) :
    (ratOpenCombInterval w k).1 < (ratOpenCombInterval w k).2 := by
  have hwQ : (0 : ℚ) < (w : ℚ) := by exact_mod_cast hw
  have hkQ : (k : ℚ) ≤ (w : ℚ) := by exact_mod_cast hk
  have hupperPos :
      (0 : ℚ) < (((k : ℚ) + (1 : ℚ) / 14) / (w : ℚ)) := by
    positivity
  have hlowerOne :
      (((k : ℚ) - (1 : ℚ) / 14) / (w : ℚ)) < 1 := by
    rw [div_lt_iff₀ hwQ]
    linarith
  have hlowerUpper :
      (((k : ℚ) - (1 : ℚ) / 14) / (w : ℚ)) <
        (((k : ℚ) + (1 : ℚ) / 14) / (w : ℚ)) := by
    apply (div_lt_div_iff_of_pos_right hwQ).2
    linarith [show (0 : ℚ) < (1 : ℚ) / 14 by norm_num]
  unfold ratOpenCombInterval
  apply max_lt
  · apply lt_min
    · norm_num
    · exact hupperPos
  · apply lt_min
    · exact hlowerOne
    · exact hlowerUpper

theorem ratOpenCombInterval_separated
    (w : ℕ) (hw : 0 < w) (k : ℕ) :
    (ratOpenCombInterval w k).2 ≤
      (ratOpenCombInterval w (k + 1)).1 := by
  have hwQ : (0 : ℚ) < (w : ℚ) := by exact_mod_cast hw
  unfold ratOpenCombInterval
  calc
    min 1 (((k : ℚ) + (1 : ℚ) / 14) / (w : ℚ))
        ≤ (((k : ℚ) + (1 : ℚ) / 14) / (w : ℚ)) := min_le_right _ _
    _ ≤ ((((k + 1 : ℕ) : ℚ) - (1 : ℚ) / 14) / (w : ℚ)) := by
      apply (div_le_div_iff_of_pos_right hwQ).2
      push_cast
      linarith [show (1 : ℚ) / 14 + (1 : ℚ) / 14 ≤ 1 by norm_num]
    _ ≤ max 0 ((((k + 1 : ℕ) : ℚ) - (1 : ℚ) / 14) / (w : ℚ)) :=
      le_max_right _ _

theorem length_ratOpenCombTail (w start count : ℕ) :
    (ratOpenCombTail w start count).length = count := by
  induction count generalizing start with
  | zero => rfl
  | succ count ih =>
      simp [ratOpenCombTail, ih]

theorem length_ratOpenCombRegion (w : ℕ) :
    (ratOpenCombRegion w).length = w + 1 := by
  exact length_ratOpenCombTail w 0 (w + 1)

theorem strictMem_ratOpenCombTail_iff
    (x : ℚ) (w start count : ℕ) :
    strictMem x (ratOpenCombTail w start count) ↔
      ∃ offset < count,
        inside x (ratOpenCombInterval w (start + offset)) := by
  induction count generalizing start with
  | zero => simp [ratOpenCombTail, strictMem]
  | succ count ih =>
      rw [ratOpenCombTail, strictMem_cons_iff, ih]
      constructor
      · rintro (hhead | ⟨offset, hoffset, hinside⟩)
        · exact ⟨0, by omega, by simpa using hhead⟩
        · refine ⟨offset + 1, by omega, ?_⟩
          simpa [Nat.add_assoc, Nat.add_comm, Nat.add_left_comm] using hinside
      · rintro ⟨offset, hoffset, hinside⟩
        cases offset with
        | zero =>
            exact Or.inl (by simpa using hinside)
        | succ offset =>
            apply Or.inr
            refine ⟨offset, by omega, ?_⟩
            simpa [Nat.add_assoc, Nat.add_comm, Nat.add_left_comm] using hinside

theorem norm_ratOpenCombTail
    (w : ℕ) (hw : 0 < w) (start count : ℕ)
    (hbound : start + count ≤ w + 1) :
    RatIntervals.Norm (ratOpenCombTail w start count) := by
  induction count generalizing start with
  | zero => trivial
  | succ count ih =>
      cases count with
      | zero =>
          simp only [ratOpenCombTail]
          exact ratOpenCombInterval_live w hw start (by omega)
      | succ count =>
          simp only [ratOpenCombTail, RatIntervals.Norm]
          refine ⟨ratOpenCombInterval_live w hw start (by omega),
            ratOpenCombInterval_separated w hw start, ?_⟩
          exact ih (start + 1) (by omega)

theorem norm_ratOpenCombRegion (w : ℕ) (hw : 0 < w) :
    RatIntervals.Norm (ratOpenCombRegion w) := by
  exact norm_ratOpenCombTail w hw 0 (w + 1) (by omega)

theorem cast_ratOpenCombInterval_left
    (w k : ℕ) (hk : k < w + 1) :
    (((ratOpenCombInterval w k).1 : ℚ) : ℝ) =
      openCombLeft w ⟨k, hk⟩ := by
  simp [ratOpenCombInterval, openCombLeft]

theorem cast_ratOpenCombInterval_right
    (w k : ℕ) (hk : k < w + 1) :
    (((ratOpenCombInterval w k).2 : ℚ) : ℝ) =
      openCombRight w ⟨k, hk⟩ := by
  simp [ratOpenCombInterval, openCombRight]

/-- Exact carrier equality between the normalized rational list and the real
finite comb, at rational test points. -/
theorem strictMem_ratOpenCombRegion_iff_mem_real
    (w : ℕ) (x : ℚ) :
    strictMem x (ratOpenCombRegion w) ↔
      (x : ℝ) ∈ openCombRegion w := by
  rw [ratOpenCombRegion, strictMem_ratOpenCombTail_iff]
  simp only [zero_add]
  constructor
  · rintro ⟨k, hk, hinside⟩
    apply Set.mem_iUnion.mpr
    refine ⟨⟨k, by omega⟩, ?_, ?_⟩
    · rw [← cast_ratOpenCombInterval_left w k (by omega)]
      exact_mod_cast hinside.1
    · rw [← cast_ratOpenCombInterval_right w k (by omega)]
      exact_mod_cast hinside.2
  · intro hreal
    obtain ⟨k, hk⟩ := Set.mem_iUnion.mp hreal
    refine ⟨(k : ℕ), by omega, ?_, ?_⟩
    · have hleftReal :
          (((ratOpenCombInterval w (k : ℕ)).1 : ℚ) : ℝ) < (x : ℝ) := by
        rw [cast_ratOpenCombInterval_left w k (by omega)]
        exact hk.1
      exact_mod_cast hleftReal
    · have hrightReal :
          (x : ℝ) < (((ratOpenCombInterval w (k : ℕ)).2 : ℚ) : ℝ) := by
        rw [cast_ratOpenCombInterval_right w k (by omega)]
        exact hk.2
      exact_mod_cast hrightReal

/-- Two-pointer strict intersection of two normalized rational danger combs. -/
def ratOpenPairRegion (w₁ w₂ : ℕ) : Region :=
  strictInterMerge (ratOpenCombRegion w₁) (ratOpenCombRegion w₂)

theorem norm_ratOpenPairRegion
    (w₁ w₂ : ℕ) (hw₁ : 0 < w₁) (hw₂ : 0 < w₂) :
    RatIntervals.Norm (ratOpenPairRegion w₁ w₂) := by
  exact norm_strictInterMerge _ _
    (norm_ratOpenCombRegion w₁ hw₁)
    (norm_ratOpenCombRegion w₂ hw₂)

theorem length_ratOpenPairRegion_le
    (w₁ w₂ : ℕ) :
    (ratOpenPairRegion w₁ w₂).length ≤ w₁ + w₂ + 1 := by
  have h₁ : ratOpenCombRegion w₁ ≠ [] := by
    intro h
    have := length_ratOpenCombRegion w₁
    rw [h] at this
    simp at this
  have h₂ : ratOpenCombRegion w₂ ≠ [] := by
    intro h
    have := length_ratOpenCombRegion w₂
    rw [h] at this
    simp at this
  have hmerge := length_strictInterMerge_le_add_sub_one
    (ratOpenCombRegion w₁) (ratOpenCombRegion w₂) h₁ h₂
  unfold ratOpenPairRegion
  rw [length_ratOpenCombRegion, length_ratOpenCombRegion] at hmerge
  omega

theorem strictMem_ratOpenPairRegion_iff
    (w₁ w₂ : ℕ) (hw₁ : 0 < w₁) (hw₂ : 0 < w₂) (x : ℚ) :
    strictMem x (ratOpenPairRegion w₁ w₂) ↔
      strictMem x (ratOpenCombRegion w₁) ∧
        strictMem x (ratOpenCombRegion w₂) := by
  exact strictMem_strictInterMerge_iff x _ _
    (norm_ratOpenCombRegion w₁ hw₁)
    (norm_ratOpenCombRegion w₂ hw₂)

theorem strictMem_ratOpenPairRegion_iff_mem_real
    (w₁ w₂ : ℕ) (hw₁ : 0 < w₁) (hw₂ : 0 < w₂) (x : ℚ) :
    strictMem x (ratOpenPairRegion w₁ w₂) ↔
      (x : ℝ) ∈ openCombRegion w₁ ∩ openCombRegion w₂ := by
  rw [strictMem_ratOpenPairRegion_iff w₁ w₂ hw₁ hw₂,
    strictMem_ratOpenCombRegion_iff_mem_real,
    strictMem_ratOpenCombRegion_iff_mem_real]
  rfl

/-- Strict membership of a rational interval at a real test point. -/
def realInside (x : ℝ) (interval : Interval) : Prop :=
  (interval.1 : ℝ) < x ∧ x < (interval.2 : ℝ)

/-- Strict membership in a rational region, interpreted over `ℝ`. -/
def realStrictMem (x : ℝ) (region : Region) : Prop :=
  ∃ interval ∈ region, realInside x interval

theorem realStrictMem_cons_iff
    {x : ℝ} {interval : Interval} {region : Region} :
    realStrictMem x (interval :: region) ↔
      realInside x interval ∨ realStrictMem x region := by
  simp [realStrictMem]

theorem realInside_clip_iff
    {x : ℝ} {left right : Interval} :
    realInside x (RatIntervals.clip left right) ↔
      realInside x left ∧ realInside x right := by
  simp [realInside, RatIntervals.clip, and_assoc, and_left_comm, and_comm]

theorem real_head_left_drop
    {x : ℝ} {left right : Interval} {leftTail rightTail : Region}
    (hrightNorm : RatIntervals.Norm (right :: rightTail))
    (hend : left.2 ≤ right.2) :
    (realStrictMem x (left :: leftTail) ∧
        realStrictMem x (right :: rightTail)) ↔
      (realInside x left ∧ realInside x right) ∨
        (realStrictMem x leftTail ∧
          realStrictMem x (right :: rightTail)) := by
  constructor
  · rintro ⟨hleft, hright⟩
    rw [realStrictMem_cons_iff] at hleft
    rcases hleft with hleftHead | hleftTail
    · rw [realStrictMem_cons_iff] at hright
      rcases hright with hrightHead | hrightTail
      · exact Or.inl ⟨hleftHead, hrightHead⟩
      · obtain ⟨other, hother, hinside⟩ := hrightTail
        have horderQ := RatIntervals.norm_head_le hrightNorm other hother
        have horderR : (right.2 : ℝ) ≤ (other.1 : ℝ) := by exact_mod_cast horderQ
        have hendR : (left.2 : ℝ) ≤ (right.2 : ℝ) := by exact_mod_cast hend
        exfalso
        linarith [hleftHead.2, hinside.1]
    · exact Or.inr ⟨hleftTail, hright⟩
  · rintro (hheads | htails)
    · exact ⟨realStrictMem_cons_iff.mpr (Or.inl hheads.1),
        realStrictMem_cons_iff.mpr (Or.inl hheads.2)⟩
    · exact ⟨realStrictMem_cons_iff.mpr (Or.inr htails.1), htails.2⟩

theorem real_head_right_drop
    {x : ℝ} {left right : Interval} {leftTail rightTail : Region}
    (hleftNorm : RatIntervals.Norm (left :: leftTail))
    (hend : right.2 < left.2) :
    (realStrictMem x (left :: leftTail) ∧
        realStrictMem x (right :: rightTail)) ↔
      (realInside x left ∧ realInside x right) ∨
        (realStrictMem x (left :: leftTail) ∧
          realStrictMem x rightTail) := by
  constructor
  · rintro ⟨hleft, hright⟩
    rw [realStrictMem_cons_iff] at hright
    rcases hright with hrightHead | hrightTail
    · rw [realStrictMem_cons_iff] at hleft
      rcases hleft with hleftHead | hleftTail
      · exact Or.inl ⟨hleftHead, hrightHead⟩
      · obtain ⟨other, hother, hinside⟩ := hleftTail
        have horderQ := RatIntervals.norm_head_le hleftNorm other hother
        have horderR : (left.2 : ℝ) ≤ (other.1 : ℝ) := by exact_mod_cast horderQ
        have hendR : (right.2 : ℝ) < (left.2 : ℝ) := by exact_mod_cast hend
        exfalso
        linarith [hrightHead.2, hinside.1]
    · exact Or.inr ⟨hleft, hrightTail⟩
  · rintro (hheads | htails)
    · exact ⟨realStrictMem_cons_iff.mpr (Or.inl hheads.1),
        realStrictMem_cons_iff.mpr (Or.inl hheads.2)⟩
    · exact ⟨htails.1, realStrictMem_cons_iff.mpr (Or.inr htails.2)⟩

/-- The strict two-pointer merge has the same exact carrier at every real
point, not only at rational test points. -/
theorem realStrictMem_strictInterMerge_iff (x : ℝ) :
    ∀ (left right : Region),
      RatIntervals.Norm left → RatIntervals.Norm right →
      (realStrictMem x (strictInterMerge left right) ↔
        realStrictMem x left ∧ realStrictMem x right)
  | [], right, _hleft, _hright => by simp [realStrictMem, strictInterMerge]
  | left, [], _hleft, _hright => by
      cases left <;> simp [realStrictMem, strictInterMerge]
  | left :: leftTail, right :: rightTail, hleftNorm, hrightNorm => by
      rw [strictInterMerge]
      by_cases hend : left.2 ≤ right.2
      · rw [if_pos hend]
        by_cases hlive : (RatIntervals.clip left right).1 <
            (RatIntervals.clip left right).2
        · rw [if_pos hlive, realStrictMem_cons_iff, realInside_clip_iff,
            realStrictMem_strictInterMerge_iff x leftTail (right :: rightTail)
              (RatIntervals.norm_tail hleftNorm) hrightNorm,
            real_head_left_drop hrightNorm hend]
        · rw [if_neg hlive,
            realStrictMem_strictInterMerge_iff x leftTail (right :: rightTail)
              (RatIntervals.norm_tail hleftNorm) hrightNorm,
            real_head_left_drop hrightNorm hend]
          constructor
          · exact Or.inr
          · rintro (hheads | htails)
            · have hoverlap : realInside x (RatIntervals.clip left right) :=
                realInside_clip_iff.mpr hheads
              have hliveR :
                  ((RatIntervals.clip left right).1 : ℝ) <
                    ((RatIntervals.clip left right).2 : ℝ) :=
                hoverlap.1.trans hoverlap.2
              have hliveQ : (RatIntervals.clip left right).1 <
                  (RatIntervals.clip left right).2 := by exact_mod_cast hliveR
              exact False.elim (hlive hliveQ)
            · exact htails
      · rw [if_neg hend]
        have hend' : right.2 < left.2 := lt_of_not_ge hend
        by_cases hlive : (RatIntervals.clip left right).1 <
            (RatIntervals.clip left right).2
        · rw [if_pos hlive, realStrictMem_cons_iff, realInside_clip_iff,
            realStrictMem_strictInterMerge_iff x (left :: leftTail) rightTail
              hleftNorm (RatIntervals.norm_tail hrightNorm),
            real_head_right_drop hleftNorm hend']
        · rw [if_neg hlive,
            realStrictMem_strictInterMerge_iff x (left :: leftTail) rightTail
              hleftNorm (RatIntervals.norm_tail hrightNorm),
            real_head_right_drop hleftNorm hend']
          constructor
          · exact Or.inr
          · rintro (hheads | htails)
            · have hoverlap : realInside x (RatIntervals.clip left right) :=
                realInside_clip_iff.mpr hheads
              have hliveR :
                  ((RatIntervals.clip left right).1 : ℝ) <
                    ((RatIntervals.clip left right).2 : ℝ) :=
                hoverlap.1.trans hoverlap.2
              have hliveQ : (RatIntervals.clip left right).1 <
                  (RatIntervals.clip left right).2 := by exact_mod_cast hliveR
              exact False.elim (hlive hliveQ)
            · exact htails
termination_by left right => left.length + right.length
decreasing_by all_goals simp_wf

theorem realStrictMem_ratOpenCombTail_iff
    (x : ℝ) (w start count : ℕ) :
    realStrictMem x (ratOpenCombTail w start count) ↔
      ∃ offset < count,
        realInside x (ratOpenCombInterval w (start + offset)) := by
  induction count generalizing start with
  | zero => simp [ratOpenCombTail, realStrictMem]
  | succ count ih =>
      rw [ratOpenCombTail, realStrictMem_cons_iff, ih]
      constructor
      · rintro (hhead | ⟨offset, hoffset, hinside⟩)
        · exact ⟨0, by omega, by simpa using hhead⟩
        · refine ⟨offset + 1, by omega, ?_⟩
          simpa [Nat.add_assoc, Nat.add_comm, Nat.add_left_comm] using hinside
      · rintro ⟨offset, hoffset, hinside⟩
        cases offset with
        | zero => exact Or.inl (by simpa using hinside)
        | succ offset =>
            apply Or.inr
            refine ⟨offset, by omega, ?_⟩
            simpa [Nat.add_assoc, Nat.add_comm, Nat.add_left_comm] using hinside

theorem realStrictMem_ratOpenCombRegion_iff
    (w : ℕ) (x : ℝ) :
    realStrictMem x (ratOpenCombRegion w) ↔ x ∈ openCombRegion w := by
  rw [ratOpenCombRegion, realStrictMem_ratOpenCombTail_iff]
  simp only [zero_add]
  constructor
  · rintro ⟨k, hk, hinside⟩
    apply Set.mem_iUnion.mpr
    refine ⟨⟨k, by omega⟩, ?_, ?_⟩
    · rw [← cast_ratOpenCombInterval_left w k (by omega)]
      exact hinside.1
    · rw [← cast_ratOpenCombInterval_right w k (by omega)]
      exact hinside.2
  · intro hreal
    obtain ⟨k, hk⟩ := Set.mem_iUnion.mp hreal
    refine ⟨(k : ℕ), by omega, ?_, ?_⟩
    · rw [cast_ratOpenCombInterval_left w k (by omega)]
      exact hk.1
    · rw [cast_ratOpenCombInterval_right w k (by omega)]
      exact hk.2

theorem realStrictMem_ratOpenPairRegion_iff
    (w₁ w₂ : ℕ) (hw₁ : 0 < w₁) (hw₂ : 0 < w₂) (x : ℝ) :
    realStrictMem x (ratOpenPairRegion w₁ w₂) ↔
      x ∈ openCombRegion w₁ ∩ openCombRegion w₂ := by
  rw [ratOpenPairRegion,
    realStrictMem_strictInterMerge_iff x _ _
      (norm_ratOpenCombRegion w₁ hw₁)
      (norm_ratOpenCombRegion w₂ hw₂),
    realStrictMem_ratOpenCombRegion_iff,
    realStrictMem_ratOpenCombRegion_iff]
  rfl

/-- Real carrier of a rational strict interval list. -/
def realRegionCarrier (region : Region) : Set ℝ :=
  {x | realStrictMem x region}

theorem realRegionCarrier_nil : realRegionCarrier [] = ∅ := by
  ext x
  simp [realRegionCarrier, realStrictMem]

theorem realRegionCarrier_cons
    (interval : Interval) (region : Region) :
    realRegionCarrier (interval :: region) =
      Set.Ioo (interval.1 : ℝ) (interval.2 : ℝ) ∪
        realRegionCarrier region := by
  ext x
  simp [realRegionCarrier, realStrictMem_cons_iff, realInside]

theorem measurableSet_realRegionCarrier (region : Region) :
    MeasurableSet (realRegionCarrier region) := by
  induction region with
  | nil => rw [realRegionCarrier_nil]; exact MeasurableSet.empty
  | cons interval region ih =>
      rw [realRegionCarrier_cons]
      exact measurableSet_Ioo.union ih

theorem head_disjoint_realRegionCarrier
    {interval : Interval} {region : Region}
    (hnorm : RatIntervals.Norm (interval :: region)) :
    Disjoint (Set.Ioo (interval.1 : ℝ) (interval.2 : ℝ))
      (realRegionCarrier region) := by
  rw [Set.disjoint_left]
  intro x hxInterval hxRegion
  obtain ⟨other, hother, hxOther⟩ := hxRegion
  have horderQ := RatIntervals.norm_head_le hnorm other hother
  have horderR : (interval.2 : ℝ) ≤ (other.1 : ℝ) := by
    exact_mod_cast horderQ
  linarith [hxInterval.2, hxOther.1]

/-- A normalized rational strict list has Lebesgue volume equal to its
`RatIntervals.length`, after casting to `ℝ`. -/
theorem volume_realRegionCarrier_eq_length
    (region : Region) (hnorm : RatIntervals.Norm region) :
    volume (realRegionCarrier region) =
      ENNReal.ofReal ((RatIntervals.length region : ℚ) : ℝ) := by
  induction region with
  | nil =>
      rw [realRegionCarrier_nil]
      simp [RatIntervals.length]
  | cons interval region ih =>
      have htailNorm := RatIntervals.norm_tail hnorm
      have hheadLiveQ := RatIntervals.norm_head_lt hnorm
      have hheadLiveR : (interval.1 : ℝ) < (interval.2 : ℝ) := by
        exact_mod_cast hheadLiveQ
      have htailNonnegQ := RatIntervals.length_nonneg region
      have htailNonnegR :
          (0 : ℝ) ≤ ((RatIntervals.length region : ℚ) : ℝ) := by
        exact_mod_cast htailNonnegQ
      rw [realRegionCarrier_cons,
        measure_union (head_disjoint_realRegionCarrier hnorm)
          (measurableSet_realRegionCarrier region),
        Real.volume_Ioo, ih htailNorm]
      rw [← ENNReal.ofReal_add (sub_nonneg.mpr hheadLiveR.le) htailNonnegR]
      apply congrArg ENNReal.ofReal
      unfold RatIntervals.length
      simp only [List.map_cons, List.sum_cons]
      have hdiffQ : (0 : ℚ) ≤ interval.2 - interval.1 := sub_nonneg.mpr hheadLiveQ.le
      rw [max_eq_right hdiffQ]
      push_cast
      ring

theorem openCombRegion_subset_Ioo (w : ℕ) :
    openCombRegion w ⊆ Set.Ioo (0 : ℝ) 1 := by
  intro x hx
  obtain ⟨k, hk⟩ := Set.mem_iUnion.mp hx
  constructor
  · exact lt_of_le_of_lt (by unfold openCombLeft; exact le_max_left _ _) hk.1
  · exact lt_of_lt_of_le hk.2 (by unfold openCombRight; exact min_le_left _ _)

/-- Circle pair event for two positive speeds. -/
def positivePairDanger (w₁ w₂ : ℕ) : Set UnitAddCircle :=
  LRCCommensuration.danger (w₁ : ℤ) 0 (1 / 14) ∩
    LRCCommensuration.danger (w₂ : ℤ) 0 (1 / 14)

theorem positivePairDanger_preimage_Ioc
    (w₁ w₂ : ℕ) (hw₁ : 0 < w₁) (hw₂ : 0 < w₂) :
    ((↑) : ℝ → UnitAddCircle) ⁻¹' positivePairDanger w₁ w₂ ∩
        Set.Ioc (0 : ℝ) 1 =
      realRegionCarrier (ratOpenPairRegion w₁ w₂) ∪ ({1} : Set ℝ) := by
  ext x
  constructor
  · rintro ⟨hpair, hx0, hx1⟩
    by_cases hxOne : x = 1
    · exact Or.inr (by simpa [hxOne])
    · apply Or.inl
      have hxIoo : x ∈ Set.Ioo (0 : ℝ) 1 := ⟨hx0, lt_of_le_of_ne hx1 hxOne⟩
      change realStrictMem x (ratOpenPairRegion w₁ w₂)
      rw [realStrictMem_ratOpenPairRegion_iff w₁ w₂ hw₁ hw₂]
      exact ⟨(mem_openCombRegion_iff_mem_danger w₁ hw₁ x hxIoo).2 hpair.1,
        (mem_openCombRegion_iff_mem_danger w₂ hw₂ x hxIoo).2 hpair.2⟩
  · rintro (hregion | hxOne)
    · have hrealPair : x ∈ openCombRegion w₁ ∩ openCombRegion w₂ :=
          (realStrictMem_ratOpenPairRegion_iff w₁ w₂ hw₁ hw₂ x).1 hregion
      have hxIoo := openCombRegion_subset_Ioo w₁ hrealPair.1
      exact ⟨⟨(mem_openCombRegion_iff_mem_danger w₁ hw₁ x hxIoo).1 hrealPair.1,
          (mem_openCombRegion_iff_mem_danger w₂ hw₂ x hxIoo).1 hrealPair.2⟩,
        hxIoo.1, hxIoo.2.le⟩
    · have hx : x = 1 := by simpa using hxOne
      subst x
      constructor
      · change ((1 : ℝ) : UnitAddCircle) ∈ positivePairDanger w₁ w₂
        rw [show ((1 : ℝ) : UnitAddCircle) = 0 by exact AddCircle.coe_period 1]
        constructor <;>
          simp [positivePairDanger, LRCCommensuration.danger,
            LRCCommensuration.runnerMap, Metric.mem_ball, dist_eq_norm]
      · norm_num

/-- Exact continuum mass of a positive-speed pair from the normalized strict
merge list. -/
theorem volume_positivePairDanger_eq_length
    (w₁ w₂ : ℕ) (hw₁ : 0 < w₁) (hw₂ : 0 < w₂) :
    volume (positivePairDanger w₁ w₂) =
      ENNReal.ofReal (((RatIntervals.length (ratOpenPairRegion w₁ w₂) : ℚ) : ℝ)) := by
  have hmeas : MeasurableSet (positivePairDanger w₁ w₂) :=
    (LRCCommensuration.measurableSet_danger _ _ _).inter
      (LRCCommensuration.measurableSet_danger _ _ _)
  have hpull := (UnitAddCircle.measurePreserving_mk 0).measure_preimage
    hmeas.nullMeasurableSet
  rw [MeasureTheory.Measure.restrict_apply' measurableSet_Ioc] at hpull
  norm_num at hpull
  rw [positivePairDanger_preimage_Ioc w₁ w₂ hw₁ hw₂,
    measure_union] at hpull
  · simpa [volume_realRegionCarrier_eq_length _
      (norm_ratOpenPairRegion w₁ w₂ hw₁ hw₂)] using hpull.symm
  · rw [Set.disjoint_left]
    intro x hxRegion hxOne
    have hrealPair :=
      (realStrictMem_ratOpenPairRegion_iff w₁ w₂ hw₁ hw₂ x).1 hxRegion
    have hxIoo := openCombRegion_subset_Ioo w₁ hrealPair.1
    have hxEq : x = 1 := by simpa using hxOne
    exact hxIoo.2.ne hxEq
  · exact measurableSet_singleton 1

#print axioms norm_ratOpenCombRegion
#print axioms strictMem_ratOpenCombRegion_iff_mem_real
#print axioms norm_ratOpenPairRegion
#print axioms length_ratOpenPairRegion_le
#print axioms strictMem_ratOpenPairRegion_iff_mem_real

end
end LRCRationalOpenComb
end LonelyRunner
