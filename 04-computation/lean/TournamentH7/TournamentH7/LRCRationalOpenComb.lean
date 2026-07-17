import TournamentH7.LRCOpenDangerComb
import TournamentH7.LRCStrictInterMerge

/-!
Normalized rational strict-open danger combs and their exact two-pointer pair
intersection.  The rational endpoints expose exact grid arithmetic while the
cast lemmas connect the lists to the real comb carrier.
-/

namespace LonelyRunner
namespace LRCRationalOpenComb

open RatIntervals LRCStrictInterMerge LRCOpenDangerComb

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

#print axioms norm_ratOpenCombRegion
#print axioms strictMem_ratOpenCombRegion_iff_mem_real
#print axioms norm_ratOpenPairRegion
#print axioms length_ratOpenPairRegion_le
#print axioms strictMem_ratOpenPairRegion_iff_mem_real

end
end LRCRationalOpenComb
end LonelyRunner
