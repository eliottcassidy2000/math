import TournamentH7.RatIntervals

/-!
Strict-open two-pointer intersection for normalized rational interval regions.
Unlike `RatIntervals.inter`, this preserves open endpoint semantics, normalization,
and the sharp linear component bound needed by the LRC(14) pair-grid ledger.
-/

namespace LonelyRunner
namespace LRCStrictInterMerge

open RatIntervals

abbrev Interval := ℚ × ℚ

def inside (x : ℚ) (interval : Interval) : Prop :=
  interval.1 < x ∧ x < interval.2

def strictMem (x : ℚ) (region : Region) : Prop :=
  ∃ interval ∈ region, inside x interval

theorem strictMem_nil (x : ℚ) : ¬ strictMem x [] := by
  simp [strictMem]

theorem strictMem_cons_iff {x : ℚ} {interval : Interval} {region : Region} :
    strictMem x (interval :: region) ↔ inside x interval ∨ strictMem x region := by
  simp [strictMem]

theorem inside_clip_iff {x : ℚ} {left right : Interval} :
    inside x (clip left right) ↔ inside x left ∧ inside x right := by
  simp [inside, clip, and_assoc, and_left_comm, and_comm]

def strictInterMerge : Region → Region → Region
  | [], _ => []
  | _, [] => []
  | left :: leftTail, right :: rightTail =>
      let overlap := clip left right
      if left.2 ≤ right.2 then
        if overlap.1 < overlap.2 then
          overlap :: strictInterMerge leftTail (right :: rightTail)
        else
          strictInterMerge leftTail (right :: rightTail)
      else
        if overlap.1 < overlap.2 then
          overlap :: strictInterMerge (left :: leftTail) rightTail
        else
          strictInterMerge (left :: leftTail) rightTail
termination_by left right => left.length + right.length
decreasing_by all_goals simp_wf

theorem head_left_drop {x : ℚ} {left right : Interval}
    {leftTail rightTail : Region}
    (hrightNorm : Norm (right :: rightTail)) (hend : left.2 ≤ right.2) :
    (strictMem x (left :: leftTail) ∧ strictMem x (right :: rightTail)) ↔
      (inside x left ∧ inside x right) ∨
        (strictMem x leftTail ∧ strictMem x (right :: rightTail)) := by
  constructor
  · rintro ⟨hleft, hright⟩
    rw [strictMem_cons_iff] at hleft
    rcases hleft with hleftHead | hleftTail
    · rw [strictMem_cons_iff] at hright
      rcases hright with hrightHead | hrightTail
      · exact Or.inl ⟨hleftHead, hrightHead⟩
      · obtain ⟨other, hother, hinside⟩ := hrightTail
        have horder := norm_head_le hrightNorm other hother
        exfalso
        linarith [hleftHead.2, hinside.1, hend, horder]
    · exact Or.inr ⟨hleftTail, hright⟩
  · rintro (hheads | htails)
    · exact ⟨strictMem_cons_iff.mpr (Or.inl hheads.1),
        strictMem_cons_iff.mpr (Or.inl hheads.2)⟩
    · exact ⟨strictMem_cons_iff.mpr (Or.inr htails.1), htails.2⟩

theorem head_right_drop {x : ℚ} {left right : Interval}
    {leftTail rightTail : Region}
    (hleftNorm : Norm (left :: leftTail)) (hend : right.2 < left.2) :
    (strictMem x (left :: leftTail) ∧ strictMem x (right :: rightTail)) ↔
      (inside x left ∧ inside x right) ∨
        (strictMem x (left :: leftTail) ∧ strictMem x rightTail) := by
  constructor
  · rintro ⟨hleft, hright⟩
    rw [strictMem_cons_iff] at hright
    rcases hright with hrightHead | hrightTail
    · rw [strictMem_cons_iff] at hleft
      rcases hleft with hleftHead | hleftTail
      · exact Or.inl ⟨hleftHead, hrightHead⟩
      · obtain ⟨other, hother, hinside⟩ := hleftTail
        have horder := norm_head_le hleftNorm other hother
        exfalso
        linarith [hrightHead.2, hinside.1, hend, horder]
    · exact Or.inr ⟨hleft, hrightTail⟩
  · rintro (hheads | htails)
    · exact ⟨strictMem_cons_iff.mpr (Or.inl hheads.1),
        strictMem_cons_iff.mpr (Or.inl hheads.2)⟩
    · exact ⟨htails.1, strictMem_cons_iff.mpr (Or.inr htails.2)⟩

theorem strictMem_strictInterMerge_iff (x : ℚ) :
    ∀ (left right : Region), Norm left → Norm right →
      (strictMem x (strictInterMerge left right) ↔
        strictMem x left ∧ strictMem x right)
  | [], right, _hleft, _hright => by simp [strictMem, strictInterMerge]
  | left, [], _hleft, _hright => by
      cases left <;> simp [strictMem, strictInterMerge]
  | left :: leftTail, right :: rightTail, hleftNorm, hrightNorm => by
      rw [strictInterMerge]
      by_cases hend : left.2 ≤ right.2
      · rw [if_pos hend]
        by_cases hlive : (clip left right).1 < (clip left right).2
        · rw [if_pos hlive, strictMem_cons_iff, inside_clip_iff,
            strictMem_strictInterMerge_iff x leftTail (right :: rightTail)
              (norm_tail hleftNorm) hrightNorm,
            head_left_drop hrightNorm hend]
        · rw [if_neg hlive,
            strictMem_strictInterMerge_iff x leftTail (right :: rightTail)
              (norm_tail hleftNorm) hrightNorm,
            head_left_drop hrightNorm hend]
          constructor
          · exact Or.inr
          · rintro (hheads | htails)
            · have hoverlap : inside x (clip left right) := inside_clip_iff.mpr hheads
              exact False.elim (hlive (hoverlap.1.trans hoverlap.2))
            · exact htails
      · rw [if_neg hend]
        have hend' : right.2 < left.2 := lt_of_not_ge hend
        by_cases hlive : (clip left right).1 < (clip left right).2
        · rw [if_pos hlive, strictMem_cons_iff, inside_clip_iff,
            strictMem_strictInterMerge_iff x (left :: leftTail) rightTail
              hleftNorm (norm_tail hrightNorm),
            head_right_drop hleftNorm hend']
        · rw [if_neg hlive,
            strictMem_strictInterMerge_iff x (left :: leftTail) rightTail
              hleftNorm (norm_tail hrightNorm),
            head_right_drop hleftNorm hend']
          constructor
          · exact Or.inr
          · rintro (hheads | htails)
            · have hoverlap : inside x (clip left right) := inside_clip_iff.mpr hheads
              exact False.elim (hlive (hoverlap.1.trans hoverlap.2))
            · exact htails
termination_by left right => left.length + right.length
decreasing_by all_goals simp_wf

theorem mem_strictInterMerge_sources {output : Interval} :
    ∀ (left right : Region), output ∈ strictInterMerge left right →
      ∃ leftSource ∈ left, ∃ rightSource ∈ right,
        output = clip leftSource rightSource
  | [], right, houtput => by simp [strictInterMerge] at houtput
  | left, [], houtput => by
      cases left <;> simp [strictInterMerge] at houtput
  | left :: leftTail, right :: rightTail, houtput => by
      rw [strictInterMerge] at houtput
      by_cases hend : left.2 ≤ right.2
      · rw [if_pos hend] at houtput
        by_cases hlive : (clip left right).1 < (clip left right).2
        · rw [if_pos hlive] at houtput
          rcases List.mem_cons.mp houtput with rfl | houtput
          · exact ⟨left, List.mem_cons_self, right, List.mem_cons_self, rfl⟩
          · obtain ⟨leftSource, hleftSource, rightSource, hrightSource, rfl⟩ :=
              mem_strictInterMerge_sources leftTail (right :: rightTail) houtput
            exact ⟨leftSource, List.mem_cons_of_mem _ hleftSource,
              rightSource, hrightSource, rfl⟩
        · rw [if_neg hlive] at houtput
          obtain ⟨leftSource, hleftSource, rightSource, hrightSource, rfl⟩ :=
            mem_strictInterMerge_sources leftTail (right :: rightTail) houtput
          exact ⟨leftSource, List.mem_cons_of_mem _ hleftSource,
            rightSource, hrightSource, rfl⟩
      · rw [if_neg hend] at houtput
        by_cases hlive : (clip left right).1 < (clip left right).2
        · rw [if_pos hlive] at houtput
          rcases List.mem_cons.mp houtput with rfl | houtput
          · exact ⟨left, List.mem_cons_self, right, List.mem_cons_self, rfl⟩
          · obtain ⟨leftSource, hleftSource, rightSource, hrightSource, rfl⟩ :=
              mem_strictInterMerge_sources (left :: leftTail) rightTail houtput
            exact ⟨leftSource, hleftSource,
              rightSource, List.mem_cons_of_mem _ hrightSource, rfl⟩
        · rw [if_neg hlive] at houtput
          obtain ⟨leftSource, hleftSource, rightSource, hrightSource, rfl⟩ :=
            mem_strictInterMerge_sources (left :: leftTail) rightTail houtput
          exact ⟨leftSource, hleftSource,
            rightSource, List.mem_cons_of_mem _ hrightSource, rfl⟩
termination_by left right => left.length + right.length
decreasing_by all_goals simp_wf

theorem norm_cons_of_live_ordered {interval : Interval} {region : Region}
    (hlive : interval.1 < interval.2) (hregion : Norm region)
    (hordered : ∀ other ∈ region, interval.2 ≤ other.1) :
    Norm (interval :: region) := by
  cases region with
  | nil => exact hlive
  | cons head tail =>
      exact ⟨hlive, hordered head List.mem_cons_self, hregion⟩

theorem norm_strictInterMerge :
    ∀ (left right : Region), Norm left → Norm right →
      Norm (strictInterMerge left right)
  | [], right, _hleft, _hright => by simp [strictInterMerge, Norm]
  | left, [], _hleft, _hright => by
      cases left <;> simp [strictInterMerge, Norm]
  | left :: leftTail, right :: rightTail, hleftNorm, hrightNorm => by
      rw [strictInterMerge]
      by_cases hend : left.2 ≤ right.2
      · rw [if_pos hend]
        have htailNorm := norm_strictInterMerge leftTail (right :: rightTail)
          (norm_tail hleftNorm) hrightNorm
        by_cases hlive : (clip left right).1 < (clip left right).2
        · rw [if_pos hlive]
          apply norm_cons_of_live_ordered hlive htailNorm
          intro output houtput
          obtain ⟨leftSource, hleftSource, rightSource, hrightSource, rfl⟩ :=
            mem_strictInterMerge_sources leftTail (right :: rightTail) houtput
          calc
            (clip left right).2 = left.2 := by simp [clip, min_eq_left hend]
            _ ≤ leftSource.1 := norm_head_le hleftNorm leftSource hleftSource
            _ ≤ (clip leftSource rightSource).1 := by
              exact le_max_left _ _
        · rw [if_neg hlive]
          exact htailNorm
      · rw [if_neg hend]
        have hend' : right.2 < left.2 := lt_of_not_ge hend
        have htailNorm := norm_strictInterMerge (left :: leftTail) rightTail
          hleftNorm (norm_tail hrightNorm)
        by_cases hlive : (clip left right).1 < (clip left right).2
        · rw [if_pos hlive]
          apply norm_cons_of_live_ordered hlive htailNorm
          intro output houtput
          obtain ⟨leftSource, hleftSource, rightSource, hrightSource, rfl⟩ :=
            mem_strictInterMerge_sources (left :: leftTail) rightTail houtput
          calc
            (clip left right).2 = right.2 := by simp [clip, min_eq_right hend'.le]
            _ ≤ rightSource.1 := norm_head_le hrightNorm rightSource hrightSource
            _ ≤ (clip leftSource rightSource).1 := by
              exact le_max_right _ _
        · rw [if_neg hlive]
          exact htailNorm
termination_by left right => left.length + right.length
decreasing_by all_goals simp_wf

theorem length_strictInterMerge_lt_add :
    ∀ (left right : Region), left ≠ [] → right ≠ [] →
      (strictInterMerge left right).length < left.length + right.length
  | [], right, hleft, _hright => by contradiction
  | left, [], _hleft, hright => by contradiction
  | left :: leftTail, right :: rightTail, _hleft, _hright => by
      rw [strictInterMerge]
      by_cases hend : left.2 ≤ right.2
      · rw [if_pos hend]
        by_cases hlive : (clip left right).1 < (clip left right).2
        · rw [if_pos hlive]
          by_cases htail : leftTail = []
          · subst leftTail
            simp [strictInterMerge]
          · have hrec := length_strictInterMerge_lt_add leftTail (right :: rightTail)
                htail (by simp)
            simp only [List.length_cons] at hrec ⊢
            omega
        · rw [if_neg hlive]
          by_cases htail : leftTail = []
          · subst leftTail
            simp [strictInterMerge]
          · have hrec := length_strictInterMerge_lt_add leftTail (right :: rightTail)
                htail (by simp)
            simp only [List.length_cons] at hrec ⊢
            omega
      · rw [if_neg hend]
        by_cases hlive : (clip left right).1 < (clip left right).2
        · rw [if_pos hlive]
          by_cases htail : rightTail = []
          · subst rightTail
            simp [strictInterMerge]
          · have hrec := length_strictInterMerge_lt_add (left :: leftTail) rightTail
                (by simp) htail
            simp only [List.length_cons] at hrec ⊢
            omega
        · rw [if_neg hlive]
          by_cases htail : rightTail = []
          · subst rightTail
            simp [strictInterMerge]
          · have hrec := length_strictInterMerge_lt_add (left :: leftTail) rightTail
                (by simp) htail
            simp only [List.length_cons] at hrec ⊢
            omega
termination_by left right => left.length + right.length
decreasing_by all_goals simp_wf

theorem length_strictInterMerge_le_add_sub_one
    (left right : Region) (hleft : left ≠ []) (hright : right ≠ []) :
    (strictInterMerge left right).length ≤ left.length + right.length - 1 := by
  have hlt := length_strictInterMerge_lt_add left right hleft hright
  omega

#print axioms strictMem_strictInterMerge_iff
#print axioms norm_strictInterMerge
#print axioms length_strictInterMerge_le_add_sub_one

end LRCStrictInterMerge
end LonelyRunner

