import TournamentH7.LRCStrictInterMerge
import TournamentH7.LRCSevenTranslate

/-!
Length comparison between the linear strict two-pointer merge and the
quadratic all-clips rational intersection.  On normalized inputs they have
the same total length, reducing the remaining pair-covariance producer to an
explicit double clip-sum reindexing.
-/

namespace LonelyRunner
namespace LRCB5MergeLength

open RatIntervals LRCStrictInterMerge

abbrev Interval := ℚ × ℚ

def clipLength (left right : Interval) : ℚ :=
  max 0 ((clip left right).2 - (clip left right).1)

theorem clipLength_eq_zero_of_end_le_start
    {left right : Interval} (h : left.2 ≤ right.1) :
    clipLength left right = 0 := by
  unfold clipLength clip
  apply max_eq_left
  have hmin : min left.2 right.2 ≤ left.2 := min_le_left _ _
  have hmax : right.1 ≤ max left.1 right.1 := le_max_right _ _
  linarith

theorem length_map_clip_eq_zero_of_end_le
    (left : Interval) :
    ∀ (right : Region),
      (∀ interval ∈ right, left.2 ≤ interval.1) →
      RatIntervals.length (right.map fun interval => clip left interval) = 0
  | [], _h => by simp [RatIntervals.length]
  | head :: tail, h => by
      have hhead := h head List.mem_cons_self
      have htail : ∀ interval ∈ tail, left.2 ≤ interval.1 :=
        fun interval hinterval => h interval (List.mem_cons_of_mem head hinterval)
      have inductionHypothesis :=
        length_map_clip_eq_zero_of_end_le left tail htail
      simp only [List.map_cons]
      change clipLength left head +
          RatIntervals.length (tail.map fun interval => clip left interval) = 0
      rw [clipLength_eq_zero_of_end_le_start hhead, inductionHypothesis, zero_add]

theorem length_inter_drop_left
    {left right : Interval} {leftTail rightTail : Region}
    (hrightNorm : Norm (right :: rightTail))
    (hend : left.2 ≤ right.2) :
    RatIntervals.length (inter (left :: leftTail) (right :: rightTail)) =
      clipLength left right +
        RatIntervals.length (inter leftTail (right :: rightTail)) := by
  have hfuture : ∀ interval ∈ rightTail, left.2 ≤ interval.1 := by
    intro interval hinterval
    exact hend.trans (norm_head_le hrightNorm interval hinterval)
  have hzero := length_map_clip_eq_zero_of_end_le left rightTail hfuture
  unfold inter
  simp only [List.flatMap_cons, List.map_cons]
  rw [RatIntervals.length_append]
  change (clipLength left right +
      RatIntervals.length (rightTail.map fun interval => clip left interval)) +
        RatIntervals.length (inter leftTail (right :: rightTail)) =
      clipLength left right +
        RatIntervals.length (inter leftTail (right :: rightTail))
  rw [hzero, add_zero]

theorem length_inter_drop_right
    {left right : Interval} {leftTail rightTail : Region}
    (hleftNorm : Norm (left :: leftTail))
    (hend : right.2 < left.2) :
    RatIntervals.length (inter (left :: leftTail) (right :: rightTail)) =
      clipLength left right +
        RatIntervals.length (inter (left :: leftTail) rightTail) := by
  have hdrop := length_inter_drop_left (left := right) (right := left)
    (leftTail := rightTail) (rightTail := leftTail) hleftNorm hend.le
  calc
    RatIntervals.length (inter (left :: leftTail) (right :: rightTail)) =
        RatIntervals.length (inter (right :: rightTail) (left :: leftTail)) :=
      RatIntervals.length_inter_comm _ _
    _ = clipLength right left +
        RatIntervals.length (inter rightTail (left :: leftTail)) := hdrop
    _ = clipLength left right +
        RatIntervals.length (inter (left :: leftTail) rightTail) := by
      rw [RatIntervals.length_inter_comm]
      congr 1
      unfold clipLength
      rw [RatIntervals.clip_comm]

/-- The linear two-pointer merge and the quadratic all-clips intersection have
exactly the same total length on normalized inputs. -/
theorem length_strictInterMerge_eq_inter :
    ∀ (left right : Region), Norm left → Norm right →
      RatIntervals.length (strictInterMerge left right) =
        RatIntervals.length (inter left right)
  | [], right, _hleft, _hright => by simp [strictInterMerge, inter, RatIntervals.length]
  | left, [], _hleft, _hright => by
      cases left with
      | nil => simp [strictInterMerge, inter, RatIntervals.length]
      | cons head tail =>
          rw [strictInterMerge, RatIntervals.length_inter_nil]
          · rfl
          · exact List.cons_ne_nil head tail
  | left :: leftTail, right :: rightTail, hleftNorm, hrightNorm => by
      rw [strictInterMerge]
      by_cases hend : left.2 ≤ right.2
      · rw [if_pos hend, length_inter_drop_left hrightNorm hend]
        by_cases hlive : (clip left right).1 < (clip left right).2
        · rw [if_pos hlive]
          have inductionHypothesis :=
            length_strictInterMerge_eq_inter leftTail (right :: rightTail)
              (norm_tail hleftNorm) hrightNorm
          change clipLength left right +
              RatIntervals.length (strictInterMerge leftTail (right :: rightTail)) =
            clipLength left right +
              RatIntervals.length (inter leftTail (right :: rightTail))
          rw [inductionHypothesis]
        · rw [if_neg hlive]
          have hclipZero : clipLength left right = 0 := by
            unfold clipLength
            exact max_eq_left (sub_nonpos.mpr (le_of_not_gt hlive))
          rw [hclipZero, zero_add]
          exact length_strictInterMerge_eq_inter leftTail (right :: rightTail)
            (norm_tail hleftNorm) hrightNorm
      · rw [if_neg hend]
        have hend' : right.2 < left.2 := lt_of_not_ge hend
        rw [length_inter_drop_right hleftNorm hend']
        by_cases hlive : (clip left right).1 < (clip left right).2
        · rw [if_pos hlive]
          have inductionHypothesis :=
            length_strictInterMerge_eq_inter (left :: leftTail) rightTail
              hleftNorm (norm_tail hrightNorm)
          change clipLength left right +
              RatIntervals.length (strictInterMerge (left :: leftTail) rightTail) =
            clipLength left right +
              RatIntervals.length (inter (left :: leftTail) rightTail)
          rw [inductionHypothesis]
        · rw [if_neg hlive]
          have hclipZero : clipLength left right = 0 := by
            unfold clipLength
            exact max_eq_left (sub_nonpos.mpr (le_of_not_gt hlive))
          rw [hclipZero, zero_add]
          exact length_strictInterMerge_eq_inter (left :: leftTail) rightTail
            hleftNorm (norm_tail hrightNorm)
termination_by left right => left.length + right.length
decreasing_by all_goals simp_wf

#print axioms length_strictInterMerge_eq_inter

end LRCB5MergeLength
end LonelyRunner
