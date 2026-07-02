/-
  TournamentH7.LRCPeelAssembly — THE PEEL AGGREGATION ASSEMBLY (kind-pasteur-2026-07-02-S14,
  HYP-3971).  Single-writer satellite.

  The peel leg of the census+peel surface, assembled quantitatively over the S11–S13
  machinery: subtracting a far speed's danger comb from a good region costs the damped
  factor plus `O(1/w)`:

      length (goodRegion (E ++ [w]) h)  ≥  (1 − 2h)·length (goodRegion E h) − err(#pieces)/w.

  Stage 1 (this part): the ITERATED PARTITION INEQUALITY —

      length (diffF L B) ≥ length L − length (inter L B)

  for a live subtrahend, with no other hypotheses.  Ingredients: the cut pieces are
  disjoint sub-pieces of their interval (`clip_cut_pieces_le`), so intersections only
  shrink under `diff1F` (`length_inter_diff1F_le`); the S12 partition identity handles
  each subtraction exactly; the S9 Fubini-swap plumbing (`length_inter_append_right`)
  splits the subtrahend.
-/
import TournamentH7.LRCGoodPipeline
import TournamentH7.LRCFarElementRate

namespace LonelyRunner
namespace RatIntervals

/-! ## Stage 1: the iterated partition inequality -/

/-- Degenerate pairs contribute nothing: filtering them preserves length. -/
theorem length_filter_live (L : Region) :
    length (L.filter fun r => decide (r.1 < r.2)) = length L := by
  induction L with
  | nil => rfl
  | cons p L ih =>
      by_cases hp : p.1 < p.2
      · rw [List.filter_cons_of_pos (by simpa using hp)]
        unfold length at ih ⊢
        simp only [List.map_cons, List.sum_cons]
        rw [ih]
      · rw [List.filter_cons_of_neg (by simpa using hp)]
        unfold length at ih ⊢
        simp only [List.map_cons, List.sum_cons]
        rw [ih, max_eq_left (by push_neg at hp; linarith)]
        ring

/-- The filtered cut has the same length as the unfiltered cut. -/
theorem length_cutF (p q : ℚ × ℚ) : length (cutF p q) = length (cut p q) :=
  length_filter_live _

/-- The filtered partition identity: subtracting one live interval loses exactly the
clip overlap. -/
theorem length_diff1F_add_inter (L : Region) {q : ℚ × ℚ} (hq : q.1 ≤ q.2) :
    length (diff1F L q) + length (inter L [q]) = length L := by
  have h := length_diff1_add_inter L hq
  have hFl : length (diff1F L q) = length (diff1 L q) := by
    unfold diff1F diff1
    induction L with
    | nil => rfl
    | cons p L ih =>
        simp only [List.flatMap_cons]
        rw [length_append, length_append, ih, length_cutF]
  rw [hFl]
  exact h

/-- The two cut pieces are disjoint sub-pieces of their interval: clipping them against
any interval `s` totals at most the whole interval's clip.  (Live `q` orders the pieces.) -/
theorem clip_cut_pieces_le (p : ℚ × ℚ) {q : ℚ × ℚ} (hq : q.1 ≤ q.2) (s : ℚ × ℚ) :
    max 0 (min (min p.2 q.1) s.2 - max p.1 s.1)
      + max 0 (min p.2 s.2 - max (max p.1 q.2) s.1)
      ≤ max 0 (min p.2 s.2 - max p.1 s.1) := by
  have hm : min p.2 q.1 ≤ max p.1 q.2 := by
    have h1 : min p.2 q.1 ≤ q.1 := min_le_right _ _
    have h2 : q.2 ≤ max p.1 q.2 := le_max_right _ _
    linarith
  rcases le_total (min (min p.2 q.1) s.2 - max p.1 s.1) 0 with h1 | h1
  · -- first piece contributes 0
    rw [max_eq_left h1, zero_add]
    -- second ≤ whole: monotone in the left endpoint
    rcases le_total (min p.2 s.2 - max (max p.1 q.2) s.1) 0 with h2 | h2
    · rw [max_eq_left h2]
      exact le_max_left _ _
    · rw [max_eq_right h2]
      have hmono : max p.1 s.1 ≤ max (max p.1 q.2) s.1 :=
        max_le_max (le_max_left _ _) (le_refl _)
      calc min p.2 s.2 - max (max p.1 q.2) s.1
          ≤ min p.2 s.2 - max p.1 s.1 := by linarith
        _ ≤ max 0 (min p.2 s.2 - max p.1 s.1) := le_max_right _ _
  · rw [max_eq_right h1]
    rcases le_total (min p.2 s.2 - max (max p.1 q.2) s.1) 0 with h2 | h2
    · rw [max_eq_left h2, add_zero]
      have hmono1 : min (min p.2 q.1) s.2 ≤ min p.2 s.2 :=
        min_le_min (min_le_left _ _) (le_refl _)
      calc min (min p.2 q.1) s.2 - max p.1 s.1
          ≤ min p.2 s.2 - max p.1 s.1 := by linarith
        _ ≤ max 0 (min p.2 s.2 - max p.1 s.1) := le_max_right _ _
    · rw [max_eq_right h2]
      -- both live: the middle gap m₁ ≤ m₂ absorbs the double count
      have key : min (min p.2 q.1) s.2 - max p.1 s.1
          + (min p.2 s.2 - max (max p.1 q.2) s.1)
          ≤ min p.2 s.2 - max p.1 s.1 := by
        have hA : min (min p.2 q.1) s.2 ≤ min p.2 q.1 := min_le_left _ _
        have hB : max p.1 q.2 ≤ max (max p.1 q.2) s.1 := le_max_left _ _
        linarith
      calc min (min p.2 q.1) s.2 - max p.1 s.1
            + (min p.2 s.2 - max (max p.1 q.2) s.1)
          ≤ min p.2 s.2 - max p.1 s.1 := key
        _ ≤ max 0 (min p.2 s.2 - max p.1 s.1) := le_max_right _ _

/-- Intersections only shrink under a live subtraction. -/
theorem length_inter_diff1F_le (L : Region) {q : ℚ × ℚ} (hq : q.1 ≤ q.2) (B : Region) :
    length (inter (diff1F L q) B) ≤ length (inter L B) := by
  induction L with
  | nil =>
      unfold diff1F inter
      simp [length]
  | cons p L ih =>
      -- diff1F (p::L) q = cutF p q ++ diff1F L q; inter distributes over ++ on the left
      have hsplit : diff1F (p :: L) q = cutF p q ++ diff1F L q := by
        unfold diff1F
        simp [List.flatMap_cons]
      have hinterL : inter (cutF p q ++ diff1F L q) B
          = inter (cutF p q) B ++ inter (diff1F L q) B := by
        unfold inter
        simp [List.flatMap_append]
      have hinterp : inter (p :: L) B = inter [p] B ++ inter L B := by
        unfold inter
        simp [List.flatMap_cons]
      rw [hsplit, hinterL, hinterp, length_append, length_append]
      have hhead : length (inter (cutF p q) B) ≤ length (inter [p] B) := by
        -- per B-interval: the two pieces clip to at most the whole
        unfold cutF cut inter
        by_cases hlive1 : p.1 < min p.2 q.1
        · by_cases hlive2 : max p.1 q.2 < p.2
          · -- both pieces live after the filter
            rw [List.filter_cons_of_pos (by simpa using hlive1),
              List.filter_cons_of_pos (by simpa using hlive2)]
            simp only [List.filter_nil, List.flatMap_cons, List.flatMap_nil,
              List.append_nil]
            rw [length_append]
            -- lengths of mapped clips: sum the per-s inequality
            unfold length
            simp only [List.map_map]
            induction B with
            | nil => simp
            | cons s B ihB =>
                simp only [List.map_cons, List.sum_cons]
                have hper := clip_cut_pieces_le p hq s
                unfold clip at hper
                simp only [Function.comp_apply]
                have := ihB
                linarith [hper, this]
          · rw [List.filter_cons_of_pos (by simpa using hlive1),
              List.filter_cons_of_neg (by simpa using hlive2)]
            simp only [List.filter_nil, List.flatMap_cons, List.flatMap_nil,
              List.append_nil]
            -- one live piece: bound by monotone endpoints
            unfold length
            simp only [List.map_map]
            induction B with
            | nil => simp
            | cons s B ihB =>
                simp only [List.map_cons, List.sum_cons]
                have hper := clip_cut_pieces_le p hq s
                unfold clip at hper
                simp only [Function.comp_apply]
                have hsecond : (0:ℚ) ≤ max 0 (min p.2 s.2 - max (max p.1 q.2) s.1) :=
                  le_max_left _ _
                linarith [hper, ihB]
        · by_cases hlive2 : max p.1 q.2 < p.2
          · rw [List.filter_cons_of_neg (by simpa using hlive1),
              List.filter_cons_of_pos (by simpa using hlive2)]
            simp only [List.filter_nil, List.flatMap_cons, List.flatMap_nil,
              List.append_nil]
            unfold length
            simp only [List.map_map]
            induction B with
            | nil => simp
            | cons s B ihB =>
                simp only [List.map_cons, List.sum_cons]
                have hper := clip_cut_pieces_le p hq s
                unfold clip at hper
                simp only [Function.comp_apply]
                have hfirst : (0:ℚ) ≤ max 0 (min (min p.2 q.1) s.2 - max p.1 s.1) :=
                  le_max_left _ _
                linarith [hper, ihB]
          · rw [List.filter_cons_of_neg (by simpa using hlive1),
              List.filter_cons_of_neg (by simpa using hlive2)]
            simp only [List.filter_nil, List.flatMap_nil]
            unfold length
            simp only [List.map_nil, List.sum_nil]
            exact length_nonneg _ |>.trans (le_of_eq (by unfold length; rfl)) |>.trans
              (le_refl _) |>.trans (by
                have : (0:ℚ) ≤ length (List.map (fun q => clip p q) B) := length_nonneg _
                unfold length at this ⊢
                simpa using this)
      linarith [ih]

/-- **THE ITERATED PARTITION INEQUALITY**: subtracting a live region loses at most the
intersection — with NO hypotheses beyond liveness of the subtrahend. -/
theorem length_diffF_ge (L : Region) :
    ∀ B : Region, (∀ q ∈ B, q.1 ≤ q.2) →
      length L - length (inter L B) ≤ length (diffF L B) := by
  intro B
  induction B generalizing L with
  | nil =>
      intro _
      rw [length_inter_nil]
      unfold diffF
      simp
  | cons q B ih =>
      intro hlive
      have hq : q.1 ≤ q.2 := hlive q (List.mem_cons_self ..)
      have hB : ∀ s ∈ B, s.1 ≤ s.2 := fun s hs => hlive s (List.mem_cons_of_mem _ hs)
      have hfold : diffF L (q :: B) = diffF (diff1F L q) B := by
        unfold diffF
        simp [List.foldl_cons]
      rw [hfold]
      have hIH := ih (diff1F L q) hB
      have hpart := length_diff1F_add_inter L hq
      have hshrink := length_inter_diff1F_le L hq B
      have hsplitB : length (inter L (q :: B)) = length (inter L [q]) + length (inter L B) := by
        have : (q :: B : Region) = [q] ++ B := rfl
        rw [this, length_inter_append_right]
      linarith

end RatIntervals
end LonelyRunner
