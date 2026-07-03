/-
  TournamentH7.LRCRealRegions — the ℝ-interval region calculus
  (kind-pasteur-2026-07-02-S23, HYP-3980, stage 1).

  The ℝ-port of the module-0/S14 machinery (the citation route lives on real
  intervals, so the Hunter ledger cannot run over ℚ-regions).  Ports: `length`,
  `inter`, `cut`/`cutF`, `diff1F`/`diffF`, the partition identity, the
  intersection-shrink lemma, the iterated partition inequality — plus the NEW
  COMMUTATION LEMMA

      length (inter (diff1F X q) B) = length (diff1F (inter X B) q),

  which is TERMWISE a max/min identity on clip/cut pieces (the set identity
  `(X∖q)∩B = (X∩B)∖q` holds piece-by-piece, exactly).  This is the key that
  unlocks the sequential Hunter ledger (stage 2).
-/
import Mathlib
import TournamentH7.LRCFatBlockChain

namespace LonelyRunner
namespace RealRegion

/-- A real region: a list of half-open intervals. -/
abbrev RRegion := List (ℝ × ℝ)

/-- Total length (degenerate pairs contribute 0). -/
noncomputable def rlength (L : RRegion) : ℝ := (L.map fun p => max 0 (p.2 - p.1)).sum

/-- Clip one interval against another. -/
noncomputable def rclip (p q : ℝ × ℝ) : ℝ × ℝ := (max p.1 q.1, min p.2 q.2)

/-- Quadratic intersection. -/
noncomputable def rinter (A B : RRegion) : RRegion :=
  A.flatMap fun p => B.map fun q => rclip p q

theorem rlength_nonneg (L : RRegion) : 0 ≤ rlength L := by
  unfold rlength
  apply List.sum_nonneg
  intro x hx
  rw [List.mem_map] at hx
  obtain ⟨p, _, rfl⟩ := hx
  exact le_max_left _ _

theorem rlength_append (A B : RRegion) : rlength (A ++ B) = rlength A + rlength B := by
  unfold rlength
  rw [List.map_append, List.sum_append]

/-- The (≤ 2) pieces of `p` minus `q`. -/
noncomputable def rcut (p q : ℝ × ℝ) : RRegion :=
  [(p.1, min p.2 q.1), (max p.1 q.2, p.2)]

/-- Region minus one interval. -/
noncomputable def rdiff1 (L : RRegion) (q : ℝ × ℝ) : RRegion := L.flatMap fun p => rcut p q

/-- Region difference: subtract each interval of `B` in turn. -/
noncomputable def rdiff (L B : RRegion) : RRegion := B.foldl rdiff1 L

theorem rlength_cut_add_clip (p : ℝ × ℝ) {q : ℝ × ℝ} (hq : q.1 ≤ q.2) :
    rlength (rcut p q) + max 0 ((rclip p q).2 - (rclip p q).1) = max 0 (p.2 - p.1) := by
  unfold rcut rclip rlength
  simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, add_zero]
  rcases le_total p.2 p.1 with hT | hT
  · -- degenerate p: everything vanishes
    rw [max_eq_left (by linarith : p.2 - p.1 ≤ 0)]
    have hL : min p.2 q.1 - p.1 ≤ 0 := by
      have := min_le_left p.2 q.1
      linarith
    have hR : p.2 - max p.1 q.2 ≤ 0 := by
      have := le_max_left p.1 q.2
      linarith
    have hC : min p.2 q.2 - max p.1 q.1 ≤ 0 := by
      have h1 := min_le_left p.2 q.2
      have h2 := le_max_left p.1 q.1
      linarith
    rw [max_eq_left hL, max_eq_left hR, max_eq_left hC]
    ring
  · rw [max_eq_right (by linarith : (0:ℝ) ≤ p.2 - p.1)]
    rcases le_total q.2 p.1 with h1 | h1
    · -- q entirely left: rcut = p, rclip dead
      have hL : min p.2 q.1 - p.1 ≤ 0 := by
        have := min_le_right p.2 q.1
        linarith
      have hC : min p.2 q.2 - max p.1 q.1 ≤ 0 := by
        have ha := min_le_right p.2 q.2
        have hb := le_max_left p.1 q.1
        linarith
      rw [max_eq_left hL, max_eq_left hC, max_eq_left h1,
        max_eq_right (by linarith : (0:ℝ) ≤ p.2 - p.1)]
      ring
    · rcases le_total p.2 q.1 with h2 | h2
      · -- q entirely right: rcut = p, rclip dead
        have hR : p.2 - max p.1 q.2 ≤ 0 := by
          have := le_max_right p.1 q.2
          linarith
        have hC : min p.2 q.2 - max p.1 q.1 ≤ 0 := by
          have ha := min_le_left p.2 q.2
          have hb := le_max_right p.1 q.1
          linarith
        rw [max_eq_left hR, max_eq_left hC, min_eq_left h2,
          max_eq_right (by linarith : (0:ℝ) ≤ p.2 - p.1)]
        ring
      · -- q meets p: rclip is live; four sub-positions of q's ends
        have hCpos : (0:ℝ) ≤ min p.2 q.2 - max p.1 q.1 := by
          have ha : max p.1 q.1 ≤ min p.2 q.2 := max_le (le_min hT h1) (le_min h2 hq)
          linarith
        rw [max_eq_right hCpos]
        rcases le_total q.1 p.1 with hq1 | hq1
        · rcases le_total p.2 q.2 with hq2 | hq2
          · -- q covers p entirely
            have hL : min p.2 q.1 - p.1 ≤ 0 := by
              have := min_le_right p.2 q.1
              linarith
            have hR : p.2 - max p.1 q.2 ≤ 0 := by
              have := le_max_right p.1 q.2
              linarith
            rw [max_eq_left hL, max_eq_left hR, min_eq_left hq2, max_eq_left hq1]
            ring
          · -- q covers the left end
            have hL : min p.2 q.1 - p.1 ≤ 0 := by
              have := min_le_right p.2 q.1
              linarith
            rw [max_eq_left hL, max_eq_right h1, min_eq_right hq2, max_eq_left hq1,
              max_eq_right (by linarith : (0:ℝ) ≤ p.2 - q.2)]
            ring
        · rcases le_total p.2 q.2 with hq2 | hq2
          · -- q covers the right end
            have hR : p.2 - max p.1 q.2 ≤ 0 := by
              have := le_max_right p.1 q.2
              linarith
            rw [max_eq_left hR, min_eq_right h2, min_eq_left hq2, max_eq_right hq1,
              max_eq_right (by linarith : (0:ℝ) ≤ q.1 - p.1)]
            ring
          · -- q strictly inside p
            rw [min_eq_right h2, max_eq_right h1, min_eq_right hq2, max_eq_right hq1,
              max_eq_right (by linarith : (0:ℝ) ≤ q.1 - p.1),
              max_eq_right (by linarith : (0:ℝ) ≤ p.2 - q.2)]
            ring

/-- RRegion-level, UNCONDITIONAL in `L`: subtracting one interval loses exactly the rclip
overlap. -/
theorem rlength_diff1_add_inter (L : RRegion) {q : ℝ × ℝ} (hq : q.1 ≤ q.2) :
    rlength (rdiff1 L q) + rlength (rinter L [q]) = rlength L := by
  induction L with
  | nil => simp [rdiff1, rinter, rlength]
  | cons p L ih =>
      unfold rdiff1 rinter at ih ⊢
      simp only [List.flatMap_cons]
      rw [rlength_append]
      have hmap : rlength (List.map (fun r => rclip p r) [q]) = max 0 ((rclip p q).2 - (rclip p q).1) := by
        unfold rlength
        simp
      have hL : rlength (List.map (fun r => rclip p r) [q] ++ L.flatMap fun p => List.map (fun r => rclip p r) [q])
          = max 0 ((rclip p q).2 - (rclip p q).1) + rlength (L.flatMap fun p => List.map (fun r => rclip p r) [q]) := by
        rw [rlength_append, hmap]
      have hlenp : rlength (p :: L) = max 0 (p.2 - p.1) + rlength L := by
        unfold rlength
        simp only [List.map_cons, List.sum_cons]
      rw [hL, hlenp, ← rlength_cut_add_clip p hq]
      linarith [ih]


/-- The two rcut pieces are disjoint sub-pieces of their interval: clipping them against
any interval `s` totals at most the whole interval's rclip.  (Live `q` orders the pieces.) -/
theorem clip_cut_pieces_le (p : ℝ × ℝ) {q : ℝ × ℝ} (hq : q.1 ≤ q.2) (s : ℝ × ℝ) :
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


end RealRegion
end LonelyRunner
