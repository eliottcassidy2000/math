/-
  TournamentH7.LRCRegionDiff — REGION DIFFERENCE (kind-pasteur-2026-07-02-S12, HYP-3969b).
  Single-writer satellite.

  The last missing module-0 operation: set difference for rational interval regions.
  Both remaining legs of the unconditional-LRC(14) surface need it:
   * GOOD-SET MATERIALIZATION — `good = diff fullWindow (danger combs)`: the bounded
     census leg's per-class base positivity becomes `0 < length (diff … …)`, a computable
     rational;
   * PEEL AGGREGATION — `Good_A(E ∪ {w}) = diff (Good_A E) (hit comb)`, with the length
     identity feeding the rate lemma's `(1 − ρ)`-damping bookkeeping.

  Design: `cut p q` = the (≤ 2) pieces of one interval minus one interval — kept possibly
  degenerate (degenerate pairs have no members and length 0, so no filtering/normalization
  is ever needed).  `diff1` maps it over a region; `diff` folds over the subtrahend.

  THE PARTITION IDENTITY (`length_cut_add_clip`, per-interval, UNCONDITIONAL):
      length (cut p q) + cliplen p q = length [p]
  which sums to `length_diff1 : length (diff1 L q) + length (inter L [q]) = length L`
  with NO hypotheses on `L` — the subtracted interval never double-counts because `cut`
  is exact at the endpoints (half-open discipline).
-/
import TournamentH7.LRCSevenTranslate

namespace LonelyRunner
namespace RatIntervals

/-! ## One interval minus one interval -/

/-- The (≤ 2) pieces of `p` minus `q`: the left piece `[p₁, min p₂ q₁)` and the right
piece `[max p₁ q₂, p₂)`.  Degenerate pieces are kept (memberless, length 0). -/
def cut (p q : ℚ × ℚ) : Region :=
  [(p.1, min p.2 q.1), (max p.1 q.2, p.2)]

/-- Region minus one interval. -/
def diff1 (L : Region) (q : ℚ × ℚ) : Region := L.flatMap fun p => cut p q

/-- Region difference: subtract each interval of `B` in turn. -/
def diff (L B : Region) : Region := B.foldl diff1 L

/-! ## Membership -/

theorem mem_cut {x : ℚ} {p q : ℚ × ℚ} :
    mem x (cut p q) ↔ (p.1 ≤ x ∧ x < p.2) ∧ ¬ (q.1 ≤ x ∧ x < q.2) := by
  unfold mem cut
  constructor
  · rintro ⟨s, hs, hx1, hx2⟩
    rcases List.mem_cons.mp hs with rfl | hs'
    · simp only at hx1 hx2
      have h1 : x < q.1 := lt_of_lt_of_le hx2 (min_le_right _ _)
      have h2 : x < p.2 := lt_of_lt_of_le hx2 (min_le_left _ _)
      exact ⟨⟨hx1, h2⟩, fun hq => absurd h1 (not_lt.mpr hq.1)⟩
    · rcases List.mem_singleton.mp hs' with rfl
      simp only at hx1 hx2
      have h1 : p.1 ≤ x := le_trans (le_max_left _ _) hx1
      have h2 : q.2 ≤ x := le_trans (le_max_right _ _) hx1
      exact ⟨⟨h1, hx2⟩, fun hq => absurd h2 (not_le.mpr hq.2)⟩
  · rintro ⟨⟨hp1, hp2⟩, hq⟩
    rcases le_or_gt q.1 x with hq1 | hq1
    · -- then ¬(x < q.2), i.e. q.2 ≤ x: the right piece
      have hq2 : q.2 ≤ x := by
        by_contra hcon
        exact hq ⟨hq1, lt_of_not_ge hcon⟩
      refine ⟨(max p.1 q.2, p.2), List.mem_cons_of_mem _ (List.mem_singleton_self _), ?_, ?_⟩
      · simp only
        exact max_le hp1 hq2
      · simpa using hp2
    · -- x < q.1: the left piece
      refine ⟨(p.1, min p.2 q.1), List.mem_cons_self .., ?_, ?_⟩
      · simpa using hp1
      · simp only
        exact lt_min hp2 hq1

theorem mem_diff1 {x : ℚ} {L : Region} {q : ℚ × ℚ} :
    mem x (diff1 L q) ↔ mem x L ∧ ¬ (q.1 ≤ x ∧ x < q.2) := by
  unfold diff1
  constructor
  · rintro ⟨s, hs, hx⟩
    simp only [List.mem_flatMap] at hs
    obtain ⟨p, hp, hsp⟩ := hs
    have := mem_cut.mp ⟨s, hsp, hx⟩
    exact ⟨⟨p, hp, this.1⟩, this.2⟩
  · rintro ⟨⟨p, hp, hxp⟩, hq⟩
    obtain ⟨s, hsp, hx⟩ := mem_cut.mpr ⟨hxp, hq⟩
    exact ⟨s, List.mem_flatMap.mpr ⟨p, hp, hsp⟩, hx⟩

/-- Membership through the full difference: in `L` and in no interval of `B`. -/
theorem mem_diff {x : ℚ} : ∀ {B : Region} {L : Region},
    mem x (diff L B) ↔ mem x L ∧ ∀ q ∈ B, ¬ (q.1 ≤ x ∧ x < q.2) := by
  intro B
  induction B with
  | nil =>
      intro L
      unfold diff
      simp
  | cons q B ih =>
      intro L
      unfold diff at ih ⊢
      simp only [List.foldl_cons]
      rw [ih, mem_diff1]
      constructor
      · rintro ⟨⟨hL, hq⟩, hB⟩
        refine ⟨hL, ?_⟩
        intro s hs
        rcases List.mem_cons.mp hs with rfl | hs'
        · exact hq
        · exact hB s hs'
      · rintro ⟨hL, hall⟩
        exact ⟨⟨hL, hall q (List.mem_cons_self ..)⟩,
          fun s hs => hall s (List.mem_cons_of_mem _ hs)⟩

/-- The difference avoids the subtrahend: `mem x (diff L B) → ¬ mem x B`. -/
theorem not_mem_of_mem_diff {x : ℚ} {L B : Region} (h : mem x (diff L B)) :
    ¬ mem x B := by
  rintro ⟨q, hq, hx⟩
  exact ((mem_diff.mp h).2 q hq) hx

/-! ## The partition identity -/

/-- Per-interval: for a LIVE subtrahend (`q.1 ≤ q.2`), the cut pieces and the clip
overlap partition the interval's length exactly.  (For degenerate `q` the two cut
pieces would overlap — liveness is essential.) -/
theorem length_cut_add_clip (p : ℚ × ℚ) {q : ℚ × ℚ} (hq : q.1 ≤ q.2) :
    length (cut p q) + max 0 ((clip p q).2 - (clip p q).1) = max 0 (p.2 - p.1) := by
  unfold cut clip length
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
  · rw [max_eq_right (by linarith : (0:ℚ) ≤ p.2 - p.1)]
    rcases le_total q.2 p.1 with h1 | h1
    · -- q entirely left: cut = p, clip dead
      have hL : min p.2 q.1 - p.1 ≤ 0 := by
        have := min_le_right p.2 q.1
        linarith
      have hC : min p.2 q.2 - max p.1 q.1 ≤ 0 := by
        have ha := min_le_right p.2 q.2
        have hb := le_max_left p.1 q.1
        linarith
      rw [max_eq_left hL, max_eq_left hC, max_eq_left h1,
        max_eq_right (by linarith : (0:ℚ) ≤ p.2 - p.1)]
      ring
    · rcases le_total p.2 q.1 with h2 | h2
      · -- q entirely right: cut = p, clip dead
        have hR : p.2 - max p.1 q.2 ≤ 0 := by
          have := le_max_right p.1 q.2
          linarith
        have hC : min p.2 q.2 - max p.1 q.1 ≤ 0 := by
          have ha := min_le_left p.2 q.2
          have hb := le_max_right p.1 q.1
          linarith
        rw [max_eq_left hR, max_eq_left hC, min_eq_left h2,
          max_eq_right (by linarith : (0:ℚ) ≤ p.2 - p.1)]
        ring
      · -- q meets p: clip is live; four sub-positions of q's ends
        have hCpos : (0:ℚ) ≤ min p.2 q.2 - max p.1 q.1 := by
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
              max_eq_right (by linarith : (0:ℚ) ≤ p.2 - q.2)]
            ring
        · rcases le_total p.2 q.2 with hq2 | hq2
          · -- q covers the right end
            have hR : p.2 - max p.1 q.2 ≤ 0 := by
              have := le_max_right p.1 q.2
              linarith
            rw [max_eq_left hR, min_eq_right h2, min_eq_left hq2, max_eq_right hq1,
              max_eq_right (by linarith : (0:ℚ) ≤ q.1 - p.1)]
            ring
          · -- q strictly inside p
            rw [min_eq_right h2, max_eq_right h1, min_eq_right hq2, max_eq_right hq1,
              max_eq_right (by linarith : (0:ℚ) ≤ q.1 - p.1),
              max_eq_right (by linarith : (0:ℚ) ≤ p.2 - q.2)]
            ring

/-- Region-level, UNCONDITIONAL in `L`: subtracting one interval loses exactly the clip
overlap. -/
theorem length_diff1_add_inter (L : Region) {q : ℚ × ℚ} (hq : q.1 ≤ q.2) :
    length (diff1 L q) + length (inter L [q]) = length L := by
  induction L with
  | nil => simp [diff1, inter, length]
  | cons p L ih =>
      unfold diff1 inter at ih ⊢
      simp only [List.flatMap_cons]
      rw [length_append]
      have hmap : length (List.map (fun r => clip p r) [q]) = max 0 ((clip p q).2 - (clip p q).1) := by
        unfold length
        simp
      have hL : length (List.map (fun r => clip p r) [q] ++ L.flatMap fun p => List.map (fun r => clip p r) [q])
          = max 0 ((clip p q).2 - (clip p q).1) + length (L.flatMap fun p => List.map (fun r => clip p r) [q]) := by
        rw [length_append, hmap]
      have hlenp : length (p :: L) = max 0 (p.2 - p.1) + length L := by
        unfold length
        simp only [List.map_cons, List.sum_cons]
      rw [hL, hlenp, ← length_cut_add_clip p hq]
      linarith [ih]

/-! ## Computational variants: degenerate pieces dropped

The plain `diff` keeps degenerate pieces, which is fine for proofs but doubles the list
per subtraction — exponential under evaluation.  The filtered variants keep the piece
count linear (a live piece per retained endpoint); membership is unchanged because
degenerate pieces are memberless. -/

/-- Computational cut: only live pieces. -/
def cutF (p q : ℚ × ℚ) : Region := (cut p q).filter fun r => decide (r.1 < r.2)

/-- Computational region-minus-interval. -/
def diff1F (L : Region) (q : ℚ × ℚ) : Region := L.flatMap fun p => cutF p q

/-- Computational region difference. -/
def diffF (L B : Region) : Region := B.foldl diff1F L

theorem mem_cutF {x : ℚ} {p q : ℚ × ℚ} :
    mem x (cutF p q) ↔ (p.1 ≤ x ∧ x < p.2) ∧ ¬ (q.1 ≤ x ∧ x < q.2) := by
  rw [← mem_cut]
  unfold cutF mem
  constructor
  · rintro ⟨s, hs, hx⟩
    exact ⟨s, (List.mem_filter.mp hs).1, hx⟩
  · rintro ⟨s, hs, hx1, hx2⟩
    refine ⟨s, List.mem_filter.mpr ⟨hs, ?_⟩, hx1, hx2⟩
    exact decide_eq_true_eq.mpr (lt_of_le_of_lt hx1 hx2)

theorem mem_diff1F {x : ℚ} {L : Region} {q : ℚ × ℚ} :
    mem x (diff1F L q) ↔ mem x L ∧ ¬ (q.1 ≤ x ∧ x < q.2) := by
  unfold diff1F
  constructor
  · rintro ⟨s, hs, hx⟩
    simp only [List.mem_flatMap] at hs
    obtain ⟨p, hp, hsp⟩ := hs
    have := mem_cutF.mp ⟨s, hsp, hx⟩
    exact ⟨⟨p, hp, this.1⟩, this.2⟩
  · rintro ⟨⟨p, hp, hxp⟩, hq⟩
    obtain ⟨s, hsp, hx⟩ := mem_cutF.mpr ⟨hxp, hq⟩
    exact ⟨s, List.mem_flatMap.mpr ⟨p, hp, hsp⟩, hx⟩

/-- Membership through the computational difference — same characterization as `mem_diff`. -/
theorem mem_diffF {x : ℚ} : ∀ {B : Region} {L : Region},
    mem x (diffF L B) ↔ mem x L ∧ ∀ q ∈ B, ¬ (q.1 ≤ x ∧ x < q.2) := by
  intro B
  induction B with
  | nil =>
      intro L
      unfold diffF
      simp
  | cons q B ih =>
      intro L
      unfold diffF at ih ⊢
      simp only [List.foldl_cons]
      rw [ih, mem_diff1F]
      constructor
      · rintro ⟨⟨hL, hq⟩, hB⟩
        refine ⟨hL, ?_⟩
        intro s hs
        rcases List.mem_cons.mp hs with rfl | hs'
        · exact hq
        · exact hB s hs'
      · rintro ⟨hL, hall⟩
        exact ⟨⟨hL, hall q (List.mem_cons_self ..)⟩,
          fun s hs => hall s (List.mem_cons_of_mem _ hs)⟩

/-! ## The witness extractor -/

/-- Positive length yields an explicit member: the left endpoint of any live interval.
With `mem_diff`, a computed `0 < length (diff window dangerCombs)` produces a rational
point avoiding every danger interval — the census leg's witness extractor. -/
theorem exists_mem_of_length_pos : ∀ {L : Region}, 0 < length L → ∃ x : ℚ, mem x L := by
  intro L
  induction L with
  | nil => intro h; simp [length] at h
  | cons p L ih =>
      intro h
      rcases le_or_gt p.2 p.1 with hdeg | hlive
      · -- degenerate head contributes 0: recurse
        have hlen : length (p :: L) = length L := by
          unfold length
          simp only [List.map_cons, List.sum_cons]
          rw [max_eq_left (by linarith : p.2 - p.1 ≤ 0)]
          ring
        rw [hlen] at h
        obtain ⟨x, q, hq, hx⟩ := ih h
        exact ⟨x, q, List.mem_cons_of_mem _ hq, hx⟩
      · exact ⟨p.1, p, List.mem_cons_self .., le_refl _, hlive⟩

end RatIntervals
end LonelyRunner
