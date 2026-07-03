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


theorem rlength_inter_nil (B : RRegion) : rlength (rinter B []) = 0 := by
  induction B with
  | nil => rfl
  | cons p B ih =>
      unfold rinter at ih ⊢
      simp only [List.flatMap_cons, List.map_nil, List.nil_append]
      exact ih

/-! ## Stage 1b: the commutation (swap) machinery -/

/-- **THE LIST-LEVEL SWAP for a single window interval**: intersecting with one
interval after a one-interval subtraction IS subtracting from the intersection —
the lists are literally equal (per piece, `min`/`max` right-commutativity). -/
theorem rinter_rdiff1_single (X : RRegion) (q d : ℝ × ℝ) :
    rinter (rdiff1 X q) [d] = rdiff1 (rinter X [d]) q := by
  induction X with
  | nil => rfl
  | cons x X ih =>
      have hL : rinter (rdiff1 (x :: X) q) [d]
          = [rclip (x.1, min x.2 q.1) d, rclip (max x.1 q.2, x.2) d]
            ++ rinter (rdiff1 X q) [d] := by
        unfold rdiff1 rinter rcut
        simp [List.flatMap_cons]
      have hR : rdiff1 (rinter (x :: X) [d]) q
          = rcut (rclip x d) q ++ rdiff1 (rinter X [d]) q := by
        unfold rdiff1 rinter
        simp [List.flatMap_cons]
      rw [hL, hR, ih]
      congr 1
      unfold rclip rcut
      have h1 : min (min x.2 q.1) d.2 = min (min x.2 d.2) q.1 := min_right_comm _ _ _
      have h2 : max (max x.1 q.2) d.1 = max (max x.1 d.1) q.2 := max_right_comm _ _ _
      rw [h1, h2]

/-- Single-window swap through a whole region subtraction (lists). -/
theorem rinter_rdiff_single (D : RRegion) (G : RRegion) (d : ℝ × ℝ) :
    rinter (rdiff G D) [d] = rdiff (rinter G [d]) D := by
  induction D generalizing G with
  | nil => rfl
  | cons q D ih =>
      have hfoldL : rdiff G (q :: D) = rdiff (rdiff1 G q) D := by
        unfold rdiff
        simp [List.foldl_cons]
      have hfoldR : rdiff (rinter G [d]) (q :: D) = rdiff (rdiff1 (rinter G [d]) q) D := by
        unfold rdiff
        simp [List.foldl_cons]
      rw [hfoldL, hfoldR, ih (rdiff1 G q), rinter_rdiff1_single]

/-- Intersection length is additive in the window list. -/
theorem rlength_inter_append_right (A B C : RRegion) :
    rlength (rinter A (B ++ C)) = rlength (rinter A B) + rlength (rinter A C) := by
  induction A with
  | nil => simp [rinter, rlength]
  | cons p A ih =>
      unfold rinter at ih ⊢
      simp only [List.flatMap_cons, List.map_append, List.append_assoc, rlength_append]
      simp only [List.map_append] at ih
      rw [ih]
      ring

/-- Cross-vanishing: clips inside `d` never meet a disjoint `q`. -/
theorem cross_vanish {d q : ℝ × ℝ} (hdis : d.2 ≤ q.1 ∨ q.2 ≤ d.1) (L : RRegion) :
    rlength (rinter (rinter L [d]) [q]) = 0 := by
  induction L with
  | nil => rfl
  | cons x L ih =>
      have hsplit : rinter (rinter (x :: L) [d]) [q]
          = [rclip (rclip x d) q] ++ rinter (rinter L [d]) [q] := by
        unfold rinter
        simp [List.flatMap_cons]
      rw [hsplit, rlength_append, ih]
      have hzero : rlength [rclip (rclip x d) q] = 0 := by
        unfold rclip rlength
        simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, add_zero]
        apply max_eq_left
        rcases hdis with h | h
        · have h1 : min (min x.2 d.2) q.2 ≤ d.2 := le_trans (min_le_left _ _) (min_le_right _ _)
          have h2 : q.1 ≤ max (max x.1 d.1) q.1 := le_max_right _ _
          linarith
        · have h1 : min (min x.2 d.2) q.2 ≤ q.2 := min_le_right _ _
          have h2 : d.1 ≤ max (max x.1 d.1) q.1 := le_trans (le_max_right _ _) (le_max_left _ _)
          linarith
      rw [hzero]
      ring

/-- Subtracting an interval `q` lying strictly before every interval of `D` does not
change the `D`-intersection mass (per-window swap + cross-vanish). -/
theorem rlength_inter_rdiff1_disjoint (L : RRegion) {q : ℝ × ℝ} (hq : q.1 ≤ q.2) :
    ∀ (D : RRegion), (∀ d ∈ D, d.2 ≤ q.1 ∨ q.2 ≤ d.1) →
    rlength (rinter (rdiff1 L q) D) = rlength (rinter L D) := by
  intro D
  induction D with
  | nil =>
      intro _
      rw [rlength_inter_nil, rlength_inter_nil]
  | cons d D ih =>
      intro hsep
      have hd := hsep d (List.mem_cons_self ..)
      have htail : ∀ r ∈ D, r.2 ≤ q.1 ∨ q.2 ≤ r.1 :=
        fun r hr => hsep r (List.mem_cons_of_mem _ hr)
      have hcons : (d :: D : RRegion) = [d] ++ D := rfl
      rw [hcons, rlength_inter_append_right, rlength_inter_append_right, ih htail]
      have hhead : rlength (rinter (rdiff1 L q) [d]) = rlength (rinter L [d]) := by
        rw [rinter_rdiff1_single]
        have hpart := rlength_diff1_add_inter (rinter L [d]) hq
        have hcross := cross_vanish hd L
        linarith
      rw [hhead]

/-- Pairwise separation of a region's intervals (sorted teeth: each interval ends
before the next begins). -/
def SortedSep : RRegion → Prop
  | [] => True
  | q :: D => (∀ d ∈ D, q.2 ≤ d.1) ∧ SortedSep D

/-- **THE EXACT REGION PARTITION** (separated subtrahend): subtracting a region of
sorted separated live intervals removes EXACTLY the intersection mass — the
sequential ledger has no loss terms. -/
theorem rlength_rdiff_partition :
    ∀ (D : RRegion) (L : RRegion), (∀ r ∈ D, r.1 ≤ r.2) → SortedSep D →
    rlength (rdiff L D) = rlength L - rlength (rinter L D) := by
  intro D
  induction D with
  | nil =>
      intro L _ _
      rw [rlength_inter_nil]
      simp [rdiff]
  | cons q D ih =>
      intro L hlive hsep
      obtain ⟨hqd, hsep'⟩ := hsep
      have hq : q.1 ≤ q.2 := hlive q (List.mem_cons_self ..)
      have hfold : rdiff L (q :: D) = rdiff (rdiff1 L q) D := by
        unfold rdiff
        rw [List.foldl_cons]
      rw [hfold, ih (rdiff1 L q) (fun r hr => hlive r (List.mem_cons_of_mem _ hr)) hsep']
      have hswap : rlength (rinter (rdiff1 L q) D) = rlength (rinter L D) :=
        rlength_inter_rdiff1_disjoint L hq D (fun d hd => Or.inr (hqd d hd))
      have hpart := rlength_diff1_add_inter L hq
      have hcons : (q :: D : RRegion) = [q] ++ D := rfl
      rw [hswap, hcons, rlength_inter_append_right]
      linarith

end RealRegion
end LonelyRunner
