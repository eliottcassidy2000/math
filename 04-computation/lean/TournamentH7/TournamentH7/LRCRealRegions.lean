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

open LRC14

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

/-! ## Stage 3 — the exact depletion identity and THE HUNTER LEDGER -/

/-- Intersection mass decomposes over the window list's singletons. -/
theorem rlength_inter_singletons (A : RRegion) :
    ∀ D : RRegion, rlength (rinter A D) = (D.map fun d => rlength (rinter A [d])).sum := by
  intro D
  induction D with
  | nil => rw [rlength_inter_nil]; rfl
  | cons d D ih =>
      have hcons : (d :: D : RRegion) = [d] ++ D := rfl
      rw [hcons, rlength_inter_append_right, ih, List.map_append, List.sum_append]
      simp

/-- Subtracting one live interval shrinks mass. -/
theorem rlength_rdiff1_le (L : RRegion) {q : ℝ × ℝ} (hq : q.1 ≤ q.2) :
    rlength (rdiff1 L q) ≤ rlength L := by
  have h := rlength_diff1_add_inter L hq
  have h2 := rlength_nonneg (rinter L [q])
  linarith

/-- **THE EXACT DEPLETION IDENTITY**: intersecting a subtracted region with a window
list loses EXACTLY the per-window joint mass with the subtrahend.  (Subtrahend `D'`
sorted-separated live; the window list `Di` arbitrary.) -/
theorem rlength_inter_rdiff_expand (D' : RRegion) (hlive : ∀ r ∈ D', r.1 ≤ r.2)
    (hsep : SortedSep D') (X : RRegion) (Di : RRegion) :
    rlength (rinter (rdiff X D') Di)
      = rlength (rinter X Di)
        - (Di.map fun d => rlength (rinter (rinter X [d]) D')).sum := by
  induction Di with
  | nil =>
      rw [rlength_inter_nil, rlength_inter_nil]
      simp only [List.map_nil, List.sum_nil]
      ring
  | cons d Di ih =>
      have hconsd : (d :: Di : RRegion) = [d] ++ Di := rfl
      rw [hconsd, rlength_inter_append_right, rlength_inter_append_right, ih]
      have hhead : rlength (rinter (rdiff X D') [d])
          = rlength (rinter X [d]) - rlength (rinter (rinter X [d]) D') := by
        rw [rinter_rdiff_single, rlength_rdiff_partition D' (rinter X [d]) hlive hsep]
      rw [hhead]
      simp only [List.map_append, List.map_cons, List.map_nil, List.sum_append,
        List.sum_cons, List.sum_nil, add_zero]
      ring

/-- Depletion monotonicity: subtraction only shrinks any window-intersection mass. -/
theorem rlength_inter_rdiff_le (D' : RRegion) (hlive : ∀ r ∈ D', r.1 ≤ r.2)
    (hsep : SortedSep D') (X : RRegion) (B : RRegion) :
    rlength (rinter (rdiff X D') B) ≤ rlength (rinter X B) := by
  rw [rlength_inter_rdiff_expand D' hlive hsep X B]
  have hnn : 0 ≤ (B.map fun d => rlength (rinter (rinter X [d]) D')).sum := by
    apply List.sum_nonneg
    intro x hx
    rw [List.mem_map] at hx
    obtain ⟨d, _, rfl⟩ := hx
    exact rlength_nonneg _
  linarith

/-- The pair credits of the sequential peel: for each CONSECUTIVE pair of danger
regions, the joint mass of the second's windows with the first, measured on the
region as depleted so far.  These are the Hunter/path-Bonferroni credit terms. -/
noncomputable def pairCredits (I : RRegion) : List RRegion → ℝ
  | [] => 0
  | [_] => 0
  | D₁ :: D₂ :: rest =>
      ((D₂.map fun d => rlength (rinter (rinter I [d]) D₁)).sum)
        + pairCredits (rdiff I D₁) (D₂ :: rest)

@[simp] theorem pairCredits_nil (I : RRegion) : pairCredits I [] = 0 := rfl

@[simp] theorem pairCredits_single (I : RRegion) (D : RRegion) :
    pairCredits I [D] = 0 := rfl

theorem pairCredits_cons₂ (I : RRegion) (D₁ D₂ : RRegion) (rest : List RRegion) :
    pairCredits I (D₁ :: D₂ :: rest)
      = ((D₂.map fun d => rlength (rinter (rinter I [d]) D₁)).sum)
        + pairCredits (rdiff I D₁) (D₂ :: rest) := rfl

/-- **THE HUNTER LEDGER** (path-tree Bonferroni, EXACT credits): the surviving mass
of a sequential peel is at least the window mass, minus each danger's FULL mass on
the ORIGINAAL window, PLUS the consecutive-pair credits.  No loss terms: the pair
credits are measured exactly where they arise. -/
theorem hunter_ledger :
    ∀ (Ds : List RRegion), (∀ D ∈ Ds, ∀ r ∈ D, r.1 ≤ r.2) → (∀ D ∈ Ds, SortedSep D) →
    ∀ (I : RRegion),
    rlength (Ds.foldl rdiff I)
      ≥ rlength I - (Ds.map fun D => rlength (rinter I D)).sum + pairCredits I Ds := by
  intro Ds
  induction Ds with
  | nil =>
      intro _ _ I
      simp
  | cons D₁ tail ih =>
      intro hlive hsep I
      have hliveD₁ : ∀ r ∈ D₁, r.1 ≤ r.2 := hlive D₁ (List.mem_cons_self ..)
      have hsepD₁ : SortedSep D₁ := hsep D₁ (List.mem_cons_self ..)
      have hlivetail : ∀ D ∈ tail, ∀ r ∈ D, r.1 ≤ r.2 :=
        fun D hD => hlive D (List.mem_cons_of_mem _ hD)
      have hseptail : ∀ D ∈ tail, SortedSep D :=
        fun D hD => hsep D (List.mem_cons_of_mem _ hD)
      have hfold : (D₁ :: tail).foldl rdiff I = tail.foldl rdiff (rdiff I D₁) :=
        List.foldl_cons ..
      have hI' : rlength (rdiff I D₁) = rlength I - rlength (rinter I D₁) :=
        rlength_rdiff_partition D₁ I hliveD₁ hsepD₁
      cases tail with
      | nil =>
          simp only [hfold, List.foldl_nil, List.map_cons, List.map_nil, List.sum_cons,
            List.sum_nil, pairCredits_single]
          linarith
      | cons D₂ rest =>
          have hih := ih hlivetail hseptail (rdiff I D₁)
          -- expand the D₂-mass on the depleted region: the head credit appears exactly
          have hexp : rlength (rinter (rdiff I D₁) D₂)
              = rlength (rinter I D₂)
                - ((D₂.map fun d => rlength (rinter (rinter I [d]) D₁)).sum) :=
            rlength_inter_rdiff_expand D₁ hliveD₁ hsepD₁ I D₂
          -- the rest of the masses only shrink under depletion
          have hrest : ((rest.map fun D => rlength (rinter (rdiff I D₁) D)).sum)
              ≤ (rest.map fun D => rlength (rinter I D)).sum := by
            apply List.sum_le_sum
            intro D hD
            exact rlength_inter_rdiff_le D₁ hliveD₁ hsepD₁ I D
          rw [hfold]
          rw [pairCredits_cons₂]
          simp only [List.map_cons, List.sum_cons] at hih ⊢
          rw [hexp] at hih
          linarith [hih, hrest, hI']

/-! ## Stage 4 — the semantic bridge: points of the peeled region are good -/

/-- A region of positive mass contains a nondegenerate interval. -/
theorem exists_pos_interval_of_rlength_pos {L : RRegion} (h : 0 < rlength L) :
    ∃ p ∈ L, p.1 < p.2 := by
  induction L with
  | nil =>
      exfalso
      unfold rlength at h
      simp at h
  | cons p L ih =>
      have hsplit : rlength (p :: L) = max 0 (p.2 - p.1) + rlength L := by
        unfold rlength
        rw [List.map_cons, List.sum_cons]
      rw [hsplit] at h
      by_cases hplt : p.1 < p.2
      · exact ⟨p, List.mem_cons_self .., hplt⟩
      · push_neg at hplt
        rw [max_eq_left (by linarith)] at h
        obtain ⟨q, hq, hqlt⟩ := ih (by linarith)
        exact ⟨q, List.mem_cons_of_mem _ hq, hqlt⟩

/-- **Fold-membership semantics**: any point of any interval of the peeled region
lies in an original interval and (weakly) avoids every subtracted interval. -/
theorem rdiff_point_good :
    ∀ (D : List (ℝ × ℝ)) (L : RRegion) (p' : ℝ × ℝ), p' ∈ rdiff L D →
    ∀ t : ℝ, p'.1 ≤ t → t ≤ p'.2 →
    (∃ p ∈ L, p.1 ≤ t ∧ t ≤ p.2) ∧ (∀ q ∈ D, t ≤ q.1 ∨ q.2 ≤ t) := by
  intro D
  induction D with
  | nil =>
      intro L p' hmem t ht1 ht2
      have hmem' : p' ∈ L := by simpa [rdiff] using hmem
      exact ⟨⟨p', hmem', ht1, ht2⟩, fun q hq => absurd hq List.not_mem_nil⟩
  | cons q D ih =>
      intro L p' hmem t ht1 ht2
      have hfold : rdiff L (q :: D) = rdiff (rdiff1 L q) D := by
        unfold rdiff
        rw [List.foldl_cons]
      rw [hfold] at hmem
      obtain ⟨⟨p'', hp''mem, hp''1, hp''2⟩, havoidD⟩ := ih (rdiff1 L q) p' hmem t ht1 ht2
      unfold rdiff1 at hp''mem
      rw [List.mem_flatMap] at hp''mem
      obtain ⟨p, hpL, hp''cut⟩ := hp''mem
      unfold rcut at hp''cut
      simp only [List.mem_cons, List.not_mem_nil, or_false] at hp''cut
      rcases hp''cut with hpc | hpc
      · -- left piece: (p.1, min p.2 q.1)
        subst hpc
        simp only at hp''1 hp''2
        have htq : t ≤ q.1 := le_trans hp''2 (min_le_right _ _)
        have htp2 : t ≤ p.2 := le_trans hp''2 (min_le_left _ _)
        refine ⟨⟨p, hpL, hp''1, htp2⟩, ?_⟩
        intro r hr
        rcases List.mem_cons.mp hr with rfl | hrD
        · exact Or.inl htq
        · exact havoidD r hrD
      · -- right piece: (max p.1 q.2, p.2)
        subst hpc
        simp only at hp''1 hp''2
        have htq : q.2 ≤ t := le_trans (le_max_right _ _) hp''1
        have htp1 : p.1 ≤ t := le_trans (le_max_left _ _) hp''1
        refine ⟨⟨p, hpL, htp1, hp''2⟩, ?_⟩
        intro r hr
        rcases List.mem_cons.mp hr with rfl | hrD
        · exact Or.inr htq
        · exact havoidD r hrD

/-- Chain-level fold semantics: a point of the fully peeled chain avoids every
interval of every danger region. -/
theorem rdiff_chain_point_good :
    ∀ (Ds : List RRegion) (L : RRegion) (p' : ℝ × ℝ), p' ∈ Ds.foldl rdiff L →
    ∀ t : ℝ, p'.1 ≤ t → t ≤ p'.2 →
    (∃ p ∈ L, p.1 ≤ t ∧ t ≤ p.2) ∧ (∀ D ∈ Ds, ∀ q ∈ D, t ≤ q.1 ∨ q.2 ≤ t) := by
  intro Ds
  induction Ds with
  | nil =>
      intro L p' hmem t ht1 ht2
      exact ⟨⟨p', by simpa using hmem, ht1, ht2⟩, fun D hD => absurd hD List.not_mem_nil⟩
  | cons D Ds ih =>
      intro L p' hmem t ht1 ht2
      rw [List.foldl_cons] at hmem
      obtain ⟨⟨p'', hp''mem, h1, h2⟩, hrest⟩ := ih (rdiff L D) p' hmem t ht1 ht2
      obtain ⟨⟨p, hpL, hp1, hp2⟩, hD⟩ := rdiff_point_good D L p'' hp''mem t h1 h2
      refine ⟨⟨p, hpL, hp1, hp2⟩, ?_⟩
      intro D' hD'
      rcases List.mem_cons.mp hD' with rfl | hD'2
      · exact hD
      · exact hrest D' hD'2

/-- `SortedSep` for a mapped range with separated consecutive images. -/
theorem sortedSep_map_range :
    ∀ (n : ℕ) (f : ℕ → ℝ × ℝ), (∀ i j : ℕ, i < j → j < n → (f i).2 ≤ (f j).1) →
    SortedSep ((List.range n).map f) := by
  intro n
  induction n with
  | zero =>
      intro f _
      simp [SortedSep]
  | succ n ih =>
      intro f hsep
      rw [List.range_succ_eq_map, List.map_cons, List.map_map]
      refine ⟨?_, ?_⟩
      · intro d hd
        rw [List.mem_map] at hd
        obtain ⟨i, hi, rfl⟩ := hd
        rw [List.mem_range] at hi
        exact hsep 0 (i + 1) (Nat.succ_pos i) (by omega)
      · exact ih (f ∘ Nat.succ) fun i j hij hj =>
          hsep (i + 1) (j + 1) (by omega) (by omega)

/-- Each tooth is a live interval. -/
theorem teeth_live {w : ℤ} (hw : 0 < w) (a b : ℝ) :
    ∀ r ∈ teeth w a b, r.1 ≤ r.2 := by
  intro r hr
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  unfold teeth at hr
  rw [List.mem_map] at hr
  obtain ⟨i, _, rfl⟩ := hr
  unfold tooth
  simp only
  gcongr
  linarith

/-- A runner's teeth are sorted-separated (gap `6/(7w)` between consecutive teeth). -/
theorem teeth_sortedSep {w : ℤ} (hw : 0 < w) (a b : ℝ) : SortedSep (teeth w a b) := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  unfold teeth
  apply sortedSep_map_range
  intro i j hij _
  unfold tooth
  simp only
  have hij' : (i : ℝ) + 1 ≤ (j : ℝ) := by exact_mod_cast hij
  gcongr
  push_cast
  linarith

/-- The window-intersection mass IS the clipped-teeth sum of `LRCBlockSix`. -/
theorem rlength_inter_window_clipsum (a b : ℝ) (D : RRegion) :
    rlength (rinter [(a, b)] D) = (D.map fun p => clipLen p a b).sum := by
  have hflat : rinter [(a, b)] D = D.map fun q => rclip (a, b) q := by
    unfold rinter
    simp [List.flatMap_cons]
  rw [hflat]
  unfold rlength
  rw [List.map_map]
  congr 1
  apply List.map_congr_left
  intro q _
  unfold rclip clipLen
  simp only [Function.comp_apply]
  rw [min_comm b q.2, max_comm a q.1]

/-- **THE HUNTER BLOCK STEP**: if the exact ledger — window mass, minus full teeth
masses, plus consecutive-pair credits — is positive, the runner block has a common
1/14-good point in the window. -/
theorem hunter_block_step (ws : List ℤ) (hpos : ∀ w ∈ ws, 0 < w) (a b : ℝ) (hab : a ≤ b)
    (hledger : 0 < (b - a)
        - ((ws.map fun (w : ℤ) => rlength (rinter [(a, b)] (teeth w a b))).sum)
        + pairCredits [(a, b)] (ws.map fun (w : ℤ) => teeth w a b)) :
    ∃ t : ℝ, a ≤ t ∧ t ≤ b ∧
      ∀ w ∈ ws, ∀ m : ℤ, (1 : ℝ) / 14 ≤ |(w : ℝ) * t - m| := by
  set Ds : List RRegion := ws.map fun (w : ℤ) => teeth w a b with hDs
  have hlive : ∀ D ∈ Ds, ∀ r ∈ D, r.1 ≤ r.2 := by
    intro D hD
    rw [hDs, List.mem_map] at hD
    obtain ⟨w, hw, rfl⟩ := hD
    exact teeth_live (hpos w hw) a b
  have hsep : ∀ D ∈ Ds, SortedSep D := by
    intro D hD
    rw [hDs, List.mem_map] at hD
    obtain ⟨w, hw, rfl⟩ := hD
    exact teeth_sortedSep (hpos w hw) a b
  have hIlen : rlength [(a, b)] = b - a := by
    unfold rlength
    simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, add_zero]
    exact max_eq_right (by linarith)
  have hled := hunter_ledger Ds hlive hsep [(a, b)]
  have hsum_eq : (Ds.map fun D => rlength (rinter [(a, b)] D)).sum
      = (ws.map fun (w : ℤ) => rlength (rinter [(a, b)] (teeth w a b))).sum := by
    rw [hDs, List.map_map]
    rfl
  have hpos' : 0 < rlength (Ds.foldl rdiff [(a, b)]) := by
    rw [hIlen, hsum_eq] at hled
    linarith
  obtain ⟨p', hp'mem, hp'lt⟩ := exists_pos_interval_of_rlength_pos hpos'
  obtain ⟨⟨p, hpI, hpt1, hpt2⟩, havoid⟩ :=
    rdiff_chain_point_good Ds [(a, b)] p' hp'mem p'.1 (le_refl _) (le_of_lt hp'lt)
  have hpab : p = (a, b) := by
    rcases List.mem_cons.mp hpI with rfl | hf
    · rfl
    · exact absurd hf List.not_mem_nil
  subst hpab
  simp only at hpt1 hpt2
  refine ⟨p'.1, hpt1, hpt2, ?_⟩
  intro w hw m
  apply good_of_avoid_teeth (hpos w hw) hpt1 hpt2
  intro r hr
  exact havoid (teeth w a b) (by rw [hDs, List.mem_map]; exact ⟨w, hw, rfl⟩) r hr

end RealRegion
end LonelyRunner
