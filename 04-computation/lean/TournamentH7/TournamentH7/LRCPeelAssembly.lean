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
  have hFl : ∀ L' : Region, length (diff1F L' q) = length (diff1 L' q) := by
    intro L'
    induction L' with
    | nil => rfl
    | cons p L' ih =>
        unfold diff1F diff1 at ih ⊢
        simp only [List.flatMap_cons]
        rw [length_append, length_append, ih, length_cutF]
  rw [hFl L]
  exact length_diff1_add_inter L hq

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

/-- Filtering degenerate pieces does not change intersection lengths (their clips are
degenerate too). -/
theorem length_inter_filter_live (X B : Region) :
    length (inter (X.filter fun r => decide (r.1 < r.2)) B) = length (inter X B) := by
  induction X with
  | nil => rfl
  | cons p X ih =>
      by_cases hp : p.1 < p.2
      · rw [List.filter_cons_of_pos (by simpa using hp)]
        unfold inter at ih ⊢
        simp only [List.flatMap_cons]
        rw [length_append, length_append, ih]
      · rw [List.filter_cons_of_neg (by simpa using hp)]
        unfold inter at ih ⊢
        simp only [List.flatMap_cons]
        rw [length_append, ih]
        have hzero : length (B.map fun r => clip p r) = 0 := by
          unfold length
          apply List.sum_eq_zero
          intro y hy
          simp only [List.map_map, List.mem_map] at hy
          obtain ⟨r, _, rfl⟩ := hy
          simp only [Function.comp_apply]
          apply max_eq_left
          unfold clip
          simp only
          push_neg at hp
          have h1 : min p.2 r.2 ≤ p.2 := min_le_left _ _
          have h2 : p.1 ≤ max p.1 r.1 := le_max_left _ _
          linarith
        rw [hzero]
        ring

/-- Clipping the two unfiltered cut pieces against a region totals at most the whole
interval's clips. -/
theorem length_inter_cut_le (p : ℚ × ℚ) {q : ℚ × ℚ} (hq : q.1 ≤ q.2) (B : Region) :
    length (inter (cut p q) B) ≤ length (inter [p] B) := by
  unfold inter cut
  simp only [List.flatMap_cons, List.flatMap_nil, List.append_nil]
  rw [length_append]
  induction B with
  | nil =>
      norm_num [length]
  | cons s B ihB =>
      unfold length at ihB ⊢
      simp only [List.map_cons, List.sum_cons]
      have hper := clip_cut_pieces_le p hq s
      simp only [clip] at ihB ⊢
      linarith [hper, ihB]

/-- Intersections only shrink under a live subtraction. -/
theorem length_inter_diff1F_le (L : Region) {q : ℚ × ℚ} (hq : q.1 ≤ q.2) (B : Region) :
    length (inter (diff1F L q) B) ≤ length (inter L B) := by
  induction L with
  | nil =>
      unfold diff1F inter
      simp [length]
  | cons p L ih =>
      have hsplit : diff1F (p :: L) q = cutF p q ++ diff1F L q := by
        unfold diff1F
        simp [List.flatMap_cons]
      have h1 : inter (cutF p q ++ diff1F L q) B
          = inter (cutF p q) B ++ inter (diff1F L q) B := by
        unfold inter
        simp [List.flatMap_append]
      have h2 : inter (p :: L) B = inter [p] B ++ inter L B := by
        unfold inter
        simp [List.flatMap_cons]
      rw [hsplit, h1, length_append, h2, length_append]
      have hcutF : length (inter (cutF p q) B) = length (inter (cut p q) B) :=
        length_inter_filter_live _ _
      have hcut := length_inter_cut_le p hq B
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

/-! ## Stage 2: the half-comb danger pair and the damped peel

The phase-0 danger comb wraps at the seam, but the danger set splits as
`frac(sx) ∈ [0,h) ∪ [1−h,1)` — TWO NON-WRAPPING combs at radius `h/2`, phases `h/2`
and `1 − h/2`.  The S11 rate lemma applies to each half verbatim; the completeness
bridge needs no wrap layer at all. -/

/-- The danger region of speed `s` at band `h`, as two non-wrapping half-combs. -/
def dangerPair (s : ℕ) (h : ℚ) : Region :=
  comb s (h / 2) (h / 2) ++ comb s (h / 2) (1 - h / 2)

/-- The half-comb completeness bridge (wrap-free): a point of `[0,1)` avoiding both
halves is `h`-far from every integer. -/
theorem not_mem_dangerPair_forall {s : ℕ} (hs : 0 < s) {h : ℚ} (hh0 : 0 < h)
    (hh1 : h ≤ 1 / 2) {x : ℚ} (hx0 : 0 ≤ x) (hx1 : x < 1)
    (hnot : ¬ mem x (dangerPair s h)) :
    ∀ m : ℤ, h ≤ |(s : ℚ) * x - m| := by
  intro m
  by_contra hcon
  apply hnot
  have habs : |(s : ℚ) * x - m| < h := lt_of_not_ge hcon
  rw [abs_lt] at habs
  have hsQ : (0 : ℚ) < (s : ℚ) := by exact_mod_cast hs
  have hsx0 : (0 : ℚ) ≤ (s : ℚ) * x := mul_nonneg hsQ.le hx0
  have hsx1 : (s : ℚ) * x < (s : ℚ) := by
    calc (s : ℚ) * x < (s : ℚ) * 1 := mul_lt_mul_of_pos_left hx1 hsQ
      _ = (s : ℚ) := mul_one _
  unfold dangerPair
  rw [mem_union]
  rcases le_or_gt (m : ℚ) ((s : ℚ) * x) with hm | hm
  · -- 0 ≤ sx − m < h: the left half, tooth k = m
    left
    rw [TournamentH7.CombPatterns.mem_comb hs]
    have hm0 : (0 : ℤ) ≤ m := by
      by_contra hneg
      have hm1 : m ≤ -1 := by omega
      have : (m : ℚ) ≤ -1 := by exact_mod_cast hm1
      linarith [habs.2]
    have hms : m < (s : ℤ) := by
      by_contra hge
      have hge2 : (s : ℤ) ≤ m := by omega
      have : (s : ℚ) ≤ (m : ℚ) := by exact_mod_cast hge2
      linarith
    refine ⟨m.toNat, ?_, ?_, ?_⟩
    · omega
    · have htn : ((m.toNat : ℕ) : ℚ) = (m : ℚ) := by exact_mod_cast Int.toNat_of_nonneg hm0
      rw [htn]
      linarith
    · have htn : ((m.toNat : ℕ) : ℚ) = (m : ℚ) := by exact_mod_cast Int.toNat_of_nonneg hm0
      rw [htn]
      linarith [habs.2]
  · -- −h < sx − m < 0: the right half, tooth k = m − 1
    right
    rw [TournamentH7.CombPatterns.mem_comb hs]
    have hm1 : (1 : ℤ) ≤ m := by
      by_contra hneg
      have hm2 : m ≤ 0 := by omega
      have : (m : ℚ) ≤ 0 := by exact_mod_cast hm2
      linarith
    have hms : m ≤ (s : ℤ) := by
      by_contra hge
      have hge2 : (s : ℤ) + 1 ≤ m := by omega
      have : (s : ℚ) + 1 ≤ (m : ℚ) := by exact_mod_cast hge2
      linarith [habs.1]
    refine ⟨(m - 1).toNat, ?_, ?_, ?_⟩
    · omega
    · have htn : (((m - 1).toNat : ℕ) : ℚ) = ((m : ℚ) - 1) := by
        have h1 : ((m - 1).toNat : ℤ) = m - 1 := Int.toNat_of_nonneg (by omega)
        have h2 : (((m - 1).toNat : ℕ) : ℚ) = (((m - 1 : ℤ)) : ℚ) := by exact_mod_cast h1
        rw [h2]
        push_cast
        ring
      rw [htn]
      linarith [habs.1]
    · have htn : (((m - 1).toNat : ℕ) : ℚ) = ((m : ℚ) - 1) := by
        have h1 : ((m - 1).toNat : ℤ) = m - 1 := Int.toNat_of_nonneg (by omega)
        have h2 : (((m - 1).toNat : ℕ) : ℚ) = (((m - 1 : ℤ)) : ℚ) := by exact_mod_cast h1
        rw [h2]
        push_cast
        ring
      rw [htn]
      linarith

/-- The wrap-free good region: the unit window minus every speed's danger pair. -/
def goodRegion2 (speeds : List ℤ) (h : ℚ) : Region :=
  diffF [((0 : ℚ), 1)] (speeds.flatMap fun s => dangerPair s.toNat h)

/-- Peeling identity: extending the speed list subtracts the new danger pair from the
good region. -/
theorem goodRegion2_append (E : List ℤ) (w : ℤ) (h : ℚ) :
    goodRegion2 (E ++ [w]) h = diffF (goodRegion2 E h) (dangerPair w.toNat h) := by
  unfold goodRegion2 diffF
  rw [List.flatMap_append, List.foldl_append]
  simp [List.flatMap_cons]

/-- Every piece of a `diffF` from the unit window is live and inside `[0,1]`. -/
theorem diffF_window_bounds :
    ∀ (B : Region) (L : Region), (∀ p ∈ L, 0 ≤ p.1 ∧ p.1 ≤ p.2 ∧ p.2 ≤ 1) →
      ∀ p ∈ diffF L B, 0 ≤ p.1 ∧ p.1 ≤ p.2 ∧ p.2 ≤ 1 := by
  intro B
  induction B with
  | nil => intro L hL; simpa [diffF] using hL
  | cons q B ih =>
      intro L hL p hp
      have hfold : diffF L (q :: B) = diffF (diff1F L q) B := by
        unfold diffF
        simp [List.foldl_cons]
      rw [hfold] at hp
      refine ih (diff1F L q) ?_ p hp
      intro r hr
      unfold diff1F at hr
      simp only [List.mem_flatMap] at hr
      obtain ⟨u, hu, hru⟩ := hr
      obtain ⟨hu0, hu12, hu1⟩ := hL u hu
      unfold cutF cut at hru
      rw [List.mem_filter] at hru
      obtain ⟨hmem, hlive⟩ := hru
      have hlive' : r.1 < r.2 := by simpa using hlive
      rcases List.mem_cons.mp hmem with rfl | hmem'
      · refine ⟨hu0, le_of_lt hlive', ?_⟩
        simp only
        exact le_trans (min_le_left _ _) hu1
      · rcases List.mem_singleton.mp hmem' with rfl
        exact ⟨le_trans hu0 (le_max_left _ _), le_of_lt hlive', hu1⟩

/-- **THE DAMPED PEEL THEOREM**: adding a far speed `w` to a family costs the good region
at most the damped factor `2h` of its length plus `4h/w` per piece:

    length (good (E ++ [w])) ≥ (1 − 2h)·length (good E) − (#pieces)·4h/w. -/
theorem damped_peel (E : List ℤ) (w : ℤ) (hw : 0 < w) {h : ℚ} (hh0 : 0 < h)
    (hh1 : h ≤ 1 / 2) :
    (1 - 2 * h) * length (goodRegion2 E h)
      - ((goodRegion2 E h).length : ℚ) * (4 * h / w.toNat)
    ≤ length (goodRegion2 (E ++ [w]) h) := by
  set G : Region := goodRegion2 E h with hG
  have hGbounds : ∀ p ∈ G, 0 ≤ p.1 ∧ p.1 ≤ p.2 ∧ p.2 ≤ 1 := by
    rw [hG]
    unfold goodRegion2
    apply diffF_window_bounds
    intro p hp
    rcases List.mem_singleton.mp hp with rfl
    norm_num
  have hwnat : 0 < w.toNat := by omega
  have hwQ : (0 : ℚ) < (w.toNat : ℚ) := by exact_mod_cast hwnat
  have hlive : ∀ q ∈ dangerPair w.toNat h, q.1 ≤ q.2 := by
    intro q hq
    unfold dangerPair at hq
    rcases List.mem_append.mp hq with hq' | hq' <;>
    · unfold comb at hq'
      simp only [List.map_map, List.mem_map, List.mem_range, Function.comp_apply] at hq'
      obtain ⟨k, _, rfl⟩ := hq'
      simp only
      apply div_le_div_of_nonneg_right ?_ hwQ.le
      linarith
  rw [goodRegion2_append, ← hG]
  have hstage1 := length_diffF_ge G (dangerPair w.toNat h) hlive
  have hsplit : length (inter G (dangerPair w.toNat h))
      = length (inter G (comb w.toNat (h / 2) (h / 2)))
        + length (inter G (comb w.toNat (h / 2) (1 - h / 2))) := by
    unfold dangerPair
    exact length_inter_append_right G _ _
  have hrate1 := length_inter_comb_near_region (w := w.toNat) hwnat
    (r := h / 2) (φ := h / 2) (by linarith) (le_refl _) (by linarith) G hGbounds
  have hrate2 := length_inter_comb_near_region (w := w.toNat) hwnat
    (r := h / 2) (φ := 1 - h / 2) (by linarith) (by linarith) (by linarith) G hGbounds
  have habs1 := (abs_le.mp hrate1).2
  have habs2 := (abs_le.mp hrate2).2
  have hbound : length (inter G (dangerPair w.toNat h))
      ≤ 2 * h * length G + (G.length : ℚ) * (4 * h / (w.toNat : ℚ)) := by
    rw [hsplit]
    have e1 : 2 * (h / 2) * length G + (G.length : ℚ) * (4 * (h / 2) / (w.toNat : ℚ))
        = h * length G + (G.length : ℚ) * (2 * h / (w.toNat : ℚ)) := by ring
    have hb1 : length (inter G (comb w.toNat (h / 2) (h / 2)))
        ≤ h * length G + (G.length : ℚ) * (2 * h / (w.toNat : ℚ)) := by
      rw [← e1]
      linarith
    have hb2 : length (inter G (comb w.toNat (h / 2) (1 - h / 2)))
        ≤ h * length G + (G.length : ℚ) * (2 * h / (w.toNat : ℚ)) := by
      rw [← e1]
      linarith
    have : (G.length : ℚ) * (2 * h / (w.toNat : ℚ)) * 2
        = (G.length : ℚ) * (4 * h / (w.toNat : ℚ)) := by ring
    linarith
  linarith

/-! ## The peel gate and the wrap-free composite gate -/

/-- **THE QUANTITATIVE PEEL GATE**: if the good region is at least `ε` long and the new
speed clears the explicit threshold, the extended good region is still positive.  This
is the step that ITERATES: peel far speeds one by one, spending `(1−2h)` damping and an
`O(1/w)` boundary fee each time. -/
theorem goodRegion2_pos_of_peel (E : List ℤ) (w : ℤ) (hw : 0 < w) {h ε : ℚ}
    (hh0 : 0 < h) (hh1 : h ≤ 1 / 2)
    (hε : ε ≤ length (goodRegion2 E h))
    (hthresh : ((goodRegion2 E h).length : ℚ) * (4 * h / w.toNat) < (1 - 2 * h) * ε) :
    0 < length (goodRegion2 (E ++ [w]) h) := by
  have hpeel := damped_peel E w hw hh0 hh1
  have h2h : (0 : ℚ) ≤ 1 - 2 * h := by linarith
  have hmono : (1 - 2 * h) * ε ≤ (1 - 2 * h) * length (goodRegion2 E h) :=
    mul_le_mul_of_nonneg_left hε h2h
  linarith

/-- The wrap-free member-safety (mirror of `good_mem_safe` over `dangerPair`). -/
theorem good2_mem_safe {speeds : List ℤ} {h : ℚ} (hh0 : 0 < h) (hh1 : h ≤ 1 / 2)
    (hpos : ∀ s ∈ speeds, 0 < s) {x : ℚ} (hx : mem x (goodRegion2 speeds h)) :
    (0 ≤ x ∧ x < 1) ∧ ∀ s ∈ speeds, ∀ m : ℤ, h ≤ |(s : ℚ) * x - m| := by
  obtain ⟨hwin, havoid⟩ := mem_diffF.mp hx
  obtain ⟨p, hp, hx1, hx2⟩ := hwin
  rcases List.mem_singleton.mp hp with rfl
  simp only at hx1 hx2
  refine ⟨⟨hx1, hx2⟩, ?_⟩
  intro s hs m
  have hspos := hpos s hs
  have hsnat : 0 < s.toNat := by omega
  have hbridge := not_mem_dangerPair_forall hsnat hh0 hh1 hx1 hx2 (fun hmem => by
    obtain ⟨q, hq, hxq⟩ := hmem
    exact havoid q (List.mem_flatMap.mpr ⟨s, hs, hq⟩) hxq) m
  have hcast : ((s.toNat : ℕ) : ℚ) = ((s : ℤ) : ℚ) := by
    exact_mod_cast Int.toNat_of_nonneg hspos.le
  rwa [hcast] at hbridge

/-- **The wrap-free composite gate**: a 13-family of positive speeds whose wrap-free good
region has positive length is lonely. -/
theorem exists_lonely_of_goodRegion2_pos (v : Fin 13 → ℤ) (hv : ∀ i, 0 < v i)
    (hlen : 0 < length (goodRegion2 (List.ofFn v) (1 / 14))) :
    ∃ t : ℝ, Lonely 14 v t := by
  obtain ⟨x, hx⟩ := exists_mem_of_length_pos hlen
  have hpos : ∀ s ∈ List.ofFn v, 0 < s := by
    intro s hs
    obtain ⟨i, rfl⟩ := List.mem_ofFn.mp hs
    exact hv i
  have hsafe := good2_mem_safe (by norm_num) (by norm_num) hpos hx
  exact ⟨((x : ℚ) : ℝ), LRC14.lonely_of_rat_forall fun i m =>
    hsafe.2 (v i) (List.mem_ofFn.mpr ⟨i, rfl⟩) m⟩

end RatIntervals
end LonelyRunner
