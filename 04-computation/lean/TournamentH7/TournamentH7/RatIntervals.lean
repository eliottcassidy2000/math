/-
  TournamentH7.RatIntervals — MODULE 0 of the LRC(14) formalization playbook (T2).
  kind-pasteur-2026-07-02-S4 (HYP-3961).  The blocking mini-library: every danger set, comb,
  window, overlap, and uncovered set in the proof is a value of this type.

  DESIGN.  A region is a plain `List (ℚ × ℚ)` of half-open intervals `[a, b)`.  The SEMANTIC
  layer needs no invariants: membership is an existential over the list, and the quadratic
  intersection satisfies `mem_inter` unconditionally (pairwise clipping distributes over the
  existential).  The MEASURE layer (`length` behaving like Lebesgue measure) holds under the
  `Norm` discipline (ordered, disjoint, nondegenerate), which the constructors provide and the
  operations preserve.  All-ℚ (playbook T1): no reals, no measure theory.
-/
import Mathlib.Data.Rat.Floor

namespace LonelyRunner
namespace RatIntervals

/-- A region: a list of half-open rational intervals `[a, b)`. -/
abbrev Region := List (ℚ × ℚ)

/-- Semantic membership: `x` lies in some listed interval. -/
def mem (x : ℚ) (L : Region) : Prop := ∃ p ∈ L, p.1 ≤ x ∧ x < p.2

instance (x : ℚ) (L : Region) : Decidable (mem x L) := by
  unfold mem; infer_instance

/-- Total length (degenerate pairs contribute 0). -/
def length (L : Region) : ℚ := (L.map fun p => max 0 (p.2 - p.1)).sum

/-- Clip one interval against another (their intersection, possibly degenerate). -/
def clip (p q : ℚ × ℚ) : ℚ × ℚ := (max p.1 q.1, min p.2 q.2)

/-- Quadratic intersection: every pairwise clip. -/
def inter (A B : Region) : Region :=
  A.flatMap fun p => B.map fun q => clip p q

/-- Translation by `t`. -/
def translate (t : ℚ) (L : Region) : Region := L.map fun p => (p.1 + t, p.2 + t)

/-- Scaling by a positive rational `s`. -/
def scale (s : ℚ) (L : Region) : Region := L.map fun p => (s * p.1, s * p.2)

/-! ## The semantic layer (no invariants needed) -/

theorem mem_clip {x : ℚ} {p q : ℚ × ℚ} :
    (clip p q).1 ≤ x ∧ x < (clip p q).2 ↔ (p.1 ≤ x ∧ x < p.2) ∧ (q.1 ≤ x ∧ x < q.2) := by
  unfold clip
  simp only [max_le_iff, lt_min_iff]
  tauto

/-- **The core semantic lemma**: membership distributes over the quadratic intersection,
with no normalization hypotheses. -/
theorem mem_inter {x : ℚ} {A B : Region} :
    mem x (inter A B) ↔ mem x A ∧ mem x B := by
  unfold mem inter
  constructor
  · rintro ⟨r, hr, hx1, hx2⟩
    simp only [List.mem_flatMap, List.mem_map] at hr
    obtain ⟨p, hp, q, hq, rfl⟩ := hr
    have := mem_clip.mp ⟨hx1, hx2⟩
    exact ⟨⟨p, hp, this.1⟩, ⟨q, hq, this.2⟩⟩
  · rintro ⟨⟨p, hp, hxp⟩, ⟨q, hq, hxq⟩⟩
    refine ⟨clip p q, ?_, ?_⟩
    · simp only [List.mem_flatMap, List.mem_map]
      exact ⟨p, hp, q, hq, rfl⟩
    · exact mem_clip.mpr ⟨hxp, hxq⟩

theorem mem_translate {x t : ℚ} {L : Region} :
    mem x (translate t L) ↔ mem (x - t) L := by
  unfold mem translate
  constructor
  · rintro ⟨r, hr, hx1, hx2⟩
    simp only [List.mem_map] at hr
    obtain ⟨p, hp, rfl⟩ := hr
    exact ⟨p, hp, by constructor <;> [linarith; linarith]⟩
  · rintro ⟨p, hp, hx1, hx2⟩
    refine ⟨(p.1 + t, p.2 + t), ?_, ?_, ?_⟩
    · simp only [List.mem_map]; exact ⟨p, hp, rfl⟩
    · simpa using by linarith
    · simpa using by linarith

theorem mem_scale {x s : ℚ} (hs : 0 < s) {L : Region} :
    mem x (scale s L) ↔ mem (x / s) L := by
  unfold mem scale
  constructor
  · rintro ⟨r, hr, hx1, hx2⟩
    simp only [List.mem_map] at hr
    obtain ⟨p, hp, rfl⟩ := hr
    refine ⟨p, hp, ?_, ?_⟩
    · rw [le_div_iff₀ hs]; linarith [hx1]
    · rw [div_lt_iff₀ hs]; linarith [hx2]
  · rintro ⟨p, hp, hx1, hx2⟩
    refine ⟨(s * p.1, s * p.2), ?_, ?_, ?_⟩
    · simp only [List.mem_map]; exact ⟨p, hp, rfl⟩
    · simp only
      rw [le_div_iff₀ hs] at hx1; linarith
    · simp only
      rw [div_lt_iff₀ hs] at hx2; linarith

/-! ## The measure layer -/

theorem length_nonneg (L : Region) : 0 ≤ length L := by
  unfold length
  induction L with
  | nil => simp
  | cons p L ih =>
      simp only [List.map_cons, List.sum_cons]
      have h1 : (0 : ℚ) ≤ max 0 (p.2 - p.1) := le_max_left _ _
      linarith

theorem length_append (A B : Region) : length (A ++ B) = length A + length B := by
  unfold length
  rw [List.map_append, List.sum_append]

theorem length_translate (t : ℚ) (L : Region) : length (translate t L) = length L := by
  unfold length translate
  rw [List.map_map]
  congr 1
  apply List.map_congr_left
  intro p _
  simp only [Function.comp_apply]
  ring_nf

theorem length_scale {s : ℚ} (hs : 0 ≤ s) (L : Region) :
    length (scale s L) = s * length L := by
  unfold length scale
  rw [List.map_map]
  induction L with
  | nil => simp
  | cons p L ih =>
      simp only [List.map_cons, List.sum_cons, Function.comp_apply] at ih ⊢
      rw [ih, mul_add]
      congr 1
      rcases le_or_gt p.2 p.1 with h | h
      · have h1 : s * p.2 - s * p.1 ≤ 0 := by nlinarith
        have h2 : p.2 - p.1 ≤ 0 := by linarith
        rw [max_eq_left h1, max_eq_left h2, mul_zero]
      · have h1 : (0:ℚ) ≤ s * p.2 - s * p.1 := by nlinarith
        have h2 : (0:ℚ) ≤ p.2 - p.1 := by linarith
        rw [max_eq_right h1, max_eq_right h2]
        ring

/-- Length of a single-pair clip is at most the first factor's span. -/
theorem length_clip_le_left (p q : ℚ × ℚ) :
    max 0 ((clip p q).2 - (clip p q).1) ≤ max 0 (p.2 - p.1) := by
  unfold clip
  simp only
  rcases le_or_gt (min p.2 q.2 - max p.1 q.1) 0 with h | h
  · rw [max_eq_left h]
    exact le_max_left _ _
  · rw [max_eq_right h.le]
    have h1 : min p.2 q.2 ≤ p.2 := min_le_left _ _
    have h2 : p.1 ≤ max p.1 q.1 := le_max_left _ _
    have h3 : min p.2 q.2 - max p.1 q.1 ≤ p.2 - p.1 := by linarith
    calc min p.2 q.2 - max p.1 q.1 ≤ p.2 - p.1 := h3
      _ ≤ max 0 (p.2 - p.1) := le_max_right _ _

/-! ## The Norm discipline: ordered, disjoint, nondegenerate -/

/-- `Norm L`: intervals are nondegenerate and strictly ordered end-to-start
(hence pairwise disjoint). -/
def Norm : Region → Prop
  | [] => True
  | [p] => p.1 < p.2
  | p :: q :: L => p.1 < p.2 ∧ p.2 ≤ q.1 ∧ Norm (q :: L)

instance : ∀ L : Region, Decidable (Norm L)
  | [] => by unfold Norm; infer_instance
  | [p] => by unfold Norm; infer_instance
  | p :: q :: L => by
      have := instDecidableNorm (q :: L)
      unfold Norm
      infer_instance

theorem norm_tail {q : ℚ × ℚ} {B : Region} (h : Norm (q :: B)) : Norm B := by
  match B with
  | [] => trivial
  | r :: B' => exact h.2.2

theorem norm_head_lt {q : ℚ × ℚ} {B : Region} (h : Norm (q :: B)) : q.1 < q.2 := by
  match B with
  | [] => exact h
  | r :: B' => exact h.1

theorem norm_head_le {q : ℚ × ℚ} {B : Region} (h : Norm (q :: B)) :
    ∀ r ∈ B, q.2 ≤ r.1 := by
  induction B generalizing q with
  | nil => intro r hr; cases hr
  | cons r B' ih =>
      intro s hs
      rcases List.mem_cons.mp hs with rfl | hs
      · exact h.2.1
      · have h1 : q.2 ≤ r.1 := h.2.1
        have h2 : r.1 < r.2 := norm_head_lt h.2.2
        have := ih h.2.2 s hs
        linarith

/-- The cursor induction: for `Norm B` with every interval starting at or after `lo`, the total
clipped length against `(p₁, p₂)` is at most `max 0 (p₂ − max p₁ lo)`. -/
theorem length_map_clip_aux {p1 p2 : ℚ} :
    ∀ {B : Region}, Norm B → ∀ lo : ℚ, (∀ q ∈ B, lo ≤ q.1) →
      length (B.map fun q => clip (p1, p2) q) ≤ max 0 (p2 - max p1 lo) := by
  intro B
  induction B with
  | nil =>
      intro _ lo _
      unfold length
      simp
  | cons q B ih =>
      intro hN lo hlo
      have hq12 : q.1 < q.2 := norm_head_lt hN
      have hloq : lo ≤ q.1 := hlo q (List.mem_cons_self ..)
      have htail : ∀ r ∈ B, q.2 ≤ r.1 := norm_head_le hN
      have hrec := ih (norm_tail hN) q.2 htail
      unfold length at hrec ⊢
      simp only [List.map_cons, List.sum_cons, List.map_map] at hrec ⊢
      -- head bound: max 0 (min p2 q.2 − max p1 q.1) ≤ max 0 (min p2 q.2 − max p1 lo)
      -- tail bound: hrec ≤ max 0 (p2 − max p1 q.2)
      -- combine by the three-case analysis
      set A : ℚ := max p1 lo with hA
      set C : ℚ := q.2 with hC
      have hCgtlo : lo < C := by rw [hC]; linarith
      have hmaxC : A ≤ max p1 C := by
        rw [hA, hC]
        rcases le_or_gt p1 lo with h | h
        · rw [max_eq_right h]
          rcases le_or_gt p1 q.2 with h2 | h2
          · rw [max_eq_right h2]; linarith
          · rw [max_eq_left h2.le]; linarith
        · rw [max_eq_left h.le]; exact le_max_left _ _
      have hhead : max 0 ((clip (p1, p2) q).2 - (clip (p1, p2) q).1)
          ≤ max 0 (min p2 C - A) := by
        unfold clip
        simp only
        have h1 : A ≤ max p1 q.1 := by
          rw [hA]
          rcases le_or_gt p1 lo with h | h
          · rw [max_eq_right h]
            calc lo ≤ q.1 := hloq
              _ ≤ max p1 q.1 := le_max_right _ _
          · rw [max_eq_left h.le]; exact le_max_left _ _
        have h2 : min p2 q.2 - max p1 q.1 ≤ min p2 C - A := by
          rw [hC]; linarith
        rcases le_or_gt (min p2 q.2 - max p1 q.1) 0 with h3 | h3
        · rw [max_eq_left h3]; exact le_max_left _ _
        · rw [max_eq_right h3.le]
          calc min p2 q.2 - max p1 q.1 ≤ min p2 C - A := h2
            _ ≤ max 0 (min p2 C - A) := le_max_right _ _
      -- the three-case combination
      have hcomb : max 0 (min p2 C - A) + max 0 (p2 - max p1 C) ≤ max 0 (p2 - A) := by
        rcases le_or_gt (min p2 C) A with hc1 | hc1
        · -- first term 0
          rw [max_eq_left (by linarith)]
          have : p2 - max p1 C ≤ p2 - A := by linarith
          rcases le_or_gt (p2 - max p1 C) 0 with hc2 | hc2
          · rw [max_eq_left hc2]; simp [le_max_iff]
          · rw [max_eq_right hc2.le]
            have : p2 - A ≥ p2 - max p1 C := by linarith
            calc 0 + (p2 - max p1 C) = p2 - max p1 C := by ring
              _ ≤ p2 - A := by linarith
              _ ≤ max 0 (p2 - A) := le_max_right _ _
        · rcases le_or_gt p2 C with hc2 | hc2
          · -- p2 ≤ C: first = p2 − A > 0; second = 0
            have hmin : min p2 C = p2 := min_eq_left hc2
            rw [hmin] at hc1 ⊢
            have hsec : p2 - max p1 C ≤ 0 := by
              have : C ≤ max p1 C := le_max_right _ _
              linarith
            rw [max_eq_left hsec, max_eq_right (by linarith : (0:ℚ) ≤ p2 - A)]
            have : max 0 (p2 - A) = p2 - A := max_eq_right (by linarith)
            linarith [this.ge]
          · -- C < p2: first = C − A; second ≤ p2 − C
            have hmin : min p2 C = C := min_eq_right hc2.le
            rw [hmin] at hc1 ⊢
            have h1 : max 0 (C - A) = C - A := max_eq_right (by linarith)
            have h2 : max 0 (p2 - max p1 C) ≤ p2 - C := by
              have hCle : C ≤ max p1 C := le_max_right _ _
              rcases le_or_gt (p2 - max p1 C) 0 with h | h
              · rw [max_eq_left h]; linarith
              · rw [max_eq_right h.le]; linarith
            have h3 : (0:ℚ) ≤ p2 - A := by linarith
            rw [h1, max_eq_right h3]
            linarith
      calc max 0 ((clip (p1, p2) q).2 - (clip (p1, p2) q).1)
            + (List.map (fun x => max 0 ((clip (p1, p2) x).2 - (clip (p1, p2) x).1)) B).sum
          ≤ max 0 (min p2 C - A) + max 0 (p2 - max p1 C) := by
            have := hrec
            simp only [Function.comp] at this
            exact add_le_add hhead this
        _ ≤ max 0 (p2 - A) := hcomb
        _ = max 0 (p2 - max p1 lo) := by rw [hA]

/-- Under `Norm B`, clipping any single interval against `B` yields total length at most the
interval's (nonnegative) span. -/
theorem length_map_clip_le (p : ℚ × ℚ) {B : Region} (hB : Norm B) :
    length (B.map fun q => clip p q) ≤ max 0 (p.2 - p.1) := by
  match B with
  | [] =>
      unfold length; simp
  | q :: B' =>
      have := length_map_clip_aux (p1 := p.1) (p2 := p.2) hB q.1
        (by
          intro r hr
          rcases List.mem_cons.mp hr with rfl | hr
          · exact le_refl _
          · have h1 : q.2 ≤ r.1 := norm_head_le hB r hr
            have h2 : q.1 < q.2 := norm_head_lt hB
            linarith)
      calc length ((q :: B').map fun r => clip p r)
          ≤ max 0 (p.2 - max p.1 q.1) := this
        _ ≤ max 0 (p.2 - p.1) := by
            apply max_le (le_max_left _ _)
            have : p.1 ≤ max p.1 q.1 := le_max_left _ _
            calc p.2 - max p.1 q.1 ≤ p.2 - p.1 := by linarith
              _ ≤ max 0 (p.2 - p.1) := le_max_right _ _

/-- **Length monotonicity**: intersecting any region with a `Norm` region does not increase
length.  (No hypothesis on `A`.) -/
theorem length_inter_le_left (A : Region) {B : Region} (hB : Norm B) :
    length (inter A B) ≤ length A := by
  unfold inter
  induction A with
  | nil => unfold length; simp
  | cons p A ih =>
      simp only [List.flatMap_cons]
      rw [length_append]
      have h1 := length_map_clip_le p hB
      unfold length at h1 ih ⊢
      simp only [List.map_cons, List.sum_cons]
      linarith

/-- Union is append; membership distributes (no invariants). -/
theorem mem_union {x : ℚ} {A B : Region} : mem x (A ++ B) ↔ mem x A ∨ mem x B := by
  unfold mem
  constructor
  · rintro ⟨p, hp, hx⟩
    rcases List.mem_append.mp hp with h | h
    · exact Or.inl ⟨p, h, hx⟩
    · exact Or.inr ⟨p, h, hx⟩
  · rintro (⟨p, hp, hx⟩ | ⟨p, hp, hx⟩)
    · exact ⟨p, List.mem_append.mpr (Or.inl hp), hx⟩
    · exact ⟨p, List.mem_append.mpr (Or.inr hp), hx⟩

/-! ## The comb constructor (ported from mac-mini-S8's concurrent module-0 build, HYP-3865;
density lemma theirs, Norm engine ours — merged kind-pasteur-2026-07-02-S5) -/

/-- The comb: `v` arcs of half-width `r/v` at centers `(k + φ)/v`, `k = 0..v−1` — the danger
set of speed `v` with phase `φ`, on one period. -/
def comb (v : ℕ) (r φ : ℚ) : Region :=
  ((List.range v).map (Nat.cast : ℕ → ℚ)).map fun k => ((k + φ - r) / v, (k + φ + r) / v)

private theorem sum_map_range_const (n : ℕ) (f : ℕ → ℚ) (c : ℚ)
    (h : ∀ k < n, f k = c) : ((List.range n).map f).sum = n * c := by
  induction n with
  | zero => simp
  | succ m ih =>
      rw [List.range_succ, List.map_append, List.sum_append]
      simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil]
      rw [ih (fun k hk => h k (Nat.lt_succ_of_lt hk)), h m (Nat.lt_succ_self m)]
      push_cast
      ring

/-- Comb density: total length `2r`, independent of the speed and phase. -/
theorem length_comb (v : ℕ) (hv : 0 < v) (r φ : ℚ) (hr : 0 ≤ r) :
    length (comb v r φ) = 2 * r := by
  unfold length comb
  have hv2 : (0 : ℚ) < (v : ℚ) := by exact_mod_cast hv
  rw [List.map_map, List.map_map]
  refine (sum_map_range_const v _ (2 * r / v) ?_).trans ?_
  · intro k _
    simp only [Function.comp_apply]
    have hdiff : ((((k : ℕ) : ℚ) + φ + r) / v) - ((((k : ℕ) : ℚ) + φ - r) / v)
        = 2 * r / v := by
      field_simp; ring
    rw [hdiff, max_eq_right (by positivity)]
  · field_simp

end RatIntervals
end LonelyRunner

/-! Compatibility shim for consumers of the concurrent mac-mini namespace (HYP-3865). -/
namespace TournamentH7.RatIntervals

abbrev RI := LonelyRunner.RatIntervals.Region
abbrev mem := LonelyRunner.RatIntervals.mem
abbrev len := LonelyRunner.RatIntervals.length
abbrev inter := LonelyRunner.RatIntervals.inter
abbrev translate := LonelyRunner.RatIntervals.translate
abbrev comb := LonelyRunner.RatIntervals.comb

theorem mem_inter {x : ℚ} {A B : RI} : mem x (inter A B) ↔ mem x A ∧ mem x B :=
  LonelyRunner.RatIntervals.mem_inter
theorem len_nonneg (l : RI) : 0 ≤ len l := LonelyRunner.RatIntervals.length_nonneg l
theorem len_append (A B : RI) : len (A ++ B) = len A + len B :=
  LonelyRunner.RatIntervals.length_append A B
theorem len_translate (t : ℚ) (l : RI) : len (translate t l) = len l :=
  LonelyRunner.RatIntervals.length_translate t l
theorem len_comb (v : ℕ) (hv : 0 < v) (r φ : ℚ) (hr : 0 ≤ r) : len (comb v r φ) = 2 * r :=
  LonelyRunner.RatIntervals.length_comb v hv r φ hr

end TournamentH7.RatIntervals
