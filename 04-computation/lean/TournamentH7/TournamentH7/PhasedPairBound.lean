/-
klein-2026-07-02-S103 — THM-604(b), the phased pair bound (the count argument, ℚ layer).

For coprime speeds P < Q and ANY phases φ, ψ, when `2r(P+Q) ≤ 1` the combs' intersection
has length at most `2r/Q` — the width of ONE deep interval. Mechanism: a positive clip
between the k-th P-interval and l-th Q-interval forces the integer `n = kQ − lP` into an
open interval of length `2r(P+Q) ≤ 1` (so `n` is unique), and coprimality makes
`(k,l) ↦ kQ − lP` injective on `range P × range Q` — hence AT MOST ONE positive clip.
This is the shift-uniform upper half of THM-604 (attainment = OriginNestQ.lean), and the
same residue-uniqueness engine as NestDecidable/LRCSevenTranslate, now with free phases.
-/
import Mathlib
import TournamentH7.RatIntervals

namespace PhasedPairBound

open LonelyRunner.RatIntervals

/-- Width-after-clip is at most the right factor's width. -/
theorem widthNN_clip_le_right (p q : ℚ × ℚ) :
    max 0 ((clip p q).2 - (clip p q).1) ≤ max 0 (q.2 - q.1) := by
  apply max_le_max (le_refl 0)
  have h1 : min p.2 q.2 ≤ q.2 := min_le_right _ _
  have h2 : q.1 ≤ max p.1 q.1 := le_max_right _ _
  unfold clip
  simp only
  linarith

/-- Exclusion sum: nonnegative entries, each ≤ c, pairwise at-most-one-nonzero ⟹ sum ≤ c. -/
theorem sum_le_of_pairwise_zero {L : List ℚ} {c : ℚ} (hc : 0 ≤ c)
    (hnn : ∀ x ∈ L, 0 ≤ x) (hle : ∀ x ∈ L, x ≤ c)
    (hpw : L.Pairwise (fun a b => a = 0 ∨ b = 0)) : L.sum ≤ c := by
  induction L with
  | nil => simpa using hc
  | cons a t ih =>
    rcases List.pairwise_cons.mp hpw with ⟨hhead, htail⟩
    rcases eq_or_lt_of_le (hnn a (by simp)) with ha | ha
    · -- a = 0
      simp only [List.sum_cons, ← ha, zero_add]
      exact ih (fun x hx => hnn x (by simp [hx])) (fun x hx => hle x (by simp [hx])) htail
    · -- a > 0 ⟹ tail all zero
      have hz : ∀ x ∈ t, x = 0 := fun x hx =>
        (hhead x hx).resolve_left (by linarith)
      have : t.sum = 0 := List.sum_eq_zero hz
      simp only [List.sum_cons, this, add_zero]
      exact hle a (by simp)

/-- The meeting inequality: the (k,l) clip is positive iff
`|(kQ − lP) + (φQ − ψP)| < r(P+Q)` (all in ℚ, P,Q > 0). -/
theorem clip_pos_iff {P Q : ℕ} (hP : 0 < (P : ℚ)) (hQ : 0 < (Q : ℚ))
    {r φ ψ : ℚ} {k l : ℚ} :
    (max ((k + φ - r) / P) ((l + ψ - r) / Q) < min ((k + φ + r) / P) ((l + ψ + r) / Q))
    ↔ |k * Q - l * P + (φ * Q - ψ * P)| < r * (P + Q) := by
  simp only [max_lt_iff, lt_min_iff, abs_lt, div_lt_div_iff₀ hP hQ,
    div_lt_div_iff₀ hQ hP, div_lt_div_iff₀ hP hP, div_lt_div_iff₀ hQ hQ]
  constructor
  · rintro ⟨⟨-, h1⟩, h2, -⟩
    constructor <;> nlinarith
  · rintro ⟨h1, h2⟩
    refine ⟨⟨?_, ?_⟩, ?_, ?_⟩ <;> nlinarith

/-- Positive-sum extraction over nonnegative lists. -/
theorem exists_pos_of_sum_pos {L : List ℚ} (hnn : ∀ x ∈ L, 0 ≤ x) (h : 0 < L.sum) :
    ∃ x ∈ L, 0 < x := by
  by_contra hcon
  push_neg at hcon
  have : L.sum ≤ 0 := by
    have hz : ∀ x ∈ L, x = 0 := fun x hx => le_antisymm (hcon x hx) (hnn x hx)
    simp [List.sum_eq_zero hz]
  linarith

/-- **THM-604(b), the phased pair bound**: for coprime `P, Q` and ANY phases, when
`2r(P+Q) ≤ 1` the comb intersection is at most one deep interval: `≤ 2r/Q`. -/
theorem phased_pair_bound {P Q : ℕ} (hcop : Nat.Coprime P Q)
    (hP : 0 < P) (hQ : 0 < Q) {r φ ψ : ℚ} (hr : 0 ≤ r)
    (hr2 : 2 * r * ((P : ℚ) + Q) ≤ 1) :
    length (inter (comb P r φ) (comb Q r ψ)) ≤ 2 * r / Q := by
  have hPQ : (0 : ℚ) < P := by exact_mod_cast hP
  have hQQ : (0 : ℚ) < Q := by exact_mod_cast hQ
  have hc : (0 : ℚ) ≤ 2 * r / Q := by positivity
  -- normalize the structure: nested range maps
  unfold comb inter length
  simp only [List.map_map, List.flatMap_map, List.map_flatMap]
  rw [List.sum_flatMap]
  -- outer exclusion
  rw [List.map_map]
  apply sum_le_of_pairwise_zero hc
  · intro x hx
    simp only [List.mem_map] at hx
    obtain ⟨k, -, rfl⟩ := hx
    apply List.sum_nonneg
    intro y hy
    simp only [List.map_map, List.mem_map] at hy
    obtain ⟨l, -, rfl⟩ := hy
    exact le_max_left 0 _
  · -- each inner sum ≤ 2r/Q
    intro x hx
    simp only [List.mem_map] at hx
    obtain ⟨k, hkmem, rfl⟩ := hx
    simp only [List.mem_range] at hkmem
    apply sum_le_of_pairwise_zero hc
    · intro y hy
      simp only [List.map_map, List.mem_map] at hy
      obtain ⟨l, -, rfl⟩ := hy
      exact le_max_left 0 _
    · intro y hy
      simp only [List.map_map, List.mem_map] at hy
      obtain ⟨l, -, rfl⟩ := hy
      refine le_trans (widthNN_clip_le_right _ _) ?_
      have hw : ((l : ℚ) + ψ + r) / Q - ((l : ℚ) + ψ - r) / Q = 2 * r / Q := by
        field_simp
        ring
      simp only [Function.comp]
      rw [hw]
      exact max_le hc (le_refl _)
    · -- inner pairwise: two distinct l cannot both be positive
      rw [List.map_map, List.pairwise_map]
      have hpr : (List.range Q).Pairwise (· < ·) := List.pairwise_lt_range
      refine hpr.imp_of_mem ?_
      intro l l' hlm hl'm hll'
      simp only [List.mem_range] at hlm hl'm
      by_contra hcon
      push_neg at hcon
      obtain ⟨h1, h2⟩ := hcon
      have hp1 : 0 < max 0 ((clip (((k:ℚ)+φ-r)/P, ((k:ℚ)+φ+r)/P)
          (((l:ℚ)+ψ-r)/Q, ((l:ℚ)+ψ+r)/Q)).2 - _.1) := lt_of_le_of_ne (le_max_left _ _) (Ne.symm h1)
      sorry
  · sorry


end PhasedPairBound
