/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-02-S38)
-/
import TournamentH7.RatIntervals

/-!
# RatIntervals: the wraparound (circle) layer — module 0's remaining spec, part 1

`RatIntervals.Region` values live on ℝ-the-line; combs and translates overhang `[0,1)`.  This
file adds the circle structure: `wrapOne`/`wrap` normalize a region into `[0,1)` (translating
each interval by an integer and splitting at `1`), and `translateCirc` is the circle translation.
Lengths are preserved (`length_wrap`, `length_translateCirc`) for regions of interval width ≤ 1
(every comb tooth, window, and clip in the project); membership gets the circle semantics
(`mem_wrap`: `x ∈ wrap L ↔ some integer translate of x lies in L`).

Design note (single-writer discipline): this is a NEW file importing the canonical module 0,
not an edit of it — after this morning's union-merge collision, Lean modules should have one
writer per session; extensions go in satellite files until the owner merges them.

Why this layer makes module 3 (Commensuration-ℚ) cleaner than its ℝ original: with half-open
intervals the seven `1/7`-translates of a `7∤P` comb tile `[0,1)` EXACTLY — no a.e., no null
sets, no Haar measure. The measure theory in the ℝ proof was an artifact of closed balls.
-/

namespace LonelyRunner
namespace RatIntervals

/-- Wrap one interval into `[0,1)`: translate so the left endpoint lands in `[0,1)`, split at
`1` if the right end overhangs. (For width > 1 the result under-covers; all project regions
have tooth width ≤ 1.) -/
def wrapOne (p : ℚ × ℚ) : Region :=
  let k : ℤ := ⌊p.1⌋
  let a := p.1 - k
  let b := p.2 - k
  if b ≤ 1 then [(a, b)] else [(a, 1), (0, b - 1)]

/-- Wrap a region into `[0,1)`. -/
def wrap (L : Region) : Region := L.flatMap wrapOne

/-- Circle translation: translate on the line, then wrap. -/
def translateCirc (t : ℚ) (L : Region) : Region := wrap (translate t L)

theorem length_wrapOne {p : ℚ × ℚ} (hw : p.2 - p.1 ≤ 1) :
    length (wrapOne p) = max 0 (p.2 - p.1) := by
  unfold wrapOne
  set k : ℤ := ⌊p.1⌋ with hk
  have h1 : (k : ℚ) ≤ p.1 := Int.floor_le p.1
  have h2 : p.1 < (k : ℚ) + 1 := Int.lt_floor_add_one p.1
  by_cases hb : p.2 - (k : ℚ) ≤ 1
  · simp only [hb, if_true, length, List.map_cons, List.map_nil, List.sum_cons,
      List.sum_nil, add_zero]
    congr 1
    ring
  · simp only [hb, if_false, length]
    simp [List.map, List.sum_cons]
    have hb' : 1 < p.2 - (k : ℚ) := lt_of_not_ge hb
    have hpos1 : 0 ≤ 1 - (p.1 - (k : ℚ)) := by linarith
    have hpos2 : 0 ≤ p.2 - (k : ℚ) - 1 := by linarith
    have hw' : 0 ≤ p.2 - p.1 := by linarith
    rw [max_eq_right hpos1, max_eq_right hpos2, max_eq_right hw']
    ring

theorem length_wrap {L : Region} (hw : ∀ p ∈ L, p.2 - p.1 ≤ 1) :
    length (wrap L) = length L := by
  induction L with
  | nil => rfl
  | cons p L ih =>
    have hp : p.2 - p.1 ≤ 1 := hw p (List.mem_cons_self ..)
    have hL : ∀ q ∈ L, q.2 - q.1 ≤ 1 := fun q hq => hw q (List.mem_cons_of_mem _ hq)
    show length (wrapOne p ++ wrap L) = length (p :: L)
    rw [length_append, ih hL, length_wrapOne hp]
    simp [length, List.map, List.sum_cons]

theorem length_translateCirc {t : ℚ} {L : Region} (hw : ∀ p ∈ L, p.2 - p.1 ≤ 1) :
    length (translateCirc t L) = length L := by
  unfold translateCirc
  rw [length_wrap, length_translate]
  intro q hq
  simp only [translate, List.mem_map] at hq
  obtain ⟨p, hp, rfl⟩ := hq
  simpa using hw p hp

/-- Circle membership through one wrapped interval: `x ∈ wrapOne p` iff some integer translate
of `x` lies in `p`.  (`x` in the fundamental domain `[0,1)`.) -/
theorem mem_wrapOne {x : ℚ} (hx0 : 0 ≤ x) (hx1 : x < 1) {p : ℚ × ℚ} (hw : p.2 - p.1 ≤ 1) :
    mem x (wrapOne p) ↔ ∃ n : ℤ, p.1 ≤ x + n ∧ x + n < p.2 := by
  unfold wrapOne
  set k : ℤ := ⌊p.1⌋ with hk
  have h1 : (k : ℚ) ≤ p.1 := Int.floor_le p.1
  have h2 : p.1 < (k : ℚ) + 1 := Int.lt_floor_add_one p.1
  by_cases hb : p.2 - (k : ℚ) ≤ 1
  · simp only [hb, if_true]
    constructor
    · rintro ⟨q, hq, hqx⟩
      rcases List.mem_singleton.mp hq with rfl
      exact ⟨k, by constructor <;> simp only at hqx <;> linarith [hqx.1, hqx.2]⟩
    · rintro ⟨n, hn1, hn2⟩
      refine ⟨_, List.mem_singleton_self _, ?_, ?_⟩ <;> simp only
      · have hnk : n = k := by
          have hlow : (k : ℚ) - 1 < (n : ℚ) := by linarith
          have hhigh : (n : ℚ) < (k : ℚ) + 1 := by linarith
          have hl : (k : ℤ) - 1 < n := by exact_mod_cast hlow
          have hh : n < (k : ℤ) + 1 := by exact_mod_cast hhigh
          omega
        subst hnk
        linarith
      · have hnk : n = k := by
          have hlow : (k : ℚ) - 1 < (n : ℚ) := by linarith
          have hhigh : (n : ℚ) < (k : ℚ) + 1 := by linarith
          have hl : (k : ℤ) - 1 < n := by exact_mod_cast hlow
          have hh : n < (k : ℤ) + 1 := by exact_mod_cast hhigh
          omega
        subst hnk
        linarith
  · simp only [hb, if_false]
    have hb' : 1 < p.2 - (k : ℚ) := lt_of_not_ge hb
    constructor
    · rintro ⟨q, hq, hqx⟩
      rcases List.mem_cons.mp hq with rfl | hq'
      · exact ⟨k, by simp only at hqx; constructor <;> linarith [hqx.1, hqx.2, hb']⟩
      · rcases List.mem_singleton.mp hq' with rfl
        exact ⟨k + 1, by simp only at hqx; push_cast; constructor <;> linarith [hqx.1, hqx.2]⟩
    · rintro ⟨n, hn1, hn2⟩
      -- n is k or k+1, decided by whether x + n < k + 1
      have hnrange : n = k ∨ n = k + 1 := by
        have hlow : (k : ℚ) - 1 < (n : ℚ) := by linarith
        have hhigh : (n : ℚ) < (k : ℚ) + 2 := by
          have : x + (n : ℚ) < p.2 := hn2
          have : p.2 ≤ p.1 + 1 := by linarith
          linarith
        have hl : (k : ℤ) - 1 < n := by exact_mod_cast hlow
        have hh : n < (k : ℤ) + 2 := by exact_mod_cast hhigh
        omega
      rcases hnrange with rfl | rfl
      · by_cases hcase : x + (k : ℚ) < (k : ℚ) + 1
        · exact ⟨_, List.mem_cons_self .., by simp only; constructor <;> linarith⟩
        · -- x ≥ 1, contradiction
          exact absurd hx1 (by linarith)
      · refine ⟨_, List.mem_cons_of_mem _ (List.mem_singleton_self _), ?_, ?_⟩ <;>
          simp only <;> push_cast at hn1 hn2 ⊢
        · linarith
        · linarith

/-- Circle membership through a wrapped region. -/
theorem mem_wrap {x : ℚ} (hx0 : 0 ≤ x) (hx1 : x < 1) {L : Region}
    (hw : ∀ p ∈ L, p.2 - p.1 ≤ 1) :
    mem x (wrap L) ↔ ∃ n : ℤ, mem (x + n) L := by
  induction L with
  | nil => simp [wrap, mem]
  | cons p L ih =>
    have hp : p.2 - p.1 ≤ 1 := hw p (List.mem_cons_self ..)
    have hL : ∀ q ∈ L, q.2 - q.1 ≤ 1 := fun q hq => hw q (List.mem_cons_of_mem _ hq)
    show mem x (wrapOne p ++ wrap L) ↔ _
    rw [mem_union, mem_wrapOne hx0 hx1 hp, ih hL]
    constructor
    · rintro (⟨n, hn⟩ | ⟨n, q, hq, hxq⟩)
      · exact ⟨n, p, List.mem_cons_self .., hn⟩
      · exact ⟨n, q, List.mem_cons_of_mem _ hq, hxq⟩
    · rintro ⟨n, q, hq, hxq⟩
      rcases List.mem_cons.mp hq with rfl | hq'
      · exact Or.inl ⟨n, hxq⟩
      · exact Or.inr ⟨n, q, hq', hxq⟩

end RatIntervals
end LonelyRunner
