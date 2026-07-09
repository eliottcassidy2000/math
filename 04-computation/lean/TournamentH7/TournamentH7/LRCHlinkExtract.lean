/-
  TournamentH7.LRCHlinkExtract — the mergeSort argmax + wrapping-gap extraction for klein-S203's
  `hlink` (opus-2026-07-09-S175).

  kps-S101 (`LRCTeethGap`) discharged the NON-WRAPPING freeness core of `hlink`.  This file supplies
  the remaining two pieces the owner named:
    (1) the `foldl max` ARGMAX extraction — `maxCircGap = G > 0` is achieved by an ADJACENT pair of the
        sorted-and-closed residue cycle `cyc = ps ++ [p0+Vmax]`;
    (2) the WRAPPING-gap case — when the max gap is the wraparound `(ps.last, p0+Vmax)`, the arc past `1`
        is still free of every integer tooth-translate (uses that `ps.last` is the max and `p0` the min).

  Built bottom-up: the `foldl max` membership helpers first (self-contained, reusable), then the
  adjacent-pair extraction, then the two freeness branches, then `hlink`.
-/
import Mathlib
import TournamentH7.LRCTeethGap

namespace LRC14

open LonelyRunner

/-! ### 1. `foldl max` argmax helpers -/

/-- The `foldl max` of a `ℕ`-list is either the initial accumulator or a member of the list. -/
theorem foldl_max_eq_or_mem : ∀ (L : List ℕ) (a : ℕ), L.foldl max a = a ∨ L.foldl max a ∈ L := by
  intro L
  induction L with
  | nil => intro a; left; rfl
  | cons x xs ih =>
    intro a
    simp only [List.foldl_cons]
    rcases ih (max a x) with h | h
    · rw [h]
      rcases max_choice a x with hx | hx
      · left; exact hx
      · right; rw [hx]; exact List.mem_cons_self
    · right; exact List.mem_cons_of_mem x h

/-- A POSITIVE `foldl max 0` is achieved by a member of the list. -/
theorem foldl_max_mem (L : List ℕ) (h : 0 < L.foldl max 0) : L.foldl max 0 ∈ L := by
  rcases foldl_max_eq_or_mem L 0 with h0 | hm
  · rw [h0] at h; exact absurd h (lt_irrefl 0)
  · exact hm

/-- Every element `d` of `List.zipWith (fun a b => b - a) c c.tail` is `y − x` for an ADJACENT pair
`x, y` of `c`: `c = l₁ ++ x :: y :: l₂` with `y − x = d`. -/
theorem mem_zipWith_sub_tail : ∀ (c : List ℕ) (d : ℕ),
    d ∈ List.zipWith (fun a b => b - a) c c.tail →
    ∃ (l₁ l₂ : List ℕ) (x y : ℕ), c = l₁ ++ x :: y :: l₂ ∧ y - x = d := by
  intro c
  induction c with
  | nil => intro d hd; simp at hd
  | cons x rest ihx =>
    cases rest with
    | nil => intro d hd; simp at hd
    | cons y rest2 =>
      intro d hd
      simp only [List.tail_cons, List.zipWith_cons_cons, List.mem_cons] at hd
      rcases hd with h | h
      · exact ⟨[], rest2, x, y, by simp, h.symm⟩
      · have h' : d ∈ List.zipWith (fun a b => b - a) (y :: rest2) (y :: rest2).tail := by
          simpa using h
        obtain ⟨l₁, l₂, x', y', heq, hd'⟩ := ihx d h'
        exact ⟨x :: l₁, l₂, x', y', by rw [heq]; simp, hd'⟩

end LRC14
