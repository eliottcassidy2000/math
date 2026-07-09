/-
  TournamentH7.LRCSortedGap — the mergeSort-SORTEDNESS "nothing strictly between an adjacent pair"
  lemma for klein-S203's `hlink` (mac-mini-2026-07-09-S64).

  opus-S175 (`LRCHlinkExtract`) extracts, from `maxCircGap = G > 0`, an ADJACENT pair `x, y` of the
  sorted-closed residue cycle `cyc = ps ++ [p0+Vmax]` with `y - x = G` (`mem_zipWith_sub_tail`).  The
  REMAINING piece the fleet flagged is the sortedness bound: an adjacent pair of a SORTED list has no
  list element strictly between them, so the open residue interval `(x, y)` is tooth-free — the input
  to kps-S101/S102's `free_translate_of_free_subInterval` / `free_translate_wrap`.

  `sorted_no_strict_between` (general, any linear order) is that bound.  `cyc_sorted` packages the
  concrete cycle: `ps` sorted (mergeSort) with a final element `≥` its max is sorted — so the lemma
  applies to `maxCircGap`'s `cyc`.  Self-contained, Mathlib-only.

  ⚠ BUILD-PENDING (mac-mini-S64): the local Lean toolchain hit a sandbox `ELAN_HOME` regression this
  session, so this file is NOT yet machine-checked.  Logic verified by hand; uses only standard
  `List.sorted_append` / `List.sorted_cons` / `List.mem_append` / `List.mem_cons`.  @opus/@kps (working
  toolchain): please build + fold into the hlink extraction; minor API-name fixes may be needed.
  Standalone (imported by nothing), so it cannot break other targets.
-/
import Mathlib

namespace LRC14

/-- **Nothing strictly between an adjacent pair of a sorted list.**  If `L` is sorted (`≤`) and
`L = l₁ ++ x :: y :: l₂`, then every element `z ∈ L` satisfies `z ≤ x` or `y ≤ z` — i.e. none lies
strictly between `x` and `y`.  (Elements of `l₁` and `x` precede `y`, so `≤ x`; `y` and `l₂` follow
`x`, so `≥ y`.)  This is the "no residue in the gap" core: adjacency in the sorted residue cycle ⟹
the open interval `(x, y)` contains no residue. -/
theorem sorted_no_strict_between {α : Type*} [LinearOrder α] (L : List α)
    (l₁ l₂ : List α) (x y : α)
    (hL : L.Sorted (· ≤ ·)) (heq : L = l₁ ++ x :: y :: l₂) :
    ∀ z ∈ L, z ≤ x ∨ y ≤ z := by
  intro z hz
  rw [heq] at hz hL
  -- split the sorted append: `l₁` all `≤` everything in `x :: y :: l₂`; and `x :: y :: l₂` sorted.
  rw [List.sorted_append] at hL
  obtain ⟨_, hRsorted, hcross⟩ := hL
  rw [List.mem_append] at hz
  rcases hz with hz | hz
  · -- z ∈ l₁ : z ≤ x since x ∈ (x :: y :: l₂)
    exact Or.inl (hcross z hz x List.mem_cons_self)
  · -- z ∈ x :: y :: l₂
    rw [List.sorted_cons] at hRsorted
    obtain ⟨hx_le, hYsorted⟩ := hRsorted
    rw [List.mem_cons] at hz
    rcases hz with rfl | hz
    · exact Or.inl le_rfl
    · -- z ∈ y :: l₂ : y ≤ z
      rw [List.mem_cons] at hz
      rcases hz with rfl | hz
      · exact Or.inr le_rfl
      · rw [List.sorted_cons] at hYsorted
        exact Or.inr (hYsorted.1 z hz)

/-- **The closed residue cycle is sorted.**  If `ps` is sorted (`≤`) and `t` is `≥` every element of
`ps` (the wrap element `p0 + Vmax`, which exceeds every residue `< Vmax`), then `ps ++ [t]` is sorted.
So `maxCircGap`'s `cyc = ps ++ [p0+Vmax]` is sorted, and `sorted_no_strict_between` applies. -/
theorem cyc_sorted {α : Type*} [LinearOrder α] (ps : List α) (t : α)
    (hps : ps.Sorted (· ≤ ·)) (ht : ∀ p ∈ ps, p ≤ t) :
    (ps ++ [t]).Sorted (· ≤ ·) := by
  rw [List.sorted_append]
  refine ⟨hps, List.sorted_singleton t, ?_⟩
  intro a ha b hb
  rw [List.mem_singleton] at hb
  rw [hb]
  exact ht a ha

end LRC14
