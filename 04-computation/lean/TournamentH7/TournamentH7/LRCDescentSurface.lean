/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-05-S93)
-/
import TournamentH7.LRCTowerLift
import TournamentH7.LRCGapDescent

/-!
# The descent surface: the l ≥ 7 leg's named formal consumer (HYP-4216)

The crutch-removal architecture (S91): the assembly supplies a CLEAR COMPONENT — an
interval on which every base runner is strictly `ρ`-far pointwise (certified later by
interval-arithmetic tables) — and a spread chain of top speeds entering at `2 ≤ w·L`.
`strict_lonely_of_clear_component` returns one point that is STRICTLY good for base and
tops simultaneously: the strict-side input hdich needs (a strictly-better-than-1/13 point
beats a tight family). Builds on `gap_interior_strict` (strictness is free in the nested
tower) and mirrors `spread_tower`'s induction in strict form.
-/

namespace LonelyRunner
namespace DescentSurface

open TowerLift

/-- **The strict spread tower**: like `GapDescent.spread_tower`, with interior points —
every dodge is strict. -/
theorem spread_tower_strict (ρ : ℝ) (hρ0 : 0 < ρ) (hρ : ρ < 1/2) (ws : List ℝ) :
    (∀ w ∈ ws, 0 < w) →
    List.IsChain (fun x y => 2 * x ≤ (1 - 2*ρ) * y) ws →
    ∀ a L : ℝ, 0 < L → (∀ w ∈ ws.head?, 2 ≤ w * L) →
    ∃ t, a < t ∧ t < a + L ∧ ∀ w ∈ ws, ∀ m : ℤ, ρ < |w * t - m| := by
  induction ws with
  | nil =>
    intro _ _ a L hL _
    exact ⟨a + L/2, by linarith, by linarith, by simp⟩
  | cons w rest ih =>
    intro hpos hchain a L hL hentry
    have hw : 0 < w := hpos w (by simp)
    have hwL : 2 ≤ w * L := hentry w (by simp)
    rw [List.isChain_cons] at hchain
    obtain ⟨hhd, hrest⟩ := hchain
    obtain ⟨c, d, hac, hda, hcd, hdc, hgood⟩ := gap_interior_strict ρ w a L hρ0 hρ hw hwL
    have hL' : (0:ℝ) < d - c := by linarith
    obtain ⟨t, hct, htd, hgoodrest⟩ :=
      ih (fun w' hw' => hpos w' (List.mem_cons_of_mem _ hw')) hrest c (d - c) hL'
        (fun w' hw' => by
          have h2 := hhd w' hw'
          have hrw : w' * (d - c) = ((1 - 2*ρ) * w') / w := by
            rw [hdc]; ring
          rw [hrw, le_div_iff₀ hw]
          linarith)
    have htd' : t < d := by linarith
    refine ⟨t, by linarith, by linarith, ?_⟩
    intro w' hw' m
    rcases List.mem_cons.mp hw' with rfl | hmem
    · exact hgood t hct htd' m
    · exact hgoodrest w' hmem m

/-- **The descent surface**: a clear component (pointwise-strict base interval) plus a
spread chain of tops yields one point strictly good for base ∪ tops — the l ≥ 7 leg's
formal consumer under the S91 architecture. -/
theorem strict_lonely_of_clear_component (ρ : ℝ) (hρ0 : 0 < ρ) (hρ : ρ < 1/2)
    (base tops : List ℝ) (a L : ℝ) (hL : 0 < L)
    (hbase : ∀ t, a < t → t < a + L → ∀ w ∈ base, ∀ m : ℤ, ρ < |w * t - m|)
    (hpos : ∀ w ∈ tops, 0 < w)
    (hchain : List.IsChain (fun x y => 2 * x ≤ (1 - 2*ρ) * y) tops)
    (hentry : ∀ w ∈ tops.head?, 2 ≤ w * L) :
    ∃ t, (∀ w ∈ base, ∀ m : ℤ, ρ < |w * t - m|) ∧
         (∀ w ∈ tops, ∀ m : ℤ, ρ < |w * t - m|) := by
  obtain ⟨t, hat, hta, htops⟩ :=
    spread_tower_strict ρ hρ0 hρ tops hpos hchain a L hL hentry
  exact ⟨t, fun w hw m => hbase t hat hta w hw m, htops⟩

#print axioms spread_tower_strict
#print axioms strict_lonely_of_clear_component

end DescentSurface
end LonelyRunner
