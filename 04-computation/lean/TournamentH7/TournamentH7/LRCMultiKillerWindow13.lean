/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-05-S79)
-/
import TournamentH7.LRCMultiKillerWindow

/-!
# The 13-band tooth geometry and sub-gap sharpening (multi-far transcription, installment 1)

The l ≥ 7 lift-rigidity closure (opus-S78/HYP-4107) consumes klein's multi-far window at band
`1/13` with base margin `β = 1/6`.  klein's S135 machinery is 14-band (`tooth` half-width
`1/14` is baked into `LRCBlockSix`).  This file transcribes the BAND-DEPENDENT GEOMETRY to
`1/13` — mirroring S135 line-by-line per the S76 lore (transcribe, don't reinvent):

  * `tooth13`/`teeth13` — half-width `1/13`, tooth width `2/(13w)`, inter-tooth gap `11/(13w)`;
  * `clipLen_tooth13_le` — each clipped tooth carries at most `2/(13w)`;
  * `tooth13_clip_pos_edges` / `tooth13_clip_disjoint` — a window with `13(b−a)w ≤ 11`
    (sub-gap) meets at most one tooth;
  * `teeth13_mass_far` — the sub-gap single-tooth mass bound, via klein's band-free
    `sum_le_single_bound`.

Installment 2 (mapped, next session): the 13-band hunter block step + `lonely13_of_window_of_fees`
/ `_multi_far`, closing the l ≥ 7 consumer end-to-end.  `clipLen`, `rlength`, `rinter`,
`pairCredits` are band-free and import as-is.
-/

namespace LonelyRunner
namespace MultiKillerWindow13

open LRC14 MultiKillerWindow

/-- One 13-band tooth of runner `w` at integer `m`: the danger zone `((m ± 1/13)/w)`. -/
noncomputable def tooth13 (w m : ℤ) : ℝ × ℝ :=
  (((m : ℝ) - 1 / 13) / w, ((m : ℝ) + 1 / 13) / w)

/-- The 13-band teeth of runner `w` that can meet the window `[a, b]`. -/
noncomputable def teeth13 (w : ℤ) (a b : ℝ) : List (ℝ × ℝ) :=
  (List.range (⌊(w : ℝ) * b⌋ + 1 - (⌈(w : ℝ) * a⌉ - 1) + 1).toNat).map fun (i : ℕ) =>
    tooth13 w (⌈(w : ℝ) * a⌉ - 1 + (i : ℤ))

/-- Each clipped 13-band tooth carries at most its width `2/(13w)`. -/
theorem clipLen_tooth13_le {w : ℤ} (hw : 0 < w) (m : ℤ) (a b : ℝ) :
    clipLen (tooth13 w m) a b ≤ 2 / (13 * (w : ℝ)) := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  unfold clipLen tooth13
  simp only
  rcases le_total (min (((m : ℝ) + 1 / 13) / w) b - max (((m : ℝ) - 1 / 13) / w) a) 0
    with h | h
  · rw [max_eq_left h]
    positivity
  · rw [max_eq_right h]
    have h1 : min (((m : ℝ) + 1 / 13) / w) b ≤ ((m : ℝ) + 1 / 13) / w := min_le_left _ _
    have h2 : ((m : ℝ) - 1 / 13) / w ≤ max (((m : ℝ) - 1 / 13) / w) a := le_max_left _ _
    have h3 : ((m : ℝ) + 1 / 13) / w - ((m : ℝ) - 1 / 13) / w = 2 / (13 * (w : ℝ)) := by
      field_simp
      ring
    linarith

/-- Positive clip forces the 13-band tooth's edges to straddle the window. -/
theorem tooth13_clip_pos_edges {w : ℤ} (hw : 0 < w) {a b : ℝ} {n : ℤ}
    (hn : 0 < clipLen (tooth13 w n) a b) :
    ((n : ℝ) - 1/13) / w < b ∧ a < ((n : ℝ) + 1/13) / w := by
  unfold clipLen tooth13 at hn
  simp only at hn
  have hn' : 0 < min (((n : ℝ) + 1/13) / w) b - max (((n : ℝ) - 1/13) / w) a := by
    by_contra hle
    push Not at hle
    rw [max_eq_left hle] at hn
    exact lt_irrefl 0 hn
  constructor
  · calc ((n : ℝ) - 1/13) / w ≤ max (((n : ℝ) - 1/13) / w) a := le_max_left _ _
      _ < min (((n : ℝ) + 1/13) / w) b := by linarith
      _ ≤ b := min_le_right _ _
  · calc a ≤ max (((n : ℝ) - 1/13) / w) a := le_max_right _ _
      _ < min (((n : ℝ) + 1/13) / w) b := by linarith
      _ ≤ ((n : ℝ) + 1/13) / w := min_le_left _ _

/-- Two different 13-band teeth cannot both clip positively into a window shorter than the
inter-tooth gap `11/(13w)`. -/
theorem tooth13_clip_disjoint {w : ℤ} (hw : 0 < w) {a b : ℝ}
    (hfar : 13 * (b - a) * (w : ℝ) ≤ 11) {m m' : ℤ} (hlt : m < m')
    (hpos : 0 < clipLen (tooth13 w m) a b) (hpos' : 0 < clipLen (tooth13 w m') a b) :
    False := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  obtain ⟨-, hra⟩ := tooth13_clip_pos_edges hw hpos
  obtain ⟨hlb, -⟩ := tooth13_clip_pos_edges hw hpos'
  rw [lt_div_iff₀ hwR] at hra
  rw [div_lt_iff₀ hwR] at hlb
  have hm1 : (m : ℝ) + 1 ≤ (m' : ℝ) := by exact_mod_cast hlt
  have hexp : 13 * (b - a) * (w : ℝ) = 13 * (b * w) - 13 * (a * w) := by ring
  rw [hexp] at hfar
  linarith

/-- **`teeth13_mass_far`**: a sub-gap window (`13(b−a)w ≤ 11`) carries at most one 13-band
tooth of runner `w` — clipped mass `≤ 2/(13w)`, no density term, no boundary fee. -/
theorem teeth13_mass_far {w : ℤ} (hw : 0 < w) (a b : ℝ)
    (hfar : 13 * (b - a) * (w : ℝ) ≤ 11) :
    ((teeth13 w a b).map fun p => clipLen p a b).sum ≤ 2 / (13 * (w : ℝ)) := by
  have hwR : (0 : ℝ) < (w : ℝ) := by exact_mod_cast hw
  apply sum_le_single_bound
  · positivity
  · intro x hx
    rw [List.mem_map] at hx
    obtain ⟨p, -, rfl⟩ := hx
    exact clipLen_nonneg p a b
  · intro x hx
    rw [List.mem_map] at hx
    obtain ⟨p, hp, rfl⟩ := hx
    unfold teeth13 at hp
    rw [List.mem_map] at hp
    obtain ⟨i, -, rfl⟩ := hp
    exact clipLen_tooth13_le hw _ a b
  · unfold teeth13
    rw [List.map_map]
    refine List.Pairwise.map _ (fun i j hij => ?_) List.pairwise_lt_range
    by_contra hcon
    push Not at hcon
    obtain ⟨h1, h2⟩ := hcon
    have hp1 : 0 < clipLen (tooth13 w (⌈(w : ℝ) * a⌉ - 1 + (i : ℤ))) a b :=
      lt_of_le_of_ne (clipLen_nonneg _ _ _) (Ne.symm h1)
    have hp2 : 0 < clipLen (tooth13 w (⌈(w : ℝ) * a⌉ - 1 + (j : ℤ))) a b :=
      lt_of_le_of_ne (clipLen_nonneg _ _ _) (Ne.symm h2)
    exact tooth13_clip_disjoint hw hfar (by omega) hp1 hp2

#print axioms teeth13_mass_far

end MultiKillerWindow13
end LonelyRunner
