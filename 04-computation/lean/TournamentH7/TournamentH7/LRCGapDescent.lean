/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-05-S81)
-/
import Mathlib

/-!
# Gap descent: the measure-free window peel (HYP-4116; corrects HYP-4107 leg 3)

## The six-top ceiling (why this file exists)

The fee criterion of `tower_step_12` (LRCTeethR) is **unsatisfiable for `|tops| ≥ 7`**
at `ρ = 1/13`, for EVERY speed configuration: a valid fee must dominate the teeth mass at
every window placement (`hfee` quantifies over all `tstar`), hence dominates its positional
MEAN, which is exactly `2ρ·L` per top (teeth of width `2ρ/w` recur with period `1/w`;
average the mass over one period of placements — Fubini).  So
`Σ fees ≥ 2ρ·l·L ≥ (14/13)·L > L` at `l = 7`, while `hcrit` demands `Σ fees < L`.
Independent of `B`, of the margin, and of the window size.  The `(13 − 2l)` denominators
in mac-mini's fee table `T_l` are therefore not an artifact of the mass bound — no
sharpened fee can cross `l = 6.5`.  (MISTAKE-105: opus-S78's HYP-4107 leg (3) claimed the
l ≥ 7 lift stratum closes through this window; it cannot.  The numeric floor 3/19 stands;
the route is replaced by THIS file plus a finite grid stratum.)

## The fix: nesting instead of accounting

`exists_gap_subinterval`: an interval `J` with `w·|J| ≥ 2` contains a **full inter-tooth
gap** — a closed subinterval of length `(1−2ρ)/w` on which `dist(w·t, ℤ) ≥ ρ` pointwise.
No measure, no fees: interval arithmetic only.

`spread_tower`: iterate the gap lemma along a list of tops whose consecutive ratios are
`≥ 2/(1−2ρ)` (`= 26/11` at `ρ = 1/13`), entry condition `2 ≤ w₁·L`.  ANY number of spread
tops is dodged inside the base window — the ceiling does not apply because the surviving
set is chosen adaptively per top (nesting), not bounded uniformly (fees).

## The corrected l ≥ 7 architecture (verified numerically, S81)

* spread tops (`w ≥ 2/L ≈ 134`, ratios `≥ 26/11`): THIS file;
* invisible lifts (`k ≡ 0 mod 13` ⟹ `v ≥ 170 > 134`): automatically descent-eligible;
* the bottom/visible cluster (which kills naive descent — tops `14..20` have teeth wider
  than the surviving window): strict `a/169` witnesses exist in 12/12 probed patterns
  with room `≥ 19/169` (`six_top_ceiling_gap_descent_opus_S81.out`, Part 3) — the S77
  kernel-row machinery one level up.  Grid rows are a finite stratum over
  `(r, k mod 13)` patterns.
-/

namespace LonelyRunner
namespace GapDescent

/-- **The gap lemma**: any interval `[a, a+L]` with `2 ≤ w·L` contains a full inter-tooth
gap of the speed-`w` comb — a closed subinterval of length `(1−2ρ)/w` on which every point
is `ρ`-far from every integer multiple.  Constructive: the gap after tooth `⌈w·a + ρ⌉`. -/
theorem exists_gap_subinterval (ρ w a L : ℝ) (hρ0 : 0 < ρ) (_hρ : ρ ≤ 1/2)
    (hw : 0 < w) (hL : 2 ≤ w * L) :
    ∃ c d : ℝ, a ≤ c ∧ d ≤ a + L ∧ d - c = (1 - 2*ρ)/w ∧
      ∀ t, c ≤ t → t ≤ d → ∀ m : ℤ, ρ ≤ |w * t - m| := by
  have hk1 : w * a + ρ ≤ ((⌈w * a + ρ⌉ : ℤ) : ℝ) := Int.le_ceil _
  have hk2 : ((⌈w * a + ρ⌉ : ℤ) : ℝ) < w * a + ρ + 1 := Int.ceil_lt_add_one _
  refine ⟨(((⌈w * a + ρ⌉ : ℤ) : ℝ) + ρ) / w, (((⌈w * a + ρ⌉ : ℤ) : ℝ) + 1 - ρ) / w,
    ?_, ?_, ?_, ?_⟩
  · rw [le_div_iff₀ hw]
    have hmul : a * w = w * a := mul_comm a w
    linarith
  · rw [div_le_iff₀ hw]
    have hmul : (a + L) * w = w * a + w * L := by ring
    linarith
  · rw [div_sub_div_same]
    congr 1
    ring
  · intro t htc htd m
    rw [div_le_iff₀ hw] at htc
    rw [le_div_iff₀ hw] at htd
    have hcomm : w * t = t * w := mul_comm w t
    by_cases hm : m ≤ ⌈w * a + ρ⌉
    · have hmR : (m : ℝ) ≤ ((⌈w * a + ρ⌉ : ℤ) : ℝ) := by exact_mod_cast hm
      exact le_abs.mpr (Or.inl (by linarith))
    · have hm1 : (⌈w * a + ρ⌉ : ℤ) + 1 ≤ m := by omega
      have hmR : ((⌈w * a + ρ⌉ : ℤ) : ℝ) + 1 ≤ (m : ℝ) := by exact_mod_cast hm1
      exact le_abs.mpr (Or.inr (by linarith))

/-- **The spread tower**: tops listed with consecutive ratios `≥ 2/(1−2ρ)` and entry
condition `2 ≤ w₁·L` are ALL dodged at a single point of the base window — nesting the
gap lemma, one top per level.  No fee accounting, hence no six-top ceiling: the list may
have any length. -/
theorem spread_tower (ρ : ℝ) (hρ0 : 0 < ρ) (hρ : ρ ≤ 1/2) (ws : List ℝ) :
    (∀ w ∈ ws, 0 < w) →
    List.IsChain (fun x y => 2 * x ≤ (1 - 2*ρ) * y) ws →
    ∀ a L : ℝ, 0 ≤ L → (∀ w ∈ ws.head?, 2 ≤ w * L) →
    ∃ t, a ≤ t ∧ t ≤ a + L ∧ ∀ w ∈ ws, ∀ m : ℤ, ρ ≤ |w * t - m| := by
  induction ws with
  | nil =>
    intro _ _ a L hL _
    exact ⟨a, le_refl a, by linarith, by simp⟩
  | cons w rest ih =>
    intro hpos hchain a L _ hentry
    have hw : 0 < w := hpos w (by simp)
    have hwL : 2 ≤ w * L := hentry w (by simp)
    rw [List.isChain_cons] at hchain
    obtain ⟨hhd, hrest⟩ := hchain
    obtain ⟨c, d, hac, hda, hdc, hgood⟩ := exists_gap_subinterval ρ w a L hρ0 hρ hw hwL
    have hL' : (0:ℝ) ≤ d - c := by
      rw [hdc]
      exact div_nonneg (by linarith) hw.le
    obtain ⟨t, hct, htd, hgoodrest⟩ :=
      ih (fun w' hw' => hpos w' (List.mem_cons_of_mem _ hw')) hrest c (d - c) hL'
        (fun w' hw' => by
          have h2 := hhd w' hw'
          have hrw : w' * (d - c) = ((1 - 2*ρ) * w') / w := by
            rw [hdc]; ring
          rw [hrw, le_div_iff₀ hw]
          linarith)
    have htd' : t ≤ d := by linarith
    refine ⟨t, le_trans hac hct, le_trans htd' hda, ?_⟩
    intro w' hw' m
    rcases List.mem_cons.mp hw' with rfl | hmem
    · exact hgood t hct htd' m
    · exact hgoodrest w' hmem m

/-- The band-`1/13` instance (the hdich lift-rigidity band): consecutive ratios `26/11`,
integer-friendly form `26·x ≤ 11·y`. -/
theorem spread_tower_13 (ws : List ℝ) (hpos : ∀ w ∈ ws, 0 < w)
    (hchain : List.IsChain (fun x y => 26 * x ≤ 11 * y) ws)
    (a L : ℝ) (hL : 0 ≤ L) (hentry : ∀ w ∈ ws.head?, 2 ≤ w * L) :
    ∃ t, a ≤ t ∧ t ≤ a + L ∧ ∀ w ∈ ws, ∀ m : ℤ, (1:ℝ)/13 ≤ |w * t - m| := by
  have h := spread_tower (1/13) (by norm_num) (by norm_num) ws hpos
    (List.IsChain.imp (fun x y hxy => by linarith) hchain) a L hL hentry
  obtain ⟨t, h1, h2, h3⟩ := h
  exact ⟨t, h1, h2, h3⟩

#print axioms exists_gap_subinterval
#print axioms spread_tower
#print axioms spread_tower_13

end GapDescent
end LonelyRunner
