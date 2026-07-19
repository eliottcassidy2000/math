/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Arithmetic consumers for THM-1235 first-lap Kakeya drift

This module kernel-checks the exact constant arithmetic, the fastest-pivot
diameter consequence, and the weighted Hamiltonian-path cut in THM-1235.
The geometric input -- freezing six moving open circle arcs on one complete
pivot lap -- remains an explicit inequality hypothesis.
-/

namespace LRC14
namespace FirstLapKakeyaDrift

/-- Six proper enlarged arcs have static length invoice
`6/7 + 2 * sum drift`; covering with strict surplus forces drift `> 1/14`. -/
theorem six_arc_invoice
    (e₁ e₂ e₃ e₄ e₅ e₆ : ℝ)
    (hcover :
      1 < (1 / 7 + 2 * |e₁|) + (1 / 7 + 2 * |e₂|) +
          (1 / 7 + 2 * |e₃|) + (1 / 7 + 2 * |e₄|) +
          (1 / 7 + 2 * |e₅|) + (1 / 7 + 2 * |e₆|)) :
    1 / 14 < |e₁| + |e₂| + |e₃| + |e₄| + |e₅| + |e₆| := by
  linarith

/-- If an enlarged radius reaches `1/2`, its drift alone exceeds the needed
`1/14` invoice. -/
theorem whole_arc_branch {e : ℝ} (he : 3 / 7 ≤ |e|) :
    1 / 14 < |e| := by
  linarith

/-- Telescoping identity behind the edge-weighted speed-path invoice. -/
theorem weighted_gap_identity (d₁ d₂ d₃ d₄ d₅ d₆ : ℝ) :
    (d₆ - d₁) + (d₆ - d₂) + (d₆ - d₃) +
        (d₆ - d₄) + (d₆ - d₅) + (d₆ - d₆) =
      (d₂ - d₁) + 2 * (d₃ - d₂) + 3 * (d₄ - d₃) +
        4 * (d₅ - d₄) + 5 * (d₆ - d₅) := by
  ring

/-- The fastest-pivot invoice forces relative diameter greater than `1/70`. -/
theorem fastest_invoice_forces_diameter
    {d₁ d₂ d₃ d₄ d₅ d₆ : ℝ}
    (horder : d₁ ≤ d₂ ∧ d₂ ≤ d₃ ∧ d₃ ≤ d₄ ∧ d₄ ≤ d₅ ∧ d₅ ≤ d₆)
    (hinvoice :
      d₆ / 14 < (d₆ - d₁) + (d₆ - d₂) + (d₆ - d₃) +
        (d₆ - d₄) + (d₆ - d₅) + (d₆ - d₆)) :
    d₆ / 70 < d₆ - d₁ := by
  rcases horder with ⟨h12, h23, h34, h45, h56⟩
  linarith

/-- Positive speeds convert the diameter form into `d₆/d₁ > 70/69`. -/
theorem fastest_invoice_forces_ratio
    {d₁ d₂ d₃ d₄ d₅ d₆ : ℝ}
    (hd₁ : 0 < d₁)
    (horder : d₁ ≤ d₂ ∧ d₂ ≤ d₃ ∧ d₃ ≤ d₄ ∧ d₄ ≤ d₅ ∧ d₅ ≤ d₆)
    (hinvoice :
      d₆ / 14 < (d₆ - d₁) + (d₆ - d₂) + (d₆ - d₃) +
        (d₆ - d₄) + (d₆ - d₅) + (d₆ - d₆)) :
    70 / 69 < d₆ / d₁ := by
  have hdiam := fastest_invoice_forces_diameter horder hinvoice
  apply (div_lt_div_iff₀ (by norm_num : (0 : ℝ) < 69) hd₁).2
  nlinarith

/-- The weighted invoice forces at least one macroscopic edge of the ordered
Hamiltonian path. -/
theorem fastest_invoice_forces_path_edge
    {d₁ d₂ d₃ d₄ d₅ d₆ : ℝ}
    (hinvoice :
      d₆ / 14 < (d₆ - d₁) + (d₆ - d₂) + (d₆ - d₃) +
        (d₆ - d₄) + (d₆ - d₅) + (d₆ - d₆)) :
    d₆ / 210 < d₂ - d₁ ∨ d₆ / 210 < d₃ - d₂ ∨
      d₆ / 210 < d₄ - d₃ ∨ d₆ / 210 < d₅ - d₄ ∨
      d₆ / 210 < d₆ - d₅ := by
  rw [weighted_gap_identity] at hinvoice
  by_contra h
  push Not at h
  rcases h with ⟨h₁, h₂, h₃, h₄, h₅⟩
  linarith

/-- Any macroscopic additive path edge gives the advertised multiplicative
`211/210` edge when its lower endpoint is positive and below `d₆`. -/
theorem additive_edge_forces_ratio {a b d₆ : ℝ}
    (ha : 0 < a) (hb6 : b ≤ d₆) (hedge : d₆ / 210 < b - a) :
    211 / 210 < b / a := by
  have hab : d₆ / 210 < b - a := hedge
  apply (div_lt_div_iff₀ (by norm_num : (0 : ℝ) < 210) ha).2
  nlinarith

/-- If only the fastest speed is guaranteed to make a full normalized lap,
its pivot invoice forces a carrier-scale `1/180` path edge. -/
theorem fastest_pivot_forces_carrier_edge
    {c d₁ d₂ d₃ d₄ d₅ d₆ : ℝ}
    (hthreshold : 7 * c / 6 ≤ d₆)
    (hinvoice :
      d₆ / 14 < (d₆ - d₁) + (d₆ - d₂) + (d₆ - d₃) +
        (d₆ - d₄) + (d₆ - d₅) + (d₆ - d₆)) :
    c / 180 < d₂ - d₁ ∨ c / 180 < d₃ - d₂ ∨
      c / 180 < d₄ - d₃ ∨ c / 180 < d₅ - d₄ ∨
      c / 180 < d₆ - d₅ := by
  rw [weighted_gap_identity] at hinvoice
  by_contra h
  push Not at h
  rcases h with ⟨h₁, h₂, h₃, h₄, h₅⟩
  linarith

/-- When the two fastest speeds make full laps, the penultimate-pivot weight
vector `(1,2,3,4,1)` improves the carrier cut to `1/132`. -/
theorem penultimate_pivot_forces_carrier_edge
    {c d₁ d₂ d₃ d₄ d₅ d₆ : ℝ}
    (hthreshold : 7 * c / 6 ≤ d₅)
    (hinvoice :
      d₅ / 14 < (d₅ - d₁) + (d₅ - d₂) + (d₅ - d₃) +
        (d₅ - d₄) + (d₅ - d₅) + (d₆ - d₅)) :
    c / 132 < d₂ - d₁ ∨ c / 132 < d₃ - d₂ ∨
      c / 132 < d₄ - d₃ ∨ c / 132 < d₅ - d₄ ∨
      c / 132 < d₆ - d₅ := by
  by_contra h
  push Not at h
  rcases h with ⟨h₁, h₂, h₃, h₄, h₅⟩
  linarith

/-- Three full-lap speeds expose the central-pivot weight vector
`(1,2,3,2,1)` and the carrier cut `1/108`. -/
theorem central_pivot_forces_carrier_edge
    {c d₁ d₂ d₃ d₄ d₅ d₆ : ℝ}
    (hthreshold : 7 * c / 6 ≤ d₄)
    (hinvoice :
      d₄ / 14 < (d₄ - d₁) + (d₄ - d₂) + (d₄ - d₃) +
        (d₄ - d₄) + (d₅ - d₄) + (d₆ - d₄)) :
    c / 108 < d₂ - d₁ ∨ c / 108 < d₃ - d₂ ∨
      c / 108 < d₄ - d₃ ∨ c / 108 < d₅ - d₄ ∨
      c / 108 < d₆ - d₅ := by
  by_contra h
  push Not at h
  rcases h with ⟨h₁, h₂, h₃, h₄, h₅⟩
  linarith

/-- Exact limiting absolute-cut constants from the full-lap threshold. -/
theorem limiting_cut_constants :
    ((7 / 6 : ℚ) / 70 = 1 / 60) ∧
      ((7 / 6 : ℚ) / 210 = 1 / 180) := by
  norm_num

end FirstLapKakeyaDrift
end LRC14
