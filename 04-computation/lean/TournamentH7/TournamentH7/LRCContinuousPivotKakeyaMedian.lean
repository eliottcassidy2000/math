/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Arithmetic consumers for THM-1236 continuous-pivot Kakeya drift

The geometric provider says that a six-comb cover pays the absolute-drift
invoice at every auxiliary scale `y >= 7*c/6`.  This module kernel-checks the
tilted `L¹` optimization: the exact constrained minimizer is
`max (7*c/6) d₄`.  No circle-cover theorem is hidden below.
-/

namespace LRC14
namespace ContinuousPivotKakeyaMedian

/-- The upper median minimizes the six-point `L¹` functional tilted by
`-y/14`.  The small tilt resolves the ordinary median interval `[d₃,d₄]`
toward its right endpoint. -/
theorem upper_median_minimizes_tilted_l1
    {d₁ d₂ d₃ d₄ d₅ d₆ y : ℝ}
    (h12 : d₁ ≤ d₂) (h23 : d₂ ≤ d₃) (h34 : d₃ ≤ d₄)
    (h45 : d₄ ≤ d₅) (h56 : d₅ ≤ d₆) :
    |d₁ - d₄| + |d₂ - d₄| + |d₃ - d₄| + |d₄ - d₄| +
          |d₅ - d₄| + |d₆ - d₄| - d₄ / 14 ≤
      |d₁ - y| + |d₂ - y| + |d₃ - y| + |d₄ - y| +
          |d₅ - y| + |d₆ - y| - y / 14 := by
  have h14 : d₁ ≤ d₄ := h12.trans (h23.trans h34)
  have h24 : d₂ ≤ d₄ := h23.trans h34
  have h15 : d₁ ≤ d₅ := h14.trans h45
  have h25 : d₂ ≤ d₅ := h24.trans h45
  have h36 : d₃ ≤ d₆ := h34.trans (h45.trans h56)
  by_cases hy : y ≤ d₄
  · have hp1 := abs_sub_le d₁ y d₄
    have hp2 := abs_sub_le d₂ y d₅
    have hp3 := abs_sub_le d₃ y d₆
    simp only [sub_self, abs_zero]
    rw [abs_of_nonpos (sub_nonpos.mpr h14),
        abs_of_nonpos (sub_nonpos.mpr h24),
        abs_of_nonpos (sub_nonpos.mpr h34),
        abs_of_nonneg (sub_nonneg.mpr h45),
        abs_of_nonneg (sub_nonneg.mpr (h45.trans h56))]
    rw [abs_of_nonpos (sub_nonpos.mpr h14)] at hp1
    rw [abs_of_nonpos (sub_nonpos.mpr h25)] at hp2
    rw [abs_of_nonpos (sub_nonpos.mpr h36)] at hp3
    rw [abs_sub_comm y d₄] at hp1
    rw [abs_sub_comm y d₅] at hp2
    rw [abs_sub_comm y d₆] at hp3
    linarith
  · have h4y : d₄ ≤ y := le_of_not_ge hy
    have h1y : d₁ ≤ y := h14.trans h4y
    have h2y : d₂ ≤ y := h24.trans h4y
    have h3y : d₃ ≤ y := h34.trans h4y
    have hp5 := abs_sub_le d₅ y d₄
    have hp6 := abs_sub_le d₆ y d₄
    simp only [sub_self, abs_zero]
    rw [abs_of_nonpos (sub_nonpos.mpr h14),
        abs_of_nonpos (sub_nonpos.mpr h24),
        abs_of_nonpos (sub_nonpos.mpr h34),
        abs_of_nonneg (sub_nonneg.mpr h45),
        abs_of_nonneg (sub_nonneg.mpr (h45.trans h56)),
        abs_of_nonpos (sub_nonpos.mpr h1y),
        abs_of_nonpos (sub_nonpos.mpr h2y),
        abs_of_nonpos (sub_nonpos.mpr h3y),
        abs_of_nonpos (sub_nonpos.mpr h4y)]
    rw [abs_of_nonneg (sub_nonneg.mpr h45)] at hp5
    rw [abs_of_nonneg (sub_nonneg.mpr (h45.trans h56))] at hp6
    rw [abs_of_nonneg (sub_nonneg.mpr h4y)] at hp5 hp6
    linarith

/-- To the right of the upper median, the tilted six-point `L¹` functional is
monotone.  Four coordinates move outwards while at most two move inwards. -/
theorem tilted_l1_monotone_right
    {d₁ d₂ d₃ d₄ d₅ d₆ a y : ℝ}
    (h14 : d₁ ≤ d₄) (h24 : d₂ ≤ d₄) (h34 : d₃ ≤ d₄)
    (h4a : d₄ ≤ a) (hay : a ≤ y) :
    |d₁ - a| + |d₂ - a| + |d₃ - a| + |d₄ - a| +
          |d₅ - a| + |d₆ - a| - a / 14 ≤
      |d₁ - y| + |d₂ - y| + |d₃ - y| + |d₄ - y| +
          |d₅ - y| + |d₆ - y| - y / 14 := by
  have h1a := h14.trans h4a
  have h2a := h24.trans h4a
  have h3a := h34.trans h4a
  have h1y := h1a.trans hay
  have h2y := h2a.trans hay
  have h3y := h3a.trans hay
  have h4y := h4a.trans hay
  have hp5 := abs_sub_le d₅ y a
  have hp6 := abs_sub_le d₆ y a
  rw [abs_of_nonpos (sub_nonpos.mpr h1a),
      abs_of_nonpos (sub_nonpos.mpr h2a),
      abs_of_nonpos (sub_nonpos.mpr h3a),
      abs_of_nonpos (sub_nonpos.mpr h4a),
      abs_of_nonpos (sub_nonpos.mpr h1y),
      abs_of_nonpos (sub_nonpos.mpr h2y),
      abs_of_nonpos (sub_nonpos.mpr h3y),
      abs_of_nonpos (sub_nonpos.mpr h4y)]
  rw [abs_of_nonneg (sub_nonneg.mpr hay)] at hp5 hp6
  linarith

/-- Exact constrained optimizer used by THM-1236. -/
theorem constrained_tilted_l1_optimizer
    {c d₁ d₂ d₃ d₄ d₅ d₆ y : ℝ}
    (h12 : d₁ ≤ d₂) (h23 : d₂ ≤ d₃) (h34 : d₃ ≤ d₄)
    (h45 : d₄ ≤ d₅) (h56 : d₅ ≤ d₆)
    (hy : 7 * c / 6 ≤ y) :
    let y₀ := max (7 * c / 6) d₄
    |d₁ - y₀| + |d₂ - y₀| + |d₃ - y₀| + |d₄ - y₀| +
          |d₅ - y₀| + |d₆ - y₀| - y₀ / 14 ≤
      |d₁ - y| + |d₂ - y| + |d₃ - y| + |d₄ - y| +
          |d₅ - y| + |d₆ - y| - y / 14 := by
  dsimp
  by_cases ha : 7 * c / 6 ≤ d₄
  · rw [max_eq_right ha]
    exact upper_median_minimizes_tilted_l1 h12 h23 h34 h45 h56
  · have h4a : d₄ ≤ 7 * c / 6 := le_of_not_ge ha
    rw [max_eq_left h4a]
    exact tilted_l1_monotone_right
      (h12.trans (h23.trans h34)) (h23.trans h34) h34 h4a hy

/-- The optimized geometric invoice follows immediately from the all-scale
continuous-pivot provider. -/
theorem optimized_invoice_of_all_pivots
    {c d₁ d₂ d₃ d₄ d₅ d₆ : ℝ}
    (hinvoice : ∀ y, 7 * c / 6 ≤ y →
      y / 14 < |d₁ - y| + |d₂ - y| + |d₃ - y| +
        |d₄ - y| + |d₅ - y| + |d₆ - y|) :
    let y₀ := max (7 * c / 6) d₄
    y₀ / 14 < |d₁ - y₀| + |d₂ - y₀| + |d₃ - y₀| +
      |d₄ - y₀| + |d₅ - y₀| + |d₆ - y₀| := by
  dsimp
  exact hinvoice _ (le_max_left _ _)

/-- At the upper median, the absolute drift is the top-three/bottom-three
half-sum imbalance. -/
theorem upper_median_half_sum_identity
    {d₁ d₂ d₃ d₄ d₅ d₆ : ℝ}
    (h14 : d₁ ≤ d₄) (h24 : d₂ ≤ d₄) (h34 : d₃ ≤ d₄)
    (h45 : d₄ ≤ d₅) (h46 : d₄ ≤ d₆) :
    |d₁ - d₄| + |d₂ - d₄| + |d₃ - d₄| + |d₄ - d₄| +
        |d₅ - d₄| + |d₆ - d₄| =
      d₄ + d₅ + d₆ - d₁ - d₂ - d₃ := by
  simp only [sub_self, abs_zero]
  rw [abs_of_nonpos (sub_nonpos.mpr h14),
      abs_of_nonpos (sub_nonpos.mpr h24),
      abs_of_nonpos (sub_nonpos.mpr h34),
      abs_of_nonneg (sub_nonneg.mpr h45),
      abs_of_nonneg (sub_nonneg.mpr h46)]
  ring

/-- Exact rational constants inherited by the three suffix-count branches. -/
theorem branch_edge_constants :
    ((7 / 6 : ℚ) / 14) / 15 = 1 / 180 ∧
      ((7 / 6 : ℚ) / 14) / 11 = 1 / 132 ∧
      ((7 / 6 : ℚ) / 14) / 9 = 1 / 108 := by
  norm_num

#print axioms upper_median_minimizes_tilted_l1
#print axioms tilted_l1_monotone_right
#print axioms constrained_tilted_l1_optimizer
#print axioms optimized_invoice_of_all_pivots
#print axioms upper_median_half_sum_identity
#print axioms branch_edge_constants

end ContinuousPivotKakeyaMedian
end LRC14
