/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Arithmetic core for THM-1264 chronological return polygons

The paper layer supplies a literal consecutive tooth subword with matching
intermediate occurrences.  This module checks its telescoping overlap law,
the resulting ratio ascent, the triangle specialization, and compact height.
-/

namespace LRC14
namespace ChronologicalReturnPolygon

noncomputable def centre (speed address : ℕ → ℝ) (q : ℕ) : ℝ :=
  address q / speed q

noncomputable def halfWidth (speed : ℕ → ℝ) (q : ℕ) : ℝ :=
  1 / (14 * speed q)

noncomputable def rawOverlap (speed address : ℕ → ℝ) (q : ℕ) : ℝ :=
  centre speed address q + halfWidth speed q -
    (centre speed address (q + 1) - halfWidth speed (q + 1))

/-- The abstract interval-chain telescope.  The endpoint half-widths combine
because the first and last teeth have the same owner speed. -/
theorem interval_chain_telescope
    (x h : ℕ → ℝ) (m : ℕ) (hh : h m = h 0) :
    Finset.sum (Finset.range m) (fun q =>
        (x q + h q) - (x (q + 1) - h (q + 1))) =
      x 0 - x m + 2 * Finset.sum (Finset.range m) h := by
  have hx := Finset.sum_range_sub x m
  have hwidth := Finset.sum_range_sub h m
  have hpointwise : ∀ q : ℕ,
      (x q + h q) - (x (q + 1) - h (q + 1)) =
        -(x (q + 1) - x q) + (h (q + 1) - h q) + 2 * h q := by
    intro q
    ring
  calc
    Finset.sum (Finset.range m) (fun q =>
        (x q + h q) - (x (q + 1) - h (q + 1))) =
        Finset.sum (Finset.range m) (fun q =>
          (-(x (q + 1) - x q) + (h (q + 1) - h q) + 2 * h q)) := by
            apply Finset.sum_congr rfl
            intro q _
            exact hpointwise q
    _ = -(Finset.sum (Finset.range m) (fun q => x (q + 1) - x q)) +
          Finset.sum (Finset.range m) (fun q => h (q + 1) - h q) +
          2 * Finset.sum (Finset.range m) h := by
            simp only [Finset.sum_add_distrib, Finset.sum_neg_distrib,
              Finset.mul_sum]
    _ = x 0 - x m + 2 * Finset.sum (Finset.range m) h := by
          rw [hx, hwidth, hh]
          ring

/-- Exact all-length return-polygon identity for a literal consecutive tooth
subword.  `address m = address 0 + R` is the positive address return. -/
theorem chronological_return_polygon_identity
    (speed address : ℕ → ℝ) (m R : ℕ)
    (hs0 : speed 0 ≠ 0)
    (hsclose : speed m = speed 0)
    (hnclose : address m = address 0 + R) :
    Finset.sum (Finset.range m) (rawOverlap speed address) =
      (1 / 7) * Finset.sum (Finset.range m) (fun q => 1 / speed q) -
        R / speed 0 := by
  have hhclose : halfWidth speed m = halfWidth speed 0 := by
    simp [halfWidth, hsclose]
  have htelescope := interval_chain_telescope
    (centre speed address) (halfWidth speed) m hhclose
  have hcentre : centre speed address 0 - centre speed address m =
      -R / speed 0 := by
    simp only [centre, hsclose, hnclose]
    field_simp [hs0]
    ring
  have hwidthsum :
      2 * Finset.sum (Finset.range m) (halfWidth speed) =
        (1 / 7) * Finset.sum (Finset.range m) (fun q => 1 / speed q) := by
    rw [Finset.mul_sum, Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro q _
    simp only [halfWidth]
    ring
  change Finset.sum (Finset.range m) (fun q =>
      (centre speed address q + halfWidth speed q) -
        (centre speed address (q + 1) - halfWidth speed (q + 1))) = _
  rw [htelescope, hcentre, hwidthsum]
  ring

/-- Positivity of a return polygon converts the exact overlap identity to the
dimensionless address-drift inequality. -/
theorem positive_return_forces_ratio_sum
    {a R overlapSum reciprocalInternal : ℝ}
    (ha : 0 < a) (hpositive : 0 < overlapSum)
    (hidentity : overlapSum =
      (1 / 7) * (1 / a + reciprocalInternal) - R / a) :
    7 * R - 1 < a * reciprocalInternal := by
  have hscale : 0 < 7 * a := by positivity
  have hmul := mul_lt_mul_of_pos_left hpositive hscale
  rw [hidentity] at hmul
  field_simp at hmul
  nlinarith

/-- A positive triangle return with integer return at least one gives a strict
factor-three ascent from its smaller internal owner to the repeated owner. -/
theorem triangle_return_forces_threefold_ascent
    {a b c R overlapSum : ℝ}
    (ha : 0 < a) (hb : 0 < b) (hc : 0 < c)
    (hR : 1 ≤ R) (hpositive : 0 < overlapSum)
    (hidentity : overlapSum =
      (1 / 7) * (1 / a + 1 / b + 1 / c) - R / a) :
    3 * min b c < a := by
  have hratio : 7 * R - 1 < a * (1 / b + 1 / c) :=
    positive_return_forces_ratio_sum ha hpositive (by
      simpa only [add_assoc] using hidentity)
  have hsum : 6 < a / b + a / c := by
    have hratio' : 7 * R - 1 < a / b + a / c := by
      calc
        7 * R - 1 < a * (1 / b + 1 / c) := hratio
        _ = a / b + a / c := by ring
    have hR6 : 6 ≤ 7 * R - 1 := by nlinarith
    linarith
  rcases le_total b c with hbc | hcb
  · rw [min_eq_left hbc]
    by_contra hnot
    have ha3 : a ≤ 3 * b := le_of_not_gt hnot
    have hinv : 1 / c ≤ 1 / b := one_div_le_one_div_of_le hb hbc
    have hinv' : c⁻¹ ≤ b⁻¹ := by simpa [one_div] using hinv
    have hac : a / c ≤ a / b := by
      simpa only [div_eq_mul_inv] using
        mul_le_mul_of_nonneg_left hinv' (le_of_lt ha)
    have hab : a / b ≤ 3 := by
      exact (div_le_iff₀ hb).2 (by nlinarith)
    linarith
  · rw [min_eq_right hcb]
    by_contra hnot
    have ha3 : a ≤ 3 * c := le_of_not_gt hnot
    have hinv : 1 / b ≤ 1 / c := one_div_le_one_div_of_le hc hcb
    have hinv' : b⁻¹ ≤ c⁻¹ := by simpa [one_div] using hinv
    have hab : a / b ≤ a / c := by
      simpa only [div_eq_mul_inv] using
        mul_le_mul_of_nonneg_left hinv' (le_of_lt ha)
    have hac : a / c ≤ 3 := by
      exact (div_le_iff₀ hc).2 (by nlinarith)
    linarith

/-- Eight strict factor-three ascents cannot fit inside the projective
`d<2345c` box; seven ascents remain numerically possible. -/
theorem factor_three_compact_height :
    3 ^ 7 < (2345 : ℕ) ∧ (2345 : ℕ) < 3 ^ 8 := by
  norm_num

/-- If one deliberately ignores the stronger fixed-six-owner DAG bound, the
weakest simple-return factor `6/5` fits 42, but not 43, projective ascents. -/
theorem six_fifths_compact_height :
    ((6 : ℚ) / 5) ^ 42 < 2345 ∧ (2345 : ℚ) < ((6 : ℚ) / 5) ^ 43 := by
  norm_num

#print axioms interval_chain_telescope
#print axioms chronological_return_polygon_identity
#print axioms positive_return_forces_ratio_sum
#print axioms triangle_return_forces_threefold_ascent
#print axioms factor_three_compact_height
#print axioms six_fifths_compact_height

end ChronologicalReturnPolygon
end LRC14
