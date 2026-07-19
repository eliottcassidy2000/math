/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Arithmetic core for THM-1252's coherent two-wall fork

The paper layer chooses every blocker from one irredundant tooth cover (with
least-speed selection as a fallback) and uses strict six-cover topology at
both walls of an ascent target tooth.  This module checks the resulting
quantum packing, exact detuned toothpick law, finite rung bound, binary-edge
reflection, and positioned-Hunter scalar invoice.
-/

namespace LRC14
namespace MinimalBlockerTwoWallFork

/-- Every blocker target tooth fits in the carrier gap around a centered
spoke, without a source/target speed-order assumption. -/
theorem target_tooth_fits_centered_carrier_margin
    {c source target : ℝ} (hc : 0 < c)
    (hsource : c < source) (htarget : c < target) :
    1 / (7 * target) + 1 / (2 * (c + source)) < 3 / (7 * c) := by
  have hsource0 : (0 : ℝ) < source := lt_trans hc hsource
  have htarget0 : (0 : ℝ) < target := lt_trans hc htarget
  have htargetDen : (0 : ℝ) < 7 * target := by positivity
  have hc7 : (0 : ℝ) < 7 * c := by positivity
  have hqDen : (0 : ℝ) < 2 * (c + source) := by positivity
  have hc4 : (0 : ℝ) < 4 * c := by positivity
  have hwidth : (1 : ℝ) / (7 * target) < 1 / (7 * c) := by
    apply (div_lt_div_iff₀ htargetDen hc7).2
    nlinarith
  have hround : (1 : ℝ) / (2 * (c + source)) < 1 / (4 * c) := by
    apply (div_lt_div_iff₀ hqDen hc4).2
    nlinarith
  have hbudget : (1 : ℝ) / (7 * c) + 1 / (4 * c) < 3 / (7 * c) := by
    field_simp
    nlinarith
  linarith

/-- Two wall credits, each no longer than one target half-tooth, fit as
disjoint inward subneedles in the target tooth. -/
theorem two_wall_quanta_fit
    {j qLeft qRight : ℝ} (hj : 0 < j)
    (hqLeft : qLeft ≤ 1 / (14 * j))
    (hqRight : qRight ≤ 1 / (14 * j)) :
    qLeft + qRight ≤ 1 / (7 * j) := by
  have hsplit : (1 : ℝ) / (7 * j) = 1 / (14 * j) + 1 / (14 * j) := by
    field_simp
    ring
  rw [hsplit]
  linarith

/-- A complete safe gap of the common wall provider strictly inside the
target danger tooth forces the first toothpick ratio. -/
theorem same_wall_label_forces_sixfold_ratio
    {j h : ℝ} (hj : 0 < j) (hh : 0 < h)
    (hgap : 6 / (7 * h) < 1 / (7 * j)) :
    6 * j < h := by
  have hcross := (div_lt_div_iff₀ (show (0 : ℝ) < 7 * h by positivity)
    (show (0 : ℝ) < 7 * j by positivity)).mp hgap
  nlinarith

/-- Exact functional form of a same-label address return: the two seam
lengths equal the detuning from the multiplier `7r-1`. -/
theorem detuned_toothpick_rung_identity
    {j h r seamSum : ℝ} (hj : j ≠ 0) (hh : h ≠ 0)
    (hseam : seamSum = 1 / (7 * j) - (7 * r - 1) / (7 * h)) :
    7 * j * h * seamSum = h - (7 * r - 1) * j := by
  rw [hseam]
  field_simp

/-- Under THM-1233's projective ratio ceiling, the integer toothpick address
return has at most 335 rungs. -/
theorem detuned_rung_index_ceiling
    {j h r : ℕ} (hj : 0 < j)
    (hcompact : h < 2345 * j)
    (hrung : (7 * r - 1) * j < h) :
    r ≤ 335 := by
  by_contra hnot
  have hr : 336 ≤ r := by omega
  have hcoef : 2351 ≤ 7 * r - 1 := by omega
  have hmul : 2351 * j ≤ (7 * r - 1) * j := by
    exact Nat.mul_le_mul_right j hcoef
  have hjump : 2345 * j < 2351 * j := by
    nlinarith
  omega

/-- A positive integer detuning on the common gcd sheet is at least that gcd. -/
theorem positive_detuning_quantum
    {g detuning : ℕ} (hpos : 0 < detuning) (hdiv : g ∣ detuning) :
    g ≤ detuning := by
  exact Nat.le_of_dvd hpos hdiv

/-- In the same-provider branch, two quanta on the wall pair plus one quantum
on a second support edge give the three-quantum positioned invoice. -/
theorem same_provider_three_quantum_invoice
    {wallPair secondEdge quantum total : ℝ}
    (hwall : 2 * quantum ≤ wallPair)
    (hsecond : quantum ≤ secondEdge)
    (htotal : wallPair + secondEdge ≤ total) :
    3 * quantum ≤ total := by
  linarith

/-- The exact detuned-rung seam formula automatically gives positivity of the
detuning whenever the two seam lengths have positive total. -/
theorem positive_seam_sum_forces_positive_detuning
    {j h r seamSum : ℝ} (hj : 0 < j) (hh : 0 < h)
    (hidentity : 7 * j * h * seamSum = h - (7 * r - 1) * j)
    (hseam : 0 < seamSum) :
    (7 * r - 1) * j < h := by
  have hscale : 0 < 7 * j * h := by positivity
  nlinarith [mul_pos hscale hseam]

/-- Exact smallest-detuning local guardrail from the paper theorem:
`(c,i,j,h)=(5,7,8,49)` has two equal wall quanta. -/
theorem smallest_detuning_guardrail :
    (49 : ℕ) - 6 * 8 = 1 ∧
      (1 : ℚ) / 5488 + 1 / 5488 = 1 / 2744 ∧
      (1 : ℚ) / (14 * (8 * 49)) = 1 / 5488 := by
  norm_num

/-- A binary THM-1248 digit makes both the original and reflected coherent
factors nonnegative.  The marked-tooth order chooses which one is used. -/
theorem binary_digit_original_reflected_nonnegative
    (c k delta : ℤ) (hk : 0 ≤ k) (hkc : k < c)
    (hdelta0 : 0 ≤ delta) (hdelta1 : delta ≤ 1) :
    0 ≤ k + delta ∧ 0 ≤ c - k - delta := by
  omega

/-- Even the terminal compact rung `r=335` is compatible with the binary
half-circle family from the paper layer. -/
theorem terminal_binary_rung_compact :
    (7 * 335 - 1) * 10004 + 1 < 2345 * 10001 := by
  norm_num

#print axioms target_tooth_fits_centered_carrier_margin
#print axioms two_wall_quanta_fit
#print axioms same_wall_label_forces_sixfold_ratio
#print axioms detuned_toothpick_rung_identity
#print axioms detuned_rung_index_ceiling
#print axioms positive_detuning_quantum
#print axioms same_provider_three_quantum_invoice
#print axioms positive_seam_sum_forces_positive_detuning
#print axioms smallest_detuning_guardrail
#print axioms binary_digit_original_reflected_nonnegative
#print axioms terminal_binary_rung_compact

end MinimalBlockerTwoWallFork
end LRC14
