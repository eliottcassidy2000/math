/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Binary centered phases land on the irredundant tooth word (THM-1256)

The paper layer combines the interval-chain topology of THM-1253 with the
coherent binary edge of THM-1254.  This module checks the exact endpoint
interpretation of the residual invoice, the local alignment/mismatch interval
lemmas, the same-provider toothpick detuning, the ABAB contradiction, and the
resulting turn-count arithmetic.
-/

namespace LRC14
namespace BinaryPhaseWordLanding

/-- Reciprocal-centre drift is exactly the determinant of the two endpoint
tooth centres. -/
theorem reciprocal_endpoint_identity
    (s₀ sᵣ n₀ nᵣ : ℚ) (hn₀ : n₀ ≠ 0) (hnᵣ : nᵣ ≠ 0) :
    n₀ * nᵣ * (s₀ / n₀ - sᵣ / nᵣ) = nᵣ * s₀ - n₀ * sᵣ := by
  field_simp [hn₀, hnᵣ]

/-- The apparent residual surplus is just the marked endpoint-address
difference. -/
theorem residual_minus_endpoint
    (P n₀ nᵣ s₀ sᵣ k delta : ℚ)
    (hP : P = nᵣ + k + delta) :
    (P * s₀ - n₀ * sᵣ) - (nᵣ * s₀ - n₀ * sᵣ) =
      s₀ * (k + delta) := by
  rw [hP]
  ring

/-- Consequently the THM-1254 invoice is address order, not a second metric
overlap credit. -/
theorem residual_invoice_iff_address_order
    {P n₀ nᵣ s₀ sᵣ : ℚ} (hs₀ : 0 < s₀) :
    nᵣ * s₀ - n₀ * sᵣ ≤ P * s₀ - n₀ * sᵣ ↔ nᵣ ≤ P := by
  constructor <;> intro h <;> nlinarith

/-- Binary digit zero puts the source centered phase after the target phase. -/
theorem digit_zero_phase_order
    {Q theta phaseDiff : ℚ}
    (hQ : 0 < Q) (htheta : 0 < theta)
    (htransport : Q * phaseDiff = theta) :
    0 < phaseDiff := by
  nlinarith

/-- Binary digit one puts the source centered phase before the target phase. -/
theorem digit_one_phase_order
    {Q theta phaseDiff : ℚ}
    (hQ : 0 < Q) (htheta : theta < 1)
    (htransport : Q * phaseDiff = theta - 1) :
    phaseDiff < 0 := by
  nlinarith

/-- If interval A precedes interval B but their marked interior phases occur
in the reverse order, A and B overlap. -/
theorem phase_position_mismatch_forces_overlap
    {aL aR bL bR x y : ℝ}
    (_hleft : aL < bL)
    (_hxA : aL < x) (hxA' : x < aR)
    (hyB : bL < y) (_hyB' : y < bR)
    (hmismatch : y < x) :
    bL < aR := by
  linarith

/-- THM-1253 separation therefore excludes a mismatch between nonconsecutive
members of the minimal interval word. -/
theorem separated_intervals_cannot_mismatch
    {aL aR bL bR x y : ℝ}
    (hleft : aL < bL)
    (hxA : aL < x) (hxA' : x < aR)
    (hyB : bL < y) (hyB' : y < bR)
    (hseparated : aR ≤ bL) :
    ¬ y < x := by
  intro hmismatch
  have := phase_position_mismatch_forces_overlap
    hleft hxA hxA' hyB hyB' hmismatch
  linarith

/-- Two adjacent phase/position mismatches would reverse the two
nonconsecutive outer marks, contradicting THM-1253 separation.  Thus the
inversion graph on marked teeth is a matching. -/
theorem adjacent_mismatches_cannot_share_a_mark
    {aR cL x y z : ℝ}
    (hxA : x < aR) (hzC : cL < z)
    (hseparated : aR ≤ cL)
    (hxy : y < x) (hyz : z < y) :
    False := by
  linarith

/-- Two overlapping ordered intervals cover every point between an interior
mark of the first and a later interior mark of the second.  Iteration gives
the paper subword statement. -/
theorem aligned_overlapping_pair_covers_segment
    {aL aR bL bR x y z : ℝ}
    (_hleft : aL < bL) (_hright : aR < bR)
    (hoverlap : bL < aR)
    (hxA : aL < x) (_hxA' : x < aR)
    (_hyB : bL < y) (hyB' : y < bR)
    (hxz : x ≤ z) (hzy : z ≤ y) :
    (aL < z ∧ z < aR) ∨ (bL < z ∧ z < bR) := by
  by_cases hzA : z < aR
  · left
    exact ⟨lt_of_lt_of_le hxA hxz, hzA⟩
  · right
    have hzA' : aR ≤ z := le_of_not_gt hzA
    exact ⟨lt_of_lt_of_le hoverlap hzA', lt_of_le_of_lt hzy hyB'⟩

/-- Exact same-provider toothpick law for a backtrack with address rung r. -/
theorem backtrack_detuning_identity
    (h j r : ℚ) (hh : h ≠ 0) (hj : j ≠ 0) :
    1 / (7 * j) - (7 * r - 1) / (7 * h) =
      (h - (7 * r - 1) * j) / (7 * j * h) := by
  field_simp [hh, hj]

/-- A positive backtrack detuning forces the repeated outer owner to be more
than six times the middle owner. -/
theorem backtrack_outer_gt_six_middle
    {h j r : ℤ}
    (hj : 0 < j) (hr : 1 ≤ r)
    (hdetune : 0 < h - (7 * r - 1) * j) :
    6 * j < h := by
  have hc : (6 : ℤ) ≤ 7 * r - 1 := by omega
  have hm : 6 * j ≤ (7 * r - 1) * j :=
    mul_le_mul_of_nonneg_right hc (le_of_lt hj)
  have hd : (7 * r - 1) * j < h := by omega
  exact lt_of_le_of_lt hm hd

/-- Oppositely oriented backtracks on one unordered owner pair are
incompatible.  This is the arithmetic core of the ABAB exclusion. -/
theorem abab_speed_contradiction
    {a b : ℤ} (_ha : 0 < a) (hb : 0 < b)
    (hab : 6 * b < a) (hba : 6 * a < b) :
    False := by
  nlinarith

/-- Five successive toothpick returns cannot fit in the compact LRC14 speed
box: each return multiplies speed by more than six, while `6^5>2345`. -/
theorem no_five_backtrack_ascent_chain
    {c x₀ x₁ x₂ x₃ x₄ x₅ : ℝ}
    (hc : 0 < c) (hx₀ : c < x₀)
    (h01 : 6 * x₀ < x₁) (h12 : 6 * x₁ < x₂)
    (h23 : 6 * x₂ < x₃) (h34 : 6 * x₃ < x₄)
    (h45 : 6 * x₄ < x₅) (hcap : x₅ < 2345 * c) :
    False := by
  nlinarith

/-- If B internal positions are backtracks and T are nonbacktracking turns,
the absence of consecutive backtracks gives the sharp half-density turn
floor. -/
theorem turn_floor_from_no_consecutive_backtracks
    {N B T : ℕ}
    (hN : 2 ≤ N)
    (hledger : B + T = N - 2)
    (hseparated : 2 * B ≤ N - 1) :
    (N - 2) / 2 ≤ T := by
  omega

#print axioms reciprocal_endpoint_identity
#print axioms residual_minus_endpoint
#print axioms residual_invoice_iff_address_order
#print axioms digit_zero_phase_order
#print axioms digit_one_phase_order
#print axioms phase_position_mismatch_forces_overlap
#print axioms separated_intervals_cannot_mismatch
#print axioms adjacent_mismatches_cannot_share_a_mark
#print axioms aligned_overlapping_pair_covers_segment
#print axioms backtrack_detuning_identity
#print axioms backtrack_outer_gt_six_middle
#print axioms abab_speed_contradiction
#print axioms no_five_backtrack_ascent_chain
#print axioms turn_floor_from_no_consecutive_backtracks

end BinaryPhaseWordLanding
end LRC14
