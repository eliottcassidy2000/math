/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Arithmetic core for THM-1250's fully located six-cover tree

The paper layer uses THM-1198's six private regions to force all six labels
into an irredundant tooth chain, then extracts a spanning tree of five actual
interior handoffs.  This module checks the gcd/lcm quantum, the exact Hunter
rearrangement, its scale-covariant form, and the finite private-stalk constants.
-/

namespace LRC14
namespace SixPrivateLocatedTree

/-- A positive handoff numerator divisible by the two speeds' gcd is at least
that gcd. -/
theorem positive_handoff_numerator_quantum
    {g numerator : ℕ} (hpos : 0 < numerator) (hdiv : g ∣ numerator) :
    g ≤ numerator := by
  exact Nat.le_of_dvd hpos hdiv

/-- The natural-number identity converting the gcd tooth quantum to the lcm
clock quantum. -/
theorem gcd_mul_lcm_identity (u v : ℕ) :
    Nat.gcd u v * Nat.lcm u v = u * v := by
  exact Nat.gcd_mul_lcm u v

/-- Six combs on a complete slow gap leave the Hunter base `invc`; five
located tree overlaps contribute with coefficient `49/6`. -/
theorem fully_located_hunter_debt
    {invc harmonic overlapSum : ℝ}
    (hcovered : 0 ≥ (6 * invc / 7) / 7 - 6 * harmonic / 49 + overlapSum) :
    invc + 49 * overlapSum / 6 ≤ harmonic := by
  linarith

/-- Replacing every actual overlap by its `1/(14*lcm)` quantum gives the
scale-covariant lcm debt. -/
theorem fully_located_lcm_debt
    {invc harmonic overlapSum lcmSum : ℝ}
    (hcovered : 0 ≥ (6 * invc / 7) / 7 - 6 * harmonic / 49 + overlapSum)
    (hquantum : lcmSum / 14 ≤ overlapSum) :
    invc + 7 * lcmSum / 12 ≤ harmonic := by
  linarith

/-- Cayley averaging on `K₆` turns all disjoint repeated handoffs into one
fixed maximum-tree Hunter certificate. -/
theorem multiplicity_averaged_debt
    {invc harmonic overlapSum maxTree totalWeight : ℝ}
    (hcovered : 0 ≥ (6 * invc / 7) / 7 - 6 * harmonic / 49 + overlapSum)
    (hquantum : maxTree / 14 ≤ overlapSum)
    (haverage : totalWeight ≤ 3 * maxTree) :
    invc + 7 * totalWeight / 36 ≤ harmonic := by
  linarith

/-- Multiplying by the positive carrier converts the normalized inequality to
`cH-1 ≥ (7c/12) sum 1/lcm`. -/
theorem scale_covariant_lcm_debt
    {c invc harmonic lcmSum : ℝ} (hc : 0 ≤ c)
    (hinv : c * invc = 1)
    (hdebt : invc + 7 * lcmSum / 12 ≤ harmonic) :
    1 + 7 * c * lcmSum / 12 ≤ c * harmonic := by
  nlinarith

/-- A small total harmonic slack forces each selected tree clock to have a
large lcm.  This is the per-edge consumer used by address compression. -/
theorem low_slack_forces_large_lcm
    {c delta ell : ℝ} (hc : 0 < c) (hd : 0 ≤ delta) (hell : 0 < ell)
    (hedge : 7 * c / (12 * ell) ≤ delta) :
    7 * c ≤ 12 * delta * ell := by
  field_simp at hedge ⊢
  nlinarith

/-- The compact ratio ceiling permits at most 2011 tooth addresses per label
on one complete carrier gap. -/
theorem tooth_address_ceiling_arithmetic :
    (6 : ℕ) * 2345 + 1 < 7 * 2011 := by
  norm_num

theorem finite_stalk_constants :
    (6 : ℕ) * 2011 = 12066 ∧
      5 * 2011 = 10055 ∧
      5 * 2011 + 1 = 10056 ∧
      49 * 2011 = 98539 ∧
      98539 * 10056 = 990908184 := by
  norm_num

/-- Pigeonhole the `1/(49c)` private mass over at most 2011 owner teeth. -/
theorem private_tooth_floor
    {invc privateMass toothMass : ℝ}
    (hprivate : invc / 49 ≤ privateMass)
    (hpartition : privateMass ≤ 2011 * toothMass) :
    invc / 98539 ≤ toothMass := by
  norm_num [show (98539 : ℝ) = 49 * 2011 by norm_num]
  linarith

/-- Cutting one private owner tooth by at most 10055 other teeth leaves a
literal private interval component of the displayed scale. -/
theorem private_interval_stalk_floor
    {invc toothMass stalkMass : ℝ}
    (htooth : invc / 98539 ≤ toothMass)
    (hcomponents : toothMass ≤ 10056 * stalkMass) :
    invc / 990908184 ≤ stalkMass := by
  norm_num [show (990908184 : ℝ) = 98539 * 10056 by norm_num]
  linarith

/-- Private mass `1/(49c)` and owner-tooth capacity `1/(7d)` force at least
`d/(7c)` owner occurrences in the chronological word. -/
theorem private_owner_recurrence
    {invc privateMass teeth d : ℝ} (hd : 0 < d)
    (hprivate : invc / 49 ≤ privateMass)
    (hcapacity : privateMass ≤ teeth / (7 * d)) :
    d * invc / 7 ≤ teeth := by
  have h := le_trans hprivate hcapacity
  have hcross := (div_le_div_iff₀ (by norm_num : (0 : ℝ) < 49)
    (show (0 : ℝ) < 7 * d by positivity)).mp h
  nlinarith

/-- If the same label exits and re-enters around an interior private stalk,
the full safe gap between its consecutive teeth fits strictly inside the
owner tooth, forcing the exact toothpick ratio branch. -/
theorem same_label_private_stalk_forces_toothpick
    {u b : ℝ} (hu : 0 < u) (hb : 0 < b)
    (hcontain : 6 / (7 * u) < 1 / (7 * b)) :
    6 * b < u := by
  have hcross := (div_lt_div_iff₀ (show (0 : ℝ) < 7 * u by positivity)
    (show (0 : ℝ) < 7 * b by positivity)).mp hcontain
  nlinarith

/-- A chronological tooth handoff splits its pair-sum exactly into reciprocal
address drift and actual overlap. -/
theorem handoff_drift_overlap_conservation
    {u v drift overlap : ℝ}
    (hoverlap : 14 * u * v * overlap = u + v - 14 * drift) :
    14 * drift + 14 * u * v * overlap = u + v := by
  linarith

/-- Exact mixed-circuit interface between a chronological path on the fast
speed clocks and one centered blocker edge. -/
theorem mixed_chronological_centered_identity
    {P N s₀ sᵣ n₀ nᵣ delta residual : ℝ}
    (hn₀ : n₀ ≠ 0) (hnᵣ : nᵣ ≠ 0)
    (hdelta : delta = s₀ / n₀ - sᵣ / nᵣ)
    (hresidual : residual = P * s₀ - N * sᵣ) :
    residual = N * nᵣ * delta + (s₀ / n₀) * (P * n₀ - N * nᵣ) := by
  rw [hdelta, hresidual]
  field_simp
  ring

/-- Shifting from fast clocks `s` to centered clocks `c+s` adds the exact
affine carrier drift `(P-N)c`; its sign cannot be discarded. -/
theorem affine_carrier_shift_identity
    {P N c s₀ sᵣ shifted residual : ℝ}
    (hshifted : shifted = P * (c + s₀) - N * (c + sᵣ))
    (hresidual : residual = P * s₀ - N * sᵣ) :
    residual = shifted - (P - N) * c := by
  rw [hshifted, hresidual]
  ring

/-- The mixed identity gives an exact sign/invoice dichotomy: either the
address product descends, or the fast residual pays the whole path drift. -/
theorem mixed_sign_or_drift_invoice
    {residual N nᵣ delta ratio addressGap : ℝ}
    (hidentity : residual = N * nᵣ * delta + ratio * addressGap)
    (hN : 0 ≤ N) (hnᵣ : 0 ≤ nᵣ) (hdelta : 0 ≤ delta)
    (hratio : 0 ≤ ratio) :
    addressGap < 0 ∨ N * nᵣ * delta ≤ residual := by
  by_cases hgap : addressGap < 0
  · exact Or.inl hgap
  · right
    have hgap0 : 0 ≤ addressGap := le_of_not_gt hgap
    have hprod : 0 ≤ ratio * addressGap := mul_nonneg hratio hgap0
    rw [hidentity]
    exact le_add_of_nonneg_right hprod

/-- Concrete guardrail: centered shifted holonomy can be positive while its
unshifted fast-speed residual is negative. -/
theorem affine_sign_reversal_guardrail :
    (9 : ℤ) * (5 + 9) - 1 * (5 + 82) = 39 ∧
      (9 : ℤ) * 9 - 1 * 82 = -1 := by
  norm_num

/-- Common-dilate packets do not make the new debt disappear: simultaneous
scaling cancels between the carrier and the lcm clock. -/
theorem common_dilate_invariance
    {scale c ell : ℝ} (hs : scale ≠ 0) :
    (scale * c) / (scale * ell) = c / ell := by
  field_simp

#print axioms positive_handoff_numerator_quantum
#print axioms gcd_mul_lcm_identity
#print axioms fully_located_hunter_debt
#print axioms fully_located_lcm_debt
#print axioms multiplicity_averaged_debt
#print axioms scale_covariant_lcm_debt
#print axioms low_slack_forces_large_lcm
#print axioms tooth_address_ceiling_arithmetic
#print axioms finite_stalk_constants
#print axioms private_tooth_floor
#print axioms private_interval_stalk_floor
#print axioms private_owner_recurrence
#print axioms same_label_private_stalk_forces_toothpick
#print axioms handoff_drift_overlap_conservation
#print axioms mixed_chronological_centered_identity
#print axioms affine_carrier_shift_identity
#print axioms mixed_sign_or_drift_invoice
#print axioms affine_sign_reversal_guardrail
#print axioms common_dilate_invariance

end SixPrivateLocatedTree
end LRC14
