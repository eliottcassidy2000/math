/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Arithmetic core for THM-1244 slowest-spoke handoff debt

The paper layer extracts a closed slowest-spoke safe component and two
gcd-quantized handoffs from an irredundant open tooth cover.  This module
checks the centered radius invoice, two-label span constants, overlap quantum,
hybrid Hunter rearrangements, and private-mass constants.
-/

namespace LRC14
namespace SlowestSpokeHandoffDebt

/-- Cleared arithmetic behind the centered-component radius floor. -/
theorem centered_radius_floor_cleared
    {c d q x r : ℝ} (hc : 0 ≤ c) (hq : q = c + d)
    (hx : 2 * q * x ≤ 1)
    (hr : 14 * d * q * r = 6 * q - 14 * c * q * x) :
    6 * d - c ≤ 14 * d * q * r := by
  have hmul := mul_le_mul_of_nonneg_left hx (show 0 ≤ 7 * c by positivity)
  rw [hq] at hr hmul ⊢
  nlinarith

/-- The component floor is strictly larger than `5/(14d)` whenever `d>c`. -/
theorem component_floor_strict_cleared {c d : ℝ} (hcd : c < d) :
    5 * (c + d) < 2 * (6 * d - c) := by
  nlinarith

/-- In the close two-label branch the total tooth span is below
`2/(7d1)`, hence below the spoke component floor. -/
theorem close_two_label_span
    {d₁ u v : ℝ} (hd : 0 < d₁) (hdu : d₁ < u) (huv : u < v) :
    1 / (7 * u) + 1 / (7 * v) < 2 / (7 * d₁) ∧
      2 / (7 * d₁) < 5 / (14 * d₁) := by
  have hu0 : 0 < u := lt_trans hd hdu
  have hv0 : 0 < v := lt_trans hu0 huv
  have huInv : (1 : ℝ) / u < 1 / d₁ :=
    one_div_lt_one_div_of_lt hd hdu
  have hvInv : (1 : ℝ) / v < 1 / d₁ :=
    one_div_lt_one_div_of_lt hd (hdu.trans huv)
  constructor
  · field_simp
    nlinarith
  · field_simp
    nlinarith

/-- In the far branch `v>6u`, two possible side extensions still have span
below `4/(21u)`, hence below the spoke component floor. -/
theorem far_two_label_span
    {d₁ u v : ℝ} (hd : 0 < d₁) (hdu : d₁ < u) (hfar : 6 * u < v) :
    1 / (7 * u) + 2 / (7 * v) < 4 / (21 * u) ∧
      4 / (21 * u) < 5 / (14 * d₁) := by
  have hu0 : 0 < u := lt_trans hd hdu
  have hv0 : 0 < v := by nlinarith
  have hvInv : (1 : ℝ) / v < 1 / (6 * u) :=
    one_div_lt_one_div_of_lt (by positivity) hfar
  have huInv : (1 : ℝ) / u < 1 / d₁ :=
    one_div_lt_one_div_of_lt hd hdu
  constructor
  · field_simp
    nlinarith
  · field_simp
    nlinarith

/-- A positive overlap numerator divisible by `g` is at least `g`. -/
theorem positive_overlap_numerator_quantum
    {g numerator : ℕ} (hpos : 0 < numerator) (hdiv : g ∣ numerator) :
    g ≤ numerator := by
  exact Nat.le_of_dvd hpos hdiv

/-- Abstract max-extension form of the five-active-wall Hunter debt. -/
theorem hybrid_hunter_max_extension
    {L harmonic Q extension : ℝ}
    (hcovered : 0 ≥ 2 * L / 7 - 6 * harmonic / 49 + Q + extension) :
    Q + extension ≤ 6 * harmonic / 49 - 2 * L / 7 := by
  linarith

/-- Taking no unlocated extension gives the basic harmonic/seam debt. -/
theorem scalar_handoff_debt
    {L harmonic Q : ℝ}
    (hcovered : 0 ≥ 2 * L / 7 - 6 * harmonic / 49 + Q) :
    7 * L / 3 + 49 * Q / 6 ≤ harmonic := by
  linarith

/-- Insert the exact component and two-seam floors into the scalar debt. -/
theorem explicit_scalar_handoff_debt
    {componentFloor quantumFloor L Q harmonic : ℝ}
    (hL : componentFloor ≤ L) (hQ : quantumFloor ≤ Q)
    (hdebt : 7 * L / 3 + 49 * Q / 6 ≤ harmonic) :
    7 * componentFloor / 3 + 49 * quantumFloor / 6 ≤ harmonic := by
  linarith

/-- Extending the located rank-two forest to a five-label spanning tree
leaves exactly the displayed two-edge gcd debt. -/
theorem spanning_tree_handoff_debt
    {L harmonic Q gcdDebt : ℝ}
    (hcovered : 0 ≥ 4 * L / 13 - 6 * harmonic / 49 + Q -
      13 * gcdDebt / 196) :
    98 * L / 39 + 49 * Q / 6 ≤ harmonic + 13 * gcdDebt / 24 := by
  linarith

theorem private_mass_constant :
    (432 : ℚ) / 1729 - 4 / 21 = 44 / 741 ∧
      (44 : ℚ) / 741 / 5 = 44 / 3705 := by
  norm_num

/-- The global multiple-cover cap leaves the claimed private mass inside K. -/
theorem private_mass_from_component
    {component multiple uniqueMass invc : ℝ}
    (hcomponent : (432 / 1729) * invc < component)
    (hmultiple : multiple ≤ (4 / 21) * invc)
    (hprivate : component - multiple ≤ uniqueMass) :
    (44 / 741) * invc < uniqueMass := by
  have hconst : (432 / 1729 : ℝ) - 4 / 21 = 44 / 741 := by norm_num
  nlinarith

/-- Five blocker labels force one private owner above one fifth of the mass. -/
theorem private_owner_pigeonhole
    {uniqueMass owner invc : ℝ}
    (hprivate : (44 / 741) * invc < uniqueMass)
    (howner : uniqueMass ≤ 5 * owner) :
    (44 / 3705) * invc < owner := by
  nlinarith

#print axioms centered_radius_floor_cleared
#print axioms component_floor_strict_cleared
#print axioms close_two_label_span
#print axioms far_two_label_span
#print axioms positive_overlap_numerator_quantum
#print axioms hybrid_hunter_max_extension
#print axioms scalar_handoff_debt
#print axioms explicit_scalar_handoff_debt
#print axioms spanning_tree_handoff_debt
#print axioms private_mass_constant
#print axioms private_mass_from_component
#print axioms private_owner_pigeonhole

end SlowestSpokeHandoffDebt
end LRC14
