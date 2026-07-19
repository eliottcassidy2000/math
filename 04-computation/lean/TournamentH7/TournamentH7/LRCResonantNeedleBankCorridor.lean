/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Arithmetic core for THM-1246 resonant needle-bank corridor

The paper theorem interprets the following rational interval identities as
circle-distance statements.  This module checks the corridor width, harmonic
offset bounds, projective-band cap, integer-family condition, and reciprocal
endpoint ladder without `native_decide` or proof placeholders.
-/

namespace LRC14
namespace ResonantNeedleBankCorridor

def corridorLower : ℚ := 15 / 196
def corridorUpper (H : ℚ) : ℚ := (14 * H + 13) / (196 * H)
def projectiveCap (H : ℚ) : ℚ := 182 * H / (14 * H + 13)

theorem corridor_width (H : ℚ) (hH : H ≠ 0) :
    corridorUpper H - corridorLower = (13 - H) / (196 * H) := by
  simp only [corridorUpper, corridorLower]
  field_simp
  ring

theorem corridor_positive {H : ℚ} (hH : 0 < H) (hH13 : H < 13) :
    0 < corridorUpper H - corridorLower := by
  rw [corridor_width H (ne_of_gt hH)]
  positivity

theorem lower_bank_boundary :
    14 * corridorLower = 1 + (1 : ℚ) / 14 := by
  norm_num [corridorLower]

theorem upper_bank_boundary (H : ℚ) (hH : H ≠ 0) :
    14 * H * corridorUpper H = H + (13 : ℚ) / 14 := by
  simp only [corridorUpper]
  field_simp
  ring

theorem bank_lower_offset (h : ℚ) :
    14 * h * corridorLower - h = h / 14 := by
  norm_num [corridorLower]
  ring

theorem bank_upper_offset (H h : ℚ) (hH : H ≠ 0) :
    14 * h * corridorUpper H - h = 13 * h / (14 * H) := by
  simp only [corridorUpper]
  field_simp
  ring

/-- Every harmonic offset lies in the closed LRC-safe band. -/
theorem harmonic_offset_safe
    {H h x : ℚ} (hH : 0 < H) (hh1 : 1 ≤ h) (hhH : h ≤ H)
    (hxL : corridorLower ≤ x) (hxU : x ≤ corridorUpper H) :
    (1 : ℚ) / 14 ≤ 14 * h * x - h ∧
      14 * h * x - h ≤ (13 : ℚ) / 14 := by
  have hh0 : 0 ≤ 14 * h := by positivity
  have hlower : h / 14 ≤ 14 * h * x - h := by
    rw [← bank_lower_offset h]
    exact sub_le_sub_right (mul_le_mul_of_nonneg_left hxL hh0) h
  have hupper : 14 * h * x - h ≤ 13 * h / (14 * H) := by
    rw [← bank_upper_offset H h (ne_of_gt hH)]
    exact sub_le_sub_right (mul_le_mul_of_nonneg_left hxU hh0) h
  constructor
  · nlinarith
  · have hden : (0 : ℚ) < 14 * H := by positivity
    have hratio : 13 * h / (14 * H) ≤ (13 : ℚ) / 14 := by
      apply (div_le_iff₀ hden).2
      nlinarith
    exact hupper.trans hratio

theorem cap_times_upper (H : ℚ) (hH : 0 < H) :
    projectiveCap H * corridorUpper H = (13 : ℚ) / 14 := by
  simp only [projectiveCap, corridorUpper]
  have hden : 14 * H + 13 ≠ 0 := by nlinarith
  field_simp
  ring

/-- Every projective speed ratio between one and the cap stays in the same
safe band throughout the corridor. -/
theorem projective_band_safe
    {H x lam : ℚ} (hH : 0 < H)
    (hxL : corridorLower ≤ x) (hxU : x ≤ corridorUpper H)
    (hlam1 : 1 ≤ lam) (hlamU : lam ≤ projectiveCap H) :
    (1 : ℚ) / 14 < lam * x ∧ lam * x ≤ (13 : ℚ) / 14 := by
  have hx0 : 0 ≤ x := by
    have : (0 : ℚ) < corridorLower := by norm_num [corridorLower]
    exact le_trans this.le hxL
  have hcap0 : 0 ≤ projectiveCap H := by
    simp only [projectiveCap]
    positivity
  have hlower : corridorLower ≤ lam * x := by
    have hxlam : x ≤ lam * x := by
      nlinarith [mul_nonneg (sub_nonneg.mpr hlam1) hx0]
    exact hxL.trans hxlam
  have hupper : lam * x ≤ projectiveCap H * corridorUpper H := by
    calc
      lam * x ≤ projectiveCap H * x :=
        mul_le_mul_of_nonneg_right hlamU hx0
      _ ≤ projectiveCap H * corridorUpper H :=
        mul_le_mul_of_nonneg_left hxU hcap0
  constructor
  · have : (1 : ℚ) / 14 < corridorLower := by norm_num [corridorLower]
    exact this.trans_le hlower
  · rwa [cap_times_upper H hH] at hupper

/-- The exact integral condition for the consecutive block is automatic for
all `1≤H≤12` once `a≥2`. -/
theorem integer_family_condition
    {H a : ℤ} (hH1 : 1 ≤ H) (hH12 : H ≤ 12) (ha : 2 ≤ a) :
    (12 - H) * (14 * H + 13) ≤ (168 * H - 13) * a := by
  have hcoef : 0 ≤ 168 * H - 13 := by omega
  have hmul : (168 * H - 13) * 2 ≤ (168 * H - 13) * a := by
    exact mul_le_mul_of_nonneg_left ha hcoef
  nlinarith [sq_nonneg (H - 1)]

theorem reciprocal_endpoint_step
    (H : ℚ) (hH : H ≠ 0) (hHm1 : H - 1 ≠ 0) :
    corridorUpper (H - 1) - corridorUpper H =
      13 / (196 * H * (H - 1)) := by
  simp only [corridorUpper]
  field_simp
  ring

theorem h6_width :
    corridorUpper 6 - corridorLower = (1 : ℚ) / 168 := by
  norm_num [corridorUpper, corridorLower]

#print axioms corridor_width
#print axioms corridor_positive
#print axioms lower_bank_boundary
#print axioms upper_bank_boundary
#print axioms bank_lower_offset
#print axioms bank_upper_offset
#print axioms harmonic_offset_safe
#print axioms cap_times_upper
#print axioms projective_band_safe
#print axioms integer_family_condition
#print axioms reciprocal_endpoint_step
#print axioms h6_width

end ResonantNeedleBankCorridor
end LRC14

