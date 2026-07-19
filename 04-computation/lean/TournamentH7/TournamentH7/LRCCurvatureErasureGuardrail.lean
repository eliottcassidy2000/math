/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Arithmetic consumers for THM-1239 curvature-erasure guardrails

The paper theorem identifies the seven crack intervals and the containing
danger teeth.  This module checks the self-similar coordinate identities,
the sharp `m=42` inequality ledger, the two-blocker margins, and the explicit
global witness arithmetic.
-/

namespace LRC14
namespace CurvatureErasureGuardrail

/-- The selected slow gap becomes exactly `[-1,11]` under
`u=14a(t-1/7)`. -/
theorem scaled_gap_endpoints (m a : ℚ) (ha : a = 7 * m + 1) (ha0 : a ≠ 0) :
    14 * a * ((14 * m + 1) / (14 * a) - 1 / 7) = -1 ∧
      14 * a * ((14 * m + 13) / (14 * a) - 1 / 7) = 11 := by
  subst a
  constructor <;> field_simp <;> ring

/-- Exact image of the `(a+r)` tooth centred at address `m+1`. -/
theorem scaled_tooth_endpoints (m a r : ℚ)
    (ha : a = 7 * m + 1) (har : a + r ≠ 0) :
    14 * a * ((14 * (m + 1) - 1) / (14 * (a + r)) - 1 / 7) =
        a * (11 - 2 * r) / (a + r) ∧
      14 * a * ((14 * (m + 1) + 1) / (14 * (a + r)) - 1 / 7) =
        a * (13 - 2 * r) / (a + r) := by
  subst a
  constructor <;> field_simp <;> ring

/-- All seven cracks fit strictly inside their odd `z=14a` needles once
`m≥42`; the controlling row is `21/(a+3)<1/14`. -/
theorem one_blocker_margin_ledger {m a : ℝ}
    (hm : 42 ≤ m) (ha : a = 7 * m + 1) :
    6 / (a + 6) < 1 / 14 ∧
      18 / (a + 2) < 1 / 14 ∧
      21 / (a + 3) < 1 / 14 ∧
      20 / (a + 4) < 1 / 14 ∧
      15 / (a + 5) < 1 / 14 ∧
      6 / (a + 6) < 1 / 14 ∧
      11 / (a + 1) < 1 / 14 := by
  subst a
  have h1 : 0 < 7 * m + 2 := by nlinarith
  have h2 : 0 < 7 * m + 3 := by nlinarith
  have h3 : 0 < 7 * m + 4 := by nlinarith
  have h4 : 0 < 7 * m + 5 := by nlinarith
  have h5 : 0 < 7 * m + 6 := by nlinarith
  have h6 : 0 < 7 * m + 7 := by nlinarith
  constructor
  · rw [show 7 * m + 1 + 6 = 7 * m + 7 by ring]
    apply (div_lt_iff₀ h6).2
    norm_num
    nlinarith
  constructor
  · rw [show 7 * m + 1 + 2 = 7 * m + 3 by ring]
    apply (div_lt_iff₀ h2).2
    norm_num
    nlinarith
  constructor
  · rw [show 7 * m + 1 + 3 = 7 * m + 4 by ring]
    apply (div_lt_iff₀ h3).2
    norm_num
    nlinarith
  constructor
  · rw [show 7 * m + 1 + 4 = 7 * m + 5 by ring]
    apply (div_lt_iff₀ h4).2
    norm_num
    nlinarith
  constructor
  · rw [show 7 * m + 1 + 5 = 7 * m + 6 by ring]
    apply (div_lt_iff₀ h5).2
    norm_num
    nlinarith
  constructor
  · rw [show 7 * m + 1 + 6 = 7 * m + 7 by ring]
    apply (div_lt_iff₀ h6).2
    norm_num
    nlinarith
  · rw [show 7 * m + 1 + 1 = 7 * m + 2 by ring]
    apply (div_lt_iff₀ h1).2
    norm_num
    nlinarith

/-- The threshold `a>291` is exactly `m≥42` for `a=7m+1`. -/
theorem sharp_one_blocker_threshold (m : ℕ) :
    42 ≤ m ↔ 291 < 7 * m + 1 := by
  omega

/-- Cleared integer margins for the two-blocker construction are all positive
from `m=8`. -/
theorem two_blocker_margin_ledger (m : ℕ) (hm : 8 ≤ m) :
    0 < 14 * m - 5 ∧ 0 < 14 * m - 75 ∧
      0 < 14 * m - 84 ∧ 0 < 14 * m - 111 ∧
      0 < 14 * m - 110 ∧ 0 < 14 * m - 81 ∧
      0 < 14 * m - 24 ∧
      0 < 27 ∧ 0 < 70 ∧ 0 < 85 ∧ 0 < 72 ∧ 0 < 31 := by
  omega

/-- Exact deletion-quartet margin on the selected gap. -/
theorem nonbad_quartet_margin (m a : ℚ)
    (ha : a = 7 * m + 1) (ha0 : a ≠ 0) :
    1 - 3 * ((14 * m + 13) / (14 * a)) - 2 / 7 =
      (28 * m - 29) / (14 * a) := by
  subst a
  have ha0' : 1 + m * 7 ≠ 0 := by
    intro h
    apply ha0
    nlinarith
  field_simp [ha0']
  ring

/-- Exact residues of the explicit thirteen-speed global witness. -/
theorem explicit_global_witness_residues :
    (1 * 44) % 199 = 44 ∧ (2 * 44) % 199 = 88 ∧
      (3 * 44) % 199 = 132 ∧ (4 * 44) % 199 = 176 ∧
      (295 * 44) % 199 = 45 ∧ (296 * 44) % 199 = 89 ∧
      (297 * 44) % 199 = 133 ∧ (298 * 44) % 199 = 177 ∧
      (299 * 44) % 199 = 22 ∧ (300 * 44) % 199 = 66 ∧
      (301 * 44) % 199 = 110 ∧ (302 * 44) % 199 = 154 ∧
      (4130 * 44) % 199 = 33 := by
  norm_num

/-- The minimum residue distance `22/199` clears radius `1/14`. -/
theorem explicit_global_witness_margin :
    (1 / 14 : ℚ) < 22 / 199 := by
  norm_num

/-- A crack-blocking address window has length below one whenever the crack
itself is shorter than one danger tooth. -/
theorem address_window_unique_scale
    {z length : ℝ} (hz : 0 ≤ z) (hlength : 0 ≤ length)
    (hshort : z * length < 1 / 7) :
    0 < 1 / 7 - z * length ∧ 1 / 7 - z * length < 1 := by
  constructor <;> nlinarith

#print axioms scaled_gap_endpoints
#print axioms scaled_tooth_endpoints
#print axioms one_blocker_margin_ledger
#print axioms sharp_one_blocker_threshold
#print axioms two_blocker_margin_ledger
#print axioms nonbad_quartet_margin
#print axioms explicit_global_witness_residues
#print axioms explicit_global_witness_margin
#print axioms address_window_unique_scale

end CurvatureErasureGuardrail
end LRC14
