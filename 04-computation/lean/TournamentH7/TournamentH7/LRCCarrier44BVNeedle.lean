import Mathlib.Tactic

/-!
# Arithmetic consumer for THM-1261's carrier-44 BV certificate

The external certificate supplies exact finite and tail load bounds on the 22
reflection-representative phase rows.  This module checks the rational margins,
the fixed-point-free reflection geometry, representative coverage, and the
six-load contradiction.  JSON completeness and analytic integration remain
explicit external providers.
-/

namespace LRC14.Carrier44BVNeedle

def finiteCap : ℚ := 3624398201278 / 21749999997651

def tailCap : ℚ := 51367660709987 / 308249999966709

theorem certified_cap_margins :
    finiteCap < 1 / 6 ∧ tailCap < 1 / 6 ∧ finiteCap ≤ tailCap := by
  norm_num [finiteCap, tailCap]

theorem reflected_gap_left (k : ℚ) :
    1 - (14 * k + 13) / (14 * 44) =
      (14 * (43 - k) + 1) / (14 * 44) := by
  ring

theorem reflected_gap_right (k : ℚ) :
    1 - (14 * k + 1) / (14 * 44) =
      (14 * (43 - k) + 13) / (14 * 44) := by
  ring

theorem phase_reflection_involution (k : ℤ) : 43 - (43 - k) = k := by
  ring

/-- The carrier-44 phase reflection has no integral fixed point. -/
theorem phase_reflection_no_fixed_point (k : ℤ) : k ≠ 43 - k := by
  omega

/-- Every phase modulo 44 is represented by `0,...,21`; phases above 21
reflect into `0,...,21`. -/
theorem phase_has_stored_representative (k : ℤ)
    (hk0 : 0 ≤ k) (hk43 : k ≤ 43) :
    (0 ≤ k ∧ k ≤ 21) ∨
      (0 ≤ 43 - k ∧ 43 - k ≤ 21) := by
  omega

theorem finite_or_tail_load_le_cap (load : ℚ)
    (hload : load ≤ finiteCap ∨ load ≤ tailCap) : load ≤ tailCap := by
  rcases certified_cap_margins with ⟨_, _, hfinite⟩
  rcases hload with hload | hload
  · exact hload.trans hfinite
  · exact hload

theorem six_load_sum_lt_one (load : Fin 6 → ℚ)
    (hload : ∀ i, load i ≤ tailCap) :
    ∑ i, load i < 1 := by
  calc
    ∑ i, load i ≤ ∑ _i : Fin 6, tailCap :=
      Finset.sum_le_sum fun i _hi ↦ hload i
    _ = 6 * tailCap := by norm_num [Fin.sum_univ_succ]
    _ < 1 := by norm_num [tailCap]

theorem six_load_cover_contradiction (load : Fin 6 → ℚ)
    (hload : ∀ i, load i ≤ tailCap)
    (hcover : 1 ≤ ∑ i, load i) : False := by
  have hsum := six_load_sum_lt_one load hload
  linarith

#print axioms certified_cap_margins
#print axioms reflected_gap_left
#print axioms reflected_gap_right
#print axioms phase_reflection_involution
#print axioms phase_reflection_no_fixed_point
#print axioms phase_has_stored_representative
#print axioms finite_or_tail_load_le_cap
#print axioms six_load_sum_lt_one
#print axioms six_load_cover_contradiction

end LRC14.Carrier44BVNeedle
