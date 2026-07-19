import Mathlib.Tactic

/-!
# Arithmetic consumer for THM-1265's carrier-45 BV certificate

The external certificate supplies exact finite and tail load bounds on the 23
reflection-representative phase rows.  This module checks the rational margins,
the fixed phase and reflection geometry, representative coverage, and the
six-load contradiction.  JSON completeness and analytic integration remain
explicit external providers.
-/

namespace LRC14.Carrier45BVNeedle

def finiteCap : ℚ := 12991498270615 / 77999999992356

def tailCap : ℚ := 294191299975319 / 1765399999814633

theorem certified_cap_margins :
    finiteCap < 1 / 6 ∧ tailCap < 1 / 6 ∧ finiteCap ≤ tailCap := by
  norm_num [finiteCap, tailCap]

theorem reflected_gap_left (k : ℚ) :
    1 - (14 * k + 13) / (14 * 45) =
      (14 * (44 - k) + 1) / (14 * 45) := by
  ring

theorem reflected_gap_right (k : ℚ) :
    1 - (14 * k + 1) / (14 * 45) =
      (14 * (44 - k) + 13) / (14 * 45) := by
  ring

theorem phase_reflection_involution (k : ℤ) : 44 - (44 - k) = k := by
  ring

/-- Carrier 45 has the unique integral fixed phase `k=22`. -/
theorem phase_reflection_fixed_iff (k : ℤ) : k = 44 - k ↔ k = 22 := by
  constructor <;> omega

/-- Every phase modulo 45 is represented by `0,...,22`; phases above 22
reflect into `0,...,21`, while phase 22 is fixed. -/
theorem phase_has_stored_representative (k : ℤ)
    (hk0 : 0 ≤ k) (hk44 : k ≤ 44) :
    (0 ≤ k ∧ k ≤ 22) ∨
      (0 ≤ 44 - k ∧ 44 - k ≤ 22) := by
  omega

theorem fixed_phase_is_represented : (0 : ℤ) ≤ 22 ∧ (22 : ℤ) ≤ 22 := by
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
#print axioms phase_reflection_fixed_iff
#print axioms phase_has_stored_representative
#print axioms fixed_phase_is_represented
#print axioms finite_or_tail_load_le_cap
#print axioms six_load_sum_lt_one
#print axioms six_load_cover_contradiction

end LRC14.Carrier45BVNeedle
