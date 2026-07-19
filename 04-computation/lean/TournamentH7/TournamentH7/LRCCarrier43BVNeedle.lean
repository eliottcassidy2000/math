import Mathlib.Tactic

/-!
# Arithmetic consumer for THM-1259's carrier-43 BV certificate

The external certificate supplies exact finite and tail load bounds on the 22
reflection-representative phase rows.  This module checks the rational margins,
the reflection geometry including the fixed phase 21, representative coverage,
and the six-load contradiction.  JSON and analytic integration remain explicit
external providers.
-/

namespace LRC14.Carrier43BVNeedle

def finiteCap : ℚ := 2331708819707 / 13999999998460

def tailCap : ℚ := 75385062494229 / 452374999952953

theorem certified_cap_margins :
    finiteCap < 1 / 6 ∧ tailCap < 1 / 6 ∧ finiteCap ≤ tailCap := by
  norm_num [finiteCap, tailCap]

theorem reflected_gap_left (k : ℚ) :
    1 - (14 * k + 13) / (14 * 43) =
      (14 * (42 - k) + 1) / (14 * 43) := by
  ring

theorem reflected_gap_right (k : ℚ) :
    1 - (14 * k + 1) / (14 * 43) =
      (14 * (42 - k) + 13) / (14 * 43) := by
  ring

theorem phase_reflection_involution (k : ℤ) : 42 - (42 - k) = k := by
  ring

/-- Phase 21 is the unique integral fixed point of `k ↦ 42-k`. -/
theorem phase_reflection_fixed_iff (k : ℤ) : k = 42 - k ↔ k = 21 := by
  omega

/-- Every phase modulo 43 is represented by `0,...,21`; phases above 21
reflect into `0,...,20`. -/
theorem phase_has_stored_representative (k : ℤ)
    (hk0 : 0 ≤ k) (hk42 : k ≤ 42) :
    (0 ≤ k ∧ k ≤ 21) ∨
      (0 ≤ 42 - k ∧ 42 - k ≤ 20) := by
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
#print axioms finite_or_tail_load_le_cap
#print axioms six_load_sum_lt_one
#print axioms six_load_cover_contradiction

end LRC14.Carrier43BVNeedle
