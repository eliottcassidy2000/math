import Mathlib.Tactic

/-!
# Arithmetic consumer for THM-1257's carrier-42 BV certificate

The external exact certificate supplies, on each of 21 reflected phase rows,
a probability density whose finite-speed loads are at most `finiteCap` and
whose BV-tail loads are at most `tailCap`.  This module checks the exact
rational margins, the full reflection-pair reduction, and the six-load
union-bound contradiction.

The JSON completeness, 200-bin integrations, and analytic Stieltjes/BV bound
remain explicit external providers.  No such statement is hidden here.
-/

namespace LRC14.Carrier42BVNeedle

/-- Largest exact finite-speed load in the 21-row certificate. -/
def finiteCap : ℚ := 28499599523095 / 170999999979651

/-- Largest exact analytic-tail load in the 21-row certificate. -/
def tailCap : ℚ := 65421208326679 / 392583333284653

/-- Both exact maxima are strictly below one sixth; the finite maximum is the
larger common cap for the abstract six-load consumer. -/
theorem certified_cap_margins :
    finiteCap < 1 / 6 ∧ tailCap < 1 / 6 ∧ tailCap ≤ finiteCap := by
  norm_num [finiteCap, tailCap]

/-- Reflection about `1/2` sends the right endpoint of `G_k(42)` to the left
endpoint of `G_(41-k)(42)`. -/
theorem reflected_gap_left (k : ℚ) :
    1 - (14 * k + 13) / (14 * 42) =
      (14 * (41 - k) + 1) / (14 * 42) := by
  ring

/-- Reflection about `1/2` sends the left endpoint of `G_k(42)` to the right
endpoint of `G_(41-k)(42)`. -/
theorem reflected_gap_right (k : ℚ) :
    1 - (14 * k + 1) / (14 * 42) =
      (14 * (41 - k) + 13) / (14 * 42) := by
  ring

/-- The phase reflection is an involution. -/
theorem phase_reflection_involution (k : ℤ) : 41 - (41 - k) = k := by
  ring

/-- Because 41 is odd, no integral phase is fixed by `k ↦ 41-k`. -/
theorem phase_reflection_no_fixed_point (k : ℤ) : k ≠ 41 - k := by
  omega

/-- Every reduced phase `0<=k<=41` is either a stored representative
`0<=k<=20` or reflects to one.  Thus there are 21, not 22, representative
rows. -/
theorem phase_has_stored_representative (k : ℤ)
    (hk0 : 0 ≤ k) (hk41 : k ≤ 41) :
    (0 ≤ k ∧ k ≤ 20) ∨
      (0 ≤ 41 - k ∧ 41 - k ≤ 20) := by
  omega

/-- Any load supplied by either exact certificate layer lies below the common
finite cap. -/
theorem finite_or_tail_load_le_cap (load : ℚ)
    (hload : load ≤ finiteCap ∨ load ≤ tailCap) : load ≤ finiteCap := by
  rcases certified_cap_margins with ⟨_, _, htail⟩
  rcases hload with hload | hload
  · exact hload
  · exact hload.trans htail

/-- Six loads bounded by the common THM-1257 cap have total mass below one. -/
theorem six_load_sum_lt_one (load : Fin 6 → ℚ)
    (hload : ∀ i, load i ≤ finiteCap) :
    ∑ i, load i < 1 := by
  calc
    ∑ i, load i ≤ ∑ _i : Fin 6, finiteCap :=
      Finset.sum_le_sum fun i _hi ↦ hload i
    _ = 6 * finiteCap := by norm_num [Fin.sum_univ_succ]
    _ < 1 := by norm_num [finiteCap]

/-- A putative six-comb cover forces dual mass at least one and contradicts
the exact carrier-42 caps. -/
theorem six_load_cover_contradiction (load : Fin 6 → ℚ)
    (hload : ∀ i, load i ≤ finiteCap)
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

end LRC14.Carrier42BVNeedle
