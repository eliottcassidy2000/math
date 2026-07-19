import Mathlib.Tactic

/-!
# Arithmetic consumer for the carrier-41 BV needle certificate (THM-1255)

The external exact certificate supplies, for each reflected phase row, a
probability density whose finite-speed loads are at most `finiteCap` and whose
analytic BV-tail loads are at most `tailCap`.  This module checks the two strict
rational margins, reflection of the carrier-41 gap, and the six-load union-bound
consumer.

The 200-bin exact integrations, completeness of the 21 certificate rows, and
the Stieltjes/BV estimate are deliberately external providers.  In particular,
this file does not pretend to formalize the JSON certificate or integration
theory.
-/

namespace LRC14.Carrier41BVNeedle

/-- Largest exact finite-speed load in the 21-row certificate. -/
def finiteCap : ℚ := 5498776467337 / 32999999996568

/-- Largest exact analytic-tail load in the 21-row certificate. -/
def tailCap : ℚ := 957696499918301 / 5746999999316107

/-- The two independently replayed maxima are both strictly below `1/6`; the
tail maximum is also the common cap used by the abstract consumer. -/
theorem certified_cap_margins :
    finiteCap < 1 / 6 ∧ tailCap < 1 / 6 ∧ finiteCap ≤ tailCap := by
  norm_num [finiteCap, tailCap]

/-- Reflection about `1/2` sends the right endpoint of `G_k(41)` to the left
endpoint of `G_(40-k)(41)`. -/
theorem reflected_gap_left (k : ℚ) :
    1 - (14 * k + 13) / (14 * 41) =
      (14 * (40 - k) + 1) / (14 * 41) := by
  ring

/-- Reflection about `1/2` sends the left endpoint of `G_k(41)` to the right
endpoint of `G_(40-k)(41)`. -/
theorem reflected_gap_right (k : ℚ) :
    1 - (14 * k + 1) / (14 * 41) =
      (14 * (40 - k) + 13) / (14 * 41) := by
  ring

/-- Any load controlled by either exact certificate layer is strictly below
one sixth. -/
theorem finite_or_tail_load_lt_one_sixth (load : ℚ)
    (hload : load ≤ finiteCap ∨ load ≤ tailCap) : load < 1 / 6 := by
  rcases certified_cap_margins with ⟨hfinite, htail, _⟩
  rcases hload with hload | hload <;> linarith

/-- Six loads bounded by the common carrier-41 cap have total mass strictly
less than one.  This is the arithmetic endpoint of the union-bound proof. -/
theorem six_load_sum_lt_one (load : Fin 6 → ℚ)
    (hload : ∀ i, load i ≤ tailCap) :
    ∑ i, load i < 1 := by
  calc
    ∑ i, load i ≤ ∑ _i : Fin 6, tailCap :=
      Finset.sum_le_sum fun i _hi ↦ hload i
    _ = 6 * tailCap := by norm_num [Fin.sum_univ_succ]
    _ < 1 := by norm_num [tailCap]

/-- A putative cover forces dual mass at least one, contradicting the six
certified load bounds. -/
theorem six_load_cover_contradiction (load : Fin 6 → ℚ)
    (hload : ∀ i, load i ≤ tailCap)
    (hcover : 1 ≤ ∑ i, load i) : False := by
  have hsum := six_load_sum_lt_one load hload
  linarith

#print axioms certified_cap_margins
#print axioms reflected_gap_left
#print axioms reflected_gap_right
#print axioms finite_or_tail_load_lt_one_sixth
#print axioms six_load_sum_lt_one
#print axioms six_load_cover_contradiction

end LRC14.Carrier41BVNeedle
