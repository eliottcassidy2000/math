import Mathlib

/-!
# Six-comb near-tiling curvature (THM-1219)

For `a = 7m`, the complete `m`-th safe gap of the `a`-comb is almost tiled
by the six danger teeth belonging to `a+1,...,a+6` and centred at `m+1`.
After replacing `14m` by `2a`, consecutive endpoint subtraction gives five
internal defects with numerators `11,9,7,5,3`, and the terminal defect has
numerator `13`.  Their total is therefore of order `a⁻²`, not `a⁻¹`.

This file checks the rational endpoint identities, their full mod-seven
form, their factorization through active-pair slack `14D-(x+y)`, the uniform
denominator sandwich, the relative slow-gap sandwich, and the threshold
beyond which the survivor is smaller than the five-comb floor `1/(49a)`.  The geometric facts
that these are the unique six teeth meeting the slow gap and that their
chronological union has exactly these complement components remain the
explicit interval-arrangement provider in THM-1219.
-/

namespace LRC14.SixCombNearTiling

/-- Exact total length of the six complement components.  Terms are ordered
from the right boundary backwards through the five adjacent handoffs. -/
def survivorLength (a : ℚ) : ℚ :=
  13 / (14 * a * (a + 1)) +
  11 / (14 * (a + 1) * (a + 2)) +
  9 / (14 * (a + 2) * (a + 3)) +
  7 / (14 * (a + 3) * (a + 4)) +
  5 / (14 * (a + 4) * (a + 5)) +
  3 / (14 * (a + 5) * (a + 6))

/-- Physical length of one complete safe gap of the slow `a`-comb. -/
def slowGapLength (a : ℚ) : ℚ :=
  6 / (7 * a)

/-- Adjacent handoff subtraction.  In the application `r=1,...,5` and
`2a+14 = 14(m+1)`. -/
theorem adjacent_gap_identity (a r : ℚ)
    (har : a + r ≠ 0) (har1 : a + r + 1 ≠ 0) :
    (2 * a + 13) / (14 * (a + r)) -
        (2 * a + 15) / (14 * (a + r + 1)) =
      (13 - 2 * r) / (14 * (a + r) * (a + r + 1)) := by
  field_simp
  ring

/-- The left endpoint of the `a+6` tooth precedes the slow-gap endpoint by
exactly `6/[14a(a+6)]`; hence it creates no left boundary survivor. -/
theorem left_boundary_identity (a : ℚ) (ha : a ≠ 0) (ha6 : a + 6 ≠ 0) :
    (2 * a + 13) / (14 * (a + 6)) -
        (2 * a + 1) / (14 * a) =
      -6 / (14 * a * (a + 6)) := by
  field_simp
  ring

/-- The other endpoint of the `a+6` tooth extends past the slow-gap endpoint
by `(2a-6)/[14a(a+6)]`, positive once `a>3`. -/
theorem left_boundary_inside_identity
    (a : ℚ) (ha : a ≠ 0) (ha6 : a + 6 ≠ 0) :
    (2 * a + 15) / (14 * (a + 6)) -
        (2 * a + 1) / (14 * a) =
      (2 * a - 6) / (14 * a * (a + 6)) := by
  field_simp
  ring

/-- Exact right-boundary defect between the slow gap and the `a+1` tooth. -/
theorem terminal_gap_identity (a : ℚ) (ha : a ≠ 0) (ha1 : a + 1 ≠ 0) :
    (2 * a + 13) / (14 * a) -
        (2 * a + 15) / (14 * (a + 1)) =
      13 / (14 * a * (a + 1)) := by
  field_simp
  ring

/-! ## Full residue spectrum and active-pair coordinates -/

/-- Adjacent handoff for `a=7m+s`; eliminating `m` leaves the residue-class
correction `-2s`. -/
theorem residue_adjacent_gap_identity (a s r : ℚ)
    (har : a + r ≠ 0) (har1 : a + r + 1 ≠ 0) :
    (2 * a - 2 * s + 13) / (14 * (a + r)) -
        (2 * a - 2 * s + 15) / (14 * (a + r + 1)) =
      (13 - 2 * s - 2 * r) / (14 * (a + r) * (a + r + 1)) := by
  field_simp
  ring

/-- Signed left-boundary handoff in residue class `s`. -/
theorem residue_left_boundary_identity (a s : ℚ)
    (ha : a ≠ 0) (ha6 : a + 6 ≠ 0) :
    (2 * a - 2 * s + 13) / (14 * (a + 6)) -
        (2 * a - 2 * s + 1) / (14 * a) =
      (12 * s - 6) / (14 * a * (a + 6)) := by
  field_simp
  ring

/-- Signed terminal handoff in residue class `s`. -/
theorem residue_terminal_gap_identity (a s : ℚ)
    (ha : a ≠ 0) (ha1 : a + 1 ≠ 0) :
    (2 * a - 2 * s + 13) / (14 * a) -
        (2 * a - 2 * s + 15) / (14 * (a + 1)) =
      (13 - 2 * s) / (14 * a * (a + 1)) := by
  field_simp
  ring

/-- A signed tooth handoff is exactly the located-maximizer pair margin
`D/(x+y)-1/14`, rescaled from distance to time. -/
theorem active_pair_handoff_identity (x y D : ℚ)
    (hx : x ≠ 0) (hy : y ≠ 0) (hsum : x + y ≠ 0) :
    (x + y) / (x * y) * (D / (x + y) - 1 / 14) =
      (14 * D - (x + y)) / (14 * x * y) := by
  field_simp

/-- The three kinds of handoff slack on `a=7m+s`: internal, terminal, and
left boundary. -/
theorem residue_active_slacks (m s r : ℚ) :
    14 * (m + 1) -
        ((7 * m + s + r) + (7 * m + s + r + 1)) = 13 - 2 * s - 2 * r ∧
      14 * (m + 1) -
        ((7 * m + s) + (7 * m + s + 1)) = 13 - 2 * s ∧
      14 * (m + s) -
        ((7 * m + s) + (7 * m + s + 6)) = 12 * s - 6 := by
  constructor
  · ring
  constructor <;> ring

/-- Sums of the positive active-pair slacks in residue classes `0,...,6`.
The minimum charge `42` is attained in residue class one. -/
theorem residue_curvature_charge_table :
    (11 + 9 + 7 + 5 + 3 + 13 : ℚ) = 48 ∧
    (6 + 9 + 7 + 5 + 3 + 1 + 11 : ℚ) = 42 ∧
    (18 + 7 + 5 + 3 + 1 + 9 : ℚ) = 43 ∧
    (30 + 5 + 3 + 1 + 7 : ℚ) = 46 ∧
    (42 + 3 + 1 + 5 : ℚ) = 51 ∧
    (54 + 1 + 3 : ℚ) = 58 ∧
    (66 + 1 : ℚ) = 67 := by
  norm_num

theorem curvature_numerator_sum :
    (3 + 5 + 7 + 9 + 11 + 13 : ℚ) = 48 := by
  norm_num

/-- Every positive curvature term with `0≤r≤5` lies between the common
denominator envelopes used in THM-1219. -/
theorem positive_term_bounds
    (a r n : ℚ) (ha : 0 < a) (hr0 : 0 ≤ r) (hr5 : r ≤ 5) (hn : 0 < n) :
    n / (14 * (a + 6) ^ 2) ≤
        n / (14 * (a + r) * (a + r + 1)) ∧
      n / (14 * (a + r) * (a + r + 1)) <
        n / (14 * a ^ 2) := by
  have har0 : 0 ≤ a + r := by linarith
  have har : 0 < a + r := lt_of_lt_of_le ha (by linarith)
  have har1 : 0 < a + r + 1 := by linarith
  have ha6 : 0 < a + 6 := by linarith
  have hleft : a + r ≤ a + 6 := by linarith
  have hright : a + r + 1 ≤ a + 6 := by linarith
  have hproduct_le : (a + r) * (a + r + 1) ≤ (a + 6) * (a + 6) :=
    mul_le_mul hleft hright (le_of_lt har1) (le_of_lt ha6)
  have hpositive_piece : 0 < a * (2 * r + 1) := by
    exact mul_pos ha (by linarith)
  have hrpiece : 0 ≤ r * (r + 1) :=
    mul_nonneg hr0 (by linarith)
  have hproduct_gt : a * a < (a + r) * (a + r + 1) := by
    nlinarith
  have hdenLower : 0 < 14 * (a + 6) ^ 2 := by positivity
  have hdenActual : 0 < 14 * (a + r) * (a + r + 1) := by positivity
  have hdenUpper : 0 < 14 * a ^ 2 := by positivity
  constructor
  · apply (div_le_div_iff₀ hdenLower hdenActual).2
    exact mul_le_mul_of_nonneg_left (by nlinarith [hproduct_le]) (le_of_lt hn)
  · apply (div_lt_div_iff₀ hdenActual hdenUpper).2
    exact mul_lt_mul_of_pos_left (by nlinarith [hproduct_gt]) hn

/-- Summing the six exact curvature terms gives the sharp denominator
sandwich `24/[7(a+6)^2] ≤ S(a) < 24/(7a^2)`. -/
theorem survivorLength_bounds (a : ℚ) (ha : 0 < a) :
    24 / (7 * (a + 6) ^ 2) ≤ survivorLength a ∧
      survivorLength a < 24 / (7 * a ^ 2) := by
  have h0 := positive_term_bounds a 0 13 ha (by norm_num) (by norm_num) (by norm_num)
  have h1 := positive_term_bounds a 1 11 ha (by norm_num) (by norm_num) (by norm_num)
  have h2 := positive_term_bounds a 2 9 ha (by norm_num) (by norm_num) (by norm_num)
  have h3 := positive_term_bounds a 3 7 ha (by norm_num) (by norm_num) (by norm_num)
  have h4 := positive_term_bounds a 4 5 ha (by norm_num) (by norm_num) (by norm_num)
  have h5 := positive_term_bounds a 5 3 ha (by norm_num) (by norm_num) (by norm_num)
  norm_num at h0
  have ha12 : a + 1 + 1 = a + 2 := by ring
  have ha23 : a + 2 + 1 = a + 3 := by ring
  have ha34 : a + 3 + 1 = a + 4 := by ring
  have ha45 : a + 4 + 1 = a + 5 := by ring
  have ha56 : a + 5 + 1 = a + 6 := by ring
  have h1lo : 11 / (14 * (a + 6) ^ 2) ≤
      11 / (14 * (a + 1) * (a + 2)) := by simpa [ha12] using h1.1
  have h2lo : 9 / (14 * (a + 6) ^ 2) ≤
      9 / (14 * (a + 2) * (a + 3)) := by simpa [ha23] using h2.1
  have h3lo : 7 / (14 * (a + 6) ^ 2) ≤
      7 / (14 * (a + 3) * (a + 4)) := by simpa [ha34] using h3.1
  have h4lo : 5 / (14 * (a + 6) ^ 2) ≤
      5 / (14 * (a + 4) * (a + 5)) := by simpa [ha45] using h4.1
  have h5lo : 3 / (14 * (a + 6) ^ 2) ≤
      3 / (14 * (a + 5) * (a + 6)) := by simpa [ha56] using h5.1
  have h1hi : 11 / (14 * (a + 1) * (a + 2)) < 11 / (14 * a ^ 2) := by
    simpa [ha12] using h1.2
  have h2hi : 9 / (14 * (a + 2) * (a + 3)) < 9 / (14 * a ^ 2) := by
    simpa [ha23] using h2.2
  have h3hi : 7 / (14 * (a + 3) * (a + 4)) < 7 / (14 * a ^ 2) := by
    simpa [ha34] using h3.2
  have h4hi : 5 / (14 * (a + 4) * (a + 5)) < 5 / (14 * a ^ 2) := by
    simpa [ha45] using h4.2
  have h5hi : 3 / (14 * (a + 5) * (a + 6)) < 3 / (14 * a ^ 2) := by
    simpa [ha56] using h5.2
  constructor
  · calc
      24 / (7 * (a + 6) ^ 2) =
          13 / (14 * (a + 6) ^ 2) +
          11 / (14 * (a + 6) ^ 2) +
          9 / (14 * (a + 6) ^ 2) +
          7 / (14 * (a + 6) ^ 2) +
          5 / (14 * (a + 6) ^ 2) +
          3 / (14 * (a + 6) ^ 2) := by
        have ha6 : a + 6 ≠ 0 := ne_of_gt (by linarith)
        field_simp
        ring
      _ ≤ survivorLength a := by
        unfold survivorLength
        linarith [h0.1, h1lo, h2lo, h3lo, h4lo, h5lo]
  · calc
      survivorLength a <
          13 / (14 * a ^ 2) + 11 / (14 * a ^ 2) +
          9 / (14 * a ^ 2) + 7 / (14 * a ^ 2) +
          5 / (14 * a ^ 2) + 3 / (14 * a ^ 2) := by
        unfold survivorLength
        linarith [h0.2, h1hi, h2hi, h3hi, h4hi, h5hi]
      _ = 24 / (7 * a ^ 2) := by
        have ha0 : a ≠ 0 := ne_of_gt ha
        field_simp
        ring

/-- Dividing by the complete slow-gap length yields the relative survivor
sandwich `4a/(a+6)^2 ≤ S/|G| < 4/a`. -/
theorem relative_survivor_bounds (a : ℚ) (ha : 0 < a) :
    4 * a / (a + 6) ^ 2 ≤ survivorLength a / slowGapLength a ∧
      survivorLength a / slowGapLength a < 4 / a := by
  have hslow : 0 < slowGapLength a := by
    unfold slowGapLength
    positivity
  have hbounds := survivorLength_bounds a ha
  constructor
  · apply (le_div_iff₀ hslow).2
    calc
      4 * a / (a + 6) ^ 2 * slowGapLength a =
          24 / (7 * (a + 6) ^ 2) := by
        unfold slowGapLength
        field_simp
        ring
      _ ≤ survivorLength a := hbounds.1
  · apply (div_lt_iff₀ hslow).2
    calc
      survivorLength a < 24 / (7 * a ^ 2) := hbounds.2
      _ = 4 / a * slowGapLength a := by
        unfold slowGapLength
        field_simp
        ring

/-- Once `a>168`, the `O(a⁻²)` survivor is already below the universal
five-comb physical floor `1/(49a)`. -/
theorem survivor_below_five_comb_floor
    (a : ℚ) (ha : 168 < a) :
    survivorLength a < 1 / (49 * a) := by
  have ha0 : 0 < a := by linarith
  have hupper := (survivorLength_bounds a ha0).2
  have hcomparison : 24 / (7 * a ^ 2) < 1 / (49 * a) := by
    have hdenLeft : 0 < 7 * a ^ 2 := by positivity
    have hdenRight : 0 < 49 * a := by positivity
    apply (div_lt_div_iff₀ hdenLeft hdenRight).2
    nlinarith [sq_pos_of_pos ha0]
  exact hupper.trans hcomparison

#print axioms adjacent_gap_identity
#print axioms left_boundary_identity
#print axioms left_boundary_inside_identity
#print axioms terminal_gap_identity
#print axioms residue_adjacent_gap_identity
#print axioms residue_left_boundary_identity
#print axioms residue_terminal_gap_identity
#print axioms active_pair_handoff_identity
#print axioms residue_active_slacks
#print axioms residue_curvature_charge_table
#print axioms curvature_numerator_sum
#print axioms positive_term_bounds
#print axioms survivorLength_bounds
#print axioms relative_survivor_bounds
#print axioms survivor_below_five_comb_floor

end LRC14.SixCombNearTiling
