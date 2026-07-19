import Mathlib.Tactic

/-!
# Centered-survivor protrusion arithmetic (THM-1267)

This module is the arithmetic consumer for the paper topology and density
argument.  It checks the exact center shift, the physical-to-normalized
protrusion formula (including `max 0`), both endpoint-density inversions, the
five-load budgets, the separated envelope candidates, both `rho` inequalities,
and the final rational and integer ratio bounds.

The external paper providers are: nearest-integer selection and identification
of the complete safe component; interval-difference topology; THM-1198's
weighted survivor inequality and exact envelope; and THM-1241's strict
`d6 / d1 > 70 / 69` separation.  No analytic or covering claim is hidden in
an axiom here.
-/

namespace LRC14.CenteredSurvivorProtrusion

def threshold : ℚ := 20 / 23
def separatedCap : ℚ := 23 / 120
def physicalProtrusion (c d rho : ℚ) : ℚ :=
  max 0 ((rho + 3 / 7) / d - 3 / (7 * c))
def protrusionExpr (r rho : ℚ) : ℚ :=
  1 / 2 + 7 * rho / 6 - r / 2
def normalizedProtrusion (c d rho : ℚ) : ℚ :=
  max 0 (1 / 2 + 7 * rho / 6 - d / (2 * c))

theorem threshold_simplifications :
    (140 : ℚ) / 161 = threshold ∧ (161 : ℚ) / 840 = separatedCap := by
  norm_num [threshold, separatedCap]

theorem normalized_separation : (6 : ℚ) / 7 * (70 / 69) = threshold := by
  norm_num [threshold]

/-- Once the safe-component integer is identified as `p-k-1`, its center is
shifted from the carrier center by exactly `epsilon/d`. -/
theorem centered_component_center (c d k p t0 epsilon : ℚ)
    (hd : d ≠ 0)
    (ht0 : c * t0 = k + 1 / 2)
    (hp : p = (c + d) * t0 + epsilon) :
    (p - k - 1 / 2) / d = t0 + epsilon / d := by
  have hnum : p - k - 1 / 2 = d * t0 + epsilon := by
    rw [hp]
    calc
      (c + d) * t0 + epsilon - k - 1 / 2 =
          c * t0 - k - 1 / 2 + d * t0 + epsilon := by ring
      _ = d * t0 + epsilon := by rw [ht0]; ring
  rw [hnum]
  field_simp

/-- The centered `d` coordinate lies strictly inside the integer strip whose
floor is one below the half-integer center. -/
theorem half_integer_strip (n x : ℚ) (hx : |x| < 1 / 4) :
    n - 3 / 4 < n - 1 / 2 - x ∧ n - 1 / 2 - x < n - 1 / 4 := by
  rw [abs_lt] at hx
  constructor <;> linarith

theorem raw_protrusion_scaling (c d rho : ℚ) (hc : c ≠ 0) (hd : d ≠ 0) :
    (7 * d / 6) * ((rho + 3 / 7) / d - 3 / (7 * c)) =
      1 / 2 + 7 * rho / 6 - d / (2 * c) := by
  field_simp
  ring

/-- Scaling the physical endpoint tail by the complete safe-component length
gives the normalized protrusion formula, including the contained `max 0`
branch. -/
theorem normalized_protrusion_formula (c d rho : ℚ) (hc : 0 < c) (hd : 0 < d) :
    (7 * d / 6) * physicalProtrusion c d rho =
      normalizedProtrusion c d rho := by
  unfold physicalProtrusion normalizedProtrusion
  let raw : ℚ := (rho + 3 / 7) / d - 3 / (7 * c)
  let scaled : ℚ := 1 / 2 + 7 * rho / 6 - d / (2 * c)
  have hscale : 0 ≤ (7 * d / 6 : ℚ) := by positivity
  have halg : (7 * d / 6) * raw = scaled := by
    dsimp [raw, scaled]
    exact raw_protrusion_scaling c d rho (ne_of_gt hc) (ne_of_gt hd)
  by_cases hraw : 0 ≤ raw
  · have hscaled : 0 ≤ scaled := by
      rw [← halg]
      exact mul_nonneg hscale hraw
    rw [max_eq_right hraw, max_eq_right hscaled, halg]
  · have hraw' : raw ≤ 0 := le_of_not_ge hraw
    have hscaled : scaled ≤ 0 := by
      rw [← halg]
      exact mul_nonpos_of_nonneg_of_nonpos hscale hraw'
    rw [max_eq_left hraw', max_eq_left hscaled, mul_zero]

theorem endpoint_quantile_identities :
    (3 : ℚ) / 4 * (1 / 27) = 1 / 36 ∧
    (3 : ℚ) / 4 * (11 / 270) = 11 / 360 ∧
    (1 : ℚ) / 27 < 1 / 6 ∧ (11 : ℚ) / 270 < 1 / 6 := by
  norm_num

/-- This is the explicit `ell >= 1/6` / `ell < 1/6` split.  In the small
branch the paper density provider supplies `mass <= (3/4) ell`. -/
theorem basic_endpoint_inverse (ell mass : ℚ)
    (hmass : 1 / 36 < mass)
    (hsmall : ell < 1 / 6 → mass ≤ (3 / 4) * ell) :
    1 / 27 < ell := by
  by_cases hell : ell < 1 / 6
  · have hupper := hsmall hell
    linarith
  · have hlarge : 1 / 6 ≤ ell := le_of_not_gt hell
    linarith

theorem separated_endpoint_inverse (ell mass : ℚ)
    (hmass : 11 / 360 < mass)
    (hsmall : ell < 1 / 6 → mass ≤ (3 / 4) * ell) :
    11 / 270 < ell := by
  by_cases hell : ell < 1 / 6
  · have hupper := hsmall hell
    linarith
  · have hlarge : 1 / 6 ≤ ell := le_of_not_gt hell
    linarith

theorem basic_five_load_sum (a b c d e : ℚ)
    (ha : a < 7 / 36) (hb : b < 7 / 36) (hc : c < 7 / 36)
    (hd : d < 7 / 36) (he : e < 7 / 36) :
    a + b + c + d + e < 35 / 36 := by
  linarith

theorem separated_five_load_sum (a b c d e : ℚ)
    (ha : a ≤ 7 / 36) (hb : b ≤ 7 / 36) (hc : c ≤ 7 / 36)
    (hd : d ≤ 7 / 36) (he : e < separatedCap) :
    a + b + c + d + e < 349 / 360 := by
  dsimp [separatedCap] at he
  linarith

theorem separated_survivor_deficit (load : ℚ) (hload : load < 349 / 360) :
    11 / 360 < 1 - load := by
  linarith

/-- Exact maxima which must be compared after the first open envelope piece.
The two closest competitors are `55/288` and the BV-tail value `4/21`. -/
theorem envelope_candidate_gaps :
    separatedCap - 3 / 16 = 1 / 240 ∧
    separatedCap - 55 / 288 = 1 / 1440 ∧
    separatedCap - 13 / 72 = 1 / 90 ∧
    separatedCap - 13 / 84 = 31 / 840 ∧
    separatedCap - 8 / 45 = 1 / 72 ∧
    separatedCap - 1 / 6 = 1 / 40 ∧
    separatedCap - 4 / 21 = 1 / 840 := by
  norm_num [separatedCap]

theorem threshold_inside_first_piece :
    (6 : ℚ) / 7 < threshold ∧ threshold < 68 / 63 := by
  norm_num [threshold]

/-- On the restricted first envelope piece, `1/(6L)` is strictly below the
excluded endpoint value `23/120` whenever `L>20/23`. -/
theorem first_piece_lt_separatedCap (L : ℚ) (hL : threshold < L) :
    1 / (6 * L) < separatedCap := by
  have hLpos : 0 < L := lt_trans (by norm_num [threshold]) hL
  rw [div_lt_iff₀ (by positivity : 0 < 6 * L)]
  dsimp [threshold] at hL
  dsimp [separatedCap]
  nlinarith

theorem five_load_budget_identities :
    4 * ((7 : ℚ) / 36) + separatedCap = 349 / 360 ∧
    1 - (349 : ℚ) / 360 = 11 / 360 := by
  norm_num [separatedCap]

theorem rho_bound_of_basic_protrusion (r rho ell : ℚ)
    (hgeom : ell = protrusionExpr r rho) (hell : 1 / 27 < ell) :
    (27 * r - 25) / 63 < rho := by
  dsimp [protrusionExpr] at hgeom
  rw [hgeom] at hell
  linarith

theorem rho_bound_of_separated_protrusion (r rho ell : ℚ)
    (hgeom : ell = protrusionExpr r rho) (hell : 11 / 270 < ell) :
    (135 * r - 124) / 315 < rho := by
  dsimp [protrusionExpr] at hgeom
  rw [hgeom] at hell
  linarith

theorem basic_ratio_bound (r rho : ℚ)
    (hrho : (27 * r - 25) / 63 < rho) (hrhoMax : rho ≤ 1 / 2) :
    r < 113 / 54 := by
  linarith

theorem separated_ratio_bound (r rho : ℚ)
    (hrho : (135 * r - 124) / 315 < rho) (hrhoMax : rho ≤ 1 / 2) :
    r < 563 / 270 := by
  linarith

theorem ratio_bound_improvements :
    (563 : ℚ) / 270 < 113 / 54 ∧ (113 : ℚ) / 54 < 89 / 42 := by
  norm_num

theorem separated_ratio_cross_multiply (c d : ℚ) (hc : 0 < c)
    (hratio : d / c < 563 / 270) :
    270 * d < 563 * c := by
  have h := (div_lt_iff₀ hc).mp hratio
  nlinarith

theorem integer_strict_rounding (c d : ℤ) (h : 270 * d < 563 * c) :
    270 * d ≤ 563 * c - 1 := by
  omega

#print axioms centered_component_center
#print axioms half_integer_strip
#print axioms normalized_protrusion_formula
#print axioms basic_endpoint_inverse
#print axioms separated_endpoint_inverse
#print axioms basic_five_load_sum
#print axioms separated_five_load_sum
#print axioms envelope_candidate_gaps
#print axioms first_piece_lt_separatedCap
#print axioms rho_bound_of_basic_protrusion
#print axioms rho_bound_of_separated_protrusion
#print axioms basic_ratio_bound
#print axioms separated_ratio_bound
#print axioms separated_ratio_cross_multiply
#print axioms integer_strict_rounding

end LRC14.CenteredSurvivorProtrusion
