import Mathlib

/-!
# GCD-period projective charge (THM-1226)

The arbitrary-speed pair formula, THM-1221 branch classification, and finite
ratio-bank maxima are external providers checked by the exact referee.  This
module kernel-checks the projective charge identity, the vertex-load consumer,
the disconnected-branch and localized constants, the counterfamily threshold
arithmetic, and the protected-needle margins.  It contains no proof
placeholders and uses no native evaluator.
-/

namespace LRC14.GCDPeriodProjectiveCharge

/-- Exact factorization of the period error through the symmetric projective
height `ab/(a+b)`, when the two speeds are `g*a` and `g*b`. -/
theorem projective_charge_factorization
    (g a b variance : ℚ)
    (hg : g ≠ 0) (ha : a ≠ 0) (hb : b ≠ 0) (hab : a + b ≠ 0) :
    variance / g =
      (variance * (a * b) / (a + b)) *
        (1 / (g * a) + 1 / (g * b)) := by
  field_simp
  ring

/-- The same error may be charged wholly to either oriented endpoint. -/
theorem oriented_charge_factorization
    (g a variance : ℚ) (hg : g ≠ 0) (ha : a ≠ 0) :
    variance / g = variance * a / (g * a) := by
  field_simp

/-- Seven-vertex weighted-load consumer.  The reciprocal speed weights are
nonnegative, so a uniform vertex-load cap bounds the whole tree error. -/
theorem seven_vertex_load_consumer
    (r0 r1 r2 r3 r4 r5 r6 l0 l1 l2 l3 l4 l5 l6 C E H : ℝ)
    (hr0 : 0 ≤ r0) (hr1 : 0 ≤ r1) (hr2 : 0 ≤ r2)
    (hr3 : 0 ≤ r3) (hr4 : 0 ≤ r4) (hr5 : 0 ≤ r5) (hr6 : 0 ≤ r6)
    (hl0 : l0 ≤ C) (hl1 : l1 ≤ C) (hl2 : l2 ≤ C)
    (hl3 : l3 ≤ C) (hl4 : l4 ≤ C) (hl5 : l5 ≤ C) (hl6 : l6 ≤ C)
    (hE : E = r0*l0 + r1*l1 + r2*l2 + r3*l3 + r4*l4 + r5*l5 + r6*l6)
    (hH : H = r0 + r1 + r2 + r3 + r4 + r5 + r6) :
    E ≤ C * H := by
  have h0 := mul_le_mul_of_nonneg_left hl0 hr0
  have h1 := mul_le_mul_of_nonneg_left hl1 hr1
  have h2 := mul_le_mul_of_nonneg_left hl2 hr2
  have h3 := mul_le_mul_of_nonneg_left hl3 hr3
  have h4 := mul_le_mul_of_nonneg_left hl4 hr4
  have h5 := mul_le_mul_of_nonneg_left hl5 hr5
  have h6 := mul_le_mul_of_nonneg_left hl6 hr6
  rw [hE, hH]
  nlinarith

/-! ## Exact disconnected-spectrum constants -/

theorem disconnected_branch_constants :
    6 * ((224458 : ℚ) / 584325) = 448916 / 194775 ∧
    6 * ((85975 : ℚ) / 342804) = 85975 / 57134 ∧
    6 * ((43774 : ℚ) / 276507) = 87548 / 92169 ∧
    (85975 : ℚ) / 57134 < 448916 / 194775 ∧
    (87548 : ℚ) / 92169 < 448916 / 194775 := by
  norm_num

theorem disconnected_crown_ratio :
    ((15 : ℚ) / 154) / (448916 / 194775 + 1 / 7) =
      417375 / 10488302 := by
  norm_num

/-- Cleared-denominator consumer for the conditional localized tree bound.
The providers supply the bulk tree term and its period-error bound. -/
theorem disconnected_localized_consumer
    (L H bulk error localMass : ℝ)
    (hbulk : 15 * L ≤ 154 * bulk)
    (herror : 194775 * error ≤ 448916 * H)
    (hlocal : bulk - error ≤ localMass) :
    15 * 194775 * L - 154 * 448916 * H ≤
      154 * 194775 * localMass := by
  nlinarith

/-! ## Counterfamily and protected-needle arithmetic -/

theorem strict_high_tail_identity :
    (1 : ℚ) / 49 - 5 / 308 = 9 / 2156 ∧
    6 * (5 / 308) = 15 / 154 ∧
    6 * ((5 / 308) * (1 - 5 / 308)) = 4545 / 47432 := by
  norm_num

/-- The cleared inequality behind `x_A>5/308`. -/
theorem strict_high_of_quadratic_cut (A : ℚ) (hA : 0 < A)
    (hcut : 539 < 9 * A^2) :
    5 / 308 < 1 / 49 - 1 / (4 * A^2) := by
  have hA2 : 0 < A^2 := sq_pos_of_pos hA
  have hden : 0 < 4 * A^2 := by positivity
  have hfrac : 1 / (4 * A^2) < (9 : ℚ) / 2156 := by
    rw [div_lt_iff₀ hden]
    nlinarith
  nlinarith

theorem base_step_arithmetic :
    27720 % 14 = 0 ∧
    (1 + 2 + 3 + 5 + 7 + 11 + 13 : ℕ) = 42 ∧
    9 * (27720 : ℕ)^2 > 539 := by
  norm_num

theorem protected_integer_data :
    let q : ℕ := 27721
    let core : List ℕ := [41582, 41583, 41584, 41585, 41586, 41587]
    q % 2 = 1 ∧ core.getLast? = some 41587 ∧
      2 * 41582 = 3 * q + 1 ∧ 2 * 41587 = 3 * q + 11 := by
  norm_num

theorem core_margin_positive (q : ℚ) (hq : 17 ≤ q) :
    0 < (5*q - 77) / (14*q) := by
  have hq0 : 0 < q := by linarith
  have hnum : 0 < 5*q - 77 := by linarith
  have hden : 0 < 14*q := by positivity
  exact div_pos hnum hden

theorem danger_polynomial_positive (q : ℚ) (hq : 521 ≤ q) :
    0 < q^2 - 517*q - 1848 := by
  have hsquare : 0 ≤ (q - 521)^2 := sq_nonneg (q - 521)
  nlinarith

theorem danger_margin_positive (q : ℚ) (hq : 521 ≤ q) :
    0 < (q^2 - 517*q - 1848) / (14*q*(3*q+11)) := by
  have hq0 : 0 < q := by linarith
  have hpoly := danger_polynomial_positive q hq
  positivity

theorem danger_threshold_exact :
    (520 : ℤ)^2 - 517*520 - 1848 = -288 ∧
    (521 : ℤ)^2 - 517*521 - 1848 = 236 := by
  norm_num

#print axioms projective_charge_factorization
#print axioms oriented_charge_factorization
#print axioms seven_vertex_load_consumer
#print axioms disconnected_branch_constants
#print axioms disconnected_localized_consumer
#print axioms strict_high_of_quadratic_cut
#print axioms protected_integer_data
#print axioms core_margin_positive
#print axioms danger_margin_positive

end LRC14.GCDPeriodProjectiveCharge
