import Mathlib
import TournamentH7.CombPatterns

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

/-- THM-1166's forest-Hunter cap `6H/49` gives the sharper crown used after
the projective-charge estimate.  The second and third identities are the
protected-needle and smallest-speed consequences. -/
theorem disconnected_forest_crown_constants :
    ((15 : ℚ) / 154) / (448916 / 194775 + 6 / 49) =
      59625 / 1485836 ∧
    ((59625 : ℚ) / 1485836) / 7 = 59625 / 10400852 ∧
    7 / ((59625 : ℚ) / 10400852) = 72805964 / 59625 := by
  norm_num

/-- Arithmetic cap obtained by composing the protected-needle crown with
THM-1233's `d₆/a<2345` separated slow-gap compactification. -/
theorem disconnected_separated_packet_cap :
    2345 * ((72805964 : ℚ) / 59625) = 34145997116 / 11925 := by
  norm_num

/-- At the canonical protected-needle length, `a<13m` makes the F7 right
side nonpositive as soon as `C≥13(15/154)/7=195/1078`. -/
theorem shallow_f7_constant :
    13 * ((15 : ℚ) / 154) / 7 = 195 / 1078 ∧
    (195 : ℚ) / 1078 < 448916 / 194775 := by
  norm_num

/-- Cleared-denominator consumer for the shallow-branch vacuity.  The
hypothesis `7L≤13H` is supplied by `|I|=1/(7m)` and `min(S)<13m`. -/
theorem shallow_f7_vacuity
    (L H C localMass : ℝ)
    (hH0 : 0 ≤ H) (hlocal : 0 ≤ localMass)
    (hneedle : 7 * L ≤ 13 * H)
    (hC : 195 ≤ 1078 * C) :
    15 * L ≤ 154 * (localMass + C * H) := by
  have hscaled : 195 * H ≤ (1078 * C) * H :=
    mul_le_mul_of_nonneg_right hC hH0
  nlinarith

/-! ## Farey-height seven phase fibers -/

/-- The exact affine relation behind the strict-high but locally empty pair
`(N+1,2N+1)`. -/
theorem two_one_relation (N : ℚ) :
    2 * (N + 1) - (2 * N + 1) = 1 := by
  ring

/-- Primitive relation coordinates are unimodular.  If `qu-pv=1`, then
`(a,b)=(pk+uh,qk+vh)` has transverse coordinates
`h=qa-pb` and `k=ub-va`.  This is the presentation-invariant replacement for
the nonprimitive clock refuted in MISTAKE-184. -/
theorem primitive_relation_coordinates
    (p q u v k h a b : ℤ)
    (hbezout : q*u - p*v = 1)
    (ha : a = p*k + u*h)
    (hb : b = q*k + v*h) :
    q*a - p*b = h ∧ u*b - v*a = k := by
  constructor
  · rw [ha, hb]
    calc
      q * (p*k + u*h) - p * (q*k + v*h) = (q*u - p*v)*h := by ring
      _ = h := by rw [hbezout]; ring
  · rw [ha, hb]
    calc
      u * (q*k + v*h) - v * (p*k + u*h) = (q*u - p*v)*k := by ring
      _ = k := by rw [hbezout]; ring

/-- A unimodular coordinate change preserves the complete common-divisor
sheet, not merely the two coordinate equations. -/
theorem primitive_relation_common_divisors
    (p q u v k h a b d : ℤ)
    (hbezout : q*u - p*v = 1)
    (ha : a = p*k + u*h)
    (hb : b = q*k + v*h) :
    (d ∣ k ∧ d ∣ h) ↔ (d ∣ a ∧ d ∣ b) := by
  have hinv := primitive_relation_coordinates p q u v k h a b hbezout ha hb
  constructor
  · rintro ⟨⟨x, hx⟩, ⟨y, hy⟩⟩
    constructor
    · refine ⟨p*x + u*y, ?_⟩
      calc
        a = p*k + u*h := ha
        _ = d * (p*x + u*y) := by rw [hx, hy]; ring
    · refine ⟨q*x + v*y, ?_⟩
      calc
        b = q*k + v*h := hb
        _ = d * (q*x + v*y) := by rw [hx, hy]; ring
  · rintro ⟨⟨x, hx⟩, ⟨y, hy⟩⟩
    constructor
    · refine ⟨u*y - v*x, ?_⟩
      calc
        k = u*b - v*a := hinv.2.symm
        _ = d * (u*y - v*x) := by rw [hx, hy]; ring
    · refine ⟨q*x - p*y, ?_⟩
      calc
        h = q*a - p*b := hinv.1.symm
        _ = d * (q*x - p*y) := by rw [hx, hy]; ring

/-- In particular the positive gcd is invariant under the primitive
unimodular relation coordinates. -/
theorem primitive_relation_gcd
    (p q u v k h a b : ℤ)
    (hbezout : q*u - p*v = 1)
    (ha : a = p*k + u*h)
    (hb : b = q*k + v*h) :
    Int.gcd k h = Int.gcd a b := by
  have hdiv := primitive_relation_common_divisors p q u v k h a b
    (Int.gcd k h : ℤ) hbezout ha hb
  have hdiv' := primitive_relation_common_divisors p q u v k h a b
    (Int.gcd a b : ℤ) hbezout ha hb
  apply Nat.dvd_antisymm
  · rw [← Int.natCast_dvd_natCast]
    exact Int.dvd_coe_gcd
      (hdiv.mp ⟨Int.gcd_dvd_left k h, Int.gcd_dvd_right k h⟩).1
      (hdiv.mp ⟨Int.gcd_dvd_left k h, Int.gcd_dvd_right k h⟩).2
  · rw [← Int.natCast_dvd_natCast]
    exact Int.dvd_coe_gcd
      (hdiv'.mpr ⟨Int.gcd_dvd_left a b, Int.gcd_dvd_right a b⟩).1
      (hdiv'.mpr ⟨Int.gcd_dvd_left a b, Int.gcd_dvd_right a b⟩).2

/-- Exact enumeration of THM-605's primitive height-seven phase channels. -/
theorem height_seven_channel_list
    (p q : ℕ) (hp : 0 < p) (hpq : p ≤ q) (hheight : p + q ≤ 7)
    (hcop : Nat.Coprime p q) :
    (p = 1 ∧ q = 1) ∨ (p = 1 ∧ q = 2) ∨ (p = 1 ∧ q = 3) ∨
    (p = 1 ∧ q = 4) ∨ (p = 2 ∧ q = 3) ∨ (p = 1 ∧ q = 5) ∨
    (p = 1 ∧ q = 6) ∨ (p = 2 ∧ q = 5) ∨ (p = 3 ∧ q = 4) := by
  have hp7 : p ≤ 7 := by omega
  have hq7 : q ≤ 7 := by omega
  interval_cases p <;> interval_cases q <;> norm_num [Nat.Coprime] at *

/-- Positivity of the elementary tangent-fiber lower bound above Farey
height seven.  The trapezoid/periodization identity is the paper provider. -/
theorem fiber_lower_bound_positive
    (p q : ℝ) (hp : 0 < p) (hq : 0 < q) (hheight : 7 < p + q) :
    0 < min (1 / (7 * q)) ((p + q - 7) / (14 * p * q)) := by
  rw [lt_min_iff]
  constructor <;> positivity

/-! The non-diagonal height-seven channels have pair masses
`1/14,1/21,1/28,1/35,1/42` according to their larger coefficient.  The
following arithmetic table checks their oriented charges; the repeated
channels have the same displayed row. -/
theorem height_seven_oriented_charge_table :
    (1 / 7 : ℚ) * (1 - 1 / 7) = 6 / 49 ∧
    2 * (1 / 14 : ℚ) * (1 - 1 / 14) = 13 / 98 ∧
    3 * (1 / 21 : ℚ) * (1 - 1 / 21) = 20 / 147 ∧
    4 * (1 / 28 : ℚ) * (1 - 1 / 28) = 27 / 196 ∧
    5 * (1 / 35 : ℚ) * (1 - 1 / 35) = 34 / 245 ∧
    6 * (1 / 42 : ℚ) * (1 - 1 / 42) = 41 / 294 := by
  norm_num

theorem exact_relation_tree_constants :
    (6 : ℚ) / 49 ≤ 41 / 294 ∧
    (13 : ℚ) / 98 ≤ 41 / 294 ∧
    (20 : ℚ) / 147 ≤ 41 / 294 ∧
    (27 : ℚ) / 196 ≤ 41 / 294 ∧
    (34 : ℚ) / 245 ≤ 41 / 294 ∧
    (41 : ℚ) / 294 < 39 / 98 ∧
    (6 : ℕ)^6 = 46656 ∧ 6 * (6 : ℕ)^5 = 46656 := by
  norm_num

/-- Three-edge instance of the primitive relation-cycle holonomy identity.
The general cycle is obtained by the same telescoping substitution. -/
theorem relation_triangle_holonomy
    (s1 s2 s3 p1 p2 p3 q1 q2 q3 h1 h2 h3 : ℤ)
    (hr1 : q1*s2 - p1*s1 = h1)
    (hr2 : q2*s3 - p2*s2 = h2)
    (hr3 : q3*s1 - p3*s3 = h3) :
    (q1*q2*q3 - p1*p2*p3)*s1 =
      p2*p3*h1 + q1*p3*h2 + q1*q2*h3 := by
  linear_combination p2*p3*hr1 + q1*p3*hr2 + q1*q2*hr3

theorem relation_cycle_crown_constants :
    7 * (6 : ℕ)^6 = 326592 ∧
    7 * (13 : ℕ)^6 = 33787663 := by
  norm_num

/-- In the balanced-holonomy branch, positive weights cannot annihilate
nonzero defects of only one sign. -/
theorem balanced_triangle_nonnegative_zero
    (c1 c2 c3 h1 h2 h3 : ℝ)
    (hc1 : 0 < c1) (hc2 : 0 < c2) (hc3 : 0 < c3)
    (hh1 : 0 ≤ h1) (hh2 : 0 ≤ h2) (hh3 : 0 ≤ h3)
    (hbalance : c1*h1 + c2*h2 + c3*h3 = 0) :
    h1 = 0 ∧ h2 = 0 ∧ h3 = 0 := by
  have h1nonneg : 0 ≤ c1*h1 := mul_nonneg (le_of_lt hc1) hh1
  have h2nonneg : 0 ≤ c2*h2 := mul_nonneg (le_of_lt hc2) hh2
  have h3nonneg : 0 ≤ c3*h3 := mul_nonneg (le_of_lt hc3) hh3
  constructor
  · nlinarith
  constructor <;> nlinarith

/-- Lifting a blocker residue from speed `dⱼ` to the centered spoke
denominator `Qⱼ=c+dⱼ` adds exactly the positive carrier term `P*c`. -/
theorem centered_blocker_relation_lift
    (P N c dj Qi Qj r : ℤ)
    (hQj : Qj = c + dj)
    (hresidue : P*dj - N*Qi = r) :
    P*Qj - N*Qi = P*c + r := by
  rw [hQj]
  calc
    P * (c + dj) - N*Qi = P*c + (P*dj - N*Qi) := by ring
    _ = P*c + r := by rw [hresidue]

/-- Abstract sign consumer used after telescoping a centered blocker cycle:
positive defect and positive base force nontrivial positive holonomy. -/
theorem positive_cycle_holonomy_sign
    (qProduct pProduct base defect : ℝ)
    (hbase : 0 < base) (hdefect : 0 < defect)
    (hidentity : (qProduct - pProduct)*base = defect) :
    pProduct < qProduct := by
  nlinarith

/-! ## Composition with THM-1237's positioned protected-needle debt -/

theorem positioned_debt_dichotomy_constants :
    24 * ((5 : ℚ) / 88) = 15 / 11 ∧
    ((30 : ℚ) / 11 - 15 / 11) / 13 = 15 / 143 ∧
    ((15 : ℚ) / 143) / 6 = 5 / 286 ∧
    7 / ((616 : ℚ) / 5) = 5 / 88 := by
  norm_num

/-- If the harmonic half of THM-1237's protected-needle debt is at most
`15/11`, the six reciprocal gcd edges must pay the remaining `15/11`. -/
theorem positioned_debt_forces_gcd_sum
    (m H G : ℝ)
    (hH : 88 * m * H ≤ 5)
    (hdebt : 30 / 11 ≤ 24 * m * H + 13 * m * G) :
    15 ≤ 143 * m * G := by
  nlinarith

/-- A six-edge reciprocal sum paying `15/143` has an edge paying at least
`5/286`. -/
theorem six_edge_reciprocal_pigeonhole
    (m r0 r1 r2 r3 r4 r5 : ℝ)
    (hsum : 15 ≤ 143 * m * (r0 + r1 + r2 + r3 + r4 + r5)) :
    5 ≤ 286*m*r0 ∨ 5 ≤ 286*m*r1 ∨ 5 ≤ 286*m*r2 ∨
    5 ≤ 286*m*r3 ∨ 5 ≤ 286*m*r4 ∨ 5 ≤ 286*m*r5 := by
  by_contra h
  push Not at h
  nlinarith

/-- Convert the reciprocal-edge conclusion back to the gcd clock. -/
theorem reciprocal_edge_gcd_bound
    (m r g : ℝ) (hg : 0 ≤ g) (hrecip : r*g = 1)
    (hedge : 5 ≤ 286*m*r) :
    5*g ≤ 286*m := by
  have hmul := mul_le_mul_of_nonneg_right hedge hg
  nlinarith

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

/-- Cleared-denominator composition of the global tree floor, projective
error bound, and THM-1166 forest cover cap. -/
theorem disconnected_forest_crown_consumer
    (L H bulk error : ℝ)
    (hbulk : 15 * L ≤ 154 * bulk)
    (herror : 194775 * error ≤ 448916 * H)
    (hforest : 49 * (bulk - error) ≤ 6 * H) :
    59625 * L ≤ 1485836 * H := by
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
#print axioms disconnected_forest_crown_consumer
#print axioms disconnected_separated_packet_cap
#print axioms shallow_f7_vacuity
#print axioms primitive_relation_coordinates
#print axioms primitive_relation_common_divisors
#print axioms primitive_relation_gcd
#print axioms height_seven_channel_list
#print axioms fiber_lower_bound_positive
#print axioms height_seven_oriented_charge_table
#print axioms exact_relation_tree_constants
#print axioms relation_triangle_holonomy
#print axioms balanced_triangle_nonnegative_zero
#print axioms centered_blocker_relation_lift
#print axioms positive_cycle_holonomy_sign
#print axioms positioned_debt_forces_gcd_sum
#print axioms six_edge_reciprocal_pigeonhole
#print axioms reciprocal_edge_gcd_bound
#print axioms strict_high_of_quadratic_cut
#print axioms protected_integer_data
#print axioms core_margin_positive
#print axioms danger_margin_positive

end LRC14.GCDPeriodProjectiveCharge
