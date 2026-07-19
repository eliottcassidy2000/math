import Mathlib

/-!
# Arithmetic consumers for the universal five-comb dual density (THM-1198)

THM-1198 supplies the analytic one-comb estimate

`A(L, phi) <= 7/36`

for every `L >= 6/7`, together with the interpretation of `A` as the mass of
a danger comb under a six-bin probability density.  This module checks the
finite rational ledger downstream of that estimate:

* the density mass and total-variation constants;
* five-load noncoverage and the `1/36` dual survivor;
* the `1/42` normalized and `1/(49*c)` physical length conversions;
* the six-load overlap, unique-provider, and physical-length conversions;
* the integer first-tooth form and the phase-envelope guardrail; and
* an abstract phase-free functional-drift consumer.

The compact arrangement exhaustiveness, its identification with the integral
`A`, the equality locus, and the Stieltjes/BV integration-by-parts estimate
are deliberately external provider theorems.  No analytic statement is
silently reconstructed from the arithmetic below.
-/

namespace LRC14.FiveCombDual

/-- The universal one-comb dual-load cap supplied by THM-1198's analytic
arrangement/BV argument. -/
def oneCombCap : ℚ := 7 / 36

/-- The six rational bin heights have average one. -/
theorem six_bin_density_mass :
    ((3 : ℚ) / 4 + 13 / 12 + 7 / 6 + 7 / 6 + 13 / 12 + 3 / 4) / 6 = 1 := by
  norm_num

/-- The interior variation of the six-bin step density is `5/6`. -/
theorem six_bin_interior_variation :
    |(13 : ℚ) / 12 - 3 / 4| + |(7 : ℚ) / 6 - 13 / 12| +
        |(7 : ℚ) / 6 - 7 / 6| + |(13 : ℚ) / 12 - 7 / 6| +
        |(3 : ℚ) / 4 - 13 / 12| = 5 / 6 := by
  norm_num [abs_of_nonneg, abs_of_nonpos]

/-- Adding the two endpoint jumps gives zero-extension variation `7/3`. -/
theorem six_bin_zero_extension_variation :
    (3 : ℚ) / 4 + 5 / 6 + 3 / 4 = 7 / 3 := by
  norm_num

/-- The primitive norm `3/49` and zero-extension variation `7/3` give the
coefficient `1/7` in the BV tail. -/
theorem bv_tail_coefficient :
    ((3 : ℚ) / 49) * (7 / 3) = 1 / 7 := by
  norm_num

/-- At the analytic cutoff `L=3`, the BV upper bound is strictly below the
compact optimum, with exact margin `1/252`. -/
theorem bv_cutoff_margin :
    (1 : ℚ) / 7 + 1 / (7 * 3) = 4 / 21 ∧
      (7 : ℚ) / 36 - 4 / 21 = 1 / 252 := by
  norm_num

/-- Five quantities individually bounded by `7/36` have total at most
`35/36`.  This is the finite consumer used after integrating five danger
indicators against the dual density. -/
theorem five_load_sum_le (load : Fin 5 → ℚ)
    (hload : ∀ i, load i ≤ oneCombCap) :
    ∑ i, load i ≤ 35 / 36 := by
  calc
    ∑ i, load i ≤ ∑ _i : Fin 5, oneCombCap :=
      Finset.sum_le_sum fun i _hi ↦ hload i
    _ = 35 / 36 := by norm_num [oneCombCap, Fin.sum_univ_succ]

/-- The abstract five-load contradiction: a cover would require dual mass at
least one, while the five analytic load caps total only `35/36`. -/
theorem five_load_contradiction (load : Fin 5 → ℚ)
    (hload : ∀ i, load i ≤ oneCombCap)
    (hcover : 1 ≤ ∑ i, load i) : False := by
  have hsum := five_load_sum_le load hload
  linarith

/-- Quantitative union-bound consumer.  If `uncoveredDualMass` is at least
one minus the sum of the five comb loads, then it is at least `1/36`. -/
theorem five_load_survivor_mass (load : Fin 5 → ℚ)
    (uncoveredDualMass : ℚ)
    (hload : ∀ i, load i ≤ oneCombCap)
    (hunion : 1 - ∑ i, load i ≤ uncoveredDualMass) :
    1 / 36 ≤ uncoveredDualMass := by
  have hsum := five_load_sum_le load hload
  linarith

/-- If the dual density is bounded above by `7/6`, dual survivor mass
`1/36` forces normalized Lebesgue length at least `1/42`. -/
theorem normalized_survivor_length
    (dualMass normalizedLength : ℚ)
    (hmass : 1 / 36 ≤ dualMass)
    (hdensity : dualMass ≤ (7 / 6) * normalizedLength) :
    1 / 42 ≤ normalizedLength := by
  linarith

/-- Rescaling a normalized `c`-slow gap by `6/(7*c)` converts the normalized
`1/42` floor into the physical `1/(49*c)` floor. -/
theorem physical_survivor_length
    (c normalizedLength : ℚ) (hc : 0 < c)
    (hlength : 1 / 42 ≤ normalizedLength) :
    1 / (49 * c) ≤ (6 / (7 * c)) * normalizedLength := by
  have hscale : 0 ≤ (6 : ℚ) / (7 * c) := by positivity
  calc
    (1 : ℚ) / (49 * c) = (6 / (7 * c)) * (1 / 42) := by
      field_simp
      ring
    _ ≤ (6 / (7 * c)) * normalizedLength :=
      mul_le_mul_of_nonneg_left hlength hscale

/-- Six loads individually bounded by `7/36` have excess over the unit cover
mass at most `1/6`. -/
theorem six_load_overlap_surplus (load : Fin 6 → ℚ)
    (hload : ∀ i, load i ≤ oneCombCap) :
    (∑ i, load i) - 1 ≤ 1 / 6 := by
  have hsum : ∑ i, load i ≤ ∑ _i : Fin 6, oneCombCap :=
    Finset.sum_le_sum fun i _hi ↦ hload i
  have hcap : (∑ _i : Fin 6, oneCombCap) = (7 : ℚ) / 6 := by
    norm_num [oneCombCap, Fin.sum_univ_succ]
  rw [hcap] at hsum
  linarith

/-- A multiply-covered region with dual mass at most `1/6`, under a density
bounded below by `3/4`, has normalized length at most `2/9`. -/
theorem multiply_covered_length_ceiling
    (multiMass multiLength : ℚ)
    (hdensity : (3 / 4) * multiLength ≤ multiMass)
    (hmass : multiMass ≤ 1 / 6) :
    multiLength ≤ 2 / 9 := by
  linarith

/-- The complementary unique-provider region therefore has normalized length
at least `7/9`. -/
theorem unique_provider_length_floor
    (multiLength uniqueLength : ℚ)
    (hpartition : uniqueLength = 1 - multiLength)
    (hmulti : multiLength ≤ 2 / 9) :
    7 / 9 ≤ uniqueLength := by
  linarith

/-- The normalized unique-provider floor rescales to physical length
`2/(3*c)` in a complete `c`-slow gap. -/
theorem physical_unique_provider_length
    (c uniqueLength : ℚ) (hc : 0 < c)
    (hlength : 7 / 9 ≤ uniqueLength) :
    2 / (3 * c) ≤ (6 / (7 * c)) * uniqueLength := by
  have hscale : 0 ≤ (6 : ℚ) / (7 * c) := by positivity
  calc
    (2 : ℚ) / (3 * c) = (6 / (7 * c)) * (7 / 9) := by
      field_simp
      ring
    _ ≤ (6 / (7 * c)) * uniqueLength :=
      mul_le_mul_of_nonneg_left hlength hscale

/-- For positive integer speeds, the strict first-tooth inequality is exactly
the subtraction form recorded in THM-1198. -/
theorem integer_first_tooth_form (c d : ℕ) (hc : 1 ≤ c) :
    6 * d < 13 * c ↔ 6 * d ≤ 13 * c - 1 := by
  omega

/-- Cast-aware version of the integer first-tooth form
`d < 13*c/6  <->  6*d <= 13*c-1`. -/
theorem rational_first_tooth_iff_integer (c d : ℕ) (hc : 1 ≤ c) :
    (d : ℚ) < 13 * (c : ℚ) / 6 ↔ 6 * d ≤ 13 * c - 1 := by
  rw [← integer_first_tooth_form c d hc]
  constructor
  · intro h
    have hq : (6 : ℚ) * d < 13 * c := by nlinarith
    exact_mod_cast hq
  · intro h
    have hq : (6 : ℚ) * d < 13 * c := by exact_mod_cast h
    nlinarith

/-- The compact phase envelope is not monotone: its value at `15/7` is
already strictly larger than `1/6`. -/
theorem phase_envelope_nonmonotone_guardrail :
    (1 : ℚ) / 6 < 8 / 45 := by
  norm_num

/-- Abstract functional-drift consumer.  A pointwise majorant for six comb
loads inherits the necessary total-mass inequality from any cover. -/
theorem phase_free_functional_drift
    (load majorant : Fin 6 → ℚ)
    (hmajorize : ∀ i, load i ≤ majorant i)
    (hcover : 1 ≤ ∑ i, load i) :
    1 ≤ ∑ i, majorant i := by
  exact hcover.trans (Finset.sum_le_sum fun i _hi ↦ hmajorize i)

#print axioms six_bin_density_mass
#print axioms six_bin_interior_variation
#print axioms six_bin_zero_extension_variation
#print axioms bv_tail_coefficient
#print axioms bv_cutoff_margin
#print axioms five_load_sum_le
#print axioms five_load_contradiction
#print axioms five_load_survivor_mass
#print axioms normalized_survivor_length
#print axioms physical_survivor_length
#print axioms six_load_overlap_surplus
#print axioms multiply_covered_length_ceiling
#print axioms unique_provider_length_floor
#print axioms physical_unique_provider_length
#print axioms integer_first_tooth_form
#print axioms rational_first_tooth_iff_integer
#print axioms phase_envelope_nonmonotone_guardrail
#print axioms phase_free_functional_drift

end LRC14.FiveCombDual
