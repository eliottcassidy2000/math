/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-18)
-/
import Mathlib.Tactic

/-!
# Cleared arithmetic for the seven-wall Fano/gcd certificate

This module formalizes the denominator-cleared arithmetic consumers in
THM-1166.  It deliberately does not formalize the analytic inputs (the
three-speed overlap floor or the periodic-remainder inequality): those enter
below only through explicit hypotheses.

Every result below is kernel-checked without proof placeholders.
-/

namespace LRC14
namespace SevenWallFanoGCD

/-! ## The quadratic multiplicity certificate -/

/-- The integer numerator of `Q(C) = C(8-C)/7`. -/
def qNumerator (C : ℤ) : ℤ := C * (8 - C)

/-- The exact rational quadratic certificate. -/
def Q (C : ℤ) : ℚ := qNumerator C / 7

/-- The complete depth table for seven combs. -/
theorem q_table :
    Q 0 = 0 ∧ Q 1 = 1 ∧ Q 2 = 12 / 7 ∧ Q 3 = 15 / 7 ∧
      Q 4 = 16 / 7 ∧ Q 5 = 15 / 7 ∧ Q 6 = 12 / 7 ∧ Q 7 = 1 := by
  norm_num [Q, qNumerator]

/-- Cleared form of `Q(C) ≥ 1` on every occupied depth `1,…,7`. -/
theorem qNumerator_ge_seven {C : ℤ} (h1 : 1 ≤ C) (h7 : C ≤ 7) :
    7 ≤ qNumerator C := by
  have hslack : 0 ≤ (C - 1) * (7 - C) :=
    mul_nonneg (sub_nonneg.mpr h1) (sub_nonneg.mpr h7)
  simp only [qNumerator]
  nlinarith

/-- Rational form of the occupied-depth lower bound. -/
theorem one_le_Q {C : ℤ} (h1 : 1 ≤ C) (h7 : C ≤ 7) : 1 ≤ Q C := by
  have hclear := qNumerator_ge_seven h1 h7
  unfold Q
  rw [le_div_iff₀ (by norm_num : (0 : ℚ) < 7)]
  norm_num
  exact_mod_cast hclear

/-- The numerator is bounded above by `16`; this is the `j=4` peak. -/
theorem qNumerator_le_sixteen (C : ℤ) : qNumerator C ≤ 16 := by
  have hsquare : 0 ≤ (C - 4) ^ 2 := sq_nonneg (C - 4)
  simp only [qNumerator]
  nlinarith

/-- The peak is unique: equality occurs exactly at depth four. -/
theorem qNumerator_eq_sixteen_iff (C : ℤ) : qNumerator C = 16 ↔ C = 4 := by
  constructor
  · intro h
    have hsquare : 0 ≤ (C - 4) ^ 2 := sq_nonneg (C - 4)
    simp only [qNumerator] at h
    nlinarith
  · rintro rfl
    norm_num [qNumerator]

/-! ## From three-speed packets to uncovered mass -/

/-- The numerical constant obtained after summing the 35 triple floors and
dividing by the five appearances of every edge. -/
theorem seven_times_three_speed_floor :
    (7 : ℚ) * (1 / 24) = 7 / 24 := by
  norm_num

/-- Denominator-cleared triple-counting step.  The hypothesis is exactly
`5R ≥ 35/24`, multiplied by `24`. -/
theorem pair_mass_of_thirty_five_triples (R : ℝ) (htriples : 35 ≤ 120 * R) :
    7 ≤ 24 * R := by
  linarith

/-- If `R ≥ 7/24` and uncovered mass `U` satisfies `U ≥ 2R/7`, then
`U ≥ 1/12`, with every denominator cleared. -/
theorem uncovered_floor_of_pair_mass (R U : ℝ)
    (hpair : 7 ≤ 24 * R) (huncovered : 2 * R ≤ 7 * U) :
    1 ≤ 12 * U := by
  linarith

/-- The corresponding exact rational identity. -/
theorem exact_uncovered_floor :
    ((2 : ℚ) / 7) * (7 / 24) = 1 / 12 := by
  norm_num

/-! ## Common-dilate composition -/

/-- The protected-needle lower bound and the covered-period upper bound imply
the exact common-dilate scale restriction.  All divisions have been cleared:

* `1 ≤ 7mL` is `L ≥ 1/(7m)`;
* `12GL ≤ 11` is `L ≤ 11/(12G)`;
* the conclusion is `12G ≤ 77m`.
-/
theorem common_dilate_scale_bound (m G L : ℝ) (hm : 0 ≤ m) (hG : 0 ≤ G)
    (hneedle : 1 ≤ 7 * m * L) (hcovered : 12 * G * L ≤ 11) :
    12 * G ≤ 77 * m := by
  calc
    12 * G = (12 * G) * 1 := by ring
    _ ≤ (12 * G) * (7 * m * L) :=
      mul_le_mul_of_nonneg_left hneedle (mul_nonneg (by norm_num) hG)
    _ = (7 * m) * (12 * G * L) := by ring
    _ ≤ (7 * m) * 11 :=
      mul_le_mul_of_nonneg_left hcovered (mul_nonneg (by norm_num) hm)
    _ = 77 * m := by ring

/-! ## Exact Fano constants and budget -/

/-- Mean of one Fano-line summand as a function of its pair mass. -/
def lineMean (Rline : ℚ) : ℚ := (1 - 2 * Rline) / 7

/-- Sharp periodic-remainder error coefficient at mean `ν`. -/
def remainderError (nu : ℚ) : ℚ := nu * (1 - 21 * nu / 8)

/-- At the three-speed floor `R_line = 1/24`, the line mean is `11/84`. -/
theorem exact_line_mean : lineMean (1 / 24) = 11 / 84 := by
  norm_num [lineMean]

/-- The corresponding maximal remainder coefficient is `11/128`. -/
theorem exact_remainder_error : remainderError (11 / 84) = 11 / 128 := by
  norm_num [remainderError]

/-- The metric-Fano budget constant. -/
theorem exact_fano_budget_constant :
    (128 : ℚ) / (84 * 11) = 32 / 231 := by
  norm_num

/-- Averaging a `32/231` budget over seven lines gives the reciprocal line
bound `1617/32`. -/
theorem exact_fano_line_constant :
    (7 : ℚ) / (32 / 231) = 1617 / 32 := by
  norm_num

/-- Clearing the two sides of the metric-Fano inequality by `896m` produces
the coefficients `256` and `77`. -/
theorem exact_metric_clearing :
    (896 : ℚ) * (2 / 7) = 256 ∧ (896 : ℚ) * (11 / 128) = 77 := by
  norm_num

/-- Pure arithmetic consumer for the Fano budget.  Here `S` denotes
`Σ_line m/G_line`.  The last hypothesis is the metric-Fano inequality after
the line errors have been bounded by `11/128` and all denominators cleared. -/
theorem fano_budget_from_cleared_inequalities (m L R S : ℝ)
    (hneedle : 1 ≤ 7 * m * L) (hpair : 7 ≤ 24 * R)
    (hmetric : 256 * m * L * R ≤ 77 * S) :
    32 ≤ 231 * S := by
  have hmL : (1 : ℝ) / 7 ≤ m * L := by
    nlinarith [hneedle]
  have hR : (7 : ℝ) / 24 ≤ R := by
    nlinarith [hpair]
  have hmL_nonneg : 0 ≤ m * L :=
    le_trans (by norm_num : (0 : ℝ) ≤ 1 / 7) hmL
  have hproduct := mul_le_mul hmL hR (by norm_num : (0 : ℝ) ≤ 7 / 24) hmL_nonneg
  have hfloor : (1 : ℝ) / 24 ≤ m * L * R := by
    calc
      (1 : ℝ) / 24 = ((1 : ℝ) / 7) * (7 / 24) := by ring
      _ ≤ (m * L) * R := hproduct
      _ = m * L * R := by ring
  nlinarith [hmetric, hfloor]

/-! ## Seven-line averaging -/

/-- An explicit, kernel-small seven-line pigeonhole step.  If the sum of the
seven nonnegative reciprocal line scales is at least `32/231`, then one line
contributes at least `32/1617`.  Nonnegativity is not needed for this purely
linear conclusion. -/
theorem seven_line_averaging
    (x0 x1 x2 x3 x4 x5 x6 : ℝ)
    (hsum : 32 ≤ 231 * (x0 + x1 + x2 + x3 + x4 + x5 + x6)) :
    32 ≤ 1617 * x0 ∨ 32 ≤ 1617 * x1 ∨ 32 ≤ 1617 * x2 ∨
      32 ≤ 1617 * x3 ∨ 32 ≤ 1617 * x4 ∨ 32 ≤ 1617 * x5 ∨
      32 ≤ 1617 * x6 := by
  by_contra h
  push Not at h
  rcases h with ⟨h0, h1, h2, h3, h4, h5, h6⟩
  nlinarith

/-! ## Axiom audit -/

#print axioms q_table
#print axioms one_le_Q
#print axioms qNumerator_eq_sixteen_iff
#print axioms pair_mass_of_thirty_five_triples
#print axioms uncovered_floor_of_pair_mass
#print axioms common_dilate_scale_bound
#print axioms exact_line_mean
#print axioms exact_remainder_error
#print axioms fano_budget_from_cleared_inequalities
#print axioms seven_line_averaging

end SevenWallFanoGCD
end LRC14
