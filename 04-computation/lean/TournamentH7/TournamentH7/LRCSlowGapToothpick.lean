/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-18)
-/
import Mathlib.Tactic

/-!
# Arithmetic consumers for the slow-gap toothpick ladder

This module formalizes the denominator-cleared arithmetic layer of THM-1176.
The geometric production of a covered slow gap, the one-comb discrepancy
estimate, and the recursive nesting of slow gaps remain explicit hypotheses.
Every result below is kernel-checked without proof placeholders.
-/

open scoped BigOperators

namespace LRC14
namespace SlowGapToothpick

/-! ## Covered-length pressure -/

/-- Clearing the denominator in
`L ≤ rL/7 + endpoint` gives the residual endpoint pressure. -/
theorem endpoint_pressure_from_cover (L endpoint r : ℝ)
    (hcover : 7 * L ≤ r * L + 7 * endpoint) :
    (7 - r) * L ≤ 7 * endpoint := by
  linarith

/-- Specialization to a slow gap of length `6/(7c)` and total endpoint term
`(6/49)S`.  The hypotheses are exactly the two denominator-cleared analytic
inputs; the conclusion is the normalized pressure `cS ≥ 7-r`. -/
theorem slow_gap_pressure_from_cleared_cover (c L S r : ℝ) (hc : 0 ≤ c)
    (hlength : 7 * c * L = 6)
    (hcover : 49 * L ≤ 7 * r * L + 6 * S) :
    7 - r ≤ c * S := by
  have hgap : 7 * (7 - r) * L ≤ 6 * S := by
    nlinarith [hcover]
  have hscaled := mul_le_mul_of_nonneg_left hgap hc
  have hlength_scaled := congrArg (fun x : ℝ ⇒ (7 - r) * x) hlength
  nlinarith [hscaled, hlength_scaled]

/-! ## Equality parity obstruction -/

/-- For every cardinality `0,…,7`, the parities of `r` and `7-r` differ. -/
theorem seven_complement_parity_mismatch (r : ℕ) (hr : r ≤ 7) :
    r % 2 ≠ (7 - r) % 2 := by
  omega

/-- Abstract cleared parity contradiction from THM-1176.  Here `A` has the
parity of a sum of `r` odd complementary products, while `P` is odd.  Thus
`3A = (7-r)P` is impossible. -/
theorem cleared_equality_parity_impossible (r A P : ℕ) (hr : r ≤ 7)
    (hA : A % 2 = r % 2) (hP : P % 2 = 1) :
    3 * A ≠ (7 - r) * P := by
  intro heq
  have hleft : (3 * A) % 2 = r % 2 := by
    rw [Nat.mul_mod, hA]
    omega
  have hright : ((7 - r) * P) % 2 = (7 - r) % 2 := by
    rw [Nat.mul_mod, hP]
    omega
  have hsides : r % 2 = (7 - r) % 2 := by
    calc
      r % 2 = (3 * A) % 2 := hleft.symm
      _ = ((7 - r) * P) % 2 := congrArg (fun x : ℕ ⇒ x % 2) heq
      _ = (7 - r) % 2 := hright
  exact seven_complement_parity_mismatch r hr hsides

/-! ## Three-or-fewer faster combs -/

/-- Every normalized faster-speed ratio is below one, so the total pressure
of a nonempty `r`-packet is strictly below `r`. -/
theorem normalized_pressure_lt_card {r : ℕ} (hr : 0 < r)
    (c : ℝ) (d : Fin r → ℝ) (hc : 0 ≤ c) (hd : ∀ i, c < d i) :
    (∑ i, c / d i) < r := by
  have hterm : ∀ i, c / d i < (1 : ℝ) := by
    intro i
    exact (div_lt_one (lt_of_le_of_lt hc (hd i))).2 (hd i)
  have hsum : (∑ i, c / d i) < ∑ _i : Fin r, (1 : ℝ) := by
    apply Finset.sum_lt_sum
    · intro i _hi
      exact (hterm i).le
    · let i0 : Fin r := ⟨0, hr⟩
      exact ⟨i0, Finset.mem_univ i0, hterm i0⟩
  simpa using hsum

/-- The strict pressure law `7-r < Σ c/d_i` is incompatible with a packet
of at most three genuinely faster speeds.  The strict pressure law itself is
an explicit input supplied by the endpoint-parity argument. -/
theorem at_most_three_impossible_of_strict_pressure {r : ℕ}
    (hrpos : 0 < r) (hr3 : r ≤ 3) (c : ℝ) (d : Fin r → ℝ)
    (hc : 0 ≤ c) (hd : ∀ i, c < d i)
    (hpressure : (7 : ℝ) - r < ∑ i, c / d i) : False := by
  have hsum := normalized_pressure_lt_card hrpos c d hc hd
  have hr_real : (r : ℝ) ≤ 3 := by
    exact_mod_cast hr3
  linarith

/-! ## The exact six-speed distinctness sharpening -/

/-- Six consecutive denominators beginning at `6a-2` have normalized
reciprocal pressure strictly below one. -/
def sixEnvelope (a : ℚ) : ℚ :=
  a / (6 * a - 2) + a / (6 * a - 1) + a / (6 * a) +
    a / (6 * a + 1) + a / (6 * a + 2) + a / (6 * a + 3)

/-- The positive numerator used after clearing the six-envelope deficit. -/
def sixEnvelopeNumerator (x : ℚ) : ℚ :=
  3 * (x - 6) ^ 4 + 62 * (x - 6) ^ 3 + 423 * (x - 6) ^ 2 +
    988 * (x - 6) + 264

/-- Exact all-scale proof of the six-consecutive envelope inequality recorded
in THM-1176(19a)--(19b). -/
theorem six_consecutive_envelope_lt_one (a : ℚ) (ha : 1 ≤ a) :
    sixEnvelope a < 1 := by
  have hm2 : 0 < 6 * a - 2 := by linarith
  have hm1 : 0 < 6 * a - 1 := by linarith
  have h0 : 0 < 6 * a := by linarith
  have hp1 : 0 < 6 * a + 1 := by linarith
  have hp2 : 0 < 6 * a + 2 := by linarith
  have hp3 : 0 < 6 * a + 3 := by linarith
  have hidentity :
      1 - sixEnvelope a =
        sixEnvelopeNumerator (6 * a) /
          (6 * (6 * a - 2) * (6 * a - 1) * (6 * a + 1) *
            (6 * a + 2) * (6 * a + 3)) := by
    field_simp [sixEnvelope, sixEnvelopeNumerator, ne_of_gt hm2, ne_of_gt hm1,
      ne_of_gt h0, ne_of_gt hp1, ne_of_gt hp2, ne_of_gt hp3]
    <;> ring
  have hy : 0 ≤ 6 * a - 6 := by linarith
  have hnum : 0 < sixEnvelopeNumerator (6 * a) := by
    unfold sixEnvelopeNumerator
    positivity
  have hden :
      0 < 6 * (6 * a - 2) * (6 * a - 1) * (6 * a + 1) *
        (6 * a + 2) * (6 * a + 3) := by
    positivity
  have hdeficit : 0 < 1 - sixEnvelope a := by
    rw [hidentity]
    exact div_pos hnum hden
  linarith

/-- If six distinct increasing integer speeds carry normalized reciprocal
pressure above one, their first speed is at most `6a-3`. -/
theorem first_speed_le_six_a_sub_three
    (a b0 b1 b2 b3 b4 b5 : ℕ) (ha : 1 ≤ a) (hab : a < b0)
    (h01 : b0 < b1) (h12 : b1 < b2) (h23 : b2 < b3)
    (h34 : b3 < b4) (h45 : b4 < b5)
    (hpressure :
      (1 : ℚ) < (a : ℚ) / b0 + (a : ℚ) / b1 + (a : ℚ) / b2 +
        (a : ℚ) / b3 + (a : ℚ) / b4 + (a : ℚ) / b5) :
    b0 ≤ 6 * a - 3 := by
  by_contra hcut
  have hn0 : 6 * a ≤ b0 + 2 := by omega
  have hn1 : 6 * a ≤ b1 + 1 := by omega
  have hn2 : 6 * a ≤ b2 := by omega
  have hn3 : 6 * a + 1 ≤ b3 := by omega
  have hn4 : 6 * a + 2 ≤ b4 := by omega
  have hn5 : 6 * a + 3 ≤ b5 := by omega
  let A : ℚ := a
  have hA1 : 1 ≤ A := by exact_mod_cast ha
  have hA0 : 0 ≤ A := le_trans (by norm_num) hA1
  have hb0 : 6 * A - 2 ≤ (b0 : ℚ) := by
    have : (6 : ℚ) * a ≤ b0 + 2 := by exact_mod_cast hn0
    dsimp [A]
    linarith
  have hb1 : 6 * A - 1 ≤ (b1 : ℚ) := by
    have : (6 : ℚ) * a ≤ b1 + 1 := by exact_mod_cast hn1
    dsimp [A]
    linarith
  have hb2 : 6 * A ≤ (b2 : ℚ) := by exact_mod_cast hn2
  have hb3 : 6 * A + 1 ≤ (b3 : ℚ) := by exact_mod_cast hn3
  have hb4 : 6 * A + 2 ≤ (b4 : ℚ) := by exact_mod_cast hn4
  have hb5 : 6 * A + 3 ≤ (b5 : ℚ) := by exact_mod_cast hn5
  have ht0 : A / (b0 : ℚ) ≤ A / (6 * A - 2) :=
    div_le_div_of_nonneg_left hA0 (by linarith) hb0
  have ht1 : A / (b1 : ℚ) ≤ A / (6 * A - 1) :=
    div_le_div_of_nonneg_left hA0 (by linarith) hb1
  have ht2 : A / (b2 : ℚ) ≤ A / (6 * A) :=
    div_le_div_of_nonneg_left hA0 (by linarith) hb2
  have ht3 : A / (b3 : ℚ) ≤ A / (6 * A + 1) :=
    div_le_div_of_nonneg_left hA0 (by linarith) hb3
  have ht4 : A / (b4 : ℚ) ≤ A / (6 * A + 2) :=
    div_le_div_of_nonneg_left hA0 (by linarith) hb4
  have ht5 : A / (b5 : ℚ) ≤ A / (6 * A + 3) :=
    div_le_div_of_nonneg_left hA0 (by linarith) hb5
  have hsum :
      A / (b0 : ℚ) + A / (b1 : ℚ) + A / (b2 : ℚ) +
          A / (b3 : ℚ) + A / (b4 : ℚ) + A / (b5 : ℚ) ≤
        sixEnvelope A := by
    dsimp [sixEnvelope]
    linarith
  have henvelope := six_consecutive_envelope_lt_one A hA1
  dsimp [A] at hsum henvelope
  linarith

/-! ## Exact recursive cutoffs and the toothpick ladder -/

/-- Denominator-cleared form of the five-speed integer cutoff. -/
theorem five_speed_cutoff_iff (c d : ℕ) (hc : 1 ≤ c) :
    d ≤ (5 * c - 4) / 2 ↔ 2 * d ≤ 5 * c - 4 := by
  omega

/-- Denominator-cleared form of the four-speed integer cutoff. -/
theorem four_speed_cutoff_iff (c d : ℕ) (hc : 2 ≤ c) :
    d ≤ (8 * c - 9) / 6 ↔ 6 * d ≤ 8 * c - 9 := by
  omega

/-- The five-speed cutoff implies the strict normalized ratio `d/c < 5/2`. -/
theorem five_speed_cutoff_ratio (c d : ℕ) (hc : 1 ≤ c)
    (hcut : 2 * d ≤ 5 * c - 4) :
    (d : ℚ) / c < 5 / 2 := by
  have hcQ : (0 : ℚ) < c := by exact_mod_cast (lt_of_lt_of_le Nat.zero_lt_one hc)
  rw [div_lt_iff₀ hcQ]
  have hcutQ : (2 : ℚ) * d ≤ 5 * c - 4 := by exact_mod_cast hcut
  linarith

/-- The four-speed cutoff implies the sharper strict ratio `d/c < 4/3`. -/
theorem four_speed_cutoff_ratio (c d : ℕ) (hc : 2 ≤ c)
    (hcut : 6 * d ≤ 8 * c - 9) :
    (d : ℚ) / c < 4 / 3 := by
  have hcQ : (0 : ℚ) < c := by
    exact_mod_cast (lt_of_lt_of_le Nat.zero_lt_one (le_trans Nat.one_le_two hc))
  rw [div_lt_iff₀ hcQ]
  have hcutQ : (6 : ℚ) * d ≤ 8 * c - 9 := by exact_mod_cast hcut
  linarith

/-- Cleared three-rung ladder.  The geometric nesting and four-speed cutoff
are explicit in `hnested`; the conclusion is the ratio alternative
`b1/a < 13/6` or `b2/b1 < 13/6` or `b3/b2 < 4/3`. -/
theorem toothpick_ladder_of_nested_four_speed_cutoff
    (a b1 b2 b3 : ℕ)
    (hnested : 13 * a ≤ 6 * b1 → 13 * b1 ≤ 6 * b2 →
      6 * b3 ≤ 8 * b2 - 9) :
    6 * b1 < 13 * a ∨ 6 * b2 < 13 * b1 ∨ 3 * b3 < 4 * b2 := by
  by_cases hfirst : 13 * a ≤ 6 * b1
  · by_cases hsecond : 13 * b1 ≤ 6 * b2
    · right
      right
      have hcut := hnested hfirst hsecond
      omega
    · right
      left
      omega
  · left
    omega

/-! ## Harmonic-drift consumer -/

/-- If the remaining pressure exceeds a positive deficit `R`, while order
bounds it by `q*a/b`, then clearing the positive denominator gives
`bR < qa`. -/
theorem harmonic_drift_cleared (a b R tail q : ℝ) (hb : 0 < b)
    (hdeficit : R < tail) (horder : b * tail ≤ q * a) :
    b * R < q * a := by
  have hscaled := mul_lt_mul_of_pos_left hdeficit hb
  exact hscaled.trans_le horder

/-! ## Axiom audit -/

#print axioms endpoint_pressure_from_cover
#print axioms slow_gap_pressure_from_cleared_cover
#print axioms cleared_equality_parity_impossible
#print axioms normalized_pressure_lt_card
#print axioms at_most_three_impossible_of_strict_pressure
#print axioms six_consecutive_envelope_lt_one
#print axioms first_speed_le_six_a_sub_three
#print axioms five_speed_cutoff_iff
#print axioms four_speed_cutoff_iff
#print axioms toothpick_ladder_of_nested_four_speed_cutoff
#print axioms harmonic_drift_cleared

end SlowGapToothpick
end LRC14
