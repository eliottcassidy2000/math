import Mathlib.Tactic

/-!
# Weak-turn two-owner collapse (THM-1278)

This module checks the arithmetic/combinatorial consumers of the paper
theorem.  Above `h/c = 399/5`, the THM-1198 load envelope forbids three
lower owners above `h/6`.  A fastest return packet with payment below `1/h`
uses only owners above `h/6`; its deletion-minimal two-owner word therefore
has type `AB`, `BA`, or `BAB`.  Strong turns and floods retain the ordinary
one-fastest-tooth payment, so failure of that tax forces a bank of these
weak cells.

The exact one-comb envelope, extraction of the deletion-minimal interval
word, and the two-comb component/minimality geometry are explicit paper
providers.  No statement here assumes or proves LRC(14).
-/

namespace LRCWeakTurnTwoOwnerCollapse

/-- The exact three-large-owner scalar envelope touches one at `x=399/5`
and is at most one on the whole ray to its right. -/
theorem threeLargeOwnerThreshold (x : ℚ) (hx : (399 : ℚ) / 5 ≤ x) :
    (121 : ℚ) / 126 + 19 / (6 * x) ≤ 1 := by
  have hxpos : (0 : ℚ) < 6 * x := by nlinarith
  have htail : (19 : ℚ) / (6 * x) ≤ 5 / 126 := by
    rw [div_le_iff₀ hxpos]
    nlinarith
  linarith

theorem thresholdIdentity :
    (121 : ℚ) / 126 + 19 / (6 * ((399 : ℚ) / 5)) = 1 := by
  norm_num

/-- Abstract load consumer.  The paper envelope is strict because each of
the three owners satisfies `d_i>h/6`, so a cover load at least one is
impossible at and above the threshold. -/
theorem threeLargeOwnerLoadContradiction
    (x total : ℚ)
    (hx : (399 : ℚ) / 5 ≤ x)
    (hcover : 1 ≤ total)
    (henvelope : total < (121 : ℚ) / 126 + 19 / (6 * x)) :
    False := by
  have hbound := threeLargeOwnerThreshold x hx
  linarith

/-- In a weak two-tooth packet, neither internal speed can be at or below
`h/6`: one ratio would contribute at least six and the other strictly more
than one. -/
theorem weakTwoTermOwnersAboveSixth
    (h j k : ℚ)
    (hj : 0 < j) (hk : 0 < k)
    (hjh : j < h) (hkh : k < h)
    (hweak : h / j + h / k < 7) :
    h < 6 * j ∧ h < 6 * k := by
  constructor
  · by_contra hnot
    have hj6 : 6 * j ≤ h := le_of_not_gt hnot
    have hratioJ : (6 : ℚ) ≤ h / j := (le_div_iff₀ hj).2 hj6
    have hratioK : (1 : ℚ) < h / k := (lt_div_iff₀ hk).2 (by nlinarith)
    linarith
  · by_contra hnot
    have hk6 : 6 * k ≤ h := le_of_not_gt hnot
    have hratioK : (6 : ℚ) ≤ h / k := (le_div_iff₀ hk).2 hk6
    have hratioJ : (1 : ℚ) < h / j := (lt_div_iff₀ hj).2 (by nlinarith)
    linarith

/-- Seven internal teeth are automatically strong because every reciprocal
ratio `h/j_q` is strictly larger than one. -/
theorem sevenInternalTermsStrong
    (u₁ u₂ u₃ u₄ u₅ u₆ u₇ : ℚ)
    (h₁ : 1 < u₁) (h₂ : 1 < u₂) (h₃ : 1 < u₃)
    (h₄ : 1 < u₄) (h₅ : 1 < u₅) (h₆ : 1 < u₆)
    (h₇ : 1 < u₇) :
    7 < u₁ + u₂ + u₃ + u₄ + u₅ + u₆ + u₇ := by
  linarith

/-- Abstract cardinality form of the same observation. -/
theorem weakPacketLengthCap (s : ℕ) (S : ℚ)
    (hcount : (s : ℚ) < S) (hweak : S < 7) :
    s ≤ 6 := by
  have hsQ : (s : ℚ) < 7 := hcount.trans hweak
  have hsN : s < 7 := by exact_mod_cast hsQ
  omega

/-- Two-comb component geometry supplies at most one slower tooth and
minimality at most two faster boundary teeth.  A turn therefore has length
two or three. -/
theorem minimalTwoOwnerLength
    (aCount bCount length : ℕ)
    (ha : aCount ≤ 1) (hb : bCount ≤ 2)
    (hlen : length = aCount + bCount)
    (hturn : 2 ≤ length) :
    length = 2 ∨ length = 3 := by
  omega

/-- The only length-three count vector is one slower tooth and two faster
teeth; adjacent-owner distinctness then fixes the word `B-A-B`. -/
theorem lengthThreeCounts
    (aCount bCount : ℕ)
    (ha : aCount ≤ 1) (hb : bCount ≤ 2)
    (hsum : aCount + bCount = 3) :
    aCount = 1 ∧ bCount = 2 := by
  omega

/-- Split THM-1275's forced exceptions into floods, strong turns, and weak
turns.  Only the weak count can remain unpaid at the one-tooth scale. -/
theorem paidExceptionFunnel
    (forced weak floods strong : ℕ)
    (hall : forced ≤ floods + strong + weak) :
    forced - weak ≤ floods + strong := by
  omega

theorem weakOnlyResidual
    (forced weak floods strong : ℕ)
    (hall : forced ≤ floods + strong + weak)
    (hunpaid : floods + strong = 0) :
    forced ≤ weak := by
  omega

/-- Exact examples and the self-similar abutment boundary used by the
referee. -/
theorem exactWeakAndAbutmentExamples :
    (6 : ℚ) < 60 / 11 + 60 / 50 ∧
    60 / 11 + 60 / 50 < 7 ∧
    (6 : ℚ) < 60 / 20 + 2 * (60 / 39) ∧
    60 / 20 + 2 * (60 / 39) < 7 ∧
    60 / 15 + 60 / 30 = 6 := by
  norm_num

#print axioms threeLargeOwnerThreshold
#print axioms thresholdIdentity
#print axioms threeLargeOwnerLoadContradiction
#print axioms weakTwoTermOwnersAboveSixth
#print axioms sevenInternalTermsStrong
#print axioms weakPacketLengthCap
#print axioms minimalTwoOwnerLength
#print axioms lengthThreeCounts
#print axioms paidExceptionFunnel
#print axioms weakOnlyResidual
#print axioms exactWeakAndAbutmentExamples

end LRCWeakTurnTwoOwnerCollapse
