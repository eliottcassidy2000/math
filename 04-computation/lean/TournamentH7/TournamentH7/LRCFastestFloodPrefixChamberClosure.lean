import Mathlib.Tactic

/-!
# Fastest-flood prefix-chamber arithmetic (THM-1272)

This module kernel-checks the arithmetic consumers of the paper/referee
theorem: tooth-count chamber endpoints, the `e=0` integer cut, the three
strict prefix-count jumps, the basic exception counts and five boundary
margins, the abstract envelope/tax contradiction and ray monotonicity, and
the functional `s=113/114` guardrail behind the improved cut `h/c<798`.

External providers are THM-1198's exact/BV envelope, THM-1233's prefix
component rows, the finite prefix dynamic-program exhaustiveness, and
THM-1275's dominated flood/turn tax.  No covering theorem is hidden in an
axiom here.
-/

namespace LRCFastestFloodPrefixChamberClosure

def ceilDiv (a b : ℕ) : ℕ := (a + b - 1) / b

/-- The real interval carrying tooth-count chamber `m` has these two affine
endpoints.  The left endpoint is open in the paper statement. -/
def chamberLeft (m : ℚ) : ℚ := (7 * m - 8) / 6
def chamberRight (m : ℚ) : ℚ := (7 * m - 1) / 6

theorem chamberEndpointWidth (m : ℚ) :
    chamberRight m - chamberLeft m = 7 / 6 := by
  simp [chamberLeft, chamberRight]
  ring

theorem firstChamberEndpoints :
    chamberLeft 2 = 1 ∧ chamberRight 2 = 13 / 6 ∧
    chamberLeft 3 = 13 / 6 ∧ chamberRight 3 = 10 / 3 ∧
    chamberLeft 4 = 10 / 3 ∧ chamberRight 4 = 9 / 2 := by
  norm_num [chamberLeft, chamberRight]

theorem chamberSupremumIdentities :
    (7 : ℚ) / 36 = 7 / 36 ∧
    (8 : ℚ) / 45 = 8 / 45 ∧
    (1 : ℚ) / 7 + 1 / 21 = 4 / 21 ∧
    (1 : ℚ) / 7 + (1 : ℚ) / 27 = (34 : ℚ) / 189 := by
  norm_num

/-- The integer form of THM-1267 plus `e=0`, i.e. `h ≤ 6d1`, gives the
strict ratio cut without division. -/
theorem eZeroIntegerCut (c d1 h : ℤ)
    (hfast : h ≤ 6 * d1)
    (hslow : 270 * d1 ≤ 563 * c - 1) :
    45 * h < 563 * c := by
  omega

/-- At the `e=3` boundary, strict component containment moves the integer
prefix sum from ten to eleven. -/
theorem eThreeCountJump (c h s : ℕ)
    (hlower : 1078 * c ≤ 5 * h)
    (hupper : 5 * h < 98 * c * (s + 1)) :
    11 ≤ s := by
  by_contra hnot
  have hs : s ≤ 10 := by omega
  have hs' : s + 1 ≤ 11 := by omega
  have hcap : 98 * c * (s + 1) ≤ 1078 * c := by
    calc
      98 * c * (s + 1) ≤ 98 * c * 11 := Nat.mul_le_mul_left (98 * c) hs'
      _ = 1078 * c := by ring
  have : 5 * h < 1078 * c := hupper.trans_le hcap
  omega

/-- The corresponding `e=4` jump is from 41 to 42. -/
theorem eFourCountJump (c h s : ℕ)
    (hlower : 1323 * c ≤ 2 * h)
    (hupper : 4 * h < 63 * c * (s + 1)) :
    42 ≤ s := by
  by_contra hnot
  have hs : s ≤ 41 := by omega
  have hs' : s + 1 ≤ 42 := by omega
  have hcap : 63 * c * (s + 1) ≤ 2646 * c := by
    calc
      63 * c * (s + 1) ≤ 63 * c * 42 := Nat.mul_le_mul_left (63 * c) hs'
      _ = 2646 * c := by ring
  have hupper' : 4 * h < 2646 * c := hupper.trans_le hcap
  have hlower' : 2646 * c ≤ 4 * h := by omega
  omega

/-- The basic `e=5` jump at `x=959`. -/
theorem eFiveBasicCountJump (c h s : ℕ)
    (hlower : 959 * c ≤ h)
    (hupper : h < 7 * c * (s + 1)) :
    137 ≤ s := by
  by_contra hnot
  have hs : s ≤ 136 := by omega
  have hs' : s + 1 ≤ 137 := by omega
  have hcap : 7 * c * (s + 1) ≤ 959 * c := by
    calc
      7 * c * (s + 1) ≤ 7 * c * 137 := Nat.mul_le_mul_left (7 * c) hs'
      _ = 959 * c := by ring
  have : h < 959 * c := hupper.trans_le hcap
  omega

/-- The functional `e=5` analysis starts at the earlier chamber `798=7*114`. -/
theorem eFiveFunctionalCountJump (c h s : ℕ)
    (hlower : 798 * c ≤ h)
    (hupper : h < 7 * c * (s + 1)) :
    114 ≤ s := by
  by_contra hnot
  have hs : s ≤ 113 := by omega
  have hs' : s + 1 ≤ 114 := by omega
  have hcap : 7 * c * (s + 1) ≤ 798 * c := by
    calc
      7 * c * (s + 1) ≤ 7 * c * 114 := Nat.mul_le_mul_left (7 * c) hs'
      _ = 798 * c := by ring
  have : h < 798 * c := hupper.trans_le hcap
  omega

theorem namedPrefixLoadSums :
    (7 : ℚ) / 36 + 4 / 21 + (1 / 7 + 1 / 27) = 61 / 108 ∧
    (7 : ℚ) / 36 + 4 / 21 + (1 / 7 + 1 / 27) +
        (1 / 7 + 1 / 209) = 112571 / 158004 ∧
    (7 : ℚ) / 36 + 4 / 21 + (1 / 7 + 1 / 27) +
        (1 / 7 + 1 / 209) + (1 / 7 + 1 / 657) =
          9882995 / 11534292 ∧
    (7 : ℚ) / 36 + 4 / 21 + 4 / 21 +
        (1 / 7 + 1 / 167) + (1 / 7 + 1 / 545) =
          2847097 / 3276540 := by
  norm_num

theorem basicExceptionCounts :
    ceilDiv 7 2 - 1 = 3 ∧
    ceilDiv 11 3 - 1 = 3 ∧
    ceilDiv 31 4 - 1 = 7 ∧
    ceilDiv 95 5 - 1 = 18 ∧
    ceilDiv 137 6 - 1 = 22 := by
  norm_num [ceilDiv]

/-- These are the five exact cross-multiplied boundary margins.  The second
row touches zero; the paper upper envelope is strict there. -/
theorem basicBoundaryMargins :
    (3 : ℚ) / 8 -
        (42 * (7 / 36 - 2 / 7) + 25 / 6) = 1 / 24 ∧
    (3 : ℚ) / 8 -
        ((1407 / 20) * (7 / 18 - 3 / 7) + 19 / 6) = 0 ∧
    (7 : ℚ) / 8 -
        ((1078 / 5) * (61 / 108 - 4 / 7) + 13 / 6) = 29 / 216 ∧
    (18 : ℚ) / 8 -
        ((1323 / 2) * (112571 / 158004 - 5 / 7) + 7 / 6) =
          11503 / 5016 ∧
    (22 : ℚ) / 8 -
        (959 * (9882995 / 11534292 - 6 / 7) + 1 / 6) =
          1185455 / 411939 := by
  norm_num

theorem basicNegativeRaySlopes :
    (7 : ℚ) / 36 - 2 / 7 < 0 ∧
    (7 : ℚ) / 18 - 3 / 7 < 0 ∧
    (61 : ℚ) / 108 - 4 / 7 < 0 ∧
    (112571 : ℚ) / 158004 - 5 / 7 < 0 ∧
    (9882995 : ℚ) / 11534292 - 6 / 7 < 0 := by
  norm_num

/-- Abstract endpoint consumer: a strict phase-free upper bound below a
forced tail tax is incompatible with a cover. -/
theorem envelopeTaxContradiction (F upper tax : ℚ)
    (hUpper : F < upper)
    (hTax : tax ≤ F)
    (hCompare : upper ≤ tax) :
    False := by
  linarith

/-- Once the cross-multiplied coefficient is negative, a verified boundary
comparison persists to the right until the next discrete improvement. -/
theorem negativeSlopeRay (x x0 slope constant target : ℚ)
    (hx : x0 ≤ x)
    (hSlope : slope < 0)
    (hBoundary : x0 * slope + constant ≤ target) :
    x * slope + constant ≤ target := by
  nlinarith

theorem functionalS114EtaIdentity :
    1 - (2847097 : ℚ) / 3276540 = 429443 / 3276540 := by
  norm_num

/-- Strict prefix load makes the private-count product lie just above this
number, so its ceiling is at least 538. -/
theorem functionalS114CountBracket :
    (537 : ℚ) < 36 * 114 * (1 - 2847097 / 3276540) ∧
    (36 : ℚ) * 114 * (1 - 2847097 / 3276540) < 538 := by
  norm_num

theorem functionalS114ExceptionCount : ceilDiv 538 6 - 1 = 89 := by
  norm_num [ceilDiv]

/-- The preceding chamber fails this exact certificate, while the `s=114`
chamber has a positive worst-endpoint margin. -/
theorem functionalGuardrailMargins :
    (86 : ℚ) / 8 -
        (798 * (868922 / 995715 - 6 / 7) + 1 / 6) =
          -113823 / 63220 ∧
    (89 : ℚ) / 8 -
        (805 * (2847097 / 3276540 - 6 / 7) + 1 / 6) =
          1921973 / 1310616 ∧
    (0 : ℚ) < 1921973 / 1310616 := by
  norm_num

#print axioms chamberEndpointWidth
#print axioms firstChamberEndpoints
#print axioms eZeroIntegerCut
#print axioms eThreeCountJump
#print axioms eFourCountJump
#print axioms eFiveBasicCountJump
#print axioms eFiveFunctionalCountJump
#print axioms namedPrefixLoadSums
#print axioms basicExceptionCounts
#print axioms basicBoundaryMargins
#print axioms basicNegativeRaySlopes
#print axioms envelopeTaxContradiction
#print axioms negativeSlopeRay
#print axioms functionalS114CountBracket
#print axioms functionalS114ExceptionCount
#print axioms functionalGuardrailMargins

end LRCFastestFloodPrefixChamberClosure
