/-
  TournamentH7.LRCQ6Contraction -- arithmetic kernel for the q6 contraction
  endpoint-period certificate.

  The analytic content lives in `lrc_q6_ratio_periodicity_macmini_0622s23.py`:
  for a bounded base `B` and far speed `f`, the all-missed atom satisfies

      q6(B ∪ {f}) / q6(B) <= 1/7 + periodmax_q6(B) / (f * q6(B)).

  This module records the exact rational arithmetic for the currently reported
  q6-ratio bounds.  It does not formalize the sawtooth identity or the endpoint
  period scan.
-/

namespace LonelyRunner
namespace Q6Contraction

/-- Common denominator for the consecutive k=9/k=10 endpoint-ratio bounds:
`1/7 + (6/49)/(15*q6(B))`. -/
def commonDen : Nat := 7 * 49 * 15

/-- Numerator for the k=9 consecutive-base bound, where `q6(B)=1/56`. -/
def consecK9BoundNum : Nat := 49 * 15 + 6 * 56 * 7

/-- Numerator for the k=10 consecutive-base bound, where `q6(B)=1/63`. -/
def consecK10BoundNum : Nat := 49 * 15 + 6 * 63 * 7

/-- Exact reduction of the k=9 consecutive-base q6-ratio bound to `3/5`. -/
theorem consecK9_bound_exact : consecK9BoundNum * 5 = 3 * commonDen := by
  native_decide

/-- Exact reduction of the k=10 consecutive-base q6-ratio bound to `23/35`. -/
theorem consecK10_bound_exact : consecK10BoundNum * 35 = 23 * commonDen := by
  native_decide

/-- The k=9 consecutive-base q6-ratio bound is a strict contraction. -/
theorem consecK9_strict_contraction : (3 : Nat) < 5 := by
  native_decide

/-- The k=10 consecutive-base q6-ratio bound is a strict contraction. -/
theorem consecK10_strict_contraction : (23 : Nat) < 35 := by
  native_decide

/-- The reported 15-base general q6 scout has worst ratio bound `33/35`, still
a strict contraction.  The exact row producing the bound remains in the Python
certificate; this theorem records the final arithmetic comparison. -/
theorem generalScout_strict_contraction : (33 : Nat) < 35 := by
  native_decide

/-! ### Axiom audit -/

#print axioms consecK9_bound_exact
#print axioms consecK10_bound_exact
#print axioms consecK9_strict_contraction
#print axioms consecK10_strict_contraction
#print axioms generalScout_strict_contraction

end Q6Contraction
end LonelyRunner
