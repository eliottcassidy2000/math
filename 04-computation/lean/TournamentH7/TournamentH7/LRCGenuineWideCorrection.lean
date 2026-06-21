/-
  TournamentH7.LRCGenuineWideCorrection -- finite arithmetic kernel for the
  HYP-2805/kpswf9 correction to the genuine-wide robust-margin claim.

  The exact search table lives in
  `05-knowledge/results/lrc14_genuine_wide_true_maximizer_kpswf9.out`.
  This module records only the small rational ledger:

    * the reported adjacent-doublet true-max rows are all below `cap_k`;
    * the k=10 row has the smallest cap margin among the five rows;
    * the stronger margin target `cap_k - p0 >= 4/25` fails at k=10.

  This is a correction kernel, not a full generalized-doublet theorem.  The
  remaining leg-C proof obligations are the generalized-doublet reduction,
  frozen-room bound, uniform R-tail, and optimized finite-window certificate.
-/

namespace LonelyRunner
namespace GenuineWideCorrection

/-- A reported exact true-max row from the kpswf9 adjacent-doublet correction
table.  The rationals are `p0 = p0Num / p0Den` and `cap = capNum / capDen`.
`basePrimitive` records the table's base-primitivity flag; the actual LRC
predicate uses primitivity of the full set. -/
structure MaxRow where
  k : Nat
  p0Num : Int
  p0Den : Int
  capNum : Int
  capDen : Int
  basePrimitive : Bool
deriving DecidableEq, Repr

def marginNum (r : MaxRow) : Int :=
  r.capNum * r.p0Den - r.p0Num * r.capDen

def marginDen (r : MaxRow) : Int :=
  r.capDen * r.p0Den

/-- Numerator of `(cap - p0) - 4/25`.  Positive means the robust `0.16`
margin holds; negative means it fails. -/
def robustNum (r : MaxRow) : Int :=
  25 * marginNum r - 4 * marginDen r

/-- Cross-multiplied margin comparison as a finite audit Boolean. -/
def marginLE (a b : MaxRow) : Bool :=
  decide (marginNum a * marginDen b <= marginNum b * marginDen a)

def k8 : MaxRow :=
  { k := 8, p0Num := 348487, p0Den := 1681680,
    capNum := 2243, capDen := 5880, basePrimitive := true }

def k9 : MaxRow :=
  { k := 9, p0Num := 321, p0Den := 980,
    capNum := 1979, capDen := 4004, basePrimitive := false }

def k10 : MaxRow :=
  { k := 10, p0Num := 265, p0Den := 588,
    capNum := 55, capDen := 91, basePrimitive := false }

def k11 : MaxRow :=
  { k := 11, p0Num := 5617, p0Den := 10780,
    capNum := 66, capDen := 91, basePrimitive := true }

def k12 : MaxRow :=
  { k := 12, p0Num := 347, p0Den := 588,
    capNum := 6, capDen := 7, basePrimitive := true }

def maxRow : Fin 5 -> MaxRow
  | ⟨0, _⟩ => k8
  | ⟨1, _⟩ => k9
  | ⟨2, _⟩ => k10
  | ⟨3, _⟩ => k11
  | _ => k12

theorem all_reported_rows_below_cap :
    ∀ r : Fin 5, 0 < marginNum (maxRow r) := by
  native_decide

/-- The k=10 correction row is the smallest margin among the reported rows. -/
theorem k10_is_smallest_reported_margin :
    marginLE k10 k8 = true ∧ marginLE k10 k9 = true ∧
      marginLE k10 k10 = true ∧ marginLE k10 k11 = true ∧
      marginLE k10 k12 = true := by
  native_decide

/-- The robust `0.16` margin target fails exactly at the k=10 correction row
among these reported rows. -/
theorem robust_margin_flags :
    0 <= robustNum k8 ∧ 0 <= robustNum k9 ∧ robustNum k10 < 0 ∧
      0 <= robustNum k11 ∧ 0 <= robustNum k12 := by
  native_decide

/-- Records the base-primitivity guardrail: k=9 and k=10 use non-primitive
bases even though the full sets are primitive in the search table. -/
theorem nonprimitive_base_guardrail :
    k9.basePrimitive = false ∧ k10.basePrimitive = false := by
  native_decide

/-! ### Axiom audit -/

#print axioms all_reported_rows_below_cap
#print axioms k10_is_smallest_reported_margin
#print axioms robust_margin_flags
#print axioms nonprimitive_base_guardrail

end GenuineWideCorrection
end LonelyRunner
