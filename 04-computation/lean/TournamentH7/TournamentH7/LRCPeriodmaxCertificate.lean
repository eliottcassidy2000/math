/-
  TournamentH7.LRCPeriodmaxCertificate -- finite arithmetic kernel for the
  completed THM-563 bounded-base period-max audit.

  The full exhaustive row validation lives in
  `05-knowledge/results/lrc_periodmax_THM563_general_check_COMPLETE_macmini_0621s7.out`.
  This Lean module records the exact arithmetic for the six per-k worst rows:

    periodmax(B) < 15 * (cap_k - Plat(B)).

  It does not reimplement the endpoint-period scan.  The point is to expose a
  small, auditable theorem boundary for the LRC14 proof DAG after HYP-2793:
  the integer single-far bounded-base leg is now finite arithmetic, while the
  remaining mathematical work is genuine-wide slack and continuous dilation
  glue.
-/

namespace LonelyRunner
namespace PeriodmaxCertificate

/-- One worst-row certificate.  The represented rationals are
`periodmax = pmNum / pmDen` and `margin = marginNum / marginDen`. -/
structure WorstRow where
  k : Nat
  pmNum : Int
  pmDen : Int
  marginNum : Int
  marginDen : Int
deriving DecidableEq, Repr

/-- Numerator of `15 * margin - periodmax`, over the positive denominator
`marginDen * pmDen`. -/
def headroomNum (r : WorstRow) : Int :=
  15 * r.marginNum * r.pmDen - r.pmNum * r.marginDen

/-- Numerator of `periodmax / margin`. -/
def ratioNum (r : WorstRow) : Int := r.pmNum * r.marginDen

/-- Denominator of `periodmax / margin`. -/
def ratioDen (r : WorstRow) : Int := r.pmDen * r.marginNum

/-- Cross-multiplied ratio comparison as a finite audit Boolean, assuming
positive denominators. -/
def ratioLE (a b : WorstRow) : Bool :=
  decide (ratioNum a * ratioDen b <= ratioNum b * ratioDen a)

def k8 : WorstRow :=
  { k := 8, pmNum := 2, pmDen := 1, marginNum := 1087, marginDen := 5880 }

def k9 : WorstRow :=
  { k := 9, pmNum := 86, pmDen := 49, marginNum := 129643, marginDen := 980980 }

def k10 : WorstRow :=
  { k := 10, pmNum := 3299, pmDen := 1470, marginNum := 7403, marginDen := 35672 }

def k11 : WorstRow :=
  { k := 11, pmNum := 6730439, pmDen := 1961960,
    marginNum := 6205307, marginDen := 17657640 }

def k12 : WorstRow :=
  { k := 12, pmNum := 536399, pmDen := 196196,
    marginNum := 5883083, marginDen := 17657640 }

def k13 : WorstRow :=
  { k := 13, pmNum := 122064, pmDen := 49049,
    marginNum := 734365, marginDen := 1765764 }

/-- The six per-k worst rows in the completed audit. -/
def worstRow : Fin 6 -> WorstRow
  | ⟨0, _⟩ => k8
  | ⟨1, _⟩ => k9
  | ⟨2, _⟩ => k10
  | ⟨3, _⟩ => k11
  | ⟨4, _⟩ => k12
  | _ => k13

theorem k8_headroom_positive : 0 < headroomNum k8 := by
  native_decide

theorem k9_headroom_positive : 0 < headroomNum k9 := by
  native_decide

theorem k10_headroom_positive : 0 < headroomNum k10 := by
  native_decide

theorem k11_headroom_positive : 0 < headroomNum k11 := by
  native_decide

theorem k12_headroom_positive : 0 < headroomNum k12 := by
  native_decide

theorem k13_headroom_positive : 0 < headroomNum k13 := by
  native_decide

/-- Every per-k worst row has positive headroom. -/
theorem all_worst_rows_headroom_positive :
    ∀ r : Fin 6, 0 < headroomNum (worstRow r) := by
  native_decide

/-- Among the six per-k worst-row ratios, the k=9 even AP row is the largest. -/
theorem k9_is_global_worst_among_worst_rows :
    ratioLE k8 k9 = true ∧ ratioLE k9 k9 = true ∧ ratioLE k10 k9 = true ∧
      ratioLE k11 k9 = true ∧ ratioLE k12 k9 = true ∧ ratioLE k13 k9 = true := by
  native_decide

/-- The completed audit's row counts. -/
structure CountRow where
  bases : Nat
  trivial : Nat
  scanned : Nat
  skipped : Nat
  passed : Nat
  failed : Nat
deriving DecidableEq, Repr

def count8 : CountRow := { bases := 3003, trivial := 2316, scanned := 687, skipped := 0, passed := 3003, failed := 0 }
def count9 : CountRow := { bases := 3432, trivial := 862, scanned := 2570, skipped := 0, passed := 3432, failed := 0 }
def count10 : CountRow := { bases := 3003, trivial := 216, scanned := 2787, skipped := 0, passed := 3003, failed := 0 }
def count11 : CountRow := { bases := 2002, trivial := 137, scanned := 1865, skipped := 0, passed := 2002, failed := 0 }
def count12 : CountRow := { bases := 1001, trivial := 204, scanned := 797, skipped := 0, passed := 1001, failed := 0 }
def count13 : CountRow := { bases := 364, trivial := 260, scanned := 104, skipped := 0, passed := 364, failed := 0 }

def countRow : Fin 6 -> CountRow
  | ⟨0, _⟩ => count8
  | ⟨1, _⟩ => count9
  | ⟨2, _⟩ => count10
  | ⟨3, _⟩ => count11
  | ⟨4, _⟩ => count12
  | _ => count13

def totalBases : Nat :=
  count8.bases + count9.bases + count10.bases + count11.bases + count12.bases + count13.bases

def totalTrivial : Nat :=
  count8.trivial + count9.trivial + count10.trivial + count11.trivial + count12.trivial + count13.trivial

def totalScanned : Nat :=
  count8.scanned + count9.scanned + count10.scanned + count11.scanned + count12.scanned + count13.scanned

def totalSkipped : Nat :=
  count8.skipped + count9.skipped + count10.skipped + count11.skipped + count12.skipped + count13.skipped

def totalPassed : Nat :=
  count8.passed + count9.passed + count10.passed + count11.passed + count12.passed + count13.passed

def totalFailed : Nat :=
  count8.failed + count9.failed + count10.failed + count11.failed + count12.failed + count13.failed

theorem count_totals :
    totalBases = 12805 ∧ totalTrivial = 3995 ∧ totalScanned = 8810 ∧
      totalSkipped = 0 ∧ totalPassed = 12805 ∧ totalFailed = 0 := by
  native_decide

theorem every_count_row_passes :
    ∀ r : Fin 6, (countRow r).passed = (countRow r).bases ∧ (countRow r).failed = 0 ∧
      (countRow r).skipped = 0 := by
  native_decide

/-- The k=8 dilated AP worst row uses period-max `2`, not `1`.  This records
the normalization guardrail found by the S77 brute P=840 recheck. -/
theorem k8_periodmax_is_two : k8.pmNum = 2 ∧ k8.pmDen = 1 := by
  native_decide

/-! ### Axiom audit -/

#print axioms k8_headroom_positive
#print axioms k9_headroom_positive
#print axioms k10_headroom_positive
#print axioms k11_headroom_positive
#print axioms k12_headroom_positive
#print axioms k13_headroom_positive
#print axioms all_worst_rows_headroom_positive
#print axioms k9_is_global_worst_among_worst_rows
#print axioms count_totals
#print axioms every_count_row_passes
#print axioms k8_periodmax_is_two

end PeriodmaxCertificate
end LonelyRunner
