import Mathlib

/-!
# The seven-wall strict-spectrum Hunter floor (THM-1221)

The analytic pair-overlap formula and its product tails, and the identification
of every speed packet with one of the finite projective cut branches, are
external providers.  The dependency-free referee for THM-1221 verifies that
classification exactly.  This module kernel-checks the arithmetic constants,
the exhaustive surviving-branch consumer, and the Hunter, common-period, and
BAD-event erosion consequences.

There are no proof placeholders and no native evaluator.
-/

namespace LRC14.SevenWallStrictSpectrum

/-! ## Exact thresholds and branch values -/

def highBar : ℚ := 1 / 63
def firstStrictHigh : ℚ := 5 / 308
def alignedFloor : ℚ := 265 / 2772
def treeFloor : ℚ := 15 / 154
def heavyCircuitCutoff : ℚ := 60 / 637

theorem first_aligned_ledger :
    firstStrictHigh + 5 * highBar = alignedFloor := by
  norm_num [firstStrictHigh, highBar, alignedFloor]

theorem strict_connected_ledger :
    6 * firstStrictHigh = treeFloor := by
  norm_num [firstStrictHigh, treeFloor]

theorem strict_two_component_ledgers :
    5 / 56 + highBar = 53 / 504 ∧
    1 / 21 + 4 * firstStrictHigh + highBar = 89 / 693 := by
  norm_num [firstStrictHigh, highBar]

theorem strong_improves_aligned :
    treeFloor - alignedFloor = 5 / 2772 ∧ alignedFloor < treeFloor := by
  norm_num [treeFloor, alignedFloor]

/-- The four globally surviving threshold-component branches.  The missing
three-component strict-high branch is impossible: the exact closed-low common
neighbor bank holds at most `2+2+2` vertices. -/
inductive SurvivingBranch
  | strictHighConnected
  | strictHighTwoComponentsOneSix
  | strictHighTwoComponentsTwoFive
  | weakHighDisconnected
  deriving DecidableEq, Repr

/-- Analytic lower values after the finite exact ratio banks have supplied
their small facts: the first strict spectral value, the `1/56` internal-Q
edge, the four large common-neighbor center ratios, and the singleton channel
multiplicities. -/
def branchFloor : SurvivingBranch → ℚ
  | .strictHighConnected => 15 / 154
  | .strictHighTwoComponentsOneSix => 53 / 504
  | .strictHighTwoComponentsTwoFive => 89 / 693
  | .weakHighDisconnected => 1284 / 13013

theorem every_branch_clears_treeFloor (branch : SurvivingBranch) :
    treeFloor ≤ branchFloor branch := by
  cases branch <;> norm_num [treeFloor, branchFloor]

/-- Honest interface to the analytic-plus-finite provider.  In applications,
`Packet` is the type of seven distinct positive speed packets.  The provider
must classify every packet into one of the four proved branches and supply its
tree-weight comparison with the exact branch minimum. -/
structure SpectrumProvider (Packet : Type*) where
  branch : Packet → SurvivingBranch
  treeWeight : Packet → ℚ
  branch_lower : ∀ packet, branchFloor (branch packet) ≤ treeWeight packet

/-- Once the exact channel/cut classification has supplied a provider, the
universal `15/154` tree floor is a kernel-checked finite consumer. -/
theorem provider_tree_floor {Packet : Type*} (provider : SpectrumProvider Packet)
    (packet : Packet) : treeFloor ≤ provider.treeWeight packet := by
  exact le_trans (every_branch_clears_treeFloor (provider.branch packet))
    (provider.branch_lower packet)

/-! ## The transparent singleton-cut subcertificate -/

def nonACredit : ℚ := 223 / 13013

def alignedSingletonLedger : Fin 5 → ℚ
  | 0 => 1284 / 13013
  | 1 => 115601 / 1171170
  | 2 => 16303 / 167310
  | 3 => 37547 / 390390
  | 4 => 103 / 1050

def strongSingletonLedger : Fin 5 → ℚ
  | 0 => 1284 / 13013
  | 1 => 115601 / 1171170
  | 2 => 28411 / 278850
  | 3 => 615257 / 5855850
  | 4 => 633883 / 5855850

theorem nonA_credit_exact :
    51 / 1183 - 2 / 77 = nonACredit ∧ highBar < nonACredit := by
  norm_num [nonACredit, highBar]

theorem singleton_ledger_clears_aligned (m : Fin 5) :
    alignedFloor < alignedSingletonLedger m := by
  fin_cases m <;> norm_num [alignedFloor, alignedSingletonLedger]

theorem strong_singleton_ledger_clears_treeFloor (m : Fin 5) :
    treeFloor < strongSingletonLedger m := by
  fin_cases m <;> norm_num [treeFloor, strongSingletonLedger]

theorem strong_singleton_uniform_minimum (m : Fin 5) :
    (1284 : ℚ) / 13013 ≤ strongSingletonLedger m := by
  fin_cases m <;> norm_num [strongSingletonLedger]

theorem strong_singleton_margin :
    (1284 : ℚ) / 13013 - treeFloor = 3 / 2366 := by
  norm_num [treeFloor]

theorem four_A_internal_credit :
    1 / 70 + 32 / 1575 + 4 / 63 = (103 : ℚ) / 1050 := by
  norm_num

/-! ## Hunter and periodic consumers -/

/-- Arithmetic endpoint of the Hunter forest inequality.  The analytic/set
provider supplies `tree ≤ safe`; the exact tree theorem is in cleared form. -/
theorem hunter_safe_floor (tree safe : ℝ)
    (htree : 15 ≤ 154 * tree) (hhunter : tree ≤ safe) :
    15 ≤ 154 * safe := by
  nlinarith

/-- A common period with safe mass `15/154` leaves covered mass at most
`139/154`.  Composing with `L ≥ 1/(7m)` gives `g/m ≤ 139/22`, with all
divisions cleared. -/
theorem common_dilate_protected_needle (m g L : ℝ)
    (hm : 0 ≤ m) (hg : 0 ≤ g)
    (hneedle : 1 ≤ 7 * m * L)
    (hcovered : 154 * g * L ≤ 139) :
    22 * g ≤ 139 * m := by
  have hscaled : 154 * g ≤ 973 * m := by
    calc
      154 * g = (154 * g) * 1 := by ring
      _ ≤ (154 * g) * (7 * m * L) :=
        mul_le_mul_of_nonneg_left hneedle (mul_nonneg (by norm_num) hg)
      _ = (7 * m) * (154 * g * L) := by ring
      _ ≤ (7 * m) * 139 :=
        mul_le_mul_of_nonneg_left hcovered (mul_nonneg (by norm_num) hm)
      _ = 973 * m := by ring
  nlinarith

theorem covered_period_constant :
    1 - treeFloor = 139 / 154 ∧ 7 * (1 - treeFloor) = 139 / 22 := by
  norm_num [treeFloor]

/-! ## THM-1203/1218 incidence and erosion budgets -/

theorem exact_bad_event_margins :
    treeFloor - 2 / 21 = 1 / 462 ∧
    treeFloor - heavyCircuitCutoff = 45 / 14014 := by
  norm_num [treeFloor, heavyCircuitCutoff]

theorem historical_aligned_margins :
    alignedFloor - 2 / 21 = 1 / 2772 ∧
    alignedFloor - heavyCircuitCutoff = 355 / 252252 := by
  norm_num [alignedFloor, heavyCircuitCutoff]

/-- If the safe set has mass at least `15/154` and an arbitrary deletion BAD
event has mass at most `2/21`, their set-difference ledger retains `1/462`.
The set-measure inequality itself is the explicit upstream provider. -/
theorem arbitrary_bad_incidence_budget (safe bad : ℝ)
    (hsafe : 15 ≤ 154 * safe) (hbad : 21 * bad ≤ 2) :
    1 ≤ 462 * (safe - bad) := by
  nlinarith

/-- In the non-AP branch THM-1218 bounds BAD mass by `60/637`, leaving the
larger exact budget `45/14014`. -/
theorem nonAP_bad_incidence_budget (safe bad : ℝ)
    (hsafe : 15 ≤ 154 * safe) (hbad : 637 * bad ≤ 60) :
    45 ≤ 14014 * (safe - bad) := by
  nlinarith

/-- A localization/discretization erosion strictly below `1/462` preserves a
positive witness even after an arbitrary BAD event is removed. -/
theorem arbitrary_bad_after_erosion (safe bad erosion remaining : ℝ)
    (hsafe : 15 ≤ 154 * safe) (hbad : 21 * bad ≤ 2)
    (herosion : 462 * erosion < 1)
    (hremaining : safe - bad - erosion ≤ remaining) :
    0 < remaining := by
  have hbudget := arbitrary_bad_incidence_budget safe bad hsafe hbad
  nlinarith

/-- The corresponding non-AP erosion budget is `45/14014`. -/
theorem nonAP_bad_after_erosion (safe bad erosion remaining : ℝ)
    (hsafe : 15 ≤ 154 * safe) (hbad : 637 * bad ≤ 60)
    (herosion : 14014 * erosion < 45)
    (hremaining : safe - bad - erosion ≤ remaining) :
    0 < remaining := by
  have hbudget := nonAP_bad_incidence_budget safe bad hsafe hbad
  nlinarith

/-! ## Axiom audit -/

#print axioms first_aligned_ledger
#print axioms strict_connected_ledger
#print axioms strict_two_component_ledgers
#print axioms every_branch_clears_treeFloor
#print axioms provider_tree_floor
#print axioms singleton_ledger_clears_aligned
#print axioms strong_singleton_ledger_clears_treeFloor
#print axioms strong_singleton_uniform_minimum
#print axioms hunter_safe_floor
#print axioms common_dilate_protected_needle
#print axioms exact_bad_event_margins
#print axioms arbitrary_bad_incidence_budget
#print axioms nonAP_bad_incidence_budget
#print axioms arbitrary_bad_after_erosion
#print axioms nonAP_bad_after_erosion

end LRC14.SevenWallStrictSpectrum
