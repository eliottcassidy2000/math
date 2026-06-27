/-
  TournamentH7.LRCObserverGluingLedger -- Lean-facing observer-gluing
  interface for the current LRC14 frontier.

  This file packages the S258/HYP-3095/HYP-3096/HYP-3097 lesson without
  asserting the still-open coverage theorem.  The new object sits between the
  concrete packet scouts and `LRCFiniteAddressBranchClosure`:

    observer charts + overlap discharge + scissors/direct-arc fields
      -> terminal discharge certificate
      -> Mreach >= 1/14

  The key guardrail is explicit: coarse mod-14 winding data is recorded as a
  non-terminal status unless a fine-scale sidecar or independent terminal
  discharge is attached.
-/

import Mathlib.Tactic
import TournamentH7.LRCFiniteAddressBranchClosure

namespace LonelyRunner
namespace LRC14

/-! ## Observer charts and coarse-tournament guardrails -/

/-- Proof observers currently used by the LRC14 gluing ledger.

The vertices are proof charts, not runners, residues, arcs, or coarse
tournament vertices. -/
inductive ObserverChart where
  | arithmetic_lift
  | direct_arc
  | normalized_arc
  | pair_scissors
  | cap_pascal
  | moment_perron
  | branch_k33
  | fine_scale_tournament
  | formal_witness
  | raw_scalar_shadow
  | coarse_mod14_winding
deriving DecidableEq, Repr, Fintype

/-- The observer-gluing frontier currently names eleven charts. -/
theorem observerChart_count : Fintype.card ObserverChart = 11 := by
  native_decide

/-- Status of the coarse mod-14 winding tournament suggested by the incoming
coverage-vs-H comparison. -/
inductive CoarseWindingStatus where
  | proper_tournament
  | degenerate_antipodal_tie
  | fine_scale_required
deriving DecidableEq, Repr, Fintype

/-- A coarse winding chart may be terminal only if it is a genuine tournament
or has already been replaced by fine-scale data. -/
def CoarseWindingStatus.chartIsTerminal :
    CoarseWindingStatus -> Prop
  | CoarseWindingStatus.proper_tournament => True
  | CoarseWindingStatus.degenerate_antipodal_tie => False
  | CoarseWindingStatus.fine_scale_required => True

/-- The antipodal-tie status is a warning, not a terminal proof chart. -/
theorem coarseWinding_degenerate_not_terminal :
    ¬ CoarseWindingStatus.chartIsTerminal
      CoarseWindingStatus.degenerate_antipodal_tie := by
  intro h
  exact h

/-! ## Legal overlap discharges -/

/-- Why an overlap between charts is legal after one chart forgets a
coordinate. -/
inductive GluingDischarge where
  | reconstructs
  | dual_annihilates
  | fiber_constant
  | descends
  | terminal_exit
  | named_debt
deriving DecidableEq, Repr, Fintype

/-- One certified chart overlap in the observer sheaf.  The boolean is kept as
a data field because the current scripts emit audit flags, but the proof field
forces it to be true before the overlap is usable. -/
structure ChartOverlapCertificate where
  left : ObserverChart
  right : ObserverChart
  destroyed : DestroyedCoordinate
  discharge : GluingDischarge
  valid : Bool
  valid_true : valid = true

/-- Certified overlaps expose their audit flag. -/
theorem chartOverlapCertificate_valid
    (cert : ChartOverlapCertificate) :
    cert.valid = true :=
  cert.valid_true

/-! ## Direct-arc denominator-net numerics -/

/-- Exact rational largest-arc numerics sufficient to claim a denominator net:
`arcNum / arcDen` is the observed direct lonely arc and `thresholdD` is a
strict reciprocal ceiling witness, encoded arithmetically as
`arcDen < thresholdD * arcNum`. -/
structure DenominatorNetNumerics where
  arcNum : Nat
  arcDen : Nat
  thresholdD : Nat
  arcNum_pos : 0 < arcNum
  arcDen_pos : 0 < arcDen
  threshold_spec : arcDen < thresholdD * arcNum

/-- Any strict reciprocal threshold has positive denominator bound. -/
theorem DenominatorNetNumerics.thresholdD_pos
    (num : DenominatorNetNumerics) :
    0 < num.thresholdD := by
  by_contra h
  have hzero : num.thresholdD = 0 := Nat.eq_zero_of_not_pos h
  have hlt : num.arcDen < num.thresholdD * num.arcNum := num.threshold_spec
  rw [hzero] at hlt
  simp at hlt

/-- S258 H7=6 boundary sample: largest direct arc `19/1372`, grid threshold
`D=73`.  This is evidence data, not a global theorem. -/
def s258H7Eq6BoundaryArc : DenominatorNetNumerics where
  arcNum := 19
  arcDen := 1372
  thresholdD := 73
  arcNum_pos := by decide
  arcDen_pos := by decide
  threshold_spec := by native_decide

/-- S258 divisor-loaded `B=8` sample: largest direct arc `1/82320`, grid
threshold `D=82321`. -/
def s258DivisorLoadedB8Arc : DenominatorNetNumerics where
  arcNum := 1
  arcDen := 82320
  thresholdD := 82321
  arcNum_pos := by decide
  arcDen_pos := by decide
  threshold_spec := by native_decide

/-! ## Observer-gluing obligations and certificates -/

/-- A proof-facing ledger row before terminal discharge is known.  It records
the S258 direct-arc fields, HYP-3097 pair-scissors sidecars, CRT/Farey lanes,
and chart-overlap audits that must be preserved by a later producer theorem. -/
structure ObserverGluingObligation where
  sourceRowId : Nat
  chart : ObserverChart
  preservedLRPredicate : Bool
  preservedLRPredicate_true : preservedLRPredicate = true
  containsMultipleOf14 : Bool
  multiple7Count : Nat
  multiple14Count : Nat
  directComponentCount : Nat
  denominatorNet : Option DenominatorNetNumerics
  h7PairShadowNum : Nat
  h7PairShadowDen : Nat
  h7PairShadowDen_pos : 0 < h7PairShadowDen
  evenPairShadowNum : Nat
  evenPairShadowDen : Nat
  evenPairShadowDen_pos : 0 < evenPairShadowDen
  mod7ScissorsSignature : Fin 7 -> Nat
  mod14ScissorsSignature : Fin 14 -> Nat
  fareyAdditiveLaneMod91 : Option Nat
  fareyProductLaneMod91 : Option Nat
  coarseWindingStatus : CoarseWindingStatus
  overlaps : List ChartOverlapCertificate
  remainingDebt : List BranchTerminalExit

/-- The obligation row carries the explicit statement that the LR predicate was
preserved through the observer chart. -/
theorem ObserverGluingObligation.preservesLR
    (obligation : ObserverGluingObligation) :
    obligation.preservedLRPredicate = true :=
  obligation.preservedLRPredicate_true

/-- A terminal observer-gluing certificate is an obligation plus the same
proof-bearing discharge used by the finite-address branch closure. -/
structure ObserverGluingCertificate (v : Fin 13 -> Int)
    extends ObserverGluingObligation where
  terminal : TerminalDischargeCertificate v

/-- Observer-gluing certificates are sound because their terminal discharge is
proof-bearing. -/
theorem observerGluingCertificate_mreach
    {v : Fin 13 -> Int}
    (cert : ObserverGluingCertificate v) :
    (1 : Real) / 14 <= Mreach v :=
  terminalDischarge_mreach cert.terminal

/-- Coverage by early gates or observer-gluing certificates.  This is the
current producer target after S258, not an asserted theorem. -/
def ObserverGluingCoverage : Prop :=
  forall v : Fin 13 -> Int, (forall i, v i ≠ 0) ->
    Nonempty (EarlyGateCertificate v) ∨ Nonempty (ObserverGluingCertificate v)

/-- Observer-gluing coverage implies the top-level LRC14 statement. -/
theorem lrc14_from_observer_gluing_coverage
    (hcoverage : ObserverGluingCoverage) :
    LRC14Statement := by
  intro v hv
  rcases hcoverage v hv with hEarly | hGlue
  · rcases hEarly with ⟨early⟩
    exact lonely_of_Mreach_ge v hv (earlyGateCertificate_mreach early)
  · rcases hGlue with ⟨cert⟩
    exact lonely_of_Mreach_ge v hv (observerGluingCertificate_mreach cert)

/-! ## Compatibility with the existing finite-address branch interface -/

/-- A finite-address branch packet can be viewed as a minimal observer-gluing
certificate.  This keeps the new S258 interface conservative over S254. -/
def ObserverGluingCertificate.ofFiniteAddressBranchPacket
    {v : Fin 13 -> Int}
    (packet : FiniteAddressBranchPacket v) :
    ObserverGluingCertificate v where
  sourceRowId := packet.sourceRowId
  chart := ObserverChart.formal_witness
  preservedLRPredicate := true
  preservedLRPredicate_true := rfl
  containsMultipleOf14 := packet.containsMultipleOf14
  multiple7Count := packet.multiple7Count
  multiple14Count := packet.multiple14Count
  directComponentCount := 0
  denominatorNet := none
  h7PairShadowNum := 0
  h7PairShadowDen := 91
  h7PairShadowDen_pos := by decide
  evenPairShadowNum := 0
  evenPairShadowDen := 91
  evenPairShadowDen_pos := by decide
  mod7ScissorsSignature := fun _ => 0
  mod14ScissorsSignature := fun _ => 0
  fareyAdditiveLaneMod91 := none
  fareyProductLaneMod91 := none
  coarseWindingStatus := CoarseWindingStatus.fine_scale_required
  overlaps := []
  remainingDebt := []
  terminal := packet.terminal

/-- The existing cutting-edge branch coverage theorem is stronger than the new
observer-gluing coverage target when branch packets are already available. -/
theorem observerGluingCoverage_of_cutting_edge
    (hcoverage : CuttingEdgeBranchCoverage) :
    ObserverGluingCoverage := by
  intro v hv
  rcases hcoverage v hv with hEarly | hPacket
  · exact Or.inl hEarly
  · rcases hPacket with ⟨packet⟩
    exact Or.inr
      ⟨ObserverGluingCertificate.ofFiniteAddressBranchPacket packet⟩

/-! ## Concrete S258 obligation row examples -/

/-- The S258 H7=6 boundary residual as a nonterminal obligation row.  It is
kept as data so later producer theorems can target the exact same fields. -/
def s258H7Eq6BoundaryObligation : ObserverGluingObligation where
  sourceRowId := 258
  chart := ObserverChart.direct_arc
  preservedLRPredicate := true
  preservedLRPredicate_true := rfl
  containsMultipleOf14 := true
  multiple7Count := 6
  multiple14Count := 3
  directComponentCount := 42
  denominatorNet := some s258H7Eq6BoundaryArc
  h7PairShadowNum := 15
  h7PairShadowDen := 91
  h7PairShadowDen_pos := by decide
  evenPairShadowNum := 15
  evenPairShadowDen := 91
  evenPairShadowDen_pos := by decide
  mod7ScissorsSignature := fun i =>
    if i.val = 0 then 6 else if i.val = 6 then 2 else 1
  mod14ScissorsSignature := fun i =>
    if i.val = 0 then 3
    else if i.val = 7 then 3
    else if i.val = 13 then 1
    else if i.val < 7 then 1
    else 0
  fareyAdditiveLaneMod91 := none
  fareyProductLaneMod91 := none
  coarseWindingStatus := CoarseWindingStatus.fine_scale_required
  overlaps := []
  remainingDebt :=
    [BranchTerminalExit.covering_moment,
      BranchTerminalExit.k33_state_lift,
      BranchTerminalExit.named_residual_debt]

/-- The S258 divisor-loaded `B=8` residual as a nonterminal obligation row. -/
def s258DivisorLoadedB8Obligation : ObserverGluingObligation where
  sourceRowId := 258
  chart := ObserverChart.normalized_arc
  preservedLRPredicate := true
  preservedLRPredicate_true := rfl
  containsMultipleOf14 := true
  multiple7Count := 2
  multiple14Count := 1
  directComponentCount := 860
  denominatorNet := some s258DivisorLoadedB8Arc
  h7PairShadowNum := 1
  h7PairShadowDen := 91
  h7PairShadowDen_pos := by decide
  evenPairShadowNum := 15
  evenPairShadowDen := 91
  evenPairShadowDen_pos := by decide
  mod7ScissorsSignature := fun i =>
    if i.val = 5 then 1 else 2
  mod14ScissorsSignature := fun i =>
    if i.val = 12 then 0 else 1
  fareyAdditiveLaneMod91 := some 14
  fareyProductLaneMod91 := some 35
  coarseWindingStatus := CoarseWindingStatus.fine_scale_required
  overlaps := []
  remainingDebt :=
    [BranchTerminalExit.covering_moment,
      BranchTerminalExit.k33_state_lift,
      BranchTerminalExit.named_residual_debt]

/-! ## Axiom audit hooks -/

#print axioms observerChart_count
#print axioms coarseWinding_degenerate_not_terminal
#print axioms chartOverlapCertificate_valid
#print axioms DenominatorNetNumerics.thresholdD_pos
#print axioms observerGluingCertificate_mreach
#print axioms lrc14_from_observer_gluing_coverage
#print axioms observerGluingCoverage_of_cutting_edge

end LRC14
end LonelyRunner
