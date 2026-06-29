/-
  TournamentH7.LRCColoredGateFormalization -- HYP-3473 Lean-facing formal
  interface for the HYP-3471 colored gate-reservoir theorem target.

  This module does not prove the finite HYP-3471 implication from geometry.
  It makes the exact formal shape of that implication explicit:

    dead component
      -> rank <= 2 survivor gate
      -> one endpoint is even (`E`) and the other is a branch wall (`B0/B1`)
      -> typed residue, branch mask, adjacency, and cover deltas are retained
      -> terminal AP/random/conductance/owner-current/two-adic/SPEC route

  The open mathematical content is a `Prop` producer obligation
  (`DeadCoverEBranchSoundness`).  The Lean contribution is the sorry-free
  packet interface, exact HYP-3471 count ledger, and conditional assembly from
  colored gate terminal packets to `LRC14Statement`.
-/

import Mathlib.Tactic
import TournamentH7.LRCFourteenSkeleton

namespace LonelyRunner
namespace LRC14

open scoped BigOperators

/-! ## Endpoint and gate words -/

/-- Endpoint colors retained by the HYP-3471 survivor-gate quotient. -/
inductive GateEndpointKind where
  | E
  | B0
  | B1
deriving DecidableEq, Repr, Fintype

/-- Branch mask of a survivor gate after the two-adic pullback. -/
inductive GateBranchMask where
  | branch0
  | branch1
  | both
deriving DecidableEq, Repr, Fintype

/-- Which adjacent bad-core block touches the survivor gate. -/
inductive GateAdjacency where
  | left_bad_edge
  | right_bad_edge
  | two_sided
deriving DecidableEq, Repr, Fintype

/-- A typed endpoint label: the wall kind plus residue modulo 14. -/
structure TypedEndpointResidue where
  kind : GateEndpointKind
  residue : Fin 14
deriving DecidableEq, Repr

/-- A survivor gate with the HYP-3471 proof payload retained. -/
structure ColoredSurvivorGate where
  leftKind : GateEndpointKind
  rightKind : GateEndpointKind
  leftResidues : List TypedEndpointResidue
  rightResidues : List TypedEndpointResidue
  branchMask : GateBranchMask
  adjacency : GateAdjacency
  b0Delta : Nat
  b1Delta : Nat
  endpointRank : Nat
deriving Repr

/-- Branch endpoints are the two odd-lift wall colors. -/
def IsBranchEndpoint : GateEndpointKind -> Prop
  | GateEndpointKind.B0 => True
  | GateEndpointKind.B1 => True
  | GateEndpointKind.E => False

/-- The HYP-3471 E/branch condition: exactly one endpoint is even and the
other is a branch wall. -/
def OneEOneBranch (gate : ColoredSurvivorGate) : Prop :=
  (gate.leftKind = GateEndpointKind.E ∧ IsBranchEndpoint gate.rightKind) ∨
  (gate.rightKind = GateEndpointKind.E ∧ IsBranchEndpoint gate.leftKind)

/-- Same-branch gates are retained as gluing sidecars, not as the basic
dead-cover exit. -/
def SameBranchGate (gate : ColoredSurvivorGate) : Prop :=
  (gate.leftKind = GateEndpointKind.B0 ∧ gate.rightKind = GateEndpointKind.B0) ∨
  (gate.leftKind = GateEndpointKind.B1 ∧ gate.rightKind = GateEndpointKind.B1)

/-- Cross-branch gates are useful diagnostics but are not E/branch gates. -/
def CrossBranchGate (gate : ColoredSurvivorGate) : Prop :=
  (gate.leftKind = GateEndpointKind.B0 ∧ gate.rightKind = GateEndpointKind.B1) ∨
  (gate.leftKind = GateEndpointKind.B1 ∧ gate.rightKind = GateEndpointKind.B0)

/-- The low-rank gate condition used by HYP-3453/HYP-3471. -/
def LowRankGate (gate : ColoredSurvivorGate) : Prop :=
  gate.endpointRank ≤ 2

/-- The sharpened gate target: low endpoint rank and an E/branch boundary. -/
def EBranchLowRankGate (gate : ColoredSurvivorGate) : Prop :=
  LowRankGate gate ∧ OneEOneBranch gate

theorem oneEOneBranch_not_sameBranch
    (gate : ColoredSurvivorGate) :
    OneEOneBranch gate -> ¬ SameBranchGate gate := by
  rcases gate with ⟨leftKind, rightKind, leftResidues, rightResidues,
    branchMask, adjacency, b0Delta, b1Delta, endpointRank⟩
  cases leftKind <;> cases rightKind <;>
    simp [OneEOneBranch, SameBranchGate, IsBranchEndpoint]

theorem oneEOneBranch_not_crossBranch
    (gate : ColoredSurvivorGate) :
    OneEOneBranch gate -> ¬ CrossBranchGate gate := by
  rcases gate with ⟨leftKind, rightKind, leftResidues, rightResidues,
    branchMask, adjacency, b0Delta, b1Delta, endpointRank⟩
  cases leftKind <;> cases rightKind <;>
    simp [OneEOneBranch, CrossBranchGate, IsBranchEndpoint]

theorem eBranchLowRank_rank_le_two
    {gate : ColoredSurvivorGate}
    (hgate : EBranchLowRankGate gate) :
    gate.endpointRank ≤ 2 :=
  hgate.1

/-! ## Formal theorem target from HYP-3471 -/

/-- Row-level finite theorem target: every dead-cover row has a low-rank
E/branch survivor gate.  The row type is abstract so Python/exact-rational
audits, hand proofs, or future Lean finite enumerators can instantiate it. -/
def DeadCoverEBranchSoundness
    (Row : Type)
    (gates : Row -> List ColoredSurvivorGate)
    (hasDeadComponent : Row -> Prop) : Prop :=
  forall row : Row, hasDeadComponent row ->
    exists gate : ColoredSurvivorGate,
      gate ∈ gates row ∧ EBranchLowRankGate gate

theorem eBranch_gate_of_dead_cover
    {Row : Type}
    {gates : Row -> List ColoredSurvivorGate}
    {hasDeadComponent : Row -> Prop}
    (sound : DeadCoverEBranchSoundness Row gates hasDeadComponent)
    {row : Row}
    (hdead : hasDeadComponent row) :
    exists gate : ColoredSurvivorGate,
      gate ∈ gates row ∧ EBranchLowRankGate gate :=
  sound row hdead

/-- Exact count ledger reported by HYP-3471.  These are audit readouts, not a
proof that the bank is universal. -/
structure ColoredGateBankCounts where
  rowsAudited : Nat
  rowsWithDeadComponents : Nat
  lowRankGates : Nat
  eBranchLowRankGates : Nat
  sameBranchLowRankGates : Nat
  crossBranchLowRankGates : Nat
  deadRowsWithEBranchLowRankGate : Nat
  deadRowsWithoutEBranchLowRankGate : Nat
  deadRowsOnlyEBranchLowRankGates : Nat
  deadRowsWithSameBranchLowRankGate : Nat
  deadRowsWithCrossBranchLowRankGate : Nat

/-- HYP-3471 aggregate counts from
`lrc14_colored_gate_reservoir_codex_20260629.py`. -/
def hyp3471BankCounts : ColoredGateBankCounts where
  rowsAudited := 135
  rowsWithDeadComponents := 130
  lowRankGates := 8666
  eBranchLowRankGates := 7002
  sameBranchLowRankGates := 1482
  crossBranchLowRankGates := 182
  deadRowsWithEBranchLowRankGate := 130
  deadRowsWithoutEBranchLowRankGate := 0
  deadRowsOnlyEBranchLowRankGates := 28
  deadRowsWithSameBranchLowRankGate := 102
  deadRowsWithCrossBranchLowRankGate := 54

theorem hyp3471_dead_rows_all_have_e_branch_gate_count :
    hyp3471BankCounts.deadRowsWithEBranchLowRankGate =
      hyp3471BankCounts.rowsWithDeadComponents ∧
    hyp3471BankCounts.deadRowsWithoutEBranchLowRankGate = 0 := by
  decide

theorem hyp3471_e_branch_gate_count_positive :
    0 < hyp3471BankCounts.eBranchLowRankGates := by
  decide

/-! ## Legal colored-gate packets and exits -/

/-- Terminal routes currently legal for a colored gate packet. -/
inductive ColoredGateTerminalExit where
  | ap84_corridor_splice
  | ap84_color_grid_placement
  | random031_gate_gluing
  | component_conductance_menger
  | owner_current_imbalance
  | two_adic_descent
  | signed_spec_rprime
  | named_residual_debt
deriving DecidableEq, Repr, Fintype

/-- A quotient is legal only if it records the sidecars it keeps. -/
structure ColoredGatePayload where
  keepsTypedResidue : Bool
  keepsEndpointKind : Bool
  keepsBranchMask : Bool
  keepsAdjacency : Bool
  keepsCoverDeltas : Bool
  keepsMirrorData : Bool
  keepsDischargeExit : Bool

/-- The full HYP-3471 colored gate word keeps exactly the proof payload needed
before any scalar/color compression is allowed. -/
def fullColoredGatePayload : ColoredGatePayload where
  keepsTypedResidue := true
  keepsEndpointKind := true
  keepsBranchMask := true
  keepsAdjacency := true
  keepsCoverDeltas := true
  keepsMirrorData := true
  keepsDischargeExit := true

theorem fullColoredGatePayload_keeps_boundary_data :
    fullColoredGatePayload.keepsTypedResidue = true ∧
    fullColoredGatePayload.keepsBranchMask = true ∧
    fullColoredGatePayload.keepsCoverDeltas = true ∧
    fullColoredGatePayload.keepsDischargeExit = true := by
  decide

/-- A terminal colored-gate packet for a concrete LRC14 row.  The `mreachFloor`
field is the final mathematical discharge supplied by the selected route. -/
structure ColoredGateTerminalPacket (v : Fin 13 -> Int) where
  exit : ColoredGateTerminalExit
  payload : ColoredGatePayload
  mreachFloor : (1 : Real) / 14 ≤ Mreach v

theorem mreach_of_coloredGateTerminalPacket
    {v : Fin 13 -> Int}
    (packet : ColoredGateTerminalPacket v) :
    (1 : Real) / 14 ≤ Mreach v :=
  packet.mreachFloor

/-- Global colored-gate coverage: every LRC14 speed family receives one of the
legal terminal colored-gate packets. -/
def ColoredGateGlobalCoverage : Prop :=
  forall v : Fin 13 -> Int, (forall i, v i ≠ 0) ->
    Nonempty (ColoredGateTerminalPacket v)

/-- Conditional top-level assembly.  This is not a proof of LRC14; it states
the exact remaining producer obligation in the colored-gate language. -/
theorem lrc14_from_colored_gate_global_coverage
    (hcoverage : ColoredGateGlobalCoverage) :
    LRC14Statement := by
  intro v hv
  rcases hcoverage v hv with ⟨packet⟩
  exact lonely_of_Mreach_ge v hv packet.mreachFloor

/-! ## Tournament-analysis carrier ledger -/

/-- Proof carriers from the HYP-3471 tournament. -/
inductive ColoredGateCarrier where
  | dead_positive_e_branch_gate
  | full_colored_gate_word
  | structural_color_sidecar
  | typed_mod14_gate_word
  | ap84_four_color_packet
  | endpoint_kind_coloring
  | numeric_14_residue
  | raw_color_count
deriving DecidableEq, Repr, Fintype

/-- Total scores from the HYP-3471 carrier tournament. -/
def coloredGateCarrierScore : ColoredGateCarrier -> Nat
  | ColoredGateCarrier.dead_positive_e_branch_gate => 66
  | ColoredGateCarrier.full_colored_gate_word => 66
  | ColoredGateCarrier.structural_color_sidecar => 61
  | ColoredGateCarrier.typed_mod14_gate_word => 58
  | ColoredGateCarrier.ap84_four_color_packet => 54
  | ColoredGateCarrier.endpoint_kind_coloring => 42
  | ColoredGateCarrier.numeric_14_residue => 31
  | ColoredGateCarrier.raw_color_count => 7

theorem coloredGateCarrier_count :
    Fintype.card ColoredGateCarrier = 8 := by
  native_decide

theorem dead_positive_e_branch_gate_score_maximal :
    forall carrier : ColoredGateCarrier,
      coloredGateCarrierScore carrier ≤
        coloredGateCarrierScore ColoredGateCarrier.dead_positive_e_branch_gate := by
  native_decide

theorem raw_color_count_score_minimal :
    forall carrier : ColoredGateCarrier,
      coloredGateCarrierScore ColoredGateCarrier.raw_color_count ≤
        coloredGateCarrierScore carrier := by
  native_decide

theorem full_word_ties_dead_positive_score :
    coloredGateCarrierScore ColoredGateCarrier.full_colored_gate_word =
      coloredGateCarrierScore ColoredGateCarrier.dead_positive_e_branch_gate := by
  decide

/-! ## Axiom audit hooks -/

#print axioms oneEOneBranch_not_sameBranch
#print axioms oneEOneBranch_not_crossBranch
#print axioms eBranch_gate_of_dead_cover
#print axioms hyp3471_dead_rows_all_have_e_branch_gate_count
#print axioms fullColoredGatePayload_keeps_boundary_data
#print axioms lrc14_from_colored_gate_global_coverage
#print axioms coloredGateCarrier_count
#print axioms dead_positive_e_branch_gate_score_maximal
#print axioms raw_color_count_score_minimal

end LRC14
end LonelyRunner
