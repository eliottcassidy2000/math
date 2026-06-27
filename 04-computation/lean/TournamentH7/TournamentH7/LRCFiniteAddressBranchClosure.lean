/-
  TournamentH7.LRCFiniteAddressBranchClosure -- S254 Lean-facing frontier for
  the current LRC14 finite-address branch-closure proof shape.

  This file intentionally does not assert the still-open analytic normalizer,
  covering-family discharge, or K33 state-lift construction.  It packages the
  proof interface now suggested by HYP-3083:

    low-apex/top-balanced residual row
      -> finite-address packet
      -> q-cusp/Hurwitz sidecar
      -> protected branch graph
      -> terminal discharge certificate
      -> concrete Mreach >= 1/14

  The contribution is a sorry-free conditional assembly theorem and a concrete
  record shape tying together the existing q-cusp, median-center, moment-dual,
  and tournament-state-lift modules.
-/

import TournamentH7.LRCMedianCenterControl
import TournamentH7.LRCModularCuspLedger
import TournamentH7.LRCMomentDual
import TournamentH7.LRCTournamentStateLift

namespace LonelyRunner
namespace LRC14

open scoped BigOperators

/-! ## Proof-route labels retained by the finite address -/

/-- Terminal exits currently allowed by the HYP-3083 branch-closure spine. -/
inductive BranchTerminalExit where
  | direct_q_witness
  | apex_majority_gamma
  | level7_lift_sieve
  | one_large_speed_peeler
  | covering_moment
  | c27_petal
  | ap_gw_boundary
  | k33_state_lift
  | thm572_state_lift_contradiction
  | named_residual_debt
deriving DecidableEq, Repr, Fintype

/-- Early gates that discharge rows before the low-apex/top-balanced residual. -/
inductive EarlyGateExit where
  | q_witness
  | apex_majority_gamma
  | level7_lift_sieve
  | one_large_speed_peeler
deriving DecidableEq, Repr, Fintype

/-- Which coordinate a quotient claims to forget, and therefore must protect. -/
inductive DestroyedCoordinate where
  | none
  | raw_runner_identity
  | exact_scale
  | endpoint_owner
  | route_label
  | q_cusp_tail
  | hurwitz_seed
  | branch_bridge
  | formal_status
deriving DecidableEq, Repr, Fintype

/-- Bridge status before and after attaching sidecars. -/
inductive BridgeStatus where
  | no_bridge
  | protected_bridge
  | raw_naked_bridge
  | residual_debt
deriving DecidableEq, Repr, Fintype

/-- Hurwitz/Markov/Pell address type retained by an arithmetic sidecar. -/
inductive HurwitzAddressKind where
  | none
  | continued_fraction_period
  | markov_depth
  | pell_wall_unit
  | vieta_mutation_word
  | zero_persistence_gate
deriving DecidableEq, Repr, Fintype

/-- Vertices of the S254 proof-carrier tournament.  These are proof
obligations and retained sidecars, not runners or q-coefficients. -/
inductive FiniteAddressBranchVertex where
  | global_packet_normalizer
  | protected_branch_graph
  | covering_moment_exit
  | k33_state_lift_exit
  | q_cusp_principal_part_guard
  | hurwitz_seed_guard
  | median_center_scheduler
  | level7_lift_sieve_gate
  | one_large_speed_gate
  | q_witness_gate
  | lean_sidecar_closure
  | raw_scalar_shadow
deriving DecidableEq, Repr, Fintype

/-- The S254 carrier tournament has the twelve vertices named in HYP-3083. -/
theorem finiteAddressBranchVertex_count :
    Fintype.card FiniteAddressBranchVertex = 12 := by
  native_decide

/-! ## q-cusp and q-Pochhammer ledgers -/

/-- A q-Pochhammer/product-tail certificate.  The product tail may be infinite
later, but its Lean-side entry has a no-negative-powers proof. -/
structure QPochhammerTailCertificate where
  coeff : Int -> Int
  noNegativePowers :
    ModularCuspLedger.HasNoNegativePowers coeff
  finiteNegativeTail :
    ModularCuspLedger.HasOnlyFiniteNegativePowers coeff

/-- No-negative-powers data is enough to build a product-tail certificate. -/
def qPochhammerTailCertificateOfNoNegative
    (coeff : Int -> Int)
    (hno : ModularCuspLedger.HasNoNegativePowers coeff) :
    QPochhammerTailCertificate where
  coeff := coeff
  noNegativePowers := hno
  finiteNegativeTail :=
    ModularCuspLedger.finite_negative_powers_of_no_negative_powers hno

/-- The S247/S254 Euler q-Pochhammer tail through q^12 is holomorphic at the
cusp and therefore legal as a nonpolar tail. -/
def qPochhammerTailS254 : QPochhammerTailCertificate :=
  qPochhammerTailCertificateOfNoNegative
    ModularCuspLedger.qPochhammerEulerCoeffTo12
    ModularCuspLedger.qPochhammerEulerCoeffTo12_no_negative_powers

/-- A q-cusp sidecar: finite principal part plus product tail and a proof that
no infinite polar tail is being smuggled into the packet. -/
structure QCuspSidecar where
  q_cusp_ledger_id : Nat
  principalPart : ModularCuspLedger.LaurentPrincipalPartPacket
  productTail : QPochhammerTailCertificate
  polarExitWord : List BranchTerminalExit
  hurwitzZeroPersistenceStatus : Bool
  illegalInfinitePolarTailFlag : Bool
  noIllegalInfinitePolarTail : illegalInfinitePolarTailFlag = false
  finitePrincipalPart :
    ModularCuspLedger.HasOnlyFiniteNegativePowers principalPart.coeff

/-- A q-cusp sidecar exposes exactly the finite-principal-part certificate. -/
theorem qCuspSidecar_finite_principal_part
    (sidecar : QCuspSidecar) :
    ModularCuspLedger.HasOnlyFiniteNegativePowers sidecar.principalPart.coeff :=
  sidecar.finitePrincipalPart

/-- The S247/S254 `j` principal part as a packet. -/
def jPrincipalPartPacketS254 :
    ModularCuspLedger.LaurentPrincipalPartPacket where
  coeff := ModularCuspLedger.jPrincipalCoeffS247
  poleOrder := 1
  coeff_vanishes_below_pole := by
    intro n hn
    have hlt : n < (-1 : Int) := by omega
    simp [ModularCuspLedger.jPrincipalCoeffS247, hlt]

/-- A concrete q-cusp sidecar for the `j`-principal-part lane. -/
def jCuspSidecarS254 : QCuspSidecar where
  q_cusp_ledger_id := 254
  principalPart := jPrincipalPartPacketS254
  productTail := qPochhammerTailS254
  polarExitWord := [BranchTerminalExit.direct_q_witness,
    BranchTerminalExit.covering_moment]
  hurwitzZeroPersistenceStatus := true
  illegalInfinitePolarTailFlag := false
  noIllegalInfinitePolarTail := rfl
  finitePrincipalPart :=
    ModularCuspLedger.packet_has_only_finite_negative_powers
      jPrincipalPartPacketS254

/-! ## Hurwitz arithmetic and protected-branch sidecars -/

/-- Hurwitz/Markov/Pell arithmetic address retained by the packet. -/
structure HurwitzArithmeticSidecar where
  kind : HurwitzAddressKind
  lagrangeMarkovDepth : Nat
  finiteSeed : List Int
  mutationWord : List Nat
  pellCassiniGap : Int
  endpointShellAddress : Nat
  destroyedCoordinate : DestroyedCoordinate
  requiredSidecarOrExit : BranchTerminalExit

/-- The empty arithmetic sidecar is still an explicit sidecar: it says no
Hurwitz scalar is being used silently. -/
def noHurwitzArithmeticSidecar : HurwitzArithmeticSidecar where
  kind := HurwitzAddressKind.none
  lagrangeMarkovDepth := 0
  finiteSeed := []
  mutationWord := []
  pellCassiniGap := 0
  endpointShellAddress := 0
  destroyedCoordinate := DestroyedCoordinate.none
  requiredSidecarOrExit := BranchTerminalExit.direct_q_witness

/-- Protected branch-graph status after attaching sidecars. -/
structure ProtectedBranchCertificate where
  protectedBranchNode : Nat
  rawBridgeStatus : BridgeStatus
  protectedBridgeStatus : BridgeStatus
  bridgeProtectionMode : List DestroyedCoordinate
  reverseVerificationPath : Bool
  noNakedProtectedBridge :
    protectedBridgeStatus ≠ BridgeStatus.raw_naked_bridge

/-- The protected sidecar blocks the exact failure mode found in the raw
scalar-star quotient. -/
theorem protectedBranch_no_raw_naked_bridge
    (cert : ProtectedBranchCertificate) :
    cert.protectedBridgeStatus ≠ BridgeStatus.raw_naked_bridge :=
  cert.noNakedProtectedBridge

/-- A minimal protected branch certificate. -/
def trivialProtectedBranchCertificate : ProtectedBranchCertificate where
  protectedBranchNode := 0
  rawBridgeStatus := BridgeStatus.no_bridge
  protectedBridgeStatus := BridgeStatus.protected_bridge
  bridgeProtectionMode := []
  reverseVerificationPath := true
  noNakedProtectedBridge := by
    intro h
    cases h

/-! ## Covering-moment dual ledger -/

/-- The moment-dual sidecar used by covering-family exits.  It records the dual
polynomial and the exact feasibility condition used by `LRCMomentDual`. -/
structure CoveringMomentDualLedger where
  E : List Int
  g : Nat -> Real
  feasible :
    forall n : Nat, n <= 6 -> (if n = 0 then (1 : Real) else 0) <= g n

/-- A covering-moment ledger gives the already-formalized p0 <= Ly reduction. -/
theorem coveringMomentDualLedger_p0_le_Ly
    (ledger : CoveringMomentDualLedger) :
    (DenseCovers.slowμ
        (DenseCovers.coverSet ledger.E)).toReal <=
      TournamentH7.MomentDual.Ly ledger.E ledger.g :=
  TournamentH7.MomentDual.p0_le_Ly ledger.E ledger.g ledger.feasible

/-! ## Terminal certificates and packet coverage -/

/-- A terminal discharge is proof-bearing: it carries a witness floor and the
handoff from that floor to concrete `Mreach`. -/
structure TerminalDischargeCertificate (v : Fin 13 -> Int) where
  exit : BranchTerminalExit
  witnessFloor : Real
  threshold_le_floor : (1 : Real) / 14 <= witnessFloor
  floor_le_mreach : witnessFloor <= Mreach v

/-- Terminal discharge soundness is just transitivity of the carried floor. -/
theorem terminalDischarge_mreach
    {v : Fin 13 -> Int}
    (cert : TerminalDischargeCertificate v) :
    (1 : Real) / 14 <= Mreach v :=
  le_trans cert.threshold_le_floor cert.floor_le_mreach

/-- Package an already-known reach bound as a terminal certificate.  This is a
conservativity check for the interface, not a new LRC14 proof. -/
noncomputable def terminalDischargeOfMreach
    (v : Fin 13 -> Int)
    (exit : BranchTerminalExit)
    (hM : (1 : Real) / 14 <= Mreach v) :
    TerminalDischargeCertificate v where
  exit := exit
  witnessFloor := (1 : Real) / 14
  threshold_le_floor := le_rfl
  floor_le_mreach := hM

/-- An early-gate certificate for rows discharged before the residual packet
normalizer is needed. -/
structure EarlyGateCertificate (v : Fin 13 -> Int) where
  exit : EarlyGateExit
  terminal : TerminalDischargeCertificate v

/-- Early-gate certificates also feed the concrete `Mreach` bridge. -/
theorem earlyGateCertificate_mreach
    {v : Fin 13 -> Int}
    (cert : EarlyGateCertificate v) :
    (1 : Real) / 14 <= Mreach v :=
  terminalDischarge_mreach cert.terminal

/-- The current low-apex/top-balanced finite-address packet target. -/
structure FiniteAddressBranchPacket (v : Fin 13 -> Int) where
  sourceRowId : Nat
  containsMultipleOf14 : Bool
  containsMultipleOf14_true : containsMultipleOf14 = true
  multiple14Count : Nat
  multiple14Count_live : 1 <= multiple14Count ∧ multiple14Count <= 6
  multiple7Count : Nat
  multiple7Count_live : 1 <= multiple7Count ∧ multiple7Count <= 6
  lowApexTopBalancedStatus : Bool
  lowApexTopBalancedStatus_true : lowApexTopBalancedStatus = true
  finiteAddressWord : List Nat
  exactScaleAddress : Nat
  destroyedCoordinate : DestroyedCoordinate
  routeCenter : CenterControlPacket v
  qCusp : QCuspSidecar
  hurwitzArithmetic : HurwitzArithmeticSidecar
  branch : ProtectedBranchCertificate
  coveringMoment : Option CoveringMomentDualLedger
  terminal : TerminalDischargeCertificate v

/-- The finite-address branch packet is sound because its terminal discharge is
proof-bearing. -/
theorem finiteAddressBranchPacket_mreach
    {v : Fin 13 -> Int}
    (packet : FiniteAddressBranchPacket v) :
    (1 : Real) / 14 <= Mreach v :=
  terminalDischarge_mreach packet.terminal

/-- Coverage by finite-address branch packets for every nonzero 13-speed
family.  This is the Lean name for the global packet theorem still missing in
the mathematics. -/
def FiniteAddressBranchCoverage : Prop :=
  forall v : Fin 13 -> Int, (forall i, v i ≠ 0) ->
    Nonempty (FiniteAddressBranchPacket v)

/-- The cutting-edge coverage shape after q-witness, level-7 lift sieve, and
one-large-speed gates: every row is either discharged early or emits the
low-apex/top-balanced finite-address packet.  THM-573 moves the live residual
from "few multiples of 14" to "at most six multiples of 7"; the packet keeps
both counts because the covering reduction still begins with a multiple of 14. -/
def CuttingEdgeBranchCoverage : Prop :=
  forall v : Fin 13 -> Int, (forall i, v i ≠ 0) ->
    Nonempty (EarlyGateCertificate v) ∨ Nonempty (FiniteAddressBranchPacket v)

/-- Finite-address branch coverage implies the top-level LRC14 statement. -/
theorem lrc14_from_finite_address_branch_coverage
    (hcoverage : FiniteAddressBranchCoverage) :
    LRC14Statement := by
  intro v hv
  rcases hcoverage v hv with ⟨packet⟩
  exact lonely_of_Mreach_ge v hv (finiteAddressBranchPacket_mreach packet)

/-- The sharper S254 frontier: early gates plus low-apex/top-balanced
finite-address packets imply LRC14. -/
theorem lrc14_from_cutting_edge_branch_coverage
    (hcoverage : CuttingEdgeBranchCoverage) :
    LRC14Statement := by
  intro v hv
  rcases hcoverage v hv with hEarly | hPacket
  · rcases hEarly with ⟨early⟩
    exact lonely_of_Mreach_ge v hv (earlyGateCertificate_mreach early)
  · rcases hPacket with ⟨packet⟩
    exact lonely_of_Mreach_ge v hv (finiteAddressBranchPacket_mreach packet)

/-- Conservative packaging check: any already-known `Mreach` lower bound can be
expressed as a finite-address branch packet. -/
noncomputable def finiteAddressBranchPacketOfMreach
    (v : Fin 13 -> Int)
    (hM : (1 : Real) / 14 <= Mreach v) :
    FiniteAddressBranchPacket v where
  sourceRowId := 254
  containsMultipleOf14 := true
  containsMultipleOf14_true := rfl
  multiple14Count := 1
  multiple14Count_live := by constructor <;> omega
  multiple7Count := 1
  multiple7Count_live := by constructor <;> omega
  lowApexTopBalancedStatus := true
  lowApexTopBalancedStatus_true := rfl
  finiteAddressWord := [254]
  exactScaleAddress := 0
  destroyedCoordinate := DestroyedCoordinate.none
  routeCenter := centerControlPacketOfMreach v hM
  qCusp := jCuspSidecarS254
  hurwitzArithmetic := noHurwitzArithmeticSidecar
  branch := trivialProtectedBranchCertificate
  coveringMoment := none
  terminal :=
    terminalDischargeOfMreach v BranchTerminalExit.direct_q_witness hM

/-- If a global reach-floor theorem is ever proved independently, it gives
finite-address coverage by the conservative packet wrapper. -/
theorem finiteAddressBranchCoverage_of_reach_floor
    (hreach : forall v : Fin 13 -> Int, (forall i, v i ≠ 0) ->
      (1 : Real) / 14 <= Mreach v) :
    FiniteAddressBranchCoverage := by
  intro v hv
  exact ⟨finiteAddressBranchPacketOfMreach v (hreach v hv)⟩

/-! ## Axiom audit hooks -/

#print axioms finiteAddressBranchVertex_count
#print axioms qCuspSidecar_finite_principal_part
#print axioms protectedBranch_no_raw_naked_bridge
#print axioms coveringMomentDualLedger_p0_le_Ly
#print axioms terminalDischarge_mreach
#print axioms finiteAddressBranchPacket_mreach
#print axioms lrc14_from_finite_address_branch_coverage
#print axioms lrc14_from_cutting_edge_branch_coverage

end LRC14
end LonelyRunner
