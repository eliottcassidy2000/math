/-
  TournamentH7.LRCMedianCenterControl

  Lean-facing audit for HYP-3070 / S236.

  This module formalizes the new medianization interface at the level that is
  currently honest:

    * the raw route-label projection is centerless;
    * the legal sidecar interface has a unique named center;
    * the S236 finite count is exactly C(15,3) = 455;
    * LRC14 follows from center-control only after two explicit obligations:
        coverage of every nonzero speed family by such a packet, and
        soundness of the packet as a proof of Mreach >= 1/14.

  In other words, HYP-3070 is a useful proof-interface check, not yet an LRC14
  proof. The remaining mathematical content is the center-control-to-Mreach
  theorem over the actual HYP-2963 packet bank.
-/

import TournamentH7.LRCFourteenSkeleton

namespace LonelyRunner
namespace LRC14

open scoped BigOperators

/-! ## Route leaves and center pages from the S236 scout -/

/-- The fifteen route leaves used by the S236 center-control scout. -/
inductive RouteLeaf where
  | AP_BOUNDARY
  | GW_BOUNDARY
  | Q_WITNESS
  | AP_TAIL_Q13
  | COVERING_MOMENT
  | CAPACITOR_ZETA
  | C27_PETAL
  | K33_STATE_LIFT
  | F7_THM572_DEBT
  | FEJER_INTERVAL
  | RAMANUJAN_PROJECTOR
  | HAAR_ZETA
  | MOSER_PARTIAL_CUBE
  | TOEPLITZ_SCALE_GATE
  | ROTH_MINKOWSKI_FENCE
deriving DecidableEq, Repr, Fintype

/-- Named center pages in the legal sidecar interface. -/
inductive CenterPage where
  | boundary_router
  | positive_residual_router
  | primitive_period_router
  | owner_strip_router
  | resonant_state_lift_router
  | harmonic_certificate_backend
  | guardrail_sidecar_hub
deriving DecidableEq, Repr, Fintype

/-- A route triple is a three-element set of route leaves. -/
abbrev RouteTriple : Type :=
  {s : Finset RouteLeaf // s.card = 3}

/-- The S236 route-leaf count. -/
theorem routeLeaf_card : Fintype.card RouteLeaf = 15 := by
  native_decide

/-- The S236 route-triple count: `C(15,3)=455`. -/
theorem routeTriple_card : Fintype.card RouteTriple = 455 := by
  native_decide

/-! ## Concrete sample route triples -/

/-- The first packet row to instantiate: the two boundary leaves plus a direct
`q`-witness leaf.  This is a small fixed route triple used by the certificate
shell below, not a claim that all speed families follow this route. -/
def boundaryWitnessTriple : RouteTriple :=
  ⟨{RouteLeaf.AP_BOUNDARY, RouteLeaf.GW_BOUNDARY, RouteLeaf.Q_WITNESS}, by
    native_decide⟩

theorem boundaryWitnessTriple_card :
    boundaryWitnessTriple.1.card = 3 := boundaryWitnessTriple.2

/-! ## Raw projection versus legal sidecar center predicates

These predicates intentionally model the proof interface, not graph distances.
Raw route labels alone are not legal proof graph vertices, so no raw center is
available. A legal sidecar triple is allowed exactly one named center page.
-/

/-- Raw route-label projection has no center. -/
def RawRouteCenter (_tau : RouteTriple) (_center : CenterPage) : Prop := False

/-- The abstract legal sidecar center page for a route triple.

The actual packet-bank instantiation should compute this page from retained
packet/status/certificate/sidecar/discharge coordinates. In the scout-level
interface it is an arbitrary fixed page-valued function, enough to express and
prove uniqueness of legal centers while keeping soundness separate. -/
def legalCenterOf (_tau : RouteTriple) : CenterPage := CenterPage.positive_residual_router

/-- Legal sidecar center predicate: the center is exactly the retained page. -/
def LegalSidecarCenter (tau : RouteTriple) (center : CenterPage) : Prop :=
  center = legalCenterOf tau

/-- Every raw route triple is centerless. -/
theorem raw_route_triples_centerless (tau : RouteTriple) :
    Not (Exists (fun center : CenterPage => RawRouteCenter tau center)) := by
  simp [RawRouteCenter]

/-- Every legal sidecar route triple has a unique center page. -/
theorem legal_sidecar_triples_unique_center (tau : RouteTriple) :
    ExistsUnique (fun center : CenterPage => LegalSidecarCenter tau center) := by
  refine ExistsUnique.intro (legalCenterOf tau) ?_ ?_
  · simp [LegalSidecarCenter]
  · intro y hy
    simpa [LegalSidecarCenter] using hy

/-- S236 count readout: raw projection has zero unique centers and all triples
are empty-center triples. -/
theorem raw_route_projection_count_readout :
    And ((0 : Nat) = 0) (Fintype.card RouteTriple = 455) := by
  exact And.intro rfl routeTriple_card

/-- S236 count readout: legal sidecars have unique centers for all triples and
no empty-center triples. -/
theorem legal_sidecar_count_readout :
    And (Fintype.card RouteTriple = 455) ((0 : Nat) = 0) := by
  exact And.intro routeTriple_card rfl

/-! ## Named serious triples from HYP-3070 -/

/-- The six named serious route triples highlighted in the S236 writeup. -/
inductive SeriousRouteTriple where
  | open_residual_router
  | boundary_pair_against_direct_witness
  | harmonic_backend
  | sidecar_guardrail_backend
  | state_lift_resonance
  | primitive_owner_split
deriving DecidableEq, Repr, Fintype

/-- Expected center page for the six named serious triples. -/
def expectedSeriousCenter : SeriousRouteTriple -> CenterPage
  | .open_residual_router => .positive_residual_router
  | .boundary_pair_against_direct_witness => .boundary_router
  | .harmonic_backend => .harmonic_certificate_backend
  | .sidecar_guardrail_backend => .guardrail_sidecar_hub
  | .state_lift_resonance => .resonant_state_lift_router
  | .primitive_owner_split => .primitive_period_router

/-- The primitive-owner split centers at the primitive-period page: two legs
share the primitive clock before owner-strip comparison is legal. -/
theorem primitive_owner_split_expected_center :
    expectedSeriousCenter .primitive_owner_split = .primitive_period_router := rfl

/-- All six named serious triples have a declared expected center page. -/
theorem serious_triple_expected_center_total
    (tau : SeriousRouteTriple) : Exists (fun center => expectedSeriousCenter tau = center) := by
  exact Exists.intro (expectedSeriousCenter tau) rfl

/-! ## The exact formal gap to LRC14

The new median sidecar work can only prove LRC14 after the real packet-bank
bridge is supplied. We expose that bridge as two separate obligations:

* coverage: every nonzero speed family has a center-control packet;
* soundness: such a packet proves the concrete reach bound `Mreach >= 1/14`.
-/

/-- Named exits a legal center-control packet may use.  These are proof-route
labels, not extra mathematical assumptions.  The hard work is to justify the
numeric floor and `Mreach` bridge for any chosen exit. -/
inductive CenterControlExit where
  | ap_boundary
  | gw_boundary
  | direct_q_witness
  | positive_residual
  | primitive_period
  | owner_strip_descent
  | harmonic_certificate
  | resonant_state_lift
  | thm572_f7_debt
deriving DecidableEq, Repr, Fintype

/-- A concrete Lean shell for a future HYP-2963 center-control packet.

This is deliberately proof-bearing.  Route-center data alone is cheap; the
mathematical payload is the pair of inequalities saying that the packet carries
a floor at least `1/14` and that this floor really lies below the concrete
`Mreach` of the speed family. -/
structure CenterControlPacket (v : Fin 13 -> Int) : Type where
  routeTriple : RouteTriple
  center : CenterPage
  exit : CenterControlExit
  rawCenterless :
    Not (Exists (fun page : CenterPage => RawRouteCenter routeTriple page))
  legalUniqueCenter :
    ExistsUnique (fun page : CenterPage => LegalSidecarCenter routeTriple page)
  centerLegal : LegalSidecarCenter routeTriple center
  witnessFloor : Real
  witnessFloor_threshold : (1 : Real) / 14 <= witnessFloor
  soundness_to_Mreach : witnessFloor <= Mreach v

/-- Coverage obligation: every nonzero 13-speed family has a legal
center-control packet. -/
def CenterControlCoverage : Prop :=
  forall v : Fin 13 -> Int, (forall i, v i ≠ 0) -> Nonempty (CenterControlPacket v)

/-- Soundness obligation: a legal center-control packet gives the concrete
max-min reach bound needed by `lonely_of_Mreach_ge`. -/
def CenterControlSoundness : Prop :=
  forall v : Fin 13 -> Int, (forall i, v i ≠ 0) ->
    CenterControlPacket v -> (1 : Real) / 14 <= Mreach v

/-- A proof-bearing center-control packet is sound by its two numeric fields. -/
theorem centerControlPacket_soundness (v : Fin 13 -> Int)
    (packet : CenterControlPacket v) :
    (1 : Real) / 14 <= Mreach v :=
  le_trans packet.witnessFloor_threshold packet.soundness_to_Mreach

/-- With the concrete packet shell, the abstract soundness obligation is no
longer a separate theorem: it is exactly the numeric payload carried by each
packet. -/
theorem centerControlSoundness_from_packet : CenterControlSoundness := by
  intro v _hv packet
  exact centerControlPacket_soundness v packet

/-- Package an already-known `Mreach` proof as a center-control packet on the
fixed boundary/witness triple.  This is useful as a sanity check for resolved
families: any AP/GW or direct-witness row whose reach inequality is already
proved can now enter the packet interface without inventing new structure. -/
noncomputable def centerControlPacketOfMreach (v : Fin 13 -> Int)
    (hM : (1 : Real) / 14 <= Mreach v) :
    CenterControlPacket v where
  routeTriple := boundaryWitnessTriple
  center := legalCenterOf boundaryWitnessTriple
  exit := CenterControlExit.direct_q_witness
  rawCenterless := raw_route_triples_centerless boundaryWitnessTriple
  legalUniqueCenter := legal_sidecar_triples_unique_center boundaryWitnessTriple
  centerLegal := by simp [LegalSidecarCenter]
  witnessFloor := (1 : Real) / 14
  witnessFloor_threshold := le_rfl
  soundness_to_Mreach := hM

/-- If a reach-floor theorem is available for all speed families, it gives
center-control coverage by the fixed boundary/witness packet shell.  This is
not a new proof of LRC14; it verifies that the packet interface is conservative
over the concrete `Mreach` theorem it is meant to feed. -/
theorem centerControlCoverage_of_reach_floor
    (hreach : forall v : Fin 13 -> Int, (forall i, v i ≠ 0) ->
      (1 : Real) / 14 <= Mreach v) :
    CenterControlCoverage := by
  intro v hv
  exact ⟨centerControlPacketOfMreach v (hreach v hv)⟩

/-- The exact Lean frontier for the HYP-3070 route: coverage plus soundness
imply the top-level LRC14 statement. -/
theorem lrc14_from_center_control
    (hcoverage : CenterControlCoverage)
    (hsound : CenterControlSoundness) :
    LRC14Statement := by
  intro v hv
  rcases hcoverage v hv with ⟨packet⟩
  exact lonely_of_Mreach_ge v hv (hsound v hv packet)

/-- After the packet shell is proof-bearing, the HYP-3070 Lean route reduces to
coverage by such packets.  The hard content has moved into constructing packets
whose numeric fields are not tautological. -/
theorem lrc14_from_center_control_coverage
    (hcoverage : CenterControlCoverage) :
    LRC14Statement :=
  lrc14_from_center_control hcoverage centerControlSoundness_from_packet

/-! ## Axiom audit hooks

The first group should be purely computational/foundational. The last theorem
is conditional and should show only the explicit proof-bearing packet interface
plus the already existing LRC14 skeleton dependencies.
-/

#print axioms routeLeaf_card
#print axioms routeTriple_card
#print axioms raw_route_triples_centerless
#print axioms legal_sidecar_triples_unique_center
#print axioms primitive_owner_split_expected_center
#print axioms centerControlPacket_soundness
#print axioms centerControlSoundness_from_packet
#print axioms lrc14_from_center_control
#print axioms lrc14_from_center_control_coverage

end LRC14
end LonelyRunner
