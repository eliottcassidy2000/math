/-
  TournamentH7.LRCSingletonCurrentLedger -- HYP-3480 finite ledger for the
  zero-edge singleton-current packet.

  This module does not prove the geometric singleton-current lemma.  It records
  the exact HYP-3480 dispatch arithmetic as a sorry-free Lean ledger:

    six small-touch/no-hard zero-edge rows
      -> mirror-unit singleton-current packet
      -> two cover-delta-minimum rows kept as unit-certificate shadows
      -> random031 routed to the HYP-3481/HYP-3482 hard-control seam packet

  The mathematical producer still needed is the mirror-unit singleton-current
  lemma over swapped singleton B0/B1 owner pairs with route sidecar R.
-/

import Mathlib.Tactic
import TournamentH7.LRCHardOrbitCurrentJoin

namespace LonelyRunner
namespace LRC14

/-! ## HYP-3480 finite rows -/

/-- The seven rows audited by HYP-3480: six target singleton-current rows plus
the random031 hard/currentless control. -/
inductive SingletonCurrentAuditRow where
  | random001
  | random039
  | random062
  | random074
  | random086
  | random101
  | random031
deriving DecidableEq, Repr, Fintype

/-- Rows in the six-row small-touch/no-hard target packet. -/
def hyp3480TargetRows : List SingletonCurrentAuditRow :=
  [SingletonCurrentAuditRow.random001,
   SingletonCurrentAuditRow.random039,
   SingletonCurrentAuditRow.random062,
   SingletonCurrentAuditRow.random074,
   SingletonCurrentAuditRow.random086,
   SingletonCurrentAuditRow.random101]

/-- The named hard/currentless control row, now routed through HYP-3481/HYP-3482. -/
def hyp3480ControlRow : SingletonCurrentAuditRow :=
  SingletonCurrentAuditRow.random031

theorem hyp3480_target_row_count :
    hyp3480TargetRows.length = 6 := by
  decide

theorem hyp3480_control_not_target :
    hyp3480ControlRow ∉ hyp3480TargetRows := by
  simp [hyp3480ControlRow, hyp3480TargetRows]

theorem hyp3480_audited_row_count :
    Fintype.card SingletonCurrentAuditRow = 7 := by
  native_decide

/-! ## Terminal packet classes -/

/-- Terminal class assigned by the HYP-3480 audit. -/
inductive SingletonCurrentTerminalClass where
  | mirror_unit_singleton_packet
  | mirror_unit_singleton_packet_cover_delta_min_shadow
  | random031_hard_currentless_control
deriving DecidableEq, Repr, Fintype

/-- HYP-3480 terminal class by row.  Rows `039` and `074` have cover-delta
absolute minima, but still carry component-level branch-unit certificates. -/
def singletonCurrentTerminalClass :
    SingletonCurrentAuditRow -> SingletonCurrentTerminalClass
  | SingletonCurrentAuditRow.random001 =>
      SingletonCurrentTerminalClass.mirror_unit_singleton_packet
  | SingletonCurrentAuditRow.random039 =>
      SingletonCurrentTerminalClass.mirror_unit_singleton_packet_cover_delta_min_shadow
  | SingletonCurrentAuditRow.random062 =>
      SingletonCurrentTerminalClass.mirror_unit_singleton_packet
  | SingletonCurrentAuditRow.random074 =>
      SingletonCurrentTerminalClass.mirror_unit_singleton_packet_cover_delta_min_shadow
  | SingletonCurrentAuditRow.random086 =>
      SingletonCurrentTerminalClass.mirror_unit_singleton_packet
  | SingletonCurrentAuditRow.random101 =>
      SingletonCurrentTerminalClass.mirror_unit_singleton_packet
  | SingletonCurrentAuditRow.random031 =>
      SingletonCurrentTerminalClass.random031_hard_currentless_control

theorem hyp3480_cover_delta_shadow_rows_are_unit_certificate_targets :
    singletonCurrentTerminalClass SingletonCurrentAuditRow.random039 =
        SingletonCurrentTerminalClass.mirror_unit_singleton_packet_cover_delta_min_shadow ∧
      singletonCurrentTerminalClass SingletonCurrentAuditRow.random074 =
        SingletonCurrentTerminalClass.mirror_unit_singleton_packet_cover_delta_min_shadow := by
  decide

theorem hyp3480_random031_is_control :
    singletonCurrentTerminalClass hyp3480ControlRow =
      SingletonCurrentTerminalClass.random031_hard_currentless_control := by
  decide

/-! ## Aggregate count ledger -/

/-- Exact count ledger from `lrc14_zero_edge_singleton_current_codex_20260629.py`. -/
structure SingletonCurrentAuditCounts where
  auditedRows : Nat
  targetRows : Nat
  controlRows : Nat
  targetProjectionZeroRows : Nat
  targetDeadComponents : Nat
  targetComponentsWithCompleteBranchUnitTouch : Nat
  targetMirrorPairs : Nat
  targetMirrorPairsWithBranchUnitGate : Nat
  coverDeltaMinShadowRowsWithUnitCertificate : Nat
  controlDeadComponents : Nat
  controlComponentsWithAnyCompleteTouch : Nat
  controlComponentsWithCompleteBranchUnitTouch : Nat
  controlMirrorPairs : Nat
  controlMirrorPairsWithBranchUnitGate : Nat
  hardCurrentlessControlRows : Nat

/-- HYP-3480 aggregate counts. -/
def hyp3480Counts : SingletonCurrentAuditCounts where
  auditedRows := 7
  targetRows := 6
  controlRows := 1
  targetProjectionZeroRows := 6
  targetDeadComponents := 14
  targetComponentsWithCompleteBranchUnitTouch := 14
  targetMirrorPairs := 7
  targetMirrorPairsWithBranchUnitGate := 7
  coverDeltaMinShadowRowsWithUnitCertificate := 2
  controlDeadComponents := 4
  controlComponentsWithAnyCompleteTouch := 2
  controlComponentsWithCompleteBranchUnitTouch := 0
  controlMirrorPairs := 2
  controlMirrorPairsWithBranchUnitGate := 0
  hardCurrentlessControlRows := 1

theorem hyp3480_audited_row_partition :
    hyp3480Counts.targetRows + hyp3480Counts.controlRows =
      hyp3480Counts.auditedRows := by
  decide

theorem hyp3480_all_targets_are_zero_edge :
    hyp3480Counts.targetProjectionZeroRows = hyp3480Counts.targetRows := by
  decide

theorem hyp3480_all_target_components_have_branch_unit_touch :
    hyp3480Counts.targetComponentsWithCompleteBranchUnitTouch =
      hyp3480Counts.targetDeadComponents := by
  decide

theorem hyp3480_all_target_mirror_pairs_have_unit_gate :
    hyp3480Counts.targetMirrorPairsWithBranchUnitGate =
      hyp3480Counts.targetMirrorPairs := by
  decide

theorem hyp3480_cover_delta_shadows_still_have_unit_certificates :
    hyp3480Counts.coverDeltaMinShadowRowsWithUnitCertificate = 2 ∧
      hyp3480Counts.coverDeltaMinShadowRowsWithUnitCertificate ≤
        hyp3480Counts.targetRows := by
  decide

theorem hyp3480_random031_has_no_complete_branch_unit_control_touch :
    hyp3480Counts.controlComponentsWithCompleteBranchUnitTouch = 0 ∧
      hyp3480Counts.controlMirrorPairsWithBranchUnitGate = 0 := by
  decide

/-! ## Dispatch packet -/

/-- HYP-3480 terminal dispatch count split. -/
structure SingletonCurrentDispatchPacket where
  mirrorUnitRows : Nat
  coverDeltaShadowRows : Nat
  random031ControlRows : Nat
  totalRows : Nat
  complete :
    mirrorUnitRows + coverDeltaShadowRows + random031ControlRows = totalRows

/-- Four direct mirror-unit rows, two cover-delta-minimum shadow rows with unit
certificates, and one random031 control row. -/
def hyp3480DispatchPacket : SingletonCurrentDispatchPacket where
  mirrorUnitRows := 4
  coverDeltaShadowRows := 2
  random031ControlRows := 1
  totalRows := 7
  complete := by decide

theorem hyp3480_dispatch_complete :
    hyp3480DispatchPacket.mirrorUnitRows +
        hyp3480DispatchPacket.coverDeltaShadowRows +
        hyp3480DispatchPacket.random031ControlRows =
      hyp3480DispatchPacket.totalRows :=
  hyp3480DispatchPacket.complete

theorem hyp3480_dispatch_matches_counts :
    hyp3480DispatchPacket.totalRows = hyp3480Counts.auditedRows ∧
      hyp3480DispatchPacket.coverDeltaShadowRows =
        hyp3480Counts.coverDeltaMinShadowRowsWithUnitCertificate ∧
      hyp3480DispatchPacket.random031ControlRows =
        hyp3480Counts.hardCurrentlessControlRows := by
  decide

/-! ## Formal producer obligation -/

/-- Data that the still-open singleton-current lemma must consume. -/
structure MirrorUnitSingletonCurrentObligation where
  swappedSingletonB0B1OwnerPair : Prop
  mirrorCompatibleCompleteBranchUnitGatePair : Prop
  routeSidecarR : Prop
  ownerResidueTwoAdicWord : Prop

/-- The proof target left by HYP-3480.  The terminal predicate is abstract so
future geometry, owner-current, two-adic, or SPEC proofs can instantiate it. -/
def MirrorUnitSingletonCurrentLemmaTarget
    (Terminal : MirrorUnitSingletonCurrentObligation -> Prop) : Prop :=
  forall packet : MirrorUnitSingletonCurrentObligation,
    packet.swappedSingletonB0B1OwnerPair ->
    packet.mirrorCompatibleCompleteBranchUnitGatePair ->
    packet.routeSidecarR ->
    packet.ownerResidueTwoAdicWord ->
    Terminal packet

/-! ## Tournament carrier ledger -/

inductive SingletonCurrentCarrier where
  | mirror_unit_singleton_current_packet
  | component_complete_touch_certificate
  | random031_hard_control_clause
  | route_sidecar_R_join
  | cover_delta_min_gate_shadow
  | owner_residue_two_adic_shadow
  | raw_zero_edge_projection_count
  | raw_row_name_list
deriving DecidableEq, Repr, Fintype

/-- HYP-3480 proof-carrier scores from the stored tournament audit. -/
def singletonCurrentCarrierScore : SingletonCurrentCarrier -> Nat
  | SingletonCurrentCarrier.mirror_unit_singleton_current_packet => 67
  | SingletonCurrentCarrier.component_complete_touch_certificate => 62
  | SingletonCurrentCarrier.random031_hard_control_clause => 59
  | SingletonCurrentCarrier.route_sidecar_R_join => 59
  | SingletonCurrentCarrier.cover_delta_min_gate_shadow => 52
  | SingletonCurrentCarrier.owner_residue_two_adic_shadow => 45
  | SingletonCurrentCarrier.raw_zero_edge_projection_count => 8
  | SingletonCurrentCarrier.raw_row_name_list => 5

theorem singletonCurrentCarrier_count :
    Fintype.card SingletonCurrentCarrier = 8 := by
  native_decide

theorem mirror_unit_singleton_current_packet_score_maximal :
    forall carrier : SingletonCurrentCarrier,
      singletonCurrentCarrierScore carrier ≤
        singletonCurrentCarrierScore
          SingletonCurrentCarrier.mirror_unit_singleton_current_packet := by
  native_decide

theorem raw_row_name_list_score_minimal :
    forall carrier : SingletonCurrentCarrier,
      singletonCurrentCarrierScore SingletonCurrentCarrier.raw_row_name_list ≤
        singletonCurrentCarrierScore carrier := by
  native_decide

theorem random031_control_ties_route_sidecar_score :
    singletonCurrentCarrierScore
        SingletonCurrentCarrier.random031_hard_control_clause =
      singletonCurrentCarrierScore
        SingletonCurrentCarrier.route_sidecar_R_join := by
  decide

/-! ## Axiom audit hooks -/

#print axioms hyp3480_target_row_count
#print axioms hyp3480_random031_is_control
#print axioms hyp3480_audited_row_partition
#print axioms hyp3480_all_target_components_have_branch_unit_touch
#print axioms hyp3480_all_target_mirror_pairs_have_unit_gate
#print axioms hyp3480_dispatch_complete
#print axioms hyp3480_dispatch_matches_counts
#print axioms singletonCurrentCarrier_count
#print axioms mirror_unit_singleton_current_packet_score_maximal
#print axioms raw_row_name_list_score_minimal

end LRC14
end LonelyRunner
