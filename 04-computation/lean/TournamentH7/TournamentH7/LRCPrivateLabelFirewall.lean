/-
  TournamentH7.LRCPrivateLabelFirewall -- HYP-3490 finite ledger for the
  private-label firewall split.

  This module does not prove the graph-theoretic multiplicity lemma.  It
  records the exact HYP-3490 dispatch arithmetic as a sorry-free Lean ledger:

    seven random projection-edge exceptions
      -> private-label firewall rows
      -> one random031 hard/private overlap
      -> six small-touch/no-hard rows routed to HYP-3480
      -> two AP currentless rows remain non-private pair-current terminals

  The mathematical producer still needed is the private-label projection
  firewall lemma: if every E/branch-touched blocker label is private to a
  single dead component, deleting any union of such labels cannot remove a
  dead-cover projection edge.
-/

import Mathlib.Tactic
import TournamentH7.LRCSingletonCurrentLedger

namespace LonelyRunner
namespace LRC14

/-! ## HYP-3490 private-firewall rows -/

/-- The seven HYP-3472 random projection-edge exceptions whose E/branch-touched
blocker labels are private in HYP-3490. -/
inductive PrivateFirewallRow where
  | random001
  | random031
  | random039
  | random062
  | random074
  | random086
  | random101
deriving DecidableEq, Repr, Fintype

/-- All private-label firewall rows. -/
def hyp3490PrivateRows : List PrivateFirewallRow :=
  [PrivateFirewallRow.random001,
   PrivateFirewallRow.random031,
   PrivateFirewallRow.random039,
   PrivateFirewallRow.random062,
   PrivateFirewallRow.random074,
   PrivateFirewallRow.random086,
   PrivateFirewallRow.random101]

/-- The unique private-firewall row that also carries hard mirror-orbit debt. -/
def hyp3490HardPrivateOverlapRow : PrivateFirewallRow :=
  PrivateFirewallRow.random031

/-- The six private-firewall rows with no hard mirror-orbit debt, routed to
HYP-3480's singleton-current packet. -/
def hyp3490SmallTouchPrivateRows : List PrivateFirewallRow :=
  [PrivateFirewallRow.random001,
   PrivateFirewallRow.random039,
   PrivateFirewallRow.random062,
   PrivateFirewallRow.random074,
   PrivateFirewallRow.random086,
   PrivateFirewallRow.random101]

theorem hyp3490_private_row_count :
    hyp3490PrivateRows.length = 7 := by
  decide

theorem hyp3490_small_touch_private_row_count :
    hyp3490SmallTouchPrivateRows.length = 6 := by
  decide

theorem hyp3490_random031_not_small_touch_private :
    hyp3490HardPrivateOverlapRow ∉ hyp3490SmallTouchPrivateRows := by
  simp [hyp3490HardPrivateOverlapRow, hyp3490SmallTouchPrivateRows]

theorem hyp3490_private_rows_card :
    Fintype.card PrivateFirewallRow = 7 := by
  native_decide

/-! ## Row details and routes -/

/-- Terminal route class for the HYP-3490 frontier split over dead rows. -/
inductive PrivateFirewallRoute where
  | ordinary_or_hard_currented
  | ap84_edge_only_pair_current_terminal
  | random031_private_firewall_hard_overlap
  | small_touch_private_firewall
deriving DecidableEq, Repr, Fintype

/-- Exact HYP-3490 row detail for the seven private-firewall rows. -/
structure PrivateFirewallRowDetail where
  deadComponents : Nat
  projectionEdgeLabels : Nat
  eBranchGates : Nat
  touchedLabels : Nat
  multiplicityOneTouchedLabels : Nat
  hardOrbits : Nat
  hardMaxDelta : Nat
  route : PrivateFirewallRoute

/-- HYP-3490 row-local firewall details. -/
def hyp3490RowDetail : PrivateFirewallRow -> PrivateFirewallRowDetail
  | PrivateFirewallRow.random001 =>
      { deadComponents := 4
        projectionEdgeLabels := 0
        eBranchGates := 56
        touchedLabels := 8
        multiplicityOneTouchedLabels := 8
        hardOrbits := 0
        hardMaxDelta := 0
        route := PrivateFirewallRoute.small_touch_private_firewall }
  | PrivateFirewallRow.random031 =>
      { deadComponents := 4
        projectionEdgeLabels := 0
        eBranchGates := 88
        touchedLabels := 8
        multiplicityOneTouchedLabels := 8
        hardOrbits := 1
        hardMaxDelta := 7
        route := PrivateFirewallRoute.random031_private_firewall_hard_overlap }
  | PrivateFirewallRow.random039 =>
      { deadComponents := 2
        projectionEdgeLabels := 0
        eBranchGates := 12
        touchedLabels := 4
        multiplicityOneTouchedLabels := 4
        hardOrbits := 0
        hardMaxDelta := 0
        route := PrivateFirewallRoute.small_touch_private_firewall }
  | PrivateFirewallRow.random062 =>
      { deadComponents := 2
        projectionEdgeLabels := 0
        eBranchGates := 26
        touchedLabels := 4
        multiplicityOneTouchedLabels := 4
        hardOrbits := 0
        hardMaxDelta := 0
        route := PrivateFirewallRoute.small_touch_private_firewall }
  | PrivateFirewallRow.random074 =>
      { deadComponents := 2
        projectionEdgeLabels := 0
        eBranchGates := 58
        touchedLabels := 4
        multiplicityOneTouchedLabels := 4
        hardOrbits := 0
        hardMaxDelta := 0
        route := PrivateFirewallRoute.small_touch_private_firewall }
  | PrivateFirewallRow.random086 =>
      { deadComponents := 2
        projectionEdgeLabels := 0
        eBranchGates := 34
        touchedLabels := 4
        multiplicityOneTouchedLabels := 4
        hardOrbits := 0
        hardMaxDelta := 0
        route := PrivateFirewallRoute.small_touch_private_firewall }
  | PrivateFirewallRow.random101 =>
      { deadComponents := 2
        projectionEdgeLabels := 0
        eBranchGates := 48
        touchedLabels := 4
        multiplicityOneTouchedLabels := 4
        hardOrbits := 0
        hardMaxDelta := 0
        route := PrivateFirewallRoute.small_touch_private_firewall }

theorem hyp3490_all_private_rows_have_no_projection_edge_labels :
    forall row : PrivateFirewallRow,
      (hyp3490RowDetail row).projectionEdgeLabels = 0 := by
  native_decide

theorem hyp3490_all_touched_labels_are_private :
    forall row : PrivateFirewallRow,
      (hyp3490RowDetail row).multiplicityOneTouchedLabels =
        (hyp3490RowDetail row).touchedLabels := by
  native_decide

theorem hyp3490_random031_is_private_hard_overlap :
    (hyp3490RowDetail hyp3490HardPrivateOverlapRow).route =
      PrivateFirewallRoute.random031_private_firewall_hard_overlap ∧
    (hyp3490RowDetail hyp3490HardPrivateOverlapRow).hardOrbits = 1 ∧
    (hyp3490RowDetail hyp3490HardPrivateOverlapRow).hardMaxDelta = 7 := by
  decide

theorem hyp3490_random031_unique_private_hard_overlap :
    forall row : PrivateFirewallRow,
      (hyp3490RowDetail row).hardOrbits = 1 ->
        row = hyp3490HardPrivateOverlapRow := by
  native_decide

theorem hyp3490_small_touch_private_rows_have_no_hard_orbits :
    forall row : PrivateFirewallRow,
      row ∈ hyp3490SmallTouchPrivateRows ->
        (hyp3490RowDetail row).hardOrbits = 0 := by
  native_decide

/-! ## Aggregate frontier ledger -/

/-- Exact aggregate count ledger from
`lrc14_private_label_firewall_codex_20260629.py`. -/
structure PrivateLabelFirewallCounts where
  rowsAudited : Nat
  deadRows : Nat
  deadRowsWithProjectionEdgeCut : Nat
  deadRowsWithSeparatingCurrent : Nat
  privateFirewallRows : Nat
  sharedTouchRows : Nat
  mismatchPrivateFirewallVsNoEdgeCut : Nat
  mismatchSharedTouchVsEdgeCut : Nat
  ordinaryOrHardCurrentedRows : Nat
  ap84PairCurrentTerminalRows : Nat
  random031PrivateHardOverlapRows : Nat
  smallTouchPrivateRows : Nat
  noSeparatingRows : Nat

/-- HYP-3490 aggregate counts. -/
def hyp3490Counts : PrivateLabelFirewallCounts where
  rowsAudited := 135
  deadRows := 130
  deadRowsWithProjectionEdgeCut := 123
  deadRowsWithSeparatingCurrent := 121
  privateFirewallRows := 7
  sharedTouchRows := 123
  mismatchPrivateFirewallVsNoEdgeCut := 0
  mismatchSharedTouchVsEdgeCut := 0
  ordinaryOrHardCurrentedRows := 121
  ap84PairCurrentTerminalRows := 2
  random031PrivateHardOverlapRows := 1
  smallTouchPrivateRows := 6
  noSeparatingRows := 9

theorem hyp3490_projection_edge_frontier_partition :
    hyp3490Counts.deadRowsWithProjectionEdgeCut +
      hyp3490Counts.privateFirewallRows =
        hyp3490Counts.deadRows := by
  decide

theorem hyp3490_shared_touch_partition_matches_edge_cuts :
    hyp3490Counts.sharedTouchRows = hyp3490Counts.deadRowsWithProjectionEdgeCut ∧
      hyp3490Counts.mismatchPrivateFirewallVsNoEdgeCut = 0 ∧
      hyp3490Counts.mismatchSharedTouchVsEdgeCut = 0 := by
  decide

theorem hyp3490_route_split_complete :
    hyp3490Counts.ordinaryOrHardCurrentedRows +
        hyp3490Counts.ap84PairCurrentTerminalRows +
        hyp3490Counts.random031PrivateHardOverlapRows +
        hyp3490Counts.smallTouchPrivateRows =
      hyp3490Counts.deadRows := by
  decide

theorem hyp3490_private_route_split_complete :
    hyp3490Counts.random031PrivateHardOverlapRows +
        hyp3490Counts.smallTouchPrivateRows =
      hyp3490Counts.privateFirewallRows := by
  decide

theorem hyp3490_no_separating_rows_split :
    hyp3490Counts.ap84PairCurrentTerminalRows +
        hyp3490Counts.privateFirewallRows =
      hyp3490Counts.noSeparatingRows := by
  decide

/-! ## Dispatch packet -/

/-- HYP-3490 dispatch packet for the currentless frontier. -/
structure PrivateLabelFirewallDispatchPacket where
  ordinaryOrHardCurrentedRows : Nat
  ap84PairCurrentTerminalRows : Nat
  random031PrivateHardOverlapRows : Nat
  smallTouchPrivateRows : Nat
  deadRows : Nat
  complete :
    ordinaryOrHardCurrentedRows +
        ap84PairCurrentTerminalRows +
        random031PrivateHardOverlapRows +
        smallTouchPrivateRows =
      deadRows

/-- The HYP-3490 split: 121 rows already currented, 2 AP pair-current terminals,
one private/hard random031 overlap, and six private/no-hard rows. -/
def hyp3490DispatchPacket : PrivateLabelFirewallDispatchPacket where
  ordinaryOrHardCurrentedRows := 121
  ap84PairCurrentTerminalRows := 2
  random031PrivateHardOverlapRows := 1
  smallTouchPrivateRows := 6
  deadRows := 130
  complete := by decide

theorem hyp3490_dispatch_complete :
    hyp3490DispatchPacket.ordinaryOrHardCurrentedRows +
        hyp3490DispatchPacket.ap84PairCurrentTerminalRows +
        hyp3490DispatchPacket.random031PrivateHardOverlapRows +
        hyp3490DispatchPacket.smallTouchPrivateRows =
      hyp3490DispatchPacket.deadRows :=
  hyp3490DispatchPacket.complete

theorem hyp3490_dispatch_matches_counts :
    hyp3490DispatchPacket.deadRows = hyp3490Counts.deadRows ∧
      hyp3490DispatchPacket.ap84PairCurrentTerminalRows =
        hyp3490Counts.ap84PairCurrentTerminalRows ∧
      hyp3490DispatchPacket.random031PrivateHardOverlapRows =
        hyp3490Counts.random031PrivateHardOverlapRows ∧
      hyp3490DispatchPacket.smallTouchPrivateRows =
        hyp3490Counts.smallTouchPrivateRows := by
  decide

/-! ## Formal producer obligation -/

/-- Data that the still-open private-label firewall lemma must consume. -/
structure PrivateLabelFirewallObligation where
  everyTouchedLabelPrivateToOneDeadComponent : Prop
  deletedLabelsAreTouchedEBranchLabels : Prop
  projectionEdgeRequiresSharedTouchedLabel : Prop
  routeSidecarPreserved : Prop

/-- The proof target left by HYP-3490.  The terminal predicate is abstract so
future graph-current, owner-current, two-adic, or SPEC proofs can instantiate
it. -/
def PrivateLabelProjectionFirewallLemmaTarget
    (NoProjectionEdgeDeleted : PrivateLabelFirewallObligation -> Prop) : Prop :=
  forall packet : PrivateLabelFirewallObligation,
    packet.everyTouchedLabelPrivateToOneDeadComponent ->
    packet.deletedLabelsAreTouchedEBranchLabels ->
    packet.projectionEdgeRequiresSharedTouchedLabel ->
    packet.routeSidecarPreserved ->
    NoProjectionEdgeDeleted packet

/-! ## Tournament carrier ledger -/

inductive PrivateLabelFirewallCarrier where
  | private_label_firewall_theorem
  | random031_private_hard_overlap
  | small_touch_private_packet
  | ap84_nonprivate_pair_current
  | edge_support_label_axis
  | raw_pair_current_count
  | raw_exception_name
deriving DecidableEq, Repr, Fintype

/-- HYP-3490 proof-carrier scores from the stored tournament audit. -/
def privateLabelFirewallCarrierScore : PrivateLabelFirewallCarrier -> Nat
  | PrivateLabelFirewallCarrier.private_label_firewall_theorem => 67
  | PrivateLabelFirewallCarrier.random031_private_hard_overlap => 67
  | PrivateLabelFirewallCarrier.small_touch_private_packet => 64
  | PrivateLabelFirewallCarrier.ap84_nonprivate_pair_current => 59
  | PrivateLabelFirewallCarrier.edge_support_label_axis => 58
  | PrivateLabelFirewallCarrier.raw_pair_current_count => 17
  | PrivateLabelFirewallCarrier.raw_exception_name => 7

theorem privateLabelFirewallCarrier_count :
    Fintype.card PrivateLabelFirewallCarrier = 7 := by
  native_decide

theorem private_label_firewall_theorem_score_maximal :
    forall carrier : PrivateLabelFirewallCarrier,
      privateLabelFirewallCarrierScore carrier ≤
        privateLabelFirewallCarrierScore
          PrivateLabelFirewallCarrier.private_label_firewall_theorem := by
  native_decide

theorem random031_private_hard_overlap_ties_firewall_theorem_score :
    privateLabelFirewallCarrierScore
        PrivateLabelFirewallCarrier.random031_private_hard_overlap =
      privateLabelFirewallCarrierScore
        PrivateLabelFirewallCarrier.private_label_firewall_theorem := by
  decide

theorem raw_exception_name_score_minimal :
    forall carrier : PrivateLabelFirewallCarrier,
      privateLabelFirewallCarrierScore
          PrivateLabelFirewallCarrier.raw_exception_name ≤
        privateLabelFirewallCarrierScore carrier := by
  native_decide

/-! ## Axiom audit hooks -/

#print axioms hyp3490_private_row_count
#print axioms hyp3490_all_touched_labels_are_private
#print axioms hyp3490_random031_is_private_hard_overlap
#print axioms hyp3490_random031_unique_private_hard_overlap
#print axioms hyp3490_projection_edge_frontier_partition
#print axioms hyp3490_route_split_complete
#print axioms hyp3490_dispatch_complete
#print axioms hyp3490_dispatch_matches_counts
#print axioms privateLabelFirewallCarrier_count
#print axioms private_label_firewall_theorem_score_maximal
#print axioms raw_exception_name_score_minimal

end LRC14
end LonelyRunner
