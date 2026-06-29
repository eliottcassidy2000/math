/-
  TournamentH7.LRCHardOrbitCurrentJoin -- HYP-3479 finite ledger for joining
  hard mirror-orbit debt to boundary-current exits.

  This module does not prove the geometric separating-current transfer. It
  records the exact HYP-3479 dispatch arithmetic as a sorry-free Lean ledger:

    hard mirror orbits
      -> separating-current transfer for seven orbits
      -> random031 named clause for the one remaining orbit

  The mathematical producer still needed is the separating-current transfer
  theorem referenced in HYP-3472/HYP-3479.
-/

import Mathlib.Tactic
import TournamentH7.LRCColoredGateFormalization

namespace LonelyRunner
namespace LRC14

/-! ## HYP-3479 finite counts -/

/-- Count ledger for the HYP-3479 hard-orbit/current join. -/
structure HardOrbitCurrentJoinCounts where
  hardOrbits : Nat
  hardRows : Nat
  sameBranchHardOrbits : Nat
  crossBranchHardOrbits : Nat
  hardRowsWithProjectionEdgeCut : Nat
  hardRowsWithoutProjectionEdgeCut : Nat
  hardRowsWithSeparatingCurrent : Nat
  hardRowsWithoutSeparatingCurrent : Nat
  hardOrbitsWithSeparatingCurrent : Nat
  hardOrbitsWithoutSeparatingCurrent : Nat
  edgeExceptionsWithoutHardOrbit : Nat
  separatingExceptionsWithoutHardOrbit : Nat
  ap84HardRows : Nat

/-- Exact aggregate counts from
`lrc14_hard_orbit_current_join_codex_20260629.py`. -/
def hyp3479JoinCounts : HardOrbitCurrentJoinCounts where
  hardOrbits := 8
  hardRows := 7
  sameBranchHardOrbits := 6
  crossBranchHardOrbits := 2
  hardRowsWithProjectionEdgeCut := 6
  hardRowsWithoutProjectionEdgeCut := 1
  hardRowsWithSeparatingCurrent := 6
  hardRowsWithoutSeparatingCurrent := 1
  hardOrbitsWithSeparatingCurrent := 7
  hardOrbitsWithoutSeparatingCurrent := 1
  edgeExceptionsWithoutHardOrbit := 6
  separatingExceptionsWithoutHardOrbit := 8
  ap84HardRows := 0

theorem hyp3479_hard_orbit_class_partition :
    hyp3479JoinCounts.sameBranchHardOrbits +
      hyp3479JoinCounts.crossBranchHardOrbits =
        hyp3479JoinCounts.hardOrbits := by
  decide

theorem hyp3479_projection_edge_row_partition :
    hyp3479JoinCounts.hardRowsWithProjectionEdgeCut +
      hyp3479JoinCounts.hardRowsWithoutProjectionEdgeCut =
        hyp3479JoinCounts.hardRows /\
    hyp3479JoinCounts.hardRowsWithoutProjectionEdgeCut = 1 := by
  decide

theorem hyp3479_separating_current_row_partition :
    hyp3479JoinCounts.hardRowsWithSeparatingCurrent +
      hyp3479JoinCounts.hardRowsWithoutSeparatingCurrent =
        hyp3479JoinCounts.hardRows /\
    hyp3479JoinCounts.hardRowsWithoutSeparatingCurrent = 1 := by
  decide

theorem hyp3479_separating_current_orbit_partition :
    hyp3479JoinCounts.hardOrbitsWithSeparatingCurrent +
      hyp3479JoinCounts.hardOrbitsWithoutSeparatingCurrent =
        hyp3479JoinCounts.hardOrbits /\
    hyp3479JoinCounts.hardOrbitsWithoutSeparatingCurrent = 1 := by
  decide

theorem hyp3479_no_ap84_hard_rows :
    hyp3479JoinCounts.ap84HardRows = 0 := by
  decide

/-! ## Dispatch packet -/

/-- Terminal labels for the HYP-3479 hard-orbit dispatch. -/
inductive HardOrbitTerminalExit where
  | boundary_current_transfer
  | random031_named_clause
deriving DecidableEq, Repr, Fintype

/-- A finite dispatch packet for the hard-orbit family. -/
structure HardOrbitDispatchPacket where
  boundaryCurrentTransferOrbits : Nat
  random031ClauseOrbits : Nat
  totalHardOrbits : Nat
  complete :
    boundaryCurrentTransferOrbits + random031ClauseOrbits = totalHardOrbits

/-- HYP-3479 dispatches seven hard orbits by boundary current and leaves one
named random031 clause. -/
def hyp3479DispatchPacket : HardOrbitDispatchPacket where
  boundaryCurrentTransferOrbits := 7
  random031ClauseOrbits := 1
  totalHardOrbits := 8
  complete := by decide

theorem hyp3479_dispatch_complete :
    hyp3479DispatchPacket.boundaryCurrentTransferOrbits +
      hyp3479DispatchPacket.random031ClauseOrbits =
        hyp3479DispatchPacket.totalHardOrbits :=
  hyp3479DispatchPacket.complete

theorem hyp3479_dispatch_matches_join_counts :
    hyp3479DispatchPacket.totalHardOrbits = hyp3479JoinCounts.hardOrbits /\
    hyp3479DispatchPacket.boundaryCurrentTransferOrbits =
      hyp3479JoinCounts.hardOrbitsWithSeparatingCurrent /\
    hyp3479DispatchPacket.random031ClauseOrbits =
      hyp3479JoinCounts.hardOrbitsWithoutSeparatingCurrent := by
  decide

/-! ## Tournament carrier ledger -/

inductive HardOrbitCurrentCarrier where
  | hard_orbit_current_join
  | singleton_intersection_ledger
  | separating_current_transfer
  | random031_named_gluing_debt
  | hard_mirror_orbit_ledger
  | dead_touch_gate_universal_lemma
  | raw_exception_set_overlap
  | raw_hard_delta_count
deriving DecidableEq, Repr, Fintype

/-- HYP-3479 proof-carrier scores from the stored tournament audit. -/
def hardOrbitCurrentCarrierScore : HardOrbitCurrentCarrier -> Nat
  | HardOrbitCurrentCarrier.hard_orbit_current_join => 60
  | HardOrbitCurrentCarrier.singleton_intersection_ledger => 58
  | HardOrbitCurrentCarrier.separating_current_transfer => 54
  | HardOrbitCurrentCarrier.random031_named_gluing_debt => 54
  | HardOrbitCurrentCarrier.hard_mirror_orbit_ledger => 47
  | HardOrbitCurrentCarrier.dead_touch_gate_universal_lemma => 42
  | HardOrbitCurrentCarrier.raw_exception_set_overlap => 29
  | HardOrbitCurrentCarrier.raw_hard_delta_count => 19

theorem hardOrbitCurrentCarrier_count :
    Fintype.card HardOrbitCurrentCarrier = 8 := by
  native_decide

theorem hard_orbit_current_join_score_maximal :
    forall carrier : HardOrbitCurrentCarrier,
      hardOrbitCurrentCarrierScore carrier <=
        hardOrbitCurrentCarrierScore
          HardOrbitCurrentCarrier.hard_orbit_current_join := by
  native_decide

theorem raw_hard_delta_count_score_minimal :
    forall carrier : HardOrbitCurrentCarrier,
      hardOrbitCurrentCarrierScore
          HardOrbitCurrentCarrier.raw_hard_delta_count <=
        hardOrbitCurrentCarrierScore carrier := by
  native_decide

theorem separating_current_ties_random031_score :
    hardOrbitCurrentCarrierScore
        HardOrbitCurrentCarrier.separating_current_transfer =
      hardOrbitCurrentCarrierScore
        HardOrbitCurrentCarrier.random031_named_gluing_debt := by
  decide

/-! ## Axiom audit hooks -/

#print axioms hyp3479_hard_orbit_class_partition
#print axioms hyp3479_separating_current_orbit_partition
#print axioms hyp3479_dispatch_complete
#print axioms hyp3479_dispatch_matches_join_counts
#print axioms hardOrbitCurrentCarrier_count
#print axioms hard_orbit_current_join_score_maximal
#print axioms raw_hard_delta_count_score_minimal

end LRC14
end LonelyRunner
