/-
  TournamentH7.LRCProofFrontier -- S259 / HYP-3107 Lean-facing frontier ledger
  for the current LRC(14) proof edge.

  This module does not prove LRC(14).  It records the current split after the
  S31ag/S63/S258/S64 updates:

    * the cap RHS is solved as a pair-Pascal / THM-576 field, with THM-577
      upgrading the k=10,11 values to symbolic closed-form overlap values;
    * `LRCObserverGluingLedger` is the current producer interface for turning
      chart overlaps into terminal certificates;
    * HYP-3099 shows cap optimality is a bounded but non-transitive exchange
      problem, and that apex-7-to-H=7 is only a numerical coincidence;
    * HYP-3100 / `LRCBleedingEdgeFrontier` packages observer, equivalence,
      Pascal, polynomial, and moment-degree sidecars around finite-address
      packets.
    * HYP-3093/HYP-3097 separate equinumerosity, equidecomposability, and
      equidistribution before a quotient is allowed to forget data.
    * the live bounded-core crux is coverage extremality at k=8,9,10;
    * the unbounded twin is the Node-3 effective Erdős-Turán peel;
    * the coarse winding-tournament H bridge degenerates at the apex-7 ties, so
      any tournament transfer must move to a fine-scale observable.

  The proofs below are sorry-free logical glue and exact rational checks.  The
  analytic/combinatorial frontier remains as named `Prop` obligations.
-/

import TournamentH7.LRCFourteenSkeleton
import TournamentH7.LRCFiniteAddressBranchClosure
import TournamentH7.LRCObserverGluingLedger
import TournamentH7.LRCBleedingEdgeFrontier

namespace LonelyRunner
namespace LRC14

open scoped BigOperators

/-! ## Frontier nodes and experiment vertices -/

/-- Status labels for the current LRC14 frontier. -/
inductive FrontierStatus where
  | solved
  | finite_import
  | open_obligation
  | false_or_degenerate_route
  | conditional_glue
deriving DecidableEq, Repr, Fintype

/-- Current proof-frontier nodes.  These are proof obligations / observers, not
runners or raw residue classes. -/
inductive FrontierNode where
  | q_witness_gate
  | level7_lift_sieve
  | pair_pascal_cap_rhs
  | bounded_coverage_extremality
  | reflection_perron_certificate
  | node3_effective_peel
  | fine_modp_winding_transfer
  | finite_ruler_glue
  | finite_address_packet_glue
  | terminal_mreach_readout
deriving DecidableEq, Repr, Fintype

/-- The S259 frontier ledger has ten named proof nodes. -/
theorem frontierNode_count : Fintype.card FrontierNode = 10 := by
  native_decide

/-- Status assignment for the current state-of-proof map. -/
def frontierStatus : FrontierNode -> FrontierStatus
  | FrontierNode.q_witness_gate => FrontierStatus.solved
  | FrontierNode.level7_lift_sieve => FrontierStatus.solved
  | FrontierNode.pair_pascal_cap_rhs => FrontierStatus.solved
  | FrontierNode.bounded_coverage_extremality => FrontierStatus.open_obligation
  | FrontierNode.reflection_perron_certificate => FrontierStatus.open_obligation
  | FrontierNode.node3_effective_peel => FrontierStatus.open_obligation
  | FrontierNode.fine_modp_winding_transfer => FrontierStatus.open_obligation
  | FrontierNode.finite_ruler_glue => FrontierStatus.open_obligation
  | FrontierNode.finite_address_packet_glue => FrontierStatus.conditional_glue
  | FrontierNode.terminal_mreach_readout => FrontierStatus.solved

/-- The solved nodes are exactly the already-proved gates/readouts in this
frontier ledger. -/
theorem solved_frontier_nodes :
    (Finset.univ.filter
      (fun v : FrontierNode => frontierStatus v == FrontierStatus.solved)).card = 4 := by
  native_decide

/-- Creative next-test families suggested by the current frontier. -/
inductive FrontierExperiment where
  | eberlein_hankel_dip
  | reflection_perron_dual
  | fine_modp_winding_H
  | observer_sheaf_overlap
  | finite_ruler_sampler
  | node3_erdos_turan_tail
  | pascal_scissors_nonequidecomp
  | equivalence_triad_audit
deriving DecidableEq, Repr, Fintype

/-- Eight concrete experiment vertices are tracked for the next frontier pass. -/
theorem frontierExperiment_count :
    Fintype.card FrontierExperiment = 8 := by
  native_decide

/-- Each experiment is attached to the proof node it is meant to refine. -/
def experimentPrimaryNode : FrontierExperiment -> FrontierNode
  | FrontierExperiment.eberlein_hankel_dip =>
      FrontierNode.pair_pascal_cap_rhs
  | FrontierExperiment.reflection_perron_dual =>
      FrontierNode.reflection_perron_certificate
  | FrontierExperiment.fine_modp_winding_H =>
      FrontierNode.fine_modp_winding_transfer
  | FrontierExperiment.observer_sheaf_overlap =>
      FrontierNode.finite_address_packet_glue
  | FrontierExperiment.finite_ruler_sampler =>
      FrontierNode.finite_ruler_glue
  | FrontierExperiment.node3_erdos_turan_tail =>
      FrontierNode.node3_effective_peel
  | FrontierExperiment.pascal_scissors_nonequidecomp =>
      FrontierNode.bounded_coverage_extremality
  | FrontierExperiment.equivalence_triad_audit =>
      FrontierNode.finite_address_packet_glue

/-! ## Solved cap RHS as pair-Pascal arithmetic -/

/-- The pair-Pascal mass `C(k+1,2) / C(14,2)`, written as
`k*(k+1)/(2*91)` for exact rational normalization. -/
def pairPascalMassRat (k : Nat) : Rat :=
  (((k + 1) * k : Nat) : Rat) / 182

/-- In the dense rows, the cap table agrees with the pair-Pascal mass.  THM-577
upgrades the k=10,11 entries from search-backed dense values to symbolic
closed-form overlap values; this Lean lemma records the exact rational ledger. -/
theorem capRat_eq_pairPascalMassRat_dense
    (k : Nat) (h10 : 10 <= k) (h13 : k <= 13) :
    capRat k = pairPascalMassRat k := by
  interval_cases k <;> norm_num [capRat, pairPascalMassRat]

/-- The k=9 binding row is one `1/4004` unit below pure pair-Pascal mass. -/
theorem capRat_k9_pairPascal_debt :
    capRat 9 = pairPascalMassRat 9 - (1 : Rat) / 4004 := by
  norm_num [capRat, pairPascalMassRat]

/-- The k=8 binding row carries the larger higher-order debt recorded in S63. -/
theorem capRat_k8_pairPascal_debt :
    capRat 8 = pairPascalMassRat 8 - (1081 : Rat) / 76440 := by
  norm_num [capRat, pairPascalMassRat]

/-- A cap ledger keeps the count shadow and the higher-order debt separate. -/
structure PairPascalCapLedger where
  k : Nat
  cap : Rat
  pairMass : Rat
  higherOrderDebt : Rat
  cap_eq_pairMass_sub_debt : cap = pairMass - higherOrderDebt

/-- Solved cap ledger for k=10: no higher-order debt. -/
def capLedgerK10 : PairPascalCapLedger where
  k := 10
  cap := capRat 10
  pairMass := pairPascalMassRat 10
  higherOrderDebt := 0
  cap_eq_pairMass_sub_debt := by
    norm_num [capRat, pairPascalMassRat]

/-- Binding cap ledger for k=9: the first higher-order defect is explicit. -/
def capLedgerK9 : PairPascalCapLedger where
  k := 9
  cap := capRat 9
  pairMass := pairPascalMassRat 9
  higherOrderDebt := (1 : Rat) / 4004
  cap_eq_pairMass_sub_debt := capRat_k9_pairPascal_debt

/-- Binding cap ledger for k=8: the order-4 debt is explicit. -/
def capLedgerK8 : PairPascalCapLedger where
  k := 8
  cap := capRat 8
  pairMass := pairPascalMassRat 8
  higherOrderDebt := (1081 : Rat) / 76440
  cap_eq_pairMass_sub_debt := capRat_k8_pairPascal_debt

/-! ## Open frontier obligations as formal predicates -/

/-- A coverage-extremality statement: consecutive/AP clusters dominate the
chosen coverage observable for k=8,9,10.  The observable is intentionally an
argument so future work can instantiate it with `p0`, a moment dual, or a
finite exact table. -/
def CoverageExtremality
    (coverageMass : List Int -> Real)
    (consecutiveCluster : Nat -> List Int) : Prop :=
  forall k : Nat, k = 8 ∨ k = 9 ∨ k = 10 ->
    forall E : List Int, E.length = k ->
      coverageMass E <= coverageMass (consecutiveCluster k)

/-- Node-3 effective peel: every row in the declared Node-3 tail has a direct
Mreach lower bound.  This is the unbounded twin of the bounded-core crux. -/
def Node3EffectivePeel
    (isNode3Tail : (Fin 13 -> Int) -> Prop) : Prop :=
  forall v : Fin 13 -> Int, (forall i, v i ≠ 0) ->
    isNode3Tail v -> (1 : Real) / 14 <= Mreach v

/-- A residual classifier separates rows already discharged by early gates from
rows that must enter the finite-address residual packet. -/
def FrontierResidualClassifier
    (isResidual : (Fin 13 -> Int) -> Prop) : Prop :=
  forall v : Fin 13 -> Int, (forall i, v i ≠ 0) ->
    Nonempty (EarlyGateCertificate v) ∨ isResidual v

/-- The residual packet theorem still missing at the frontier. -/
def ResidualToFiniteAddressPackets
    (isResidual : (Fin 13 -> Int) -> Prop) : Prop :=
  forall v : Fin 13 -> Int, (forall i, v i ≠ 0) ->
    isResidual v -> Nonempty (FiniteAddressBranchPacket v)

/-- A broader residual producer target: residual rows may emit proof-bearing
observer-gluing certificates before being compressed to finite-address branch
packets. -/
def ResidualToObserverGluingCertificates
    (isResidual : (Fin 13 -> Int) -> Prop) : Prop :=
  forall v : Fin 13 -> Int, (forall i, v i ≠ 0) ->
    isResidual v -> Nonempty (ObserverGluingCertificate v)

/-- Coarse winding-tournament transfer is allowed only after a no-tie guard.
S31ag found the apex-7 coarse tournament degenerates at k>=8; this predicate is
the formal replacement hook for a future fine-scale mod-p observable. -/
def FineScaleWindingTransfer
    (windingObservable : List Int -> Nat)
    (coverageMass : List Int -> Real) : Prop :=
  forall E E' : List Int,
    windingObservable E <= windingObservable E' ->
      coverageMass E <= coverageMass E'

/-- The full proof-frontier packet: fields marked as open are retained as data,
while the final theorem below uses only the already-formal packet classifier
and terminal readout. -/
structure BleedingEdgeFrontier where
  residual : (Fin 13 -> Int) -> Prop
  coverageMass : List Int -> Real
  consecutiveCluster : Nat -> List Int
  isNode3Tail : (Fin 13 -> Int) -> Prop
  windingObservable : List Int -> Nat
  residualClassifier : FrontierResidualClassifier residual
  residualToPacket : ResidualToFiniteAddressPackets residual
  residualToObserverGluing : ResidualToObserverGluingCertificates residual
  coverageExtremality :
    CoverageExtremality coverageMass consecutiveCluster
  node3Peel : Node3EffectivePeel isNode3Tail
  fineScaleWinding :
    FineScaleWindingTransfer windingObservable coverageMass

/-! ## Conditional assembly to the existing LRC14 skeleton -/

/-- If all rows are either early-gate discharged or emit a finite-address packet,
then the existing cutting-edge branch theorem proves LRC14. -/
theorem lrc14_from_residual_packet_frontier
    {isResidual : (Fin 13 -> Int) -> Prop}
    (hclass : FrontierResidualClassifier isResidual)
    (hpacket : ResidualToFiniteAddressPackets isResidual) :
    LRC14Statement := by
  apply lrc14_from_cutting_edge_branch_coverage
  intro v hv
  rcases hclass v hv with hEarly | hResidual
  · exact Or.inl hEarly
  · exact Or.inr (hpacket v hv hResidual)

/-- If every residual row emits an observer-gluing certificate, the S259
observer-gluing interface proves LRC14 directly. -/
theorem lrc14_from_residual_observer_gluing_frontier
    {isResidual : (Fin 13 -> Int) -> Prop}
    (hclass : FrontierResidualClassifier isResidual)
    (hglue : ResidualToObserverGluingCertificates isResidual) :
    LRC14Statement := by
  apply lrc14_from_observer_gluing_coverage
  intro v hv
  rcases hclass v hv with hEarly | hResidual
  · exact Or.inl hEarly
  · exact Or.inr (hglue v hv hResidual)

/-- The packaged S259 frontier implies LRC14 exactly when its residual packet
field is filled.  No analytic fact is asserted here; the theorem is only the
formal wiring from the frontier interface to the existing `Mreach` readout. -/
theorem lrc14_from_bleeding_edge_frontier
    (frontier : BleedingEdgeFrontier) :
    LRC14Statement :=
  lrc14_from_residual_packet_frontier
    frontier.residualClassifier frontier.residualToPacket

/-- The same frontier can also be closed through observer-gluing terminal
certificates.  This is the formal bridge from the S258/S259 observer ledger to
the top-level LRC14 statement. -/
theorem lrc14_from_bleeding_edge_observer_gluing_frontier
    (frontier : BleedingEdgeFrontier) :
    LRC14Statement :=
  lrc14_from_residual_observer_gluing_frontier
    frontier.residualClassifier frontier.residualToObserverGluing

/-- The upstream bleeding-edge packet wrapper is another conservative terminal
route for this frontier: once every row has a decorated packet, the top-level
statement follows through the embedded finite-address certificate. -/
theorem lrc14_from_bleeding_edge_packet_wrapper_frontier
    (hcoverage : BleedingEdgeFrontierCoverage) :
    LRC14Statement :=
  lrc14_from_bleeding_edge_frontier_coverage hcoverage

/-! ## Axiom audit hooks -/

#print axioms frontierNode_count
#print axioms solved_frontier_nodes
#print axioms frontierExperiment_count
#print axioms capRat_eq_pairPascalMassRat_dense
#print axioms capRat_k9_pairPascal_debt
#print axioms capRat_k8_pairPascal_debt
#print axioms lrc14_from_residual_packet_frontier
#print axioms lrc14_from_residual_observer_gluing_frontier
#print axioms lrc14_from_bleeding_edge_frontier
#print axioms lrc14_from_bleeding_edge_observer_gluing_frontier
#print axioms lrc14_from_bleeding_edge_packet_wrapper_frontier

end LRC14
end LonelyRunner
