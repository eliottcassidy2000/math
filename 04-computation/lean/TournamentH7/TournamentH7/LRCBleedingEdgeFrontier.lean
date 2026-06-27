/-
  TournamentH7.LRCBleedingEdgeFrontier -- Lean-facing wrapper for the current
  LRC14 observer-gluing frontier.

  This file does not assert the open coverage/extremality theorem.  It packages
  the invariants now appearing in HYP-3096/HYP-3097 and the S258 observer
  ledger as proof-carrying sidecars around the existing finite-address branch
  packet:

    finite-address packet
      + observer charts
      + equivalence triad
      + polynomial/witness-route ledger
      + Pascal pair-mass ledger
      + moment-degree ledger
      -> concrete Mreach >= 1/14

  The soundness theorem is intentionally conservative: the new frontier packet
  proves LRC14 only through the embedded `FiniteAddressBranchPacket`.
-/

import TournamentH7.LRCFiniteAddressBranchClosure

namespace LonelyRunner
namespace LRC14

/-! ## Observer charts and equivalence shadows -/

/-- Proof-frontier observer charts.  These are invariant carriers and proof
obligations, not runners, residues, or raw scalar caps. -/
inductive ObserverChartKind where
  | finite_address_sheaf
  | arithmetic_crt_lift
  | normalized_arc_ruler
  | direct_lonely_components
  | moment_perron
  | branch_k33_state
  | q_cusp_principal_part
  | equidistribution_measure
  | equinumerosity_pascal
  | equidecomposition_scissors
  | raw_scalar_shadow
deriving DecidableEq, Repr, Fintype

/-- The current observer-chart carrier has eleven named vertices. -/
theorem observerChartKind_count :
    Fintype.card ObserverChartKind = 11 := by
  decide

/-- Legal exits for a quotient that forgets a coordinate.  The last constructor
is deliberately illegal unless a chart proves it is not being used. -/
inductive ObserverChartExit where
  | reconstructible
  | dual_annihilated
  | fiber_constant
  | certificate_descent
  | terminal_atom
  | named_residual_debt
  | raw_unprotected
deriving DecidableEq, Repr, Fintype

/-- The three sameness shadows tracked by HYP-3097. -/
inductive EquivalenceShadow where
  | equidistribution
  | equinumerosity
  | equidecomposability
deriving DecidableEq, Repr, Fintype

/-- The LRC14 sameness audit has exactly the three shadows named in HYP-3097. -/
theorem equivalenceShadow_count :
    Fintype.card EquivalenceShadow = 3 := by
  decide

/-- A proof-frontier observer chart records what it preserves, what it forgets,
where it must overlap next, and why it is not an unprotected raw quotient. -/
structure ObserverChartLedger where
  chart : ObserverChartKind
  preservedLRPredicate : Bool
  destroyedCoordinate : DestroyedCoordinate
  overlapTarget : ObserverChartKind
  exit : ObserverChartExit
  retainedSidecarCount : Nat
  noRawUnprotected : exit ≠ ObserverChartExit.raw_unprotected

/-- Every observer chart exposes its "no raw unprotected quotient" guard. -/
theorem observerChartLedger_no_raw_unprotected
    (chart : ObserverChartLedger) :
    chart.exit ≠ ObserverChartExit.raw_unprotected :=
  chart.noRawUnprotected

/-- The equivalence triad is legal only when equidistribution is not used as a
bare scalar quotient. -/
structure EquivalenceTriadLedger where
  distributionLawRetained : Bool
  cardinalShadowRetained : Bool
  scissorsFiberRetained : Bool
  observerCutOrbitRetained : Bool
  interactionOrderDefect : Nat
  namedResidualDebt : Bool
  noBareEquidistribution :
    distributionLawRetained = true ->
      cardinalShadowRetained = true ∨
      scissorsFiberRetained = true ∨
      namedResidualDebt = true

/-- The triad guard says a distribution shadow needs either cardinal/scissors
payload or named residual debt. -/
theorem equivalenceTriadLedger_no_bare_distribution
    (ledger : EquivalenceTriadLedger)
    (h : ledger.distributionLawRetained = true) :
    ledger.cardinalShadowRetained = true ∨
      ledger.scissorsFiberRetained = true ∨
      ledger.namedResidualDebt = true :=
  ledger.noBareEquidistribution h

/-! ## Polynomial-route, Pascal, and moment-degree ledgers -/

/-- Exact nonnegative rational readout stored as numerator/denominator. -/
structure NonnegativeRationalReadout where
  num : Nat
  den : Nat
  den_pos : 0 < den

def zeroReadout : NonnegativeRationalReadout where
  num := 0
  den := 1
  den_pos := by decide

/-- HYP-3096 witness-route sidecar.  The component and denominator fields are
frontier obligations, not new terminal proofs. -/
structure PolynomialWitnessLedger where
  crtC7LiftStatus : Bool
  crtC2DyadicLiftStatus : Bool
  discreteGridStatus : Bool
  observedDirectComponentCount : Nat
  directLonelyMeasureFloor : NonnegativeRationalReadout
  largestLonelyArcFloor : NonnegativeRationalReadout
  denominatorNetThreshold : Nat
  finiteBadDenominatorBudget : Nat
  destroyedCoordinate : DestroyedCoordinate
  terminalExit : BranchTerminalExit

/-- The Pascal/pair-mass sidecar keeps the row-14 arithmetic as a proof-bearing
checksum instead of a floating numerical analogy. -/
structure PascalPairMassLedger where
  pairApex : Nat
  binom14_4 : Nat
  binom14_5 : Nat
  binom14_6 : Nat
  affineCompletion4004 : Nat
  cap9DefectDenominator : Nat
  cap8DefectNumerator : Nat
  cap8DefectDenominator : Nat
  pascal1001Is11Apex : binom14_4 = 11 * pairApex
  pascal2002Is22Apex : binom14_5 = 22 * pairApex
  pascal3003Is33Apex : binom14_6 = 33 * pairApex
  completion4004Is44Apex : affineCompletion4004 = 44 * pairApex
  completion4004IsPascalSum : affineCompletion4004 = binom14_4 + binom14_6

/-- The row-14 pair-mass checksum behind `1001,2002,3003,4004`. -/
def row14PascalPairMassLedger : PascalPairMassLedger where
  pairApex := 91
  binom14_4 := 1001
  binom14_5 := 2002
  binom14_6 := 3003
  affineCompletion4004 := 4004
  cap9DefectDenominator := 4004
  cap8DefectNumerator := 1081
  cap8DefectDenominator := 76440
  pascal1001Is11Apex := rfl
  pascal2002Is22Apex := rfl
  pascal3003Is33Apex := rfl
  completion4004Is44Apex := rfl
  completion4004IsPascalSum := rfl

/-- Readout tying the prompt constants to the pair-normalized Pascal ledger. -/
theorem row14PascalPairMassLedger_readout :
    row14PascalPairMassLedger.pairApex = 91 ∧
      row14PascalPairMassLedger.binom14_4 = 1001 ∧
      row14PascalPairMassLedger.binom14_5 = 2002 ∧
      row14PascalPairMassLedger.binom14_6 = 3003 ∧
      row14PascalPairMassLedger.affineCompletion4004 = 4004 ∧
      row14PascalPairMassLedger.cap9DefectDenominator = 4004 := by
  exact ⟨rfl, rfl, rfl, rfl, rfl, rfl⟩

/-- Moment-degree status for the cap RHS versus cover-bound LHS. -/
inductive MomentDegreeStatus where
  | degree2_pairwise_suffices
  | degree3_needed
  | degree4_needed
  | unresolved_frontier
deriving DecidableEq, Repr, Fintype

/-- S31ag/S258 warning: pairwise cap mass and cover-bound extremality can have
different moment-degree thresholds. -/
structure CoverBoundMomentDegreeLedger where
  bindingRowK : Nat
  capShadowDegree : Nat
  coverExtremalityDegree : Nat
  capStatus : MomentDegreeStatus
  coverStatus : MomentDegreeStatus
  rhsPairwiseCapVerified : Bool
  lhsDegree2Closes : Bool
  higherMomentDebt : Bool
  terminalExit : BranchTerminalExit

/-- The current k=10 frontier: cap side remains pairwise, while cover
extremality needs a higher-degree certificate. -/
def k10CoverBoundMomentDegreeLedger : CoverBoundMomentDegreeLedger where
  bindingRowK := 10
  capShadowDegree := 2
  coverExtremalityDegree := 3
  capStatus := MomentDegreeStatus.degree2_pairwise_suffices
  coverStatus := MomentDegreeStatus.degree3_needed
  rhsPairwiseCapVerified := true
  lhsDegree2Closes := false
  higherMomentDebt := true
  terminalExit := BranchTerminalExit.named_residual_debt

theorem k10CoverBoundMomentDegreeLedger_readout :
    k10CoverBoundMomentDegreeLedger.bindingRowK = 10 ∧
      k10CoverBoundMomentDegreeLedger.capShadowDegree = 2 ∧
      k10CoverBoundMomentDegreeLedger.coverExtremalityDegree = 3 ∧
      k10CoverBoundMomentDegreeLedger.rhsPairwiseCapVerified = true ∧
      k10CoverBoundMomentDegreeLedger.lhsDegree2Closes = false := by
  exact ⟨rfl, rfl, rfl, rfl, rfl⟩

/-! ## Default S258/S31ag frontier sidecars -/

def finiteAddressObserverChart : ObserverChartLedger where
  chart := ObserverChartKind.finite_address_sheaf
  preservedLRPredicate := true
  destroyedCoordinate := DestroyedCoordinate.none
  overlapTarget := ObserverChartKind.direct_lonely_components
  exit := ObserverChartExit.certificate_descent
  retainedSidecarCount := 6
  noRawUnprotected := by
    intro h
    cases h

def directArcDebtObserverChart : ObserverChartLedger where
  chart := ObserverChartKind.direct_lonely_components
  preservedLRPredicate := true
  destroyedCoordinate := DestroyedCoordinate.endpoint_owner
  overlapTarget := ObserverChartKind.normalized_arc_ruler
  exit := ObserverChartExit.named_residual_debt
  retainedSidecarCount := 5
  noRawUnprotected := by
    intro h
    cases h

def pascalScissorsObserverChart : ObserverChartLedger where
  chart := ObserverChartKind.equinumerosity_pascal
  preservedLRPredicate := false
  destroyedCoordinate := DestroyedCoordinate.raw_runner_identity
  overlapTarget := ObserverChartKind.equidecomposition_scissors
  exit := ObserverChartExit.fiber_constant
  retainedSidecarCount := 4
  noRawUnprotected := by
    intro h
    cases h

def rawScalarDebtObserverChart : ObserverChartLedger where
  chart := ObserverChartKind.raw_scalar_shadow
  preservedLRPredicate := false
  destroyedCoordinate := DestroyedCoordinate.formal_status
  overlapTarget := ObserverChartKind.finite_address_sheaf
  exit := ObserverChartExit.named_residual_debt
  retainedSidecarCount := 0
  noRawUnprotected := by
    intro h
    cases h

/-- A minimal observer-gluing spine matching the S258 ledger warning. -/
def defaultObserverChartsS258 : List ObserverChartLedger :=
  [finiteAddressObserverChart,
    directArcDebtObserverChart,
    pascalScissorsObserverChart,
    rawScalarDebtObserverChart]

theorem defaultObserverChartsS258_length :
    defaultObserverChartsS258.length = 4 := rfl

/-- Default triad sidecar: all three shadows are attached, with the first
higher-order defect named rather than forgotten. -/
def defaultEquivalenceTriadLedger : EquivalenceTriadLedger where
  distributionLawRetained := true
  cardinalShadowRetained := true
  scissorsFiberRetained := true
  observerCutOrbitRetained := true
  interactionOrderDefect := 3
  namedResidualDebt := true
  noBareEquidistribution := by
    intro _
    exact Or.inl rfl

/-- S258 divisor-loaded sample sidecar: the direct chart has positive computable
arcs, but the component count and denominator threshold remain proof debt. -/
def s258DirectArcDebtLedger : PolynomialWitnessLedger where
  crtC7LiftStatus := true
  crtC2DyadicLiftStatus := false
  discreteGridStatus := false
  observedDirectComponentCount := 860
  directLonelyMeasureFloor :=
    { num := 122687, den := 11771760, den_pos := by decide }
  largestLonelyArcFloor :=
    { num := 1, den := 82320, den_pos := by decide }
  denominatorNetThreshold := 82321
  finiteBadDenominatorBudget := 82320
  destroyedCoordinate := DestroyedCoordinate.endpoint_owner
  terminalExit := BranchTerminalExit.named_residual_debt

theorem s258DirectArcDebtLedger_readout :
    s258DirectArcDebtLedger.observedDirectComponentCount = 860 ∧
      s258DirectArcDebtLedger.largestLonelyArcFloor.num = 1 ∧
      s258DirectArcDebtLedger.largestLonelyArcFloor.den = 82320 ∧
      s258DirectArcDebtLedger.denominatorNetThreshold = 82321 := by
  exact ⟨rfl, rfl, rfl, rfl⟩

/-! ## Bleeding-edge packet closure -/

/-- The bleeding-edge frontier packet is a conservative wrapper around the
already-sound finite-address branch packet plus the new observer sidecars. -/
structure BleedingEdgeFrontierPacket (v : Fin 13 -> Int) where
  finiteAddress : FiniteAddressBranchPacket v
  observerCharts : List ObserverChartLedger
  equivalenceTriad : EquivalenceTriadLedger
  polynomialWitness : PolynomialWitnessLedger
  pascalPairMass : PascalPairMassLedger
  momentDegree : CoverBoundMomentDegreeLedger

/-- Soundness is inherited from the embedded finite-address terminal
certificate. -/
theorem bleedingEdgeFrontierPacket_mreach
    {v : Fin 13 -> Int}
    (packet : BleedingEdgeFrontierPacket v) :
    (1 : Real) / 14 <= Mreach v :=
  finiteAddressBranchPacket_mreach packet.finiteAddress

/-- Coverage by bleeding-edge frontier packets for every nonzero 13-speed row.
This is still a global theorem-to-prove, not an axiom introduced here. -/
def BleedingEdgeFrontierCoverage : Prop :=
  forall v : Fin 13 -> Int, (forall i, v i ≠ 0) ->
    Nonempty (BleedingEdgeFrontierPacket v)

/-- Bleeding-edge frontier coverage implies the top-level LRC14 statement. -/
theorem lrc14_from_bleeding_edge_frontier_coverage
    (hcoverage : BleedingEdgeFrontierCoverage) :
    LRC14Statement := by
  intro v hv
  rcases hcoverage v hv with ⟨packet⟩
  exact lonely_of_Mreach_ge v hv (bleedingEdgeFrontierPacket_mreach packet)

/-- Conservative packaging check: an already-known reach bound can be decorated
with the current frontier sidecars. -/
noncomputable def bleedingEdgeFrontierPacketOfMreach
    (v : Fin 13 -> Int)
    (hM : (1 : Real) / 14 <= Mreach v) :
    BleedingEdgeFrontierPacket v where
  finiteAddress := finiteAddressBranchPacketOfMreach v hM
  observerCharts := defaultObserverChartsS258
  equivalenceTriad := defaultEquivalenceTriadLedger
  polynomialWitness := s258DirectArcDebtLedger
  pascalPairMass := row14PascalPairMassLedger
  momentDegree := k10CoverBoundMomentDegreeLedger

/-- Any independent global reach-floor theorem gives bleeding-edge coverage by
the conservative packet wrapper. -/
theorem bleedingEdgeFrontierCoverage_of_reach_floor
    (hreach : forall v : Fin 13 -> Int, (forall i, v i ≠ 0) ->
      (1 : Real) / 14 <= Mreach v) :
    BleedingEdgeFrontierCoverage := by
  intro v hv
  exact ⟨bleedingEdgeFrontierPacketOfMreach v (hreach v hv)⟩

/-! ## Axiom audit hooks -/

#print axioms observerChartKind_count
#print axioms equivalenceShadow_count
#print axioms observerChartLedger_no_raw_unprotected
#print axioms equivalenceTriadLedger_no_bare_distribution
#print axioms row14PascalPairMassLedger_readout
#print axioms k10CoverBoundMomentDegreeLedger_readout
#print axioms s258DirectArcDebtLedger_readout
#print axioms bleedingEdgeFrontierPacket_mreach
#print axioms lrc14_from_bleeding_edge_frontier_coverage

end LRC14
end LonelyRunner
