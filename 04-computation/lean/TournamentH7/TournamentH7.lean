/-
  TournamentH7 — Lean 4 formalization of forbidden H-values in tournaments.

  ## Theorems formalised

  - `Tournament.H_ne_seven`         (THM-343)
  - `Tournament.H_ne_twentyone`     (HYP-1753)
  - `Tournament.H_ne_sixtythree_le_seven` (finite n ≤ 7 verification)
  - `Tournament.redei_existence`    (Rédei 1934, axiom)
  - `Tournament.redei_parity`       (Rédei 1934, axiom)
  - `Tournament.H_pos`              (corollary: H ≥ 1)
  - `Tournament.H_ne_two`           (corollary of parity)
  - `Tournament.H_ne_even`          (corollary of parity)
  - `Tournament.H_not_in_forbidden_pair`
  - `Tournament.H_not_in_forbidden_trio_le_seven`

  ## Module layout

    · TournamentH7.Basic        — Tournament structure.
    · TournamentH7.RootSigns    — Type-A root-sign atoms and walk telescoping.
    · TournamentH7.RootPackets  — Open-walk boundaries and closed root packets.
    · TournamentH7.NaturalOperationDigraphs — natural-number operation shadows.
    · TournamentH7.Cycles       — DirectedCycle, isOdd.
    · TournamentH7.SCC          — Reachability, IsSCC, HamPath, H(T).
    · TournamentH7.OCF          — External axioms (OCF, Moon, etc.).
    · TournamentH7.H7           — Theorem H_ne_seven (THM-343).
    · TournamentH7.Forbidden    — Chain step, binomial bound, extended OCF.
    · TournamentH7.H21          — Theorem H_ne_twentyone (HYP-1753).
    · TournamentH7.H63          — Finite theorem H_ne_sixtythree_le_seven.
    · TournamentH7.Redei        — Rédei 1934 axioms + corollaries.
    · TournamentH7.HSpectrum    — Universal pair + finite trio bundles.
    · TournamentH7.Tilings      — Base path, tile-complement T̃, score formula
                                  (project-novel; oracle-2026-05-11-S1).
    · TournamentH7.GridReflection — THM-280 (grid reflection ↔ complement;
                                  project-novel; opus-2026-04-03-S27).
    · TournamentH7.GoodCuts     — Good-cut buckets and exact spectrum.
    · TournamentH7.StaircaseConnectivity — Top good-cut bucket ↔ SC.
    · TournamentH7.BucketBalance — Abstract bucket balances and Boolean-cube masks.
    · TournamentH7.StaircaseBucketTransport — Concrete staircase transport checksums.
    · TournamentH7.ProductSum   — Product-sum defect normal form.
    · TournamentH7.LRCDeathChain — Finite LRC death-chain/live-depth quotient.
    · TournamentH7.LRCFactorialAtom — Finite factorial atom identities for Q₀.
    · TournamentH7.LRCBooleanTypeCut — finite arithmetic for the HYP-2791
                                      Boolean/type signed cut.
    · TournamentH7.LRCPeriodmaxCertificate — finite arithmetic kernel for the
                                      completed THM-563 period-max audit.
    · TournamentH7.LRCGenuineWideCorrection — finite arithmetic kernel for the
                                      HYP-2805 genuine-wide margin correction.
    · TournamentH7.LRCQ6Contraction — finite arithmetic kernel for the q6
                                      contraction endpoint-period certificate.
    · TournamentH7.LRCGk8SingleFar — finite arithmetic kernel for the HYP-2829
                                      gK8 single-far margin.
    · TournamentH7.LRCDoubletWitnessFloor — arithmetic checksum for the
                                      genuine-wide doublet rho*/witness scout.
    · TournamentH7.LRCLowerThresholdNeighborhood — arithmetic checksum for the
                                      HYP-2847 lower-threshold neighborhood
                                      width scout.
    · TournamentH7.LRCMreachConcrete — concrete Mreach compactness bridge
                                      for the LRC14 reach-to-lonely step.
    · TournamentH7.LRCApex7Floor — apex-7 denominator-14 obstruction.
    · TournamentH7.LRCUnitGrid14 — exact denominator-14 unit-grid sieve.
    · TournamentH7.LRCBindingPair — denominator-14 binding-pair arithmetic.
    · TournamentH7.LRCTournamentStateLift — formal closure endpoint for
                                      HYP-2908's tournament-state lift route.
    · TournamentH7.LRCWitnessAttainment — general `distZ`/margin witness
                                      attainment on a compact period.
    · TournamentH7.LRCWitnessAttainmentBridge — equivalence between the
                                      margin and concrete Mreach interfaces.
    · TournamentH7.LRCMaxGapPigeonhole — finite max-gap pigeonhole and the
                                      seven-gap equality boundary for `hnu1`.
    · TournamentH7.LRCDenseCovers — pointwise dense-cover inclusion
                                      behind `D(E) <= p0(E)`.
    · TournamentH7.LRCCoverBound — elementary cover-bound cores:
                                      coverSet monotonicity and small-k vanishing.
    · TournamentH7.LRCMarginalUniform — marginal one-speed sector measure
                                      atom for HYP-2840 decorrelation boxes.
    · TournamentH7.LRCCoverBoxes — assignment-box union skeleton for
                                      coverSet and the loose union bound.
    · TournamentH7.LRCArcComplexity — disjoint-cell arc-complexity support
                                      for the THM-546/HYP-2840 hp0cap route.
    · TournamentH7.LRCGoodSet — concrete measurable `GOOD` event carrier.
    · TournamentH7.LRCBonferroniMeasure — probability-measure Bonferroni
                                      inequality for `GOOD(E) ∩ G_P`.
    · TournamentH7.LRCEventMeasureBridge — generic event-to-shape measure
                                      handoffs for Bonferroni and `D <= p0`.
    · TournamentH7.LRCWitnessFloorConcrete — concrete slow-time witness-floor
                                      lower bound from `safeSet` and `coverSet`.
    · TournamentH7.LRCFourteenSkeleton — sorry-free conditional LRC14 proof DAG
                                      with named analytic obligations.
    · TournamentH7.LRCWitnessBonferroni — sorry-free Bonferroni/p0 reduction
                                      of the global witness floor to measure
                                      nodes plus exact arithmetic.
    · TournamentH7.LRCWitnessPartA — finite-Vmax error-budget glue for the
                                      LRC14 direct witness implication.
    · TournamentH7.LRCMedianCenterControl — HYP-3070 route-triple
                                      center-control interface and exact
                                      formal gap to LRC14.
    · TournamentH7.LRCL7Discrepancy — Finite integer core of the L7 discrepancy.
    · TournamentH7.LRCModularCuspLedger — S246 q-Pochhammer/modular-cusp
                                      finite-principal-part ledger and
                                      sixth-power collision face map.
    · TournamentH7.LRCFiniteAddressBranchClosure — S254 finite-address
                                      branch-packet interface joining q-cusp,
                                      Hurwitz, protected-branch, moment-dual,
                                      median-center, and terminal exits.
    · TournamentH7.LRCObserverGluingLedger — S259 observer-chart gluing
                                      ledger joining direct-arc, normalized
                                      arc, pair-scissors, CRT/Farey, and
                                      finite-address terminal certificates.
    · TournamentH7.LRCBleedingEdgeFrontier — observer-gluing frontier wrapper
                                      joining equivalence-triad, Pascal
                                      pair-mass, polynomial-route, and
                                      moment-degree sidecars to the
                                      finite-address packet.
    · TournamentH7.LRCProofFrontier — S259/HYP-3107 proof-frontier ledger
                                      joining solved gates, cap arithmetic,
                                      observer-gluing, and residual packet
                                      producer obligations.
    · TournamentH7.LRCColoredGateFormalization — HYP-3473 colored-gate
                                      interface for HYP-3471's low-rank
                                      E/branch producer and terminal packets.
    · TournamentH7.LRCHardOrbitCurrentJoin — HYP-3479 finite ledger joining
                                      hard mirror-orbit debt to
                                      boundary-current dispatch.
    · TournamentH7.LRCSingletonCurrentLedger — HYP-3480 finite ledger for
                                      zero-edge singleton-current dispatch.
    · TournamentH7.LRCPrivateLabelFirewall — HYP-3490 finite ledger for
                                      private-label firewall dispatch.
    · TournamentH7.Verify       — Axiom audit (#print axioms).
-/

import TournamentH7.Basic
import TournamentH7.RootSigns
import TournamentH7.RootPackets
import TournamentH7.NaturalOperationDigraphs
import TournamentH7.Cycles
import TournamentH7.SCC
import TournamentH7.OCF
import TournamentH7.H7
import TournamentH7.Forbidden
import TournamentH7.H21
import TournamentH7.H63
import TournamentH7.Redei
import TournamentH7.HSpectrum
import TournamentH7.Tilings
import TournamentH7.GridReflection
import TournamentH7.StaircaseModel
import TournamentH7.SelfComplementary
import TournamentH7.AntiAutomorphism
import TournamentH7.HPIPIdentity
import TournamentH7.Iso
import TournamentH7.IsoProperties
import TournamentH7.IsomorphismClasses
import TournamentH7.SCCounts
import TournamentH7.SCFraction
import TournamentH7.SmallTournaments
import TournamentH7.ForbiddenHCounting
import TournamentH7.StaircaseTileModel
import TournamentH7.GoodCuts
import TournamentH7.StaircaseConnectivity
import TournamentH7.BucketBalance
import TournamentH7.Examples
import TournamentH7.PaleyAxiomatic
import TournamentH7.HSpectrumExtended
import TournamentH7.ExamplesExtended
import TournamentH7.ApexBridge
import TournamentH7.StaircaseBucketTransport
import TournamentH7.RedeiFromOCF
import TournamentH7.HSpectrumClean
import TournamentH7.HSpectrumSmallN
import TournamentH7.SmallHEnumerations
import TournamentH7.BasePathSink
import TournamentH7.IsoCharacterizations
import TournamentH7.ScoreSequence
import TournamentH7.Paley3
import TournamentH7.TransitiveH
import TournamentH7.ProductSum
import TournamentH7.OpSymmetry
import TournamentH7.LonelyRunner
import TournamentH7.LonelyRunnerMathlib
import TournamentH7.LRCGapSurplusLedger
import TournamentH7.LRCWitnessCert
import TournamentH7.LRCCertTable
import TournamentH7.LRCLadderPack
import TournamentH7.LRCSevenTranslate
import TournamentH7.LRC14CertRoute
import TournamentH7.LRC14Dispatch
import TournamentH7.LRCFarElementRate
import TournamentH7.LRCRatWitness
import TournamentH7.LRCRegionDiff
import TournamentH7.LRCGoodPipeline
import TournamentH7.LRCPeelAssembly
import TournamentH7.LRC14WindowWiring
import TournamentH7.LRC14ConcreteSurface
import TournamentH7.LRCWindowData
import TournamentH7.LRCTopRatioPeel
import TournamentH7.WindowData20
import TournamentH7.WindowDispatch
import TournamentH7.LRC14CoveringFarSurface
import TournamentH7.LRCWindowData
import TournamentH7.LRCTopRatioPeel22
import TournamentH7.LRCIntervalTransport
import TournamentH7.RatIntervals
import TournamentH7.LRCDeathChain
import TournamentH7.LRCFactorialAtom
import TournamentH7.LRCBooleanTypeCut
import TournamentH7.LRCPeriodmaxCertificate
import TournamentH7.LRCGenuineWideCorrection
import TournamentH7.LRCQ6Contraction
import TournamentH7.LRCGk8SingleFar
import TournamentH7.LRCDoubletWitnessFloor
import TournamentH7.LRCLowerThresholdNeighborhood
import TournamentH7.LRCMreachConcrete
import TournamentH7.LRCApex7Floor
import TournamentH7.LRCUnitGrid14
import TournamentH7.LRCBindingPair
import TournamentH7.LRCTournamentStateLift
import TournamentH7.LRCGapReach
import TournamentH7.LRCWitnessAttainment
import TournamentH7.LRCWitnessAttainmentBridge
import TournamentH7.LRCMaxGapPigeonhole
import TournamentH7.LRCDenseCovers
import TournamentH7.LRCCoverBound
import TournamentH7.LRCMarginalUniform
import TournamentH7.LRCCoverBoxes
import TournamentH7.LRCArcComplexity
import TournamentH7.LRCGoodSet
import TournamentH7.LRCBonferroniMeasure
import TournamentH7.LRCEventMeasureBridge
import TournamentH7.LRCWitnessFloorConcrete
import TournamentH7.LRCFourteenSkeleton
import TournamentH7.LRCWitnessBonferroni
import TournamentH7.LRCWitnessPartA
import TournamentH7.LRCMedianCenterControl
import TournamentH7.LRCL7Discrepancy
import TournamentH7.LRCApexShell
import TournamentH7.LRCModularCuspLedger
import TournamentH7.LRCFiniteAddressBranchClosure
import TournamentH7.LRCObserverGluingLedger
import TournamentH7.LRCBleedingEdgeFrontier
import TournamentH7.LRCProofFrontier
import TournamentH7.LRCColoredGateFormalization
import TournamentH7.LRCHardOrbitCurrentJoin
import TournamentH7.LRCSingletonCurrentLedger
import TournamentH7.LRCPrivateLabelFirewall
import TournamentH7.LRCCrystallographicThetaFrontier
import TournamentH7.Verify
import TournamentH7.LRCUnitResidue
import TournamentH7.PolygonMirskyNewman
import TournamentH7.DangerousPatterns
import TournamentH7.BonferroniTruncation
import TournamentH7.CombPatterns
import TournamentH7.LRC14Assembly
import TournamentH7.LRCCommensuration
import TournamentH7.RatIntervalsWrap
import TournamentH7.LRCWitnessWindow
import TournamentH7.LRCLadderPerLevel
import TournamentH7.CommensurationQ
import TournamentH7.OriginNestQ
import TournamentH7.LRC14CompletenessSurface
import TournamentH7.FarElementRate
import TournamentH7.RateLemma
-- import TournamentH7.LRCWindowPack1  -- WIP: 0xC0000005 crash at file scale (both decide flavors); see HYP-3916 forensics
import TournamentH7.LRCKernelGate
import TournamentH7.LRC14AxiomAudit  -- klein-S113: #print axioms footprint of the LRC(14) endgame surface
