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
import TournamentH7.LRCChainPeel
import TournamentH7.LRCPairBlock
import TournamentH7.LRCBlockSix
import TournamentH7.LRCFatBlockChain
import TournamentH7.WindowData20
import TournamentH7.WindowDispatch
import TournamentH7.LRC14CoveringFarSurface
import TournamentH7.LRCWindowData22
import TournamentH7.LRCRealRegions
import TournamentH7.LRCIntervalTransport
import TournamentH7.LRCLadderFattening
import TournamentH7.LRCAliasingBound
import TournamentH7.LRCPLFourier
import TournamentH7.LRCGridPort
import TournamentH7.LRCP2Eval
import TournamentH7.LRCChainAnchor
import TournamentH7.LRCPairUpgrade
import TournamentH7.LRCChainMeasures
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
import TournamentH7.JointRateCore
import TournamentH7.LRCFarPeelCore  -- opus-S49 (HYP-4042): far-peel measure core (kernel-pure) on the landed region rate lemma
import TournamentH7.LRCScaleSeparation  -- opus-S50 (THM-608): scale-separation/cluster-absorption (kernel-pure) -- the renormalization core, large-magnitude side
import TournamentH7.LRCSimulPeel
import TournamentH7.LRCTrapArea
import TournamentH7.LRCSpreadPairFloor  -- klein-S118 (kps-S25 mathlib-v4.30.0 repair): spread pair-floor Stages 1-3
import TournamentH7.LRCStarSafe  -- kps-S26: star-safe positivity capstone (real danger sets, c<=7 err-free)
import TournamentH7.LRCSpread13  -- kps-S28: sharp bounded-ratio window (max <= 13 min => lonely, ratio 7->13)
import TournamentH7.LRCHlargeRoute  -- opus-S54: hlarge case-routing skeleton (route 1 = ratio<=13 => spread13; reduce to the gap obligation), kernel-pure
import TournamentH7.LRCDenominatorRoute  -- kps-S29: LRC(14) as a bounded-denominator search (single-obligation framing)
import TournamentH7.LRCFarPeelGood  -- kps-S30: far-peel positivity on goodRegion2 (opus steps 2+3 fused, no comb bridge)
import TournamentH7.LRCFarPeelDeepWell
import TournamentH7.LRCResidueLiar  -- kps-S3: residue-liar family {1..11,13,12K} lonely at (5K+2)/(12K+5), kernel-pure infinite cert (HYP-4078)
import TournamentH7.LRCOneSwapLadders  -- kps-S5: deep one-swap covering ladders drop-9/11/13 (the covering-min quartet), kernel-pure residue-table certs (HYP-4082)
import TournamentH7.LRCTwoKillerLadder  -- klein-S130: first MULTI-killer CoveringFarLonely class {1..11,14,156K} lonely at 13K/(156K+1), kernel-pure (HYP-4092)
import TournamentH7.LRCDominantPeel  -- kps-S5: SHARP dominant peel -- hdom closed at linear 13x threshold (vs quadratic far_peel), LRC13+Lipschitz+covering (HYP-4086)
import TournamentH7.LRCEndgameAssembly
import TournamentH7.LRCGridValue  -- klein-S138: the VALUE form (marginQ evaluator + margin_le_of_grid + AP tightness M({1..12})=1/13 EXACT) (HYP-4121)
import TournamentH7.LRCMergeGridAttainment  -- klein-S137: THM-592 merge-grid attainment (M attained at m/(vi+vj)) (HYP-4114)
import TournamentH7.LRCTeethR  -- klein-S136: radius-parametric window stack + level-12 tower step (HYP-4107)
import TournamentH7.LRCMultiKillerWindow  -- klein-S134: multi-killer window (<=6 tops via zero-credit Hunter) + cite_margin (HYP-4099)
import TournamentH7.LRCHcompSurface  -- klein-S133: hcomp <= TightLooseDichotomy + CornerLonely; LRC(14) <= citation + the two (HYP-4096)
import TournamentH7.LRCTightAPFreeRider  -- klein-S132: CRT tight-AP free-rider (extremal half of hcomp, mac-mini S47 handoff)
import TournamentH7.LRCOneWindowPeel  -- klein-S132: one-window peel (killer threshold B/(14(beta-1/14))) + S131 primitivity split (HYP-4095)  -- kps-S6: LRC(14) <= LRC(13) + compressed (dominant peel discharges hdom); pins the open leaf (HYP-4091)
import TournamentH7.LRCBandMargin  -- kps-S1(2026-07-05): reciprocal-window margin a/(a+b) + spread12_lonely13; tight-locus rigidity's 1st necessary condition tight=>w_max>=12 w_min (HYP-4096)
-- import TournamentH7.LRCWindowPack1  -- WIP: 0xC0000005 crash at file scale (both decide flavors); see HYP-3916 forensics
import TournamentH7.LRCKernelGate
import TournamentH7.LRC14AxiomAudit  -- klein-S113: #print axioms footprint of the LRC(14) endgame surface
import TournamentH7.LRCTopRatioPeel13  -- klein-S114: 13-ratio peel (sharpens kps 91)
import TournamentH7.LRCFarCutSplit  -- klein-S115: far-count-7 dispatch (integrates the endgame residuals)
import TournamentH7.LRCHunterLedger  -- klein-S116: path-Hunter (Bonferroni) inequality, the 7-wall crossing heart
import TournamentH7.LRCLedgerAssembly  -- klein-S117: ledger-positivity assembly (singles + pair-floor -> hledger, c<=7)
import TournamentH7.LRCDeepWellWitness  -- klein-S119: covering-min extremizer's witness, general n (deep well lonely at t*=n/Phi6(n) via zeta6 rotation)
import TournamentH7.LRCDeepWellLonely  -- klein-S119b: wires the witness to Lonely 14 (deepWell14 {1..12,182} Lonely at 14/183)
import TournamentH7.LRCDeepWellReach  -- kps-2026-07-11-S127 (cont.60): the single-killer balance bound in the reach/margin sSup API -- reach_ge_of_covering13 (Fin 13 covering-reach atom) + deepWell_reach (reach({1..12,182}) >= 14/183, the covering-min VALUE, keeping the full margin where klein's Lonely keeps only 1/14) + deepWell_reach_gt_tight (>1/14, minimizer is loose). kernel-pure [propext,Classical.choice,Quot.sound].
import TournamentH7.LRCSingleKillerLadder  -- kps-2026-07-11-S127 (cont.61): the GENERAL single-killer ladder S_c={1..12,182c}, c>=1 -- singleKiller_cover (parametrized clearing cert at q=182c+1, rotation 14c, band [14c,168c+1], killer 182c==-1 mod q) + singleKiller_reach (reach(S_c) >= 14c/(182c+1)) + singleKiller_reach_ge_floor (reach(S_c) >= 14/183 for ALL c>=1, closing the entire single-killer class in the reach API, tight only at deep well c=1). kernel-pure [propext,Classical.choice,Quot.sound].
import TournamentH7.LRCMarginMeasure  -- klein-S120: THM-613 margin->measure bridge formalized (meas lonely >= 2(M-1/14)/vmax; deep well >= 13/233142)
import TournamentH7.LRCEvenDescent  -- klein-S121: THM-612 tower step formalized (even part of a lonely family descends to U lonely at m*t)
import TournamentH7.LRCFolding  -- klein-S123: THM-615 folding-identity engine (reach half-shift for odd, even (+1/2)-periodicity, lattice fold, Psi/extremity)
import TournamentH7.LRCLargeTightener  -- klein-S125: THM-615 Lemma 3 discrepancy core (large multiplier hits reach=1/6 => folded value >= 1/12; the loose-end confinement half)
import TournamentH7.LRCSevenGap  -- kps-S24: seven-gap deficit + cluster sweep step (the CLUSTERED-block lane)

import TournamentH7.LRCDilation  -- mac-mini-S24: WLOG gcd=1 (dilation), closes HYP-4043 far-peel gap

import TournamentH7.LRCBaseFloor  -- mac-mini-S25: step-1 real->rational strict-good bridge (THM-609)
import TournamentH7.IntervalEscape
import TournamentH7.LRCResiduePinning
import TournamentH7.LRCWindowMargin13
import TournamentH7.LRCKernelGate13
import TournamentH7.LRCLiftRigidityRows
import TournamentH7.LRCHarmonicGate  -- kps-S2(2026-07-05): rational-point margin certificate (t=k/s, margin mu/s, decidable integer hypothesis) + ratio-gate reduction of TightLooseDichotomy to spread bases; 2nd harmonic closes {1..11,24}; at beta=1/13 the gate subsumes the tight branch (HYP-4102)
import TournamentH7.LRCMidrangeWitness  -- klein-S141(2026-07-05): the MIDRANGE WITNESS (general) -- at t=1/(vmin+vmax) every speed clears by vmin/(vmin+vmax); compressed corollary vmax<=13vmin => >=1/14; the >=-half (THM-526/opus-S85 combs), kernel-pure (HYP-4161)
import TournamentH7.LRCLiftFloorRows  -- mac-mini-S53(2026-07-05): the multi-lift floor certificate rows through the atom -- block46 (the {4,6}->{17,19} floor family, margin >= 2/25 at t=6/25) + the six 14r-ladder rows (margins >= 14/(13(r+1)), r=7..12; r=12 = the n=13 deep well 14/169); all decide-discharged (HYP-4109, certificates from S52/HYP-4103)
import TournamentH7.LRCLiftFloorAssembly  -- mac-mini-S54b(2026-07-05): the CORRECTED lift-floor package -- HYP-4177's 14/169 domination REFUTED (block lift {4,6}->{17,19} at 2/25, unique height-1 double below); ladder_law (six rows quantified into the parametric law), ladder_law_floor (uniform >= 14/169), lift_floor_beta_ladder (>= 2/25 = the dichotomy beta, the level the assembly consumes) (HYP-4212)
import TournamentH7.LRCKStratification  -- mac-mini-2026-07-06-S1: the k-stratification core for the 3/38 cell (dInt_scale quotient descent + binder_dvd + grid_div_38 + pair_sum_dvd k=1 anchor); composes with kps's cluster-gcd ladder for the k-reduction (HYP-4232)
import TournamentH7.LRCTorusRate  -- mac-mini-2026-07-06-S10: the elementary lift-convergence rate (Lipschitz + net) making the (A)-subsumption preprint-free -- lipschitz_ge_at_near/two_point_rate/exists_net_ge => M(v^N) >= M(U) - L/(2N), no Jain-Kravitz (HYP-4342)
import TournamentH7.LRCUniformCell  -- mac-mini-2026-07-06-S3: the UNIFORM CELL LEMMA core -- the stratification apparatus parametric over every Farey cell (c,q): dInt_scale_cell, binder_dvd_cell, grid_div_cell, pair_sum_dvd_cell + NEW binder_unit_cell (quotient binders are units mod q) + witness_determined_cell (coprime cancellation) (HYP-4252)
import TournamentH7.LRCCertCompleteness  -- kps-S3(2026-07-05): certificate completeness (margin => atom cert, modulus bound s <= B/(2(beta'-beta))+1) + gap filters (not_loose covering/near-unit/mod-13,14 pinning INTO the gap, +-pair covering) (HYP-4105)
import TournamentH7.LRCGridAttainment  -- kps-S4(2026-07-05): THM-592 grid attainment (margin max at m/(|vi|+|vj|), perturbation proof) + grid_margin_domination + loose_branch_cert_2B (modulus <= 2B, zero slack) (HYP-4108)
import TournamentH7.LRCMergeExclusion  -- kps-S5(2026-07-05): the merge exclusion formal (margin value = k/s on the grid; gap (1/13,2/25) forces |vi|+|vj| >= 38 via omega incl. the k=2 parity kill) (HYP-4110)
import TournamentH7.LRCPeelCompression  -- kps-S6(2026-07-05): the 24B top-compression of gap violators (citation floor 1/12 beats gap ceiling 2/25 by 1/300; real-dialect interval escape + peel; gap_compressed_24) (HYP-4112)
import TournamentH7.LRCGapLadder  -- kps-S7(2026-07-05): the l>=2 peel ladder at rho=2/25 on klein's S136 stack (gap_tower_step + gap_ladder_rung: order-statistic compression C_l = 150l(13-l)/((2l-1)(25-4l)), l<=6) (HYP-4115)
import TournamentH7.LRCPairWalk  -- kps-S8(2026-07-05): the pair walk (3-step tooth-boundary walk, linarith kill) => gap_pair_rung: every pair of a gap violator has min <= 22B -- 3x sharper than the S7 density C_2 (HYP-4117)
import TournamentH7.LRCTripleWalk  -- kps-S9(2026-07-05): the l=3 walk by recursion into the pair walk (revisit gap covered by the other two = walk_core; balance 4max<=7min saturates at 21/25 exactly) => gap_triple_rung: balanced triples min <= 12B (HYP-4118)
import TournamentH7.LRCTemplateSurface  -- kps-S10(2026-07-05): the Q50 crux pinned formally -- TemplateWitness (decidable) + TemplateDichotomy + lrc14_of_template_and_corner: LRC(14) <= cite + TemplateDichotomy + CornerLonely; loose side = fixed finite template family (HYP-4127)
import TournamentH7.LRCMultiFoldRows  -- kps-S17(2026-07-05): consecutive multi-fold tower floor rows D_2..D_6 >= 2/25 (corrected: the closed-form law is FALSE at l=4,5 -- witnesses 17/155, 19/155; the floor survives) -- the kps side of HYP-4212 (HYP-4217)
import TournamentH7.LRCClusterGcd  -- kps-S18(2026-07-05): THE CLUSTER-GCD LADDER kernel-pure (gap_gcd_rung: (25-8|S|)*gcd(complement) <= 25(Sum_S+|S|) for |S|<=3; citation + 1/d-periodicity + tooth-visit pigeonhole; the absolute-height mechanism for gcd-clusters) (HYP-4227)
import TournamentH7.LRCClusterGcdSharp  -- kps-S19(2026-07-06): THE SHARP VISIT COUNT ((4/25)D+3w via block split + coprime permutation + wrapped-arc two-piece clip) => gap_gcd_rung_sharp: (25-4|S|)*gcd(complement) <= 75*Sum_S for |S|<=6 -- the FULL ladder with the pole at 25/4; critical-path for mac-mini's HYP-4232 k-reduction (HYP-4237)
import TournamentH7.LRCTorusSplit  -- kps-S20(2026-07-06): THE COUPLED-TORUS SPLIT RUNG, rho-parametric (torus_split_rung: rho-covered coupled 2-torus systems need 2*rho*(#lifted) >= 1 via base citation + rho-parametric sharp count) => torus_clear_gap: every coupled 2-torus (lift-limit) system with <= 6 lifted has a 2/25-clear point -- the (A) window of mac-mini's HYP-4262 J-K reduction is EMPTY of <= 6-lifted coupled values; + the 3/38 density-wall dialect (HYP-4247)
import TournamentH7.LRCCircleCover  -- kps-S20d(2026-07-06): THE (A)-WINDOW REDUCTION (circle_clear_of_density: 2*rho*|S| < 1 => a clear theta on the circle, grid-point form; CircleClearFloor = the NAMED covering obligation, PROVED l<=6, Newman-shaped + numerically confirmed l=7..11; torus_A_window_empty: base citation + the floor => a proper coupled 2-torus has a 2/25-clear point -- the whole (A) window reduces to distinct-freq covering-impossibility at 7<=l<=11) (HYP-4247)
import TournamentH7.LRCLadderLoose  -- kps-S21(2026-07-06): THE m/(12m+1) FAREY LADDER (ladder_family_loose: for m>=2, {1..11,12m} clears 2/25 at t=m/(12m+1) via rational_point_margin, residues i*m and 11m+1 in [m,11m+1]) -- the resonant 12|v case of the {1..11,v} slice, complements opus HYP-4356 (12-nmid-v loose); ladder rungs m=1 (1/13) m=2 (2/25) ARE the gap edges; broad census confirms gap (1/13,2/25) empty (HYP-4357)
import TournamentH7.LRCBandFloor  -- kps-S51(2026-07-06): THE GENERAL AVOID-BAND 2/25-FLOOR (loose_of_band: forall i, mu <= (v_i c)%q <= q-mu with 2q<=25mu => M>=2/25 at t=c/q) -- the single tool under EVERY covering layer (q=25 mu=2 mod-25; q<=12 mu=1 small-mod; 13<=q<=32 moat, loose_at_17/loose_at_32 instances); rational_point_margin (HYP-4667)
import TournamentH7.LRCSmallModFloor  -- kps-S44(2026-07-06): THE SMALL-MODULUS 1/q FLOOR (zero_avoid_floor / no_multiple_floor / loose_of_no_multiple_12: for q<=12 the 2q/25 band is <1 so clearing = avoid residue 0 = no multiple of q => M>=1/q>2/25) -- the q<=12 layer of the case-2 covering system, closes every family missing a small multiple; rational_point_margin at mu=1 (HYP-4597)
import TournamentH7.LRCMod25Floor  -- kps-S41(2026-07-06): THE mod-25 2/25-FLOOR CERTIFICATE (mod25_covering_floor / loose_of_mod25_covering: a rotation c off {0,+-1} mod 25 forces M>=2/25 at t=c/25, direct rational_point_margin instance s=25 mu=2) -- the lower-bound half of the Freiman crux: >=3-defect 12-families (no mult of 25) are mod-25-clearable => loose; pins (G) residual to a finite mod-25 covering fact (HYP-4567)
import TournamentH7.LRCMod25Transversal  -- mac-mini-S33(2026-07-06): the CLEARING ROTATION EXISTS from a MISSED PAIR (covering_of_misses_pair / loose_of_misses_pair) -- supplies the EXISTENCE half kps's mod25_covering_floor assumed: miss the antipodal pair {a,-a} mod 25 (no mult of 25) => c:=a^-1 clears every speed off {0,+-1} => M>=2/25. PROVES the non-transversal (miss-a-pair) branch of the crux (G); residual = full transversals (pair-blocker => AP) (HYP-4642)
import TournamentH7.LRCGapCandidate  -- kps-S24(2026-07-06): GAP CANDIDATES ARE DIVISIBILITY-RICH (gap_candidate_has_multiple: a 12-family not loose at 2/25 contains a multiple of every k in {2..12}, via opus int_far_of_not_dvd_k + 1/k > 2/25; gap_candidate_prime_powers: multiples of 5,7,8,9,11,12 forced) -- the covering-system structure of gap candidates, AP = minimal; sharpens the additive/collision side of the S23 split (HYP-4417)
import TournamentH7.LRCSingleLift  -- kps-S28(2026-07-06): THE SINGLE-LIFT RIGIDITY (ap_single_lift_loose: every single-13-lift of the AP {1..12} is loose at 2/25, via 12 rational_point_margin certificates; tightest = runner 6->19 at 2/23) -- the formal BASE CASE of the density floor / near-AP fiber rigidity, the residue-pinned side of the S23 split (HYP-4447)
import TournamentH7.LRCHarmonicAP  -- kps-S30(2026-07-06): THE AP IS THE HARMONIC FAMILY (second_diff_zero_iff_ap: vanishing second differences v(i+2)-2v(i+1)+v(i)=0 <=> arithmetic progression) -- the elementary algebraic heart of "AP = the flat/extremal family" behind the density floor; the (1,-2,1) harmonic relations characterize L(AP) (opus HYP-4446), the shortest/dominant safe theta-sum terms = the additive-energy heart of safe(AP)=0 (HYP-4457)
import TournamentH7.LRCMultiKillerWindow13
import TournamentH7.LRCLiftPigeonhole
import TournamentH7.LRCGapDescent
import TournamentH7.LRCLiftRowsL7
import TournamentH7.LRCTowerLift
import TournamentH7.LRCTelescope
import TournamentH7.LRCDescentSurface
import TournamentH7.LRCClearCert
import TournamentH7.LRCClearRows
import TournamentH7.LRCClearRowsB5
import TournamentH7.LRCTwoBand
import TournamentH7.LRCRayTransport
import TournamentH7.LRCTorusProjection
import TournamentH7.LRCCircleUnshifted
import TournamentH7.LRCTorusReduction
import TournamentH7.LRCRankRigidity
import TournamentH7.LRCAPProtection
import TournamentH7.LRCDivisorProtection
import TournamentH7.LRCPinnedFloor
import TournamentH7.LRCCovererDichotomy
import TournamentH7.LRCWitnessDenominator
import TournamentH7.LRCDilationInvariance
import TournamentH7.LRCLeaveOneOut
import TournamentH7.LRCRelationLattice
import TournamentH7.LRCMinimalSumset  -- mac-mini-2026-07-06-S20: minimal-sumset bound |S+S| >= 2|S|-1 (the Freiman base of the density-floor AP-uniqueness; opus-S112 theta-sum frame) (HYP-4482)
import TournamentH7.LRCFareyGap
import TournamentH7.LRCLonelyOpen
import TournamentH7.LRCSubfamilyCap
import TournamentH7.LRCSpectrumWindow
import TournamentH7.LRCBinderInfeasible
import TournamentH7.LRCConsecutiveBlock
import TournamentH7.LRCEvenBranchWitness  -- mac-mini-2026-07-06-S29: even-branch clearance, N=12 canonical mediant family {1..10,12,33} reach >= 3/35 > 2/25 => NOT a gap member (parity, not 38=2*19 compositeness) (THM-632, HYP-4572)
import TournamentH7.LRCLadderD1  -- mac-mini-2026-07-06-S33: the d=1 ladder bound (opus-S123's remaining piece) -- canonical d=1 family {1..11,x} has reach >= 2/25 for all x!=12, two-witness proof (t=1/12 generic; t=k/(12k+1) resonant), kernel-pure (THM-633, HYP-4632)
import TournamentH7.LRCDecorrelation  -- mac-mini-2026-07-06-S38: decorrelation atom reach_decorr (M(V)>=M(K)-B/L for v_i=b_i+L*k_i) + escape_loose_of_lift_floor (r<=11 escape families loose, M>2/25 via LRC<=12+decorrelation); handles the S36 covering-escape obstruction (THM-636)
import TournamentH7.LRCDecorrelation13  -- kps-2026-07-11-S127 (cont.48): the 13-runner decorrelation atom at the DIRECT 1/14 threshold -- reach_decorr13 (reach(V)>=reach(K)-B/L, Fin 13) + escape_loose13_le12 (<=12 distinct lifts => LRC(<=13) reach(K)>=1/13, L>2366 => reach(V)>1/14) + escape_loose13_le6 (<=6 lifts => LRC(7) reach(K)>=1/7, L>182, the DC even-heavy sharp threshold, mac-mini cont.49). The large-diameter-half atom, machine-checked. kernel-pure [propext,Classical.choice,Quot.sound].
import TournamentH7.LRCReachTransport  -- kps-2026-07-11-S127 (cont.50): DILATION PRESERVES LOOSENESS -- reach_dilate_ge (reach(c*v)>=reach(v) for c>=1, via the scaled witness t0/c; no periodicity) + loose_dilate (reach v > 1/14 => reach(c*v) > 1/14). The formal underpinning of MISTAKE-140 (boxeph-S19: "min M grows with diameter" is a sampling artifact; the DC class stratifies by STRUCTURE not diameter) -- machine-checks the arrow reducing the unbounded near-dilate slice to bounded loose cores. kernel-pure [propext,Classical.choice,Quot.sound].
import TournamentH7.LRCCoveringReach  -- mac-mini-2026-07-06-S34: covering-reach atom (rational_point_margin -> reach>=mu/q) + d=2 generic (11 nmid x,y => reach>=1/11>2/25); the uniform Lean shape of the q<=39 covering system (HYP-4587)
import TournamentH7.LRCRoute2Assembly  -- opus-2026-07-06-S129: (A)<=(C) wired end-to-end -- torus_loose_of_rank2 composes GREEN projection floor + GREEN rigidity wrapper; rank-2 torus is loose given the C-bridge; + JKReduction obligation + lrc14_via_route2 top-level conditional (HYP-4652)
import TournamentH7.LRCCoveringStrata  -- kps-S51(2026-07-06): DISCHARGING STRATA of opus's CoveringComplete -- smallmod_hasWitness/no_multiple_hasWitness (a family missing a multiple of some q<=12 has HasCoveringWitness (q,1,1)); produces the covering witness for the Fan-Sun gcd / small-mod layer (klein-S147), feeding opus-S129's (C)<=CoveringComplete (HYP-4677)
import TournamentH7.LRCCoveringDisjunction  -- opus-2026-07-06-S129: branch-2 endpoint -- loose_of_covering_set packages mac-mini's per-q reach_ge_of_covering atom into the OR-over-q disjunction the residue check emits (with mu/q>=2/25 folded in); reduces branch 2 to a purely arithmetic witness-existence (HYP-4652)
import TournamentH7.LRCTenSwapLadder  -- kps-S54(2026-07-07): the k/(10k+7) one-swap ladder {1..9,11,12,13,10k}=AP[10->10k] lonely (tenSwap_lonely, first rung 2/27 at k=2) + gw_lonely (Goddyn-Wong {1..11,13,24}=AP[12->24] tight M=1/14, the SECOND tight family besides AP -- corrects S53 MISTAKE-100 recurrence); kernel-pure, closed-form near-tight census families (HYP-4717)
import TournamentH7.LRCSaturatedReduction  -- kps-S56(2026-07-07): the SIEVE reduction (leg 1 of the sieve+coarse+decorrelation decomposition) -- lrc14_iff_saturated_lonely: LRC(14) <=> every SATURATED 13-family (mult of every q in {2..14}) is lonely, via counterexample_needs_all_divisors; kernel-pure (HYP-4737)
import TournamentH7.LRCCoarseReduction  -- kps-S53(2026-07-06): THE COARSE/SCALE REDUCTION as a Route-1 tool (formalizes kps-S52 M(v)>=M(K)-A/L) -- reach_transfer_coarse (k lonely at s0 margin mu, |a_i|<=A => v_i=a_i+L*k_i lonely at s0/L margin>=mu-A/L) + lonely14_of_coarse_le12 (v clusters into <=12 groups at scale L, A/L<=1/13-1/14 => Lonely 14 via SETTLED LRC<=13); grounds the MULTI-SCALE branch of Route 1, leaving the density node only single-scale families; kernel-pure (HYP-4707)
import TournamentH7.LRCTailDiameter
import TournamentH7.LRCGoodDilation  -- mac-mini-S42(2026-07-07): THE TAIL-DIAMETER THEOREM (monad-S2 HYP-4817, kps-S59 monotonicity) -- good_anti (subset antitonicity of the empty-arc good set) + good_translate (rotation invariance) + muGood_ge_AP76 (diam<=75 chain) + hlarge_floor_of_diam_le (mu_{1/7}(E)>=m_P for every family of diameter<=75, conditional ONLY on the AP76 Farey-cell ledger certificate); measure enters only via measure_mono
import TournamentH7.LRCArcCount  -- opus-S169: the ARC-COUNT PIGEONHOLE for good-period existence -- good_period_of_arccount (#arcs < rho*.V => a grid point j/V is Good = a good period; the finite-Vmax glue of the covering capstone, dissociated branch). exists_gridpoint_Ico + exists_long_Ico elementary/unconditional; rho*>=D3 (density floor) and #arcs<=c.spread (bounded-arc-count) enter as the two inputs.
import TournamentH7.LRCRieszCertificate  -- opus-S173: soundness of the RIESZ-PRODUCT loneliness certificate (singular-series/lonely-measure route, THM-515/HYP-2540) -- riesz_certificate (R>=0 & int M*R < int R => M<1 on positive measure => S loose) + no_certificate_of_ae_covered (validity: tight S => no certificate). The POSITIVE-DEFINITE route that sidesteps the covering-W Mertens wall (opus-S172); analytic core (a uniform R => inf L>0) stays open.
import TournamentH7.LRCAP20Certificate  -- death-star-S2: AP20 density-floor certificate PROVEN (unconditional, 4.5x margin); makes diameter<=19 (census) floor unconditional
import TournamentH7.LRCGoodSetBridge  -- death-star-S3: Good ⊆ goodSet bridge -- transfers the muGood diameter floors to the concrete goodSet event that witnessG2 measures (goodSet_measure_ge_mP: bounded-diameter concrete-event floor >= m_P)
import TournamentH7.LRCRandom031ProofPackets  -- death-star-S3 housekeeping: wire in kps-2026-07-01-S28 module (committed but manifest-wiring never committed; builds clean, sorry-free, flagged by kps-S64/S65)
import TournamentH7.LRCFareyRoof  -- opus-S135(2026-07-07): THE FAREY ROOF pointwise theorems (THM-637) -- lemmaA/lemmaB (first return from above/below via the divisibility contradiction q(ip'-aq')=mq'+i), lemmaC (gap bound via three exhibited indices), zero_gap_empty + fract_q_mul: maxgap(AP_k,x) = gap-at-0 = q(x-p/q)+q'(p'/q'-x) on each open Farey-k cell; the pointwise engine of the diameter floor (kps-S59 HYP-4797 subset lemma + monad-S2 HYP-4817 exact crossing n=76/77); kernel-pure, no sorry/native_decide (HYP-4782)
import TournamentH7.LRCFareyRoofBridge  -- boxeph-2026-07-07: ROOF->GOOD BRIDGE -- good_of_roof_gt (on an open Farey-k cell, roof>theta => x in Good theta (AP_k), witness arc a=(q'x-p')+(roof-theta)/2 sits strictly inside the empty roof-interval, contradiction via zero_gap_empty) + roof_superlevel_subset_good + muGood_ge_of_cell_interval + muGood_ge_sum_intervals (measure atoms: roof-superlevel intervals lower-bound muGood); the reusable pointwise link that reduces AP76Certificate to a PURE real-superlevel (Farey-sum) measure computation; kernel-pure, no sorry/native_decide
import TournamentH7.LRCAP30Floor  -- boxeph-2026-07-07: AP30 CERTIFICATE (unconditional) -- mu_{1/7}(AP30)>=m_P via the Farey roof (bridge) + the two endpoint roof-superlevel intervals (0,6/203),(197/203,1); DERIVES death-star's hand-built AP20 arcs from opus's roof and extends the unconditional diameter floor from <=19 to <=29; ap30_certificate + ap30_certificate_icc0 (feeds TailDiameter.muGood_ge_APD); kernel-pure
import TournamentH7.LRCAP44Floor  -- boxeph-2026-07-07: AP44 CERTIFICATE (unconditional, MULTI-NODE) -- mu_{1/7}(AP44)>=m_P via the Farey roof (bridge) + 4 roof intervals across nodes 0,1/2,1; extends the unconditional diameter floor to <=43 (beyond the 2-interval endpoint method D<=29/30); ap44_certificate + ap44_certificate_icc0; the 1/2 node demonstrates the multi-node aggregation that scales to the tight AP76; kernel-pure
import TournamentH7.LRCAP76Certificate  -- opus-2026-07-07-S145: THE TIGHT AP76 CERTIFICATE (unconditional, kernel-pure) -- muGood(1/7)(AP76) >= 2314528732/40290957525 (the exact value) via boxeph's roof bridge + all 24 q<=6 Farey-76 roof-superlevel intervals (executes boxeph-S2's handoff); DISCHARGES TailDiameter.AP76Certificate: the k=13 diameter<=75 hlarge leg is now UNCONDITIONAL (hlarge_floor_diam75_unconditional); axioms [propext, Classical.choice, Quot.sound]
import TournamentH7.LRCComposedRealization  -- death-star-2026-07-09-S2 (HYP-5750, boxeph H1 handoff): the LEAN CORE of LEM-014 (P-separated composed realization) -- nearInt_ge_of_close (slow-leg 1-Lipschitz transport = G_P^eps erosion pointwise) + nearInt_clear_of_driftGap_single (cluster-leg per-runner clearance factored from klein-S205's all-13 fold) + minReach/Mreach_ge_of_composed_realization (per-runner [cluster-bound at grid j] OR [eroded slow safety at x, |tau-x|<=Delta] => lonely at tau=(j+a+g/2)/Vmax); the ratio>13 multi-scale instrument (P=empty needs only LRCPureClusterCorner); measure-side inputs (robust floor x*) = LEM-014/kps-S113/THM-669-670 hypotheses; S3 adds round_time_close + Mreach_ge_of_composed_realization_round (j=round(Vmax*x), Delta=3/(2Vmax) pre-instantiated -- plug-and-play for the robust-floor consumers); kernel-pure
import TournamentH7.LRCFourierCompletion  -- death-star-2026-07-09-S13 (HYP-5890): LEM-022's FOURIER COMPLETION STAGE A, kernel-pure -- two_mul_le_sin_pi_mul (Jordan packaged: 2t <= sin(pi t) on [0,1/2]) + norm_exp_I_sub_one (the factoring ||exp(i theta) - 1|| = 2|sin(theta/2)| via E*(E - E'), E E' = 1) + norm_expSum_le (THE INTERVAL EXPONENTIAL-SUM BOUND: ||Sum_{r<len} exp(2 pi i h(lo+r)/q)|| <= q/(2d) under the sine witness 2(d/q) <= |sin(pi h/q)|; geometric sum + ||z^len - 1|| <= 2 + the 4d/q denominator floor); Stage B remaining = self-derived orthogonality + the completion identity + assembly with the S10 box count and S12 dyadic assembly => the fully formal |C_w - b^2/q| <= 5 q log2^2(q+..)/P(w); the C wall is breached
import TournamentH7.LRCHyperbolaBox  -- death-star-2026-07-09-S10 (HYP-5870): LEM-022's COMBINATORIAL HEART in Lean -- the circle-metric API on ZMod q (cdist, subadditivity via val_add + omega, negation-invariance, signed representative sgnRep with |sgnRep| = cdist and cdist <= natAbs of any Z-rep) + card_mulsep_in_Icc (the MULTIPLICATIVE-separation pigeonhole: pairwise S <= |x-y|*D in [a,b] => (|T|-1)*S <= (b-a)*D; division only inside the fiber map, third instance of the clamped-fiber pattern) + hyperbola_box_count (THE SEPARATION COUNT: forall h != 0, P <= cdist h * cdist(wh) => (#{k != 0 : cdist k <= K, cdist(wk) <= M} - 1)*P <= 4KM, division-free); S12 adds harmonic_ratio_sum_mul_le (THE DYADIC ASSEMBLY, LEM-022 Step 3, kernel-pure, constant 20 < paper's 24 via the per-fiber dichotomy: empty fiber => 0, nonempty => witness forces P < 2^{i+j+2}; Nat.log dyadic classes + sum_fiberwise_of_maps_to; ZERO complex analysis); LEM-022 Lean surface now = ONLY the Fourier completion (AddChar/ZMod + sine bound); kernel-pure
import TournamentH7.LRCLem012NearAP  -- death-star-2026-07-09-S7 (HYP-5810): LEM-012 IN LEAN, the coherent-branch closer -- lem012_nearAP_free_gap (kernel-pure): a cluster = AP(a,d,L) + <=5 strays with 7L+mQ+m < 6Q+13 (= 7(L-1) < (6-m)(Q+1)) and Q < V (MISTAKE-129) has a period j in [1,Q] leaving an open arc g > 1/7 avoiding EVERY translate of EVERY tooth e*j/V; mechanism: Mathlib Dirichlet (kps-S111 pattern) + one-sided AP-orbit interval + complement window (translate_notMem_window, translate_unique_in_window) + Fin(m+1) pigeonhole (exists_free_piece); output = the hfree shape of klein-S205 driftGap / LRCComposedRealization + kps-S99 dispatch's near-AP branch input; THE BRANCH THAT KILLS EVERY MODULAR CERTIFICATE (THM-675: B5 = 0/260 on the 40->41 near-dilation) is now formalized closed; no fract, no circular sort
import TournamentH7.LRCDiscreteBonferroni  -- death-star-2026-07-09-S5 (HYP-5780, klein-S210 handoff): THM-671'S HISTOGRAM BOUND in Lean -- bandCount/liveCount/momentS/B5 decidable integer defs over multipliers p in (0,q); B5_le_liveCount (sum mac-mini-S101 BonferroniTruncation pointwise engine over multipliers + sum swap = the quintic moment form); exists_live_of_B5_pos; mreach_ge_of_B5_pos (B5>0 at ANY modulus => live p => t=p/q witness via kps-S114 pairsum dispatch => Mreach>=1/14 -- the decidable loneliness certificate, decide-shaped no analysis); DEMO: b5_ap_fourteen_pos (B5(AP{1..13},14)=5>0 by kernel decide, histogram {0:6,1:6,6:1}) + mreach_AP_via_histogram (the equality extremal certified end-to-end through the THM-671->THM-668 pipe); kernel-pure, NO native_decide; S6 (HYP-5800) adds THE DEPTH-GENERAL LADDER: bonf D (odd D) + bonf_le_liveCount + mreach_ge_of_bonf_pos (depth = analysis dial: 5 generic, 7+ near-coherent per THM-675's escalation) + bonf13_eq_liveCount (THE EXACTNESS ENDPOINT: coverage <= 13 => depth-13 truncation IS the indicator => complete decidable LM test, no Bonferroni loss) + the part-6 sockets mreach_ge_of_bonf_sum_pos (THM-673's aggregated target: sum over moduli positive => lonely) + dvd_Ioc_card_le (THM-673(A) dispersal core: nonzero relation value <= K*V has <= K divisors in (V,2V], cofactor injection)
import TournamentH7.LRCB5Periodic  -- kps-2026-07-11-S127 (cont.34): B5 is RESIDUE-PERIODIC (B5_congr_mod: B5 v q depends only on v mod q) => the bounded-window ruler search (hasWindowRuler_congr) is DIAMETER-FREE. Answers: a fixed window [8,Q] (empirically Q~43) contains a B5>0 ruler at ANY diameter. kernel-pure [propext, Classical.choice, Quot.sound].
import TournamentH7.LRCCleanRuler  -- kps-2026-07-11-S127 (cont.28, THM-707): the CLEAN-RULER reduction of the single Lean obligation hB5 -- trunc5_of_le_five (partial alternating binomial sum: coverage<=5 => depth-5 truncation = live indicator) + b5_pos_of_clean (shallow coverage + a live multiplier => B5 = liveCount > 0, no signed penalty) + exists_B5_pos_of_cleanRuler (a clean ruler supplies the per-family hB5 witness). Sharpens THM-671 B5<=liveCount to the EXACT B5 = liveCount - sum_{bandCount>=6} C(bandCount-1,5); reduces hB5 to the transparent seven-sector condition (live multiplier + no >=6-covered multiplier). kernel-pure [propext, Classical.choice, Quot.sound].
import TournamentH7.LRCPrimeCleanRuler  -- kps-2026-07-11-S127 (cont.31): the GENERAL prime clean ruler -- b5_pos_of_prime_ndvd: for ANY prime P<=13 with no speed divisible by P, q=P is a perfectly clean ruler (bandCount=0, liveCount=P-1) => 0<B5 v P. Generalizes opus-S227 THM-709 (q=13) to all primes {2,3,5,7,11,13} via kps b5_pos_of_clean. kernel-pure [propext, Classical.choice, Quot.sound].
import TournamentH7.LRCCompositeCleanRuler  -- kps-2026-07-11-S127 (cont.38): the GENERAL CLEAN q<=14 divisibility ruler (b5_pos_of_div_clean) -- for q<=14 the band = all nonzero residues so bandCount=#{i:q|v_i p}; q-nmid-all (live) AND every divisor d>=2 of q divides <=5 speeds (clean) => 0<B5 v q. Generalizes THM-712 (primes) to composites q=8,9,10,12,14 = tier-1 of the bounded-window covering (HYP-6035). kernel-pure [propext,Classical.choice,Quot.sound].
import TournamentH7.LRCThreeGapConsecutive  -- kps-2026-07-11-S127 (cont.44): the THREE-GAP consecutive-reduction clearing lemma -- bandCount_zero_of_consecutive_block: for 16<=q<=28, if the 13 residues form a consecutive block {r..r+12} in the safe band [2,q-2] (i.e. 2<=r, r+12<=q-2) then bandCount=0 (family clears); => liveCount_pos_of_consecutive_block feeds opus-S235 band-edge => M>1/14. The AP corner clears via p===+-d^{-1} (consecutive), window [16,27] (HYP-6090). kernel-pure [propext,Classical.choice,Quot.sound].
import TournamentH7.LRCPrimeRuler  -- opus-2026-07-11-S227 (THM-709): the PRIME clean ruler -- cleanRuler_of_not_dvd_13 (no speed divisible by 13 => q=13 clean, bandCount=0 everywhere, since the safe band [1,12] at q=13 = the nonzero residues). Discharges hB5 for the sub-class avoiding 13|v_i via exists_B5_pos_of_cleanRuler; hard residuals are exactly those WITH a speed =0 mod 13 (the 13-is-prime phenomenon). Elementary residue arithmetic, no decide. kernel-pure [propext, Classical.choice, Quot.sound].
import TournamentH7.LRCHighTailIdentity  -- klein-2026-07-09-S219: THM-676 IN LEAN -- the HIGH-TAIL IDENTITY bonf D = liveCount - penalty D (odd D, EXACT: the entire Bonferroni error is the high-coverage mass, penalty D = sum over multipliers of (bandCount-1).choose D, Nat.choose auto-killing coverage <= D = the discrete apex-7 form); engine = partial_alternating_choose (Pascal telescope: sum_{d<=D}(-1)^d C(c,d) = (-1)^D C(c-1,D)) + odd_truncation_closed_form (uniform per-p form riding N-truncated subtraction); THE DEPTH LADDER bonf_le_bonf_next (odd 5<=D<=11: B5<=B7<=B9<=B11<=B13=LM formal end-to-end with bonf13_eq_liveCount, via choose_ladder_dom kernel-decide on the 13-runner domain -- THM-675's escalation); penalty_lt_liveCount_of_bonf_pos (the certificate reading). Kernel-pure, no native_decide
import TournamentH7.LRCAggregatedDispersal  -- klein-2026-07-09-S220: THM-673(A) IN LEAN -- aggregated_dispersal (the V-INDEPENDENT dispersal: sum over q in (V,2V] of resonanceMass <= K * total weight, via dvd_Ioc_card_le cofactor injection per relation value -- the m!=0 half of the aggregated modular supply, klein-S211's one-line proof formal) + exists_thin_modulus (THE THIN-MODULUS PIGEONHOLE: some window modulus carries at most the V-free average -- the a-priori existence half of THM-671 part 6(i)) + mreach_ge_of_thin_transfer (the composition socket: thin modulus + the off-line/off-peak TRANSFER hypothesis => Mreach >= 1/14 via mreach_ge_of_bonf_pos; the transfer = THE named remaining mathematical content, THM-677 Add.3-4 / THM-680/681). Kernel-pure, no native_decide
import TournamentH7.LRCParityPairing  -- klein-2026-07-09-S227: THE PARITY PAIRING LAW in Lean (klein-S222/HYP-5825, the Redei-involution transplant) -- inBand_mirror (band mirror symmetry via residues-sum-to-0-mod-q) + live_mirror (p -> q-p is a live-preserving involution) + half_live_iff (the half-point q/2 is live iff EVERY speed is odd) + liveCount_parity (LM % 2 = half-point liveness: the +- involution pairs everything else; partition + card_bij') + liveCount_even_of_even_speed (COVERING => LM EVEN at every modulus: live multipliers come in exact +- pairs -- every certificate exhibiting p yields q-p free, halving search; an odd LM at a covering set = a computation-bug detector for the decide pipelines). Kernel-pure, no native_decide
import TournamentH7.LRCParityDemo  -- klein-2026-07-09-S228: the parity law's END-TO-END DEMO on klein-S206's worst covering set {1,2,3,4,7,8,9,10,11,12,13,14,17} at q=21 (live multipliers EXACTLY {4,17} = one +- pair, LM=2) -- demo_live_4 (kernel decide) + demo_live_17_by_law (THE TWIN via live_mirror, NO recomputation -- the certificate-search halving, formal) + demo_LM_even_by_decide/by_law (the validation invariant, two routes agreeing = the bug-detector pattern for the enumeration banks) + demo_mreach/demo_mreach_twin (both witnesses through mreach_ge_of_pairsum_band to Mreach >= 1/14: one decide, two witnesses). Kernel-pure
import TournamentH7.LRCParityBank  -- klein-2026-07-09-S229: THE PARITY LAW WIRED INTO THE 966 BANK (zero new native_decide -- the law does the work): bandOK_mirror (the Z-band mirror, kernel-pure, cast-free sum-to-zero route) + coveringWitnesses_twin_band (ALL 966 twin certificates (l, q-p, q) FREE from the stored witnesses) + coveringWitnesses_lonely_twin (966 SECOND loneliness proofs through mreach_ge_of_pairsum_band) + even_speed_of_coveringPrim (the q=2 covering clause supplies the even speed) + coveringWitnesses_parity (THE STANDING INVARIANT: liveCount even at EVERY modulus for every bank family -- any future evaluation finding an odd live count contradicts a kernel-pure theorem: the bug-detector, wired). New axioms: none beyond the bank's own native_decide
import TournamentH7.LRCWitnessG2Discharge  -- death-star-2026-07-09-S4 (HYP-5770): THE DE-OPAQUING PAYOFF -- witnessG2 and shapeOf are now CONCRETE in LRCFourteenSkeleton (witnessG2 := slowμ(goodSet E ∩ safeSet P); shapeOf := speeds<=13 -> P, rest -> co-offsets vs Vmax; opus-S180's "one substantive remaining Lean obligation" DONE) and the FIRST hfloor instances are DISCHARGED against the skeleton's own vocabulary: safeSet_nil + witnessG2_pure_cluster (P=[] reduces to slowμ(goodSet E)) + hfloor_pure_cluster_diam75 (witnessMP <= witnessG2([],E), via GoodSetBridge <- AP76 <- Farey roof) + hfloor_of_large_speeds (all |v|>13, pairwise spread<=75 => witnessMP <= witnessG2(shapeOf v) -- the exact hypothesis shape lrc14_from_witness_floor consumes); kernel-pure
import TournamentH7.LRCPureClusterCorner  -- death-star-2026-07-09-S1 (HYP-5710): the P=∅ (pure-cluster) corner of the realization node is CLOSED BY SUBSUMPTION -- co-offsets e=Vmax-v in [0,D] with 13D<=12Vmax => speeds in band [Vmax-D,Vmax] ratio<=13 => spread13_lonely (kps-S28) at explicit t=1/(2Vmax-D) => Mreach>=1/14 (pure_cluster_lonely + mreach_ge_of_pure_cluster + mreach_ge_of_ratio_band); the realization node's true residual is speeds-ratio>13 (multi-scale P∪L: boxeph LEM-014 / monad HYP-5717); kernel-pure
import TournamentH7.LRCGoodPeriodJ1
import TournamentH7.LRCGoodPeriodMaxgap  -- klein-S200: the MAXGAP good-period certificate (LEM-013). Decidable IsGoodPeriod/HasGoodPeriod (integer maxCircGap); native_decide witnesses that the exact 7-structured clusters where the ARC-COUNT route fails (MISTAKE-128, c>D3) still have a good period -- the maxgap is the right invariant, not #arcs. Template for the finite-check nodes.
import TournamentH7.LRCHlinkExtract  -- opus-S175: DISCHARGES klein-S203's hlink (free-gap extraction) sorry-free/kernel-pure -- the mergeSort argmax (foldl_max_mem + mem_zipWith_sub_tail + exists_gap_decomp) + BOTH freeness branches (interior via kps-S101, WRAPPING gap directly). mreach_ge_of_goodPeriod_of_hembed: good-period -> Mreach>=1/14 now needs ONLY hembed (THM-527 Part A ruler embedding, the shared blocker).
import TournamentH7.LRCSchurTriples  -- opus-S183: LEM-015 (upper bound, renamed from LEM-014 to avoid boxeph collision) -- the interval maximizes SCHUR TRIPLES E3=#{(a,b) in SxS : a+b in S} <= C(k,2); injects (a,b)->{a,a+b} into powersetCard 2 S (b>0 => a<a+b => card 2 + injective); schurTriple_interval_13: E3({1..13})=78 by decide. The additive-triple (order-3) analogue of Freiman |S+S|>=2k-1; governs the density-floor resonance sum at leading order (S182). Kernel-pure.
import TournamentH7.LRCWitnessMomentFloor  -- opus-S186: ROUTE hfloor THROUGH THE PROVED MOMENT FLOOR -- witness_floor_from_momentfloor_nodes replaces the OPEN Lemma A (nuConsec<=nu) with THM-661s moment floor (nu=mu>=minD3>=momentBar k=witnessMP+1-capRat k); the Bonferroni arithmetic momentBar+capRat-1=witnessMP is DEFINITIONAL (THM-661s (A) bar); hfloor now depends on the fleets PROVED density floor, not the never-proved standalone Lemma A. Kernel-pure.
import TournamentH7.LRC14GrandAssembly  -- monad-S6: THE TOP-LEVEL SURFACE (grand assembly: LRC14Statement from cite + ONE residual obligation, five branches discharged; pure variant kernel-pure)
import TournamentH7.LRCDetunedDispatch
import TournamentH7.LRCCommonResidue  -- monad-S7: THM-668 FORMALIZED (detuned-harmonic dispatch, kernel-pure: quarter-window + Bezout branch + triangle pigeonhole; grand-assembly branch 6)
import TournamentH7.LRCFreimanBurden  -- opus-S188: the FREIMAN DESCENT-BURDEN LOWER BOUND (THM-675 (ii), floor of the burden<=12 finite check) -- restrictedSum_card_ge : 2<=|s| => 2|s|-3 <= |{x+y : x<y in s}| via increasing-chain (min/max two disjoint chains); burden_ge_eleven : |s|=7 => burden>=11. Kernel-pure.
import TournamentH7.LRCMomentFloorConcrete  -- opus-S190: DISCHARGE the moment-floor legs against death-star-S4 concrete witnessG2 -- hbonf (toReal_bonferroni on slowμ) + hsize (clusterSize<=13) discharged; lrc14_from_momentfloor_concrete : LRC14 from {hMoment density floor, hB Lemma B, hsmall pigeonhole, hpartA reach} on concrete measures. Kernel-pure.
import TournamentH7.LRCMomentLP  -- opus-S191: the one-sided MOMENT-LP CORE (THM-661 analytic link) -- integral_le_measure_pos : prob measure, measurable W, p(w)<=1_{w>0} => integral p(W) <= (mu{W>0}).toReal; measure_pos_ge_of_moment_ge : b<=integral p(W) => b<=mu(GOOD). The missing link between LRCD3FloorCert (rational D3>=bar) and mu(GOOD). Kernel-pure.
import TournamentH7.LRCFareyCellMoment  -- opus-S193: the SINGLE-CELL FAREY MOMENT IDENTITY -- cell_moment : B!=0 => integral_a^b (A+Bx)^i = ((A+Bb)^(i+1)-(A+Ba)^(i+1))/(B(i+1)) via FTC (antiderivative (A+Bx)^(i+1)/(B(i+1))); cell_moment_const : integral_a^b A^i = (b-a)A^i. The atomic per-cell contribution of the Farey moment identity m_i=integral W^i (feeds opus-S192 momentLP_from_coeffs). Kernel-pure.
import TournamentH7.LRCFreimanAP  -- opus-S195: the FREIMAN AP STEP (n>=5) -- ap_of_min_burden : StrictMono a + |Rset a n|=2n-3 => a is an AP (all consecutive diffs = a1-a0). Interleaved-chain route (Rset = {a_i+a_{i+1}} u {a_i+a_{i+2}}); diff_two (d_{i+2}=d_i) + sum04 (a0+a4=a1+a3, needs n>=5) => AP. n>=5 ESSENTIAL (MISTAKE-133: n<=4 admits bi-arithmetic non-AP minima). The near-AP char THM-675 needs for 7-classes. Kernel-pure.
import TournamentH7.LRCFreimanAPBridge  -- opus-S196: FINSET BRIDGE for the AP step -- finset_min_burden_isAP : |s|>=5 & |restrictedSum s|=2|s|-3 => s = {c+k*d : 0<=k<|s|}, d>0 (s is LITERALLY an AP). Via enum (sorted orderEmbOfFin extended StrictMono to N) + Rset_enum_eq (Rset(enum s)|s| = restrictedSum s) + ap_of_min_burden. Directly citable by THM-675/THM-681. Kernel-pure.
import TournamentH7.LRCBurdenGap  -- opus-S198: DISJOINT-BLOCK burden bound (LEM-016(i) dominant-gap case) -- burden_ge_of_dominant_gap : ANY k-set split L|R at a gap dominating both spans (span L + span R < gap, i.e. gap > D/2) => burden >= 3k-7; k=7 corollary burden_ge_fourteen_of_dominant_gap (>=14, LEM-016(i)); k=13 => >=32 (monad THM-682 core B>=32, dominant-gap branch). three_block_card (L+^L, L+R, R+^R disjoint ranges) + restrictedSum_card_ge (within) + cauchy_davenport_add (cross). Kernel-pure.
import TournamentH7.LRCDepth4Dispatch  -- mac-mini-S65cont13: LEM-021 depth-4 dyadic dispatch -- mreach_ge_of_depth4: no 16-multiple + odd speeds miss a unit +-class mod 16 (c in {1,3,5,7}, m = c^{-1}) => band at q=16, p=m => Mreach >= 1/14, clearance 1/8. The even free pass (unique layer: 16 < 28 < 32). Producer for mreach_ge_of_pairsum_band.
import TournamentH7.LRC14CompletionAudit  -- opus-S199: the CURRENT finish-line completion audit -- lrc14_from_B5 : LRCUpTo13 + hB5 (single B5>0 residual obligation) => LRC14Statement; #print axioms = [propext, Classical.choice, Quot.sound, 2x winData22 native_decide], NO sorryAx. Certifies LRC(14) is sorry-free-sound, reduced to the LRC(<=13) citation + ONE open analytic obligation.
import TournamentH7.LRCWitnessFloorRepair  -- mac-mini-S65cont17: THE REPAIRED WITNESS-FLOOR ASSEMBLY (implements the cont.16 soundness-flag fix): lrc14_from_repaired_nodes = LRC14Statement from SATISFIABLE legs -- k=0 via the q=14 sieve (all speeds <= 13, internal), k in {1,2} positivity (engine floors 7/858, 313/9702), 3<=k<=7 and 8<=k<=13 the m_P floor (THM-530 admissible range), + hpartA. Replaces the unsatisfiable original hfloor quantification.
import TournamentH7.LRCIntervalBridge  -- mac-mini-S65cont18: the interval bridge -- witnessG2_pos_of_anchor: an explicit rational interval inside goodSet with checkable band bounds for the slow part => 0 < witnessG2 (the hk12 leg shape + the positivity core of hsmall3/hlarge); slowmu positivity from Ico-subset (no measurability needed); the engine (cont.16) emits the anchors
import TournamentH7.LRCWindow22Census  -- opus-S202: LEM-024 the 6-WITNESS PIGEONHOLE (kernel-pure, NO native_decide) -- window22_lonely : EVERY covering 13-tuple in [1,22] (distinct, positive) => lonely; min=1 via far at one of {12/25,9/26,7/27,11/28,4/23,11/26} (danger sets {2},{3},{4},{5},{6,17},{7,19}), fail-all-6 + covering + min1 = 14 forced distinct > 13 slots; min>=2 via spread13. Foundational-axioms-only REPLACEMENT for the C(22,13)=497420 winData22 native_decide census. Axioms [propext,Classical.choice,Quot.sound].
import TournamentH7.LRCGoodSetSmall  -- mac-mini-S65cont19: brick (i) -- goodSet_univ_of_length_le_two: clusters of length <= 2 have goodSet = univ (fract(-y) = 1 - fract(y) forbids double arc-occupancy); with the interval bridge + the cont.19 anchor table, the hk12 leg is anchors-in-proof-out
import TournamentH7.LRCFiniteUnionVolume  -- opus-S204: BRICK (iii) -- the finite-union volume identity mac-mini deferred. slowmu_toReal_ge_sum_of_disjoint_Ico : pairwise-disjoint anchors Ico a_i b_i in [0,1) all inside S => Sum (b_i - a_i) <= (slowmu S).toReal. Corollaries witnessG2_ge_sum_of_disjoint_anchors / witnessG2_ge_of_anchor_floor give the exact m_P floor shape for hsmall3/hlarge. Critical path: klein THM-685 (Kronecker transfer) => remaining analytic content = measure floors. Kernel-pure.
import TournamentH7.LRCPrimitivePeel  -- opus-S205: THE PRIMITIVITY PEEL -- lrc14_of_primitive : LRC(14) on primitive families (tupleGcd=1) => LRC14Statement, via lonely_of_dilate (lonely_scale). MOTIVATION: the ResidualObligation admits NON-PRIMITIVE DILATES (witness 2*[1..9,11,12,13,20], gcd 2, mu = 1/980, core Vmax=20 window-censused); since mu(c*w)=mu(w) exactly, inf mu = 0 over the residual => NO UNIFORM MEASURE FLOOR EXISTS. Peel first, then the floor is well-posed (primitive residual min mu ~ 0.0094). Kernel-pure.
import TournamentH7.LRCUnionVolume  -- mac-mini-S65cont20: brick (iii) COMPANION form (wire-priority CEDED to opus-S204 LRCFiniteUnionVolume, landed first) -- slowmu_ge_sum_of_sorted_Ico_subset: SORTED-adjacent list form at ENNReal level (ofReal(sum lengths) <= slowmu S), the shape the cont.16 engine merged lists emit natively (linear adjacent checks); opus form has the witnessG2/toReal corollaries -- use theirs for the legs, this for list-shaped data
import TournamentH7.LRCGoodSetBand  -- mac-mini-S65cont21: the GOODSET BAND LEMMA -- Ico_subset_goodSet_of_bounds (tooth + per-difference up/down band checks => Ico inside goodSet E; up: e*x in (j+1/7, j+1], down via fract-reflection: e*x in [j, j+6/7)) + Ico_subset_good_inter_safe + witnessG2_ge_of_sorted_bands (sorted engine list with per-interval tooth+band certs => witnessG2 >= sum of lengths, consuming brick (iii)). The discharge pipeline is now COMPLETE: hsmall3/hlarge = engine-lists-in, proof-out. Kernel-pure.
import TournamentH7.LRCPrimitiveAssembly  -- opus-S206: GRAND ASSEMBLY + PRIMITIVITY PEEL -- lrc14_grand_assembly_primitive : LRCUpTo13 + ResidualObligationPrimitive (residual clauses + tupleGcd v = 1) => LRC14Statement. Removes the non-primitive DILATES that made inf mu = 0 (opus-S205), so the measure-floor target becomes well-posed. residualPrimitive_of_residual + lrc14_grand_assembly_of_residual show the new surface strictly subsumes the old. lrc14_from_B5_primitive = the finish line with tupleGcd=1. Same 8 branches; kernel-pure, foundational axioms only.
import TournamentH7.LRCLegDemo  -- mac-mini-S65cont21: the FIRST MACHINE-GENERATED leg instance -- witnessG2 (([1..10], [0,1,2])) >= 1217/8820 EXACT (consecutive large cluster, 20 safeSet components, 260 rational band certificates all norm_num). Emitted by lrc14_leg_codegen_macmini_S65cont21.py from the cont.16 engine: ENGINE-LISTS-IN, PROOF-OUT is now literal. Template for all hsmall3/hlarge families. Kernel-pure.
import TournamentH7.LRCResidualMeasureFloorPrimitive  -- opus-S207: the measure floor on the PRIMITIVE residual (sharpest reduction) -- SafeMeasureFloorPrimitive (0 < volume safePeriod for residual v with tupleGcd v = 1) + lrc14_of_measureFloor_primitive : LRCUpTo13 + SafeMeasureFloorPrimitive => LRC14Statement, kernel-pure foundational-only. Composes kps SafeMeasureFloor bridge with opus-S206 primitivity peel. S205/S207: inf mu on the primitive residual ~ 0.0085 > 0 (uniform floor TRUE; minimizers at small Vmax, mu -> (6/7)^13 as Vmax grows).
import TournamentH7.LRCDissociatedAssembly  -- opus-S209: the d=2,3 DETUNED PEEL -- lrc14_grand_assembly_dissoc / lrc14_from_B5_dissoc : LRCUpTo13 + MultiDetunedDispatch (THM-678 citation: some g>=2 at detuning level 2 or 3 => lonely) + ResidualObligationDissoc (primitive residual + no g at level 2/3) => LRC14Statement. Peels the near-dilate d=2,3 detuned families (the S208 small-mu minimizers, mu~0.0085) via THM-678, shrinking the floor obligation to the DISSOCIATED residual (every g>=2 at level >=4) where mu decorrelates -> (6/7)^13. Kernel-pure, foundational-only.
import TournamentH7.LRCDetunedD2  -- opus-S210: THM-678 d=2 GENERIC, reduced to the counting -- lonely14_of_two_detuned : LRCUpTo13 + TwoDetunedClearing + (v = g*H u {i1,i2} detuned, (q1,q2)!=(2,2)) => lonely. Fully formalizes the LRC(11) harmonic reduction (via orderEmbOfFin of the complement; clearance 1/12 at every branch) + the construction core (lonely14_of_two_detuned_good). TwoDetunedClearing (Sum N_j/q_j < 1 counting: a branch c clears both detuned phases) = the SOLE remaining analytic piece. Kernel-pure.
import TournamentH7.LRCMomentFloorRepair  -- mac-mini-S65cont22: MISTAKE-136 repair -- hB of lrc14_from_momentfloor_concrete is UNSATISFIABLE over all Shape (([0],[0x13]): capRat 13 = 1 > measGP = 0); lrc14_from_momentfloor_concrete_shapes narrows hB to reachable shapes; + safeSet_anti (dedup/superset bridge) + measGP_ge_of_sorted_bands (safe-only certificate consumer, brick iii + cont.18 bands)
import TournamentH7.LRCSafeCertSize1  -- mac-mini-S65cont22: hB certificate table |S|=1 (13 families, 91 comps, cap 6/7 EXACT uniform) -- machine-generated from the cont.16 engine; capRat ladder = per-|S| minima exactly (engine-verified all 6 rows)
import TournamentH7.LRCSafeCertSize2  -- mac-mini-S65cont22: hB table |S|=2 (78 families, 917 comps, cap 66/91, argmin {1,13})
import TournamentH7.LRCSafeCertSize3  -- mac-mini-S65cont22: hB table |S|=3 (286 families, 4397 comps, cap 55/91, argmin {1,12,13})
import TournamentH7.LRCSafeCertSize4_g1to2  -- mac-mini-S65cont22: hB table |S|=4 first elem 1-2 (385 families, 6268 comps, cap 1979/4004, argmin {1,11,12,13}); g3to10 + Size5 files emitted, building
import TournamentH7.LRCSafeCertSize4_g3to10  -- mac-mini-S65cont22: hB table |S|=4 first elem 3-10 (330 families) -- TABLE COMPLETE with g1to2
import TournamentH7.LRCSafeCertSize5_c1  -- mac-mini-S65cont22: hB table |S|=5 first elem 1 (495 families incl. the SHARP argmin {1,5,7,8,9} = capRat(8) exactly)
import TournamentH7.LRCSafeCertSize5_g2to3  -- mac-mini-S65cont22: hB table |S|=5 first elem 2-3 (540 families)
import TournamentH7.LRCSafeCertSize5_g4to9  -- mac-mini-S65cont22: hB table |S|=5 first elem 4-9 (252 families) -- ALL 2380 hB certificates GREEN: every subset S of {1..13}, |S| <= 5, floored at its exact capRat ladder constant, kernel-pure
import TournamentH7.LRCTightRigidity  -- kps-S127: the EXTREMAL-UNIQUENESS reduction -- TightRigidity (mu=0 => DilatedFamily c*{1..13}, open, >= LRC-hard) => SafeMeasureFloor => LRC14; not_dilated_of_gapFamily (dilated ratio EXACTLY 13, residual >13 => never dilated). Difference-primitive COLLAPSE: dilated_primitive_eq_range (Primitive + Dilated => image|v| = {1..13} exactly, c=1) + PrimitiveTightRigidity (named, open) + primitiveTightRigidity_of_tightRigidity. Kernel-pure, foundational-only.
import TournamentH7.LRCDissociatedRigidity  -- kps-S127: WIRES the primitive collapse into opus-S209's dissociated peel -- primitive_of_tupleGcd_one (opus tupleGcd=1 == my Primitive) + gapFamily_image_ne_range (GapFamily forbids the AP {1..13}) => safeMeasureFloorPrimitive_of_primitiveTightRigidity (PrimitiveTightRigidity => opus SafeMeasureFloorPrimitive) => residualObligationDissoc_of_primitiveTightRigidity + lrc14_of_primitiveTightRigidity_dissoc (through opus dissoc assembly) + lrc14_of_primitiveTightRigidity (hMD-FREE: rigidity floors the WHOLE primitive residual, THM-678 dispatch redundant on this route => rigidity strictly stronger than the floor opus needs). Kernel-pure, foundational-only.
import TournamentH7.LRCSafeCertLeafTrees  -- mac-mini-S65cont23: the hB dispatcher LEAF TREES -- safe_floor_sorted_len1..5: every sorted tuple 1 <= p1 < ... < pl <= 13 dispatches to its cont.22 table certificate via explicit interval_cases trees (2379 leaves, lengths 3-5 split by head element). Kernel-pure.
import TournamentH7.LRCSafeCertDispatch  -- mac-mini-S65cont23: the hB SHAPEOF DISPATCHER -- hB_certs (dedup + insertionSort canonicalize any reachable small part; safeSet_congr + capRat_mono + cap_le_of_canonical => capRat k <= measGP, LEMMA B IS A THEOREM) + lrc14_from_momentfloor_certs (the moment route with hB DISCHARGED: LRC14 from exactly hMoment (THM-661) + hsmall + hpartA, foundational-axioms-only)
import TournamentH7.LRCIntervalCount  -- opus-S211: counting lemmas -- card_le_of_mem_Ioo (integers in a width-W real interval <= floor(W)+1) + bad_count_le (|bad_j over [0,g)| <= gcd*(floor(q/7)+1) via the de-circled psi=p*c-q*round injection) + sum_lt (the two bounds sum < g iff (q1,q2)!=(2,2)). Kernel-pure.
import TournamentH7.LRCTwoDetunedClearing  -- opus-S211: THM-678 d=2 COUNTING PROVED -- twoDetunedClearing : TwoDetunedClearing (union bound: two bad-branch sets cannot cover [0,g) when (q1,q2)!=(2,2), so an uncovered c clears both detuned phases). lonely14_of_two_detuned' makes generic d=2 detuned UNCONDITIONAL from LRC(<=13). Kernel-pure, foundational-only.
import TournamentH7.LRCDetunedD3  -- kps-S127: THM-678 d=3 GENERALIZATION (generic), PROVED -- v = g*H u {d1,d2,d3}, three detuned. threeDetunedClearing (THREE-set union bound reusing opus's per-coordinate LRCIntervalCount.bad_count_le for the 3rd delta; ThreeDetunedClearing holds when Sum_j badCount(d_j,g) < g) + lonely14_of_three_detuned_good (construction core) + lonely14_of_three_detuned' (UNCONDITIONAL from LRC(<=13), via LRC(10) harmonic clearance 1/11>=1/14). Generic d=3 closed; exceptional set = (2,2,*) [infinite double-half-harmonic, mod-2g lift] + finite small-q (q1 in {2,3}), analogue of opus's d=2 (2,2) residual. Kernel-pure, foundational-only.
import TournamentH7.LRCDetunedDispatchReduce  -- kps-S127: WIRES the d=2/d=3 clearings into opus's MultiDetunedDispatch + SHRINKS the citation -- lonely14_of_nonMultCard_three (nonMultCard v g = 3 + generic Sum badCount < g => lonely, extracting the 3 detuned coords) + lonely14_of_nonMultCard_two (d=2 analog, bridges the sum-count to opus's (q1,q2)!=(2,2) via badCount_of_q_two) + multiDetunedDispatch_of_exceptional (MultiDetunedDispatch <= ExceptionalDetunedDispatch: the generic bulk is PROVED, only the non-generic half-harmonic residual (2,2)/(2,2,*)+finite stays cited) + lrc14_grand_assembly_dissoc_exceptional. Kernel-pure, foundational-only.
import TournamentH7.LRCTwoDetunedLift  -- kps-S127: the (2,2) mod-2g lift TERMINATING BASE CASE, PROVED -- lonely14_of_two_detuned_lift2g : a d=2 detuned family whose ENTIRE harmonic part is divisible by 2g (no odd multipliers) is lonely, via opus's lonely14_of_two_detuned' at the DOUBLED modulus 2g (both detuned sit at q=4 there, generic count 1/2<1 fires). Closes the (2,2) congruent-half-harmonic residual exactly when the 2-adic descent terminates in one step. HONEST: the general (2,2) case (some odd multiplier m: g*m becomes a new half-harmonic q=2 at 2g) DESCENDS the 2-adic tower (monad THM-678) to the open core -- NOT closed here. Kernel-pure, foundational-only.
import TournamentH7.LRCMeasureTransfer  -- klein-S242: THM-685's transfer pipeline -- rational safe-interval certificates (SafeIvStrict) => strictly-live rulers at EVERY q > D/(y-x) (exists_grid_point + scaled strict rounding identity) => kps StrictWitness; demo: the deep well certified once at [93/1274, 96/1274], strictly live at every q >= 425
import TournamentH7.LRCMomentCitation  -- mac-mini-S65cont24: the THM-661 CITATION NODE -- THM661MomentFloor (faithful uniform covering-moment floor over NODUP clusters k=8..13, transcription validated to the rational against canon: bars AND block-nu table six-for-six) + SmallClusterFull (<=7-arcs pigeonhole companion, 3<=l<=7; l<=2 PROVED via brick i) + goodSet_dedup + momentBar_le_one/anti + hMoment_of_citations (dedup 3-way dispatch handles duplicate co-offsets) + lrc14_from_thm661_certs: LRC14 from TWO NAMED CITATIONS + hsmall + hpartA, foundational-axioms-only. hMoment GONE from the obligation surface.
import TournamentH7.LRCTestPointCore  -- klein-S243: the test-point theorems' cores (THM-690/691A/692) -- net_value_strictly_in_band + qstar_p_nonzero + residue_in_range (the P-side), pigeonhole_missed_class (the E-side), middle_chain_rigidity (THM-692's both-sides impossibility), and THE FATTENING LEMMA qstar_cert: a q*-test-point generates SafeIvStrict certificates on [(4732a-q)/4732q, (4732a+q)/4732q] -- a certificate factory feeding LRCMeasureTransfer
import TournamentH7.LRCTwoScaleWitness  -- klein-S244: THM-693 THE FINITE-V TWO-SCALE WITNESS -- the explicit multiplier w = aV + (c - aV) % q* is STRICTLY LIVE at modulus q*V for every two-scale family P u {V-e} (P-side residues + missed class c + V > 2184, V > 168e): the test-point program made constructive, no limits/transfer; demo (130000, 10001) for {1..5} u {V-E} at V = 10000
import TournamentH7.LRCMPCertSize10  -- mac-mini-S65cont25: the m_P CERTIFICATE TABLE |S|=10 -- 286 families, m_P = 14249/252252 <= mu(safeSet S); min attained EXACTLY at the THM-530 argmin {1,2,3,5,7,8,9,11,12,13}
import TournamentH7.LRCAnchor12Certs  -- mac-mini-S65cont25: 13 positivity anchors for the twelve-element families (hk12 dispatch; anchors cross-validate the cont.19 table)
import TournamentH7.LRCMPLeafTrees  -- mac-mini-S65cont25: m_P dispatch trees sizes 6..10 (5720 leaves; sizes 6-9 complete to size-10 supersets via measGP_anti + decide) + anchor_pos_of_missing (13-way missing-element dispatch)
import TournamentH7.LRCSmallDischarge  -- mac-mini-S65cont25: hsmall DISCHARGED VIA THE REPAIRED LEGS -- hk12 (goodSet univ + pigeonhole + anchors), hsmall3 (nu=1 + Bonferroni + m_P layer), hlarge (citations + hB + exact identity) => lrc14_from_two_citations: LRC(14) FROM EXACTLY [THM-661] + [<=7-arcs pigeonhole] + [hpartA]. THE MOMENT ROUTE IS CLOSED; every other node is a theorem, foundational-axioms-only
import TournamentH7.LRCReachCitation  -- mac-mini-S65cont26: the THM-527-A REACH CITATION + THE FULLY CITATION-CLOSED ASSEMBLY -- PROVED glue: nonzero_of_witnessG2_pos (the hv guard is DERIVABLE: zero speed empties safeSet, so hpartA's guard-free quantification is sound -- MISTAKE-136 probe passes), exists_config_of_witnessG2_pos (positivity => explicit slow time in good∩safe), Mreach_ge_of_minReach (the compactness half is glue). CITED: THM527ARulerEmbedding (the slow-fast ruler embedding = canon THM-527-A limit + O(1/Vmax) + the census/banks finite-V closure). PRIZE: lrc14_from_citations_only -- LRC(14) from THREE NAMED CITATIONS (THM-661 + <=7-arcs pigeonhole + THM-527-A reach) and NOTHING ELSE, foundational-axioms-only
import TournamentH7.LRCMultiScaleWitness  -- klein-S245: THM-694 THE MULTI-SCALE WITNESS -- band_lift (in-band at (q,a) survives to (qV, aV+Delta) for ANY Delta, 14Bq < V) + cluster_join + miss_at_next (coarse missed class lifts band-interior) => threeScale_strictlyLive: P u C2 u C1 strictly live at q*V2V1 with the mixed-radix multiplier w1 = w2V1 + (c1V2 - w2V1) % (q*V2); demo (78*10^12, 6010000020000); general r = the same two lemmas per scale
import TournamentH7.LRCReachDecitation  -- mac-mini-S65cont27: the ruler embedding DE-CITED to arithmetic via klein transfer -- PROVED: minReach_ge_of_strictWitness (strict-margin time => minReach >= 1/14: floor/floor+1 margin instantiation), strictWitness_abs (|v| normalization, reflection-invariant margins), rulerEmbedding_of_certificateSupply (klein strictWitness_of_cert consumes 13 SafeIvStrict certs => StrictWitness => the embedding). NARROWED CITATION: THM527ACertificateSupply = pure integer certificate existence (no measure/reach/compactness content). lrc14_from_certificate_citations: LRC(14) from THM-661 + pigeonhole + certificate existence -- all reach geometry, strictness transfer, and compactness are now THEOREMS
import TournamentH7.LRCRayWitness  -- klein-S247: THM-695 THE CONSTRUCTIVE RAY WITNESS -- ray-coupled families explicit: w = b*w_693 at q = q*bV (scaled_band inherits small/top bands; ray_speed_band: (u w) % (q*bV) = rho'_f V - delta(c+bf) under the decidable digit hypothesis); the certificate supply now covers two-scale + multi-scale + ray shapes; demo (260000, 20002)
import TournamentH7.LRCThreeDetunedCoarse  -- opus-S212: THM-678 d=3 clean coarse corollary -- lonely14_of_three_detuned_coarse: all q_i>=8 => Sum N_i/q_i <= 45/56 < 1 => lonely, no per-instance sum check (the all-fine case kps-S127 LRCDetunedD3 left generic). Builds on kps threeDetunedClearing + reuses opus bad_count_le. Kernel-pure.
import TournamentH7.LRCTestDataSupply  -- klein-S248: THM-696 TEST-DATA EXISTENCE + SUPPLY DISCHARGE -- supply_of_strictlyLive (THE FATTENING BRIDGE: any strictly-live ruler => the citation's SafeIvStrict interval data, D = 28BQ, x/y = 28Bw -+ 1: ALL witness theorems feed THM527ACertificateSupply) + qstar_exists (6 > 5 pigeonhole) + testData_exists (data from shape alone, a = 1) + twoScale_supply (THE CITATION'S CONCLUSION proved wholesale on every non-top-block two-scale family, V > 2184)
import TournamentH7.LRCFourierCompletionB  -- opus-S213: LEM-022 Fourier completion STAGE B (the t2 per-cell bound, taking over death-star dormant) -- sum_exp_orthogonality (Sum_{x<q} e_q(hx) = q*1{q|h}) + sine_cdist_witness (|sin(pi h/q)| >= 2 cdist(h)/q, via |sin| q-periodicity + Jordan) + norm_bandCoeff_le (band Fourier coeff ||B_hat(h)|| <= q/(2 cdist h), Stage A composed). Kernel-pure. Feeds the OffLine<=f(E3) gate (THM-680) / discharges mcorr<=M.
import TournamentH7.LRCFirstWindowWitness  -- klein-S249: THM-697 THE FIRST-WINDOW WITNESS -- t = (7j+6)/(7V) puts the cluster phase at EXACTLY 6/7 (residue 6V - ew, ONE inequality per cluster speed); small speeds ride the open first window (nondegenerate for EVERY P incl. all top-blocks: no test point, no missed class); leftmost-j constructive; demo: the k=8 TOP-BLOCK {9..13} with the CONSECUTIVE APEX cluster strictly live at (70000, 559) + its supply data via the fattening bridge -- S248's packed top-block gap CLOSED
import TournamentH7.LRCSevenGapRigidity  -- klein-S251/S252: THE <=7-ARCS PIGEONHOLE PROVED (goodSet_ae_full: slowmu(goodSet E) = 1 for every nodup 3 <= |E| <= 7) -- sorted-phase cyclic gaps telescope to 1 (cgap_sum); ANY gap > 1/7 hands the goodSet witness (cgap_witness: interior gap at s_i forces s_i < 6/7 so wrap terms clear 1/7 free); k <= 6 pigeonhole / k = 7 all-gaps-exactly-1/7 PERFECT NET (gap_dichotomy); net times = countable fibers x = (n+1/7)/d (pair_fiber_countable) => bad set null. Kernel-pure x6. The SmallClusterFull naming + two-citation assembly: LRCTwoCitationAssembly.
import TournamentH7.LRCPackedSupply  -- klein-S250/S251: THE PACKED-SUPPLY DICHOTOMY (THM-698's Lean half) -- packed_dichotomy (no-q*-room => P = [|E|+1,13] EXACTLY by interval containment + card equality => pmin >= 9) + firstWindow_supply (the window route wrapped to the citation's data for pmin >= 9 shapes, e <= 90, V > 2366) + packed_supply (the combined dispatcher: EVERY packed two-scale shape satisfies THM527ACertificateSupply's conclusion). Kernel-pure.
import TournamentH7.LRCFourierAggregation  -- kps-S127: LEM-022 Fourier completion STAGE B.3 off-diagonal aggregation (the missing middle brick) -- offDiag_bandSum_le (||Sum_{h!=0} bc(h) conj(bc(w h))|| <= (q^2/4) Sum 1/(cdist h cdist(w h)) from opus's per-coeff bound B.2, via triangle + termwise) + offDiag_bandSum_le_closed (compose with death-star's harmonic_ratio_sum_mul_le => closed form 5 q^2 (log2 q+1)^2 / P). REDUCES the whole Fourier-completion node to JUST the Parseval completion identity C_w = b^2/q + (1/q) Sum_{h!=0} B_hat(h) conj(B_hat(w h)). Kernel-pure, foundational-only.
import TournamentH7.LRCFourierCompletionC  -- opus-S214: LEM-022 Fourier completion STAGE B.3 -- completion_identity (C_w = (1/q) Sum_{k<q} B_hat(k) conj(B_hat(wk)), the analytic crux, via eInt char infra + expand/swap/collapse(B.1)/count) + completion_diff_bound (||C_w - b^2/q|| <= (1/q) Sum_{k!=0} ||B_hat(k)|| ||B_hat(wk)||, the k=0 split + triangle). Kernel-pure. The completion half of the t2 hyperbola bound feeding OffLine (THM-680).
import TournamentH7.LRCFourierClosed  -- kps-S127: LEM-022 Fourier completion CLOSED-FORM CAPSTONE -- completion_closed_of_coeffBound : from opus's completion_diff_bound + opus's per-coeff bound B.2 (as hyp hcoeff) + death-star's harmonic_ratio_sum_mul_le, ||C_w - b^2/q|| <= 5 q (log2 q+1)^2 / P. Two bridges: P>0 + ratio floor => w*z!=0 (B.2 hits both factors, no unit hyp); range-q sum -> ZMod q sum via natCast/val bijection (Finset.sum_bij'). CLOSES the LEM-022 Fourier node modulo its proved inputs (identity opus-C, coeff opus-B.2, harmonic death-star). Kernel-pure, foundational-only.
import TournamentH7.LRCZcorrEnergy  -- kps-S127: the ZMod q OFF-DIAGONAL ENERGY aggregation (wires zcorr_percell into offdiag) -- zcorr_one (zcorr A 1 = |A|) + sum_zcorr (units: Sum_w zcorr A w = |A|^2) + offdiag_zcorr_sq_le (per-cell zcorr A w <= M => Sum_{w!=1} zcorr^2 <= M(|A|^2-|A|), the ZMod q realization of offdiag_mcorr_sq_le) + zcorr_energy_of_hyperbola (compose death-star hyperbola_box_count via zcorr_percell => Sum_{w!=1} zcorr^2 <= (1 + 4K^2/P)(|A|^2-|A|)). Completes hyperbola_box_count -> zcorr_percell -> t2 off-diagonal energy. Kernel-pure, foundational-only.
