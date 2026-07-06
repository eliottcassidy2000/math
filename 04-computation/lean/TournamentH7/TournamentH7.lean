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
import TournamentH7.LRCGapCandidate  -- kps-S24(2026-07-06): GAP CANDIDATES ARE DIVISIBILITY-RICH (gap_candidate_has_multiple: a 12-family not loose at 2/25 contains a multiple of every k in {2..12}, via opus int_far_of_not_dvd_k + 1/k > 2/25; gap_candidate_prime_powers: multiples of 5,7,8,9,11,12 forced) -- the covering-system structure of gap candidates, AP = minimal; sharpens the additive/collision side of the S23 split (HYP-4417)
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
