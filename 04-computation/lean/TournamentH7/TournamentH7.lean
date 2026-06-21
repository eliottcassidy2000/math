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
    · TournamentH7.LRCL7Discrepancy — Finite integer core of the L7 discrepancy.
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
import TournamentH7.LRCDeathChain
import TournamentH7.LRCFactorialAtom
import TournamentH7.LRCBooleanTypeCut
import TournamentH7.LRCPeriodmaxCertificate
import TournamentH7.LRCGenuineWideCorrection
import TournamentH7.LRCQ6Contraction
import TournamentH7.LRCL7Discrepancy
import TournamentH7.Verify
