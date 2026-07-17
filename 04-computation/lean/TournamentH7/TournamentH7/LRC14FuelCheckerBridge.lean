/-
  TournamentH7.LRC14FuelCheckerBridge — certificate coverage feeds the legacy
  `hpartA` endgame node.

  This module closes a bookkeeping gap, not the global LRC(14) crux.  The
  certificate layer already proves two soundness routes:

  * `dispatch v` (the finite census pages) gives a lonely time;
  * a passing `LadderPack` containing every speed of `v` gives a lonely time
    through the module-6 fuel checker.

  The old assembly surface nevertheless asked directly for the guard-free node

      0 < witnessG2 (shapeOf v) -> 1/14 <= Mreach v.

  `EndgameCertificateComplete` is the exact remaining coverage statement saying
  that every positive-witness family is handled by one of those two certificate
  routes.  The theorems below prove, without further analytic assumptions, that
  this named coverage statement supplies `hpartA` and hence the old endgame.

  Honest global audit (2026-07-16):

  * PROVED here/imported: census-page soundness, ladder fuel-checker soundness,
    lonely-time-to-`Mreach`, and their composition into `hpartA`;
  * LEGACY WARNING: the original all-shapes `hfloor` is known false
    (`LRCWitnessFloorRepair`); its corollary below is compatibility-only;
  * REMAINING on the repaired route: its three satisfiable floor legs
    (`hk12`, `hsmall3`, `hlarge`) plus `EndgameCertificateComplete` (the actual
    certificate coverage census); existing downstream modules reduce those
    floor legs further to named canon citations;
  * SHARPER CURRENT ROUTE: `LRC14CompletionAudit` reduces LRC(14) to the
    sanctioned `cite : LRCUpTo13` plus the residual positive-`B5` supply;
  * SORRIES introduced here: none.  The `#print axioms` audit at the bottom is
    the authoritative footprint.
-/
import TournamentH7.LRC14Assembly
import TournamentH7.LRC14Dispatch
import TournamentH7.LRCWitnessFloorRepair

namespace LonelyRunner
namespace LRC14

/-- A family is certified either by a finite census-page row or by a passing
module-6 ladder pack containing all thirteen speeds. -/
def EndgameCertified (v : Fin 13 → ℤ) : Prop :=
  dispatch v ∨ ∃ D : LadderPack, ladderOK D ∧ ∀ i, v i ∈ D.speeds

/-- The exact coverage premise left after the two existing certificate
soundness routes are merged: every positive witness shape is certified. -/
def EndgameCertificateComplete : Prop :=
  ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) → EndgameCertified v

/-- Direct module-6-to-skeleton bridge: a passing ladder pack containing the
family clears the concrete `Mreach` threshold. -/
theorem Mreach_ge_of_ladder_mem (D : LadderPack) (hOK : ladderOK D)
    (v : Fin 13 → ℤ) (hmem : ∀ i, v i ∈ D.speeds) :
    (1 : ℝ) / 14 ≤ LRC14Concrete.Mreach v :=
  Mreach_ge_of_lonely (lonely_of_ladder_mem D hOK v hmem)

/-- Every merged endgame certificate clears the skeleton's reach node. -/
theorem Mreach_ge_of_endgameCertified (v : Fin 13 → ℤ)
    (hcert : EndgameCertified v) :
    (1 : ℝ) / 14 ≤ LRC14Concrete.Mreach v := by
  rcases hcert with hdispatch | ⟨D, hOK, hmem⟩
  · exact Mreach_ge_of_dispatch v hdispatch
  · exact Mreach_ge_of_ladder_mem D hOK v hmem

/-- The named certificate-coverage census supplies the legacy guard-free
`hpartA` node. -/
theorem hpartA_of_endgameCertificateComplete
    (hcomplete : EndgameCertificateComplete) :
    ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v := by
  intro v hpos
  exact Mreach_ge_of_endgameCertified v (hcomplete v hpos)

/-- Legacy two-node endgame with the formerly implicit fuel/census bookkeeping
made explicit as one certificate-coverage premise. -/
theorem lrc14_endgame_of_certificateCoverage
    (hfloor : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hcomplete : EndgameCertificateComplete) :
    LRC14Statement :=
  TournamentH7.LRC14Assembly.lrc14_endgame hfloor
    (hpartA_of_endgameCertificateComplete hcomplete)

/-- Usable repaired endgame: the three satisfiable witness-floor legs and the
explicit certificate-coverage premise imply LRC(14).  Unlike the compatibility
theorem above, this theorem never assumes the known-false all-shapes `hfloor`. -/
theorem lrc14_from_repaired_certificateCoverage
    (hk12 : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      1 ≤ clusterSize (shapeOf v) → clusterSize (shapeOf v) ≤ 2 →
      0 < witnessG2 (shapeOf v))
    (hsmall3 : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      3 ≤ clusterSize (shapeOf v) → clusterSize (shapeOf v) ≤ 7 →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hlarge : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      8 ≤ clusterSize (shapeOf v) → clusterSize (shapeOf v) ≤ 13 →
      (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hcomplete : EndgameCertificateComplete) :
    LRC14Statement :=
  Repair.lrc14_from_repaired_nodes hk12 hsmall3 hlarge
    (hpartA_of_endgameCertificateComplete hcomplete)

/-! ## Axiom audit -/

#print axioms Mreach_ge_of_ladder_mem
#print axioms Mreach_ge_of_endgameCertified
#print axioms hpartA_of_endgameCertificateComplete
#print axioms lrc14_endgame_of_certificateCoverage
#print axioms lrc14_from_repaired_certificateCoverage

end LRC14
end LonelyRunner
