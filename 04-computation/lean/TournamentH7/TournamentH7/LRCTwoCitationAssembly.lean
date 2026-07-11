/-
  TournamentH7.LRCTwoCitationAssembly -- LRC(14) FROM TWO CITATIONS
  (klein-2026-07-11-S251, HYP-5985): the ≤7-arcs pigeonhole citation
  `SmallClusterFull` is DISCHARGED by LRCSevenGapRigidity's core theorem
  (`goodSet_ae_full`: sorted-phase gap bridge + k = 7 perfect-net rigidity +
  countable-fiber null set), and the grand assembly drops to TWO citations:

      LRC14Statement  ⟸  [THM-661 moment floor] + [certificate supply].

  Kernel-pure: no native_decide, no sorry.
-/
import Mathlib
import TournamentH7.LRCSevenGapRigidity
import TournamentH7.LRCReachDecitation

namespace LonelyRunner
namespace LRC14
namespace TwoCitationAssembly

open MomentCitation SevenGapRigidity

/-- **The ≤7-arcs pigeonhole citation, DISCHARGED**: `SmallClusterFull`
is a theorem. -/
theorem smallClusterFull_proved : SmallClusterFull :=
  fun E hnd h3 h7 => goodSet_ae_full E hnd h3 h7

/-- **LRC(14) = [THM-661 moment floor] + [certificate supply].**  Down from
three citations to two. -/
theorem lrc14_from_moment_and_supply
    (h661 : THM661MomentFloor)
    (hcerts : ReachDecitation.THM527ACertificateSupply) :
    LRC14Statement :=
  ReachDecitation.lrc14_from_certificate_citations h661
    smallClusterFull_proved hcerts

end TwoCitationAssembly
end LRC14
end LonelyRunner

-- kernel-purity audit (fleet convention)
#print axioms LonelyRunner.LRC14.TwoCitationAssembly.smallClusterFull_proved
#print axioms LonelyRunner.LRC14.TwoCitationAssembly.lrc14_from_moment_and_supply
