/-
  TournamentH7.LRCParityDemo -- the parity pairing law, end-to-end on a concrete
  covering instance (klein-2026-07-09-S228; the S227 law's operational demo).

  Instance: klein-S206's worst covering 13-set {1,2,3,4,7,8,9,10,11,12,13,14,17}
  at modulus q = 21.  Its live multipliers are EXACTLY {4, 17} = {p, q − p}:
  one ± pair, LM = 2 — the smallest visible instance of the pairing structure.

  The demo shows the law working as machinery, not just as a statement:
   * ONE kernel-decide certificate (p = 4 live);
   * THE TWIN EXTRACTED BY THE LAW (live_mirror gives 17 = 21 − 4 live with NO
     recomputation — the certificate-search halving, formal);
   * the LM-EVEN validation invariant, twice: by kernel decide AND by
     liveCount_even_of_even_speed (the two routes agreeing is the bug-detector
     pattern recommended for the enumeration banks);
   * BOTH witnesses driven through kps's dispatch to Mreach ≥ 1/14.

  Kernel-pure: decide only (no native_decide), no sorry.
-/

import Mathlib
import TournamentH7.LRCParityPairing
import TournamentH7.LRCPairSumDispatch

namespace LonelyRunner
namespace LRC14Concrete

/-- klein-S206's worst covering 13-set (min good-period margin over the 966). -/
def demoSpeeds : Fin 13 → ℤ := ![1, 2, 3, 4, 7, 8, 9, 10, 11, 12, 13, 14, 17]

/-- The kernel-decide certificate: `p = 4` is live at `q = 21`. -/
theorem demo_live_4 : bandCount demoSpeeds 21 4 = 0 := by decide

/-- **THE TWIN, BY THE LAW**: `17 = 21 − 4` is live with no recomputation —
`live_mirror` extracts it from `demo_live_4`. -/
theorem demo_live_17_by_law : bandCount demoSpeeds 21 17 = 0 :=
  live_mirror demoSpeeds 21 4 (by norm_num) (by norm_num) demo_live_4

/-- Cross-check: the twin also verifies by kernel decide (the two routes agree). -/
theorem demo_live_17_by_decide : bandCount demoSpeeds 21 17 = 0 := by decide

/-- The validation invariant by kernel decide: the live count is even. -/
theorem demo_LM_even_by_decide : liveCount demoSpeeds 21 % 2 = 0 := by decide

/-- The validation invariant BY THE LAW: the family has an even speed
(`demoSpeeds 1 = 2`), so the live count is even at EVERY modulus — here `21`. -/
theorem demo_LM_even_by_law : liveCount demoSpeeds 21 % 2 = 0 :=
  liveCount_even_of_even_speed demoSpeeds 21 ⟨1, by decide⟩

/-- The certificate drives to loneliness: `Mreach ≥ 1/14` via the dispatch. -/
theorem demo_mreach : (1 : ℝ) / 14 ≤ Mreach demoSpeeds :=
  mreach_ge_of_pairsum_band demoSpeeds 4 21 (by norm_num) (by decide)

/-- The FREE twin drives to the same conclusion — one decide, two witnesses. -/
theorem demo_mreach_twin : (1 : ℝ) / 14 ≤ Mreach demoSpeeds :=
  mreach_ge_of_pairsum_band demoSpeeds 17 21 (by norm_num) (by decide)

end LRC14Concrete
end LonelyRunner

-- kernel-purity audit (fleet convention)
#print axioms LonelyRunner.LRC14Concrete.demo_live_4
#print axioms LonelyRunner.LRC14Concrete.demo_live_17_by_law
#print axioms LonelyRunner.LRC14Concrete.demo_LM_even_by_law
#print axioms LonelyRunner.LRC14Concrete.demo_mreach
#print axioms LonelyRunner.LRC14Concrete.demo_mreach_twin
