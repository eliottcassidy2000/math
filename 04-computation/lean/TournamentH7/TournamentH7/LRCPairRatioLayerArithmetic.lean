/-
  Kernel arithmetic for the weighted negative-pair ratio-layer certificate.
  Both the sharp planar-grid cap and the simpler two-path cap remain strictly
  below 13/30. The graph and continuum-transfer premises live in later rungs.
  No `sorry`; no `native_decide`.
-/

import Mathlib

namespace LonelyRunner
namespace LRCPairRatioLayerArithmetic

def negativePairTierBound : ℚ :=
    78 * (1 / 14112)
    + 73 * (1 / 9996 - 1 / 14112)
    + 72 * (5 / 37632 - 1 / 9996)
    + 70 * (1 / 3136 - 5 / 37632)
    + 67 * (1 / 1764 - 1 / 3136)
    + 63 * (5 / 2646 - 1 / 1764)
    + 56 * (2 / 441 - 5 / 2646)
    + 42 * (4 / 539 - 2 / 441)
    + 22 * (5 / 588 - 4 / 539)
    + 12 * (6 / 637 - 5 / 588)

theorem negativePairTierBound_eq :
    negativePairTierBound = 175847381 / 411675264 := by
  norm_num [negativePairTierBound]

theorem negativePairTierBound_lt_target :
    negativePairTierBound < 13 / 30 := by
  norm_num [negativePairTierBound]

theorem negativePairTier_margin :
    13 / 30 - negativePairTierBound = 12726167 / 2058376320 := by
  norm_num [negativePairTierBound]

/-- The simplest formal tier bound uses only that each of the two top ratio
colors (`12` and `13`) is a path forest, hence contributes at most `12` edges.
It avoids the optional planar-grid sharpening from `24` to `22`. -/
def negativePairTierBoundPathOnly : ℚ :=
    78 * (1 / 14112)
    + 73 * (1 / 9996 - 1 / 14112)
    + 72 * (5 / 37632 - 1 / 9996)
    + 70 * (1 / 3136 - 5 / 37632)
    + 67 * (1 / 1764 - 1 / 3136)
    + 63 * (5 / 2646 - 1 / 1764)
    + 56 * (2 / 441 - 5 / 2646)
    + 42 * (4 / 539 - 2 / 441)
    + 24 * (5 / 588 - 4 / 539)
    + 12 * (6 / 637 - 5 / 588)

theorem negativePairTierBoundPathOnly_eq :
    negativePairTierBoundPathOnly = 176738453 / 411675264 := by
  norm_num [negativePairTierBoundPathOnly]

theorem negativePairTierBoundPathOnly_margin :
    13 / 30 - negativePairTierBoundPathOnly = 8270807 / 2058376320 := by
  norm_num [negativePairTierBoundPathOnly]

theorem negativePairTierBoundPathOnly_lt_target :
    negativePairTierBoundPathOnly < 13 / 30 := by
  norm_num [negativePairTierBoundPathOnly]

#print axioms negativePairTierBound_eq
#print axioms negativePairTierBound_lt_target
#print axioms negativePairTier_margin
#print axioms negativePairTierBoundPathOnly_lt_target

end LRCPairRatioLayerArithmetic
end LonelyRunner

