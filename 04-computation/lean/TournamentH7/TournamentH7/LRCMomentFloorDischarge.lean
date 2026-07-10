/-
  TournamentH7.LRCMomentFloorDischarge — witnessG2-concretization plumbing, first leg
  (mac-mini-2026-07-09-S65 cont.14).

  death-star-S4 de-opaqued the skeleton (`witnessG2` and `shapeOf` are now concrete `def`s);
  opus-S186cont2 recorded that under the concrete definitions the four legs
  `hbonf, hB, hsmall, hsize` of `lrc14_from_momentfloor_nodes` become proof terms.  This file
  executes the plumbing for the first of them:

  * `clusterSize_shapeOf_le` — **hsize DISCHARGED**: the cluster of `shapeOf v` is the
    filter-map of a 13-entry list, so its length is ≤ 13.  Pure list arithmetic.
  * `lrc14_from_momentfloor_nodes'` — the assembly re-exported with `hsize` supplied
    internally: LRC(14) from FIVE parameters {hbonf, hMoment, hB, hsmall, hpartA}.

  Remaining legs (their target proof terms, per opus-S186cont2, now unblocked):
  `hbonf` ← `LRCBonferroniMeasure.toReal_bonferroni` (needs the nuShape/measGP measure
  identities against `slowμ (GOOD ∩ safe)`); `hsmall` ← the k ≤ 7 pigeonhole (`GOOD = univ`,
  so `witnessG2 = measGP ≥ m_P` — the `LRCWitnessG2Discharge` instances are the base);
  `hB` ← Lemma B on the concrete `safeSet`.  Each is a measure-identity bridge, not analysis.
-/
import Mathlib
import TournamentH7.LRCWitnessMomentFloor

namespace LonelyRunner
namespace LRC14

open MomentFloor

/-- **hsize, discharged against the concrete `shapeOf`.**  The cluster component of
`shapeOf v` is `(speeds.filter (13 < ·)).map (Vmax − ·)` with `speeds` a 13-entry list, so its
length is at most 13. -/
theorem clusterSize_shapeOf_le (v : Fin 13 → ℤ) : clusterSize (shapeOf v) ≤ 13 := by
  simp only [clusterSize, shapeOf]
  calc (((List.ofFn fun i => |v i|).filter
          (fun a => decide (13 < a))).map
          (fun a => (List.ofFn fun i => |v i|).foldr max 0 - a)).length
      = ((List.ofFn fun i => |v i|).filter (fun a => decide (13 < a))).length :=
        List.length_map ..
    _ ≤ (List.ofFn fun i => |v i|).length := List.length_filter_le _ _
    _ = 13 := by simp

/-- The moment-floor assembly with `hsize` supplied internally: **LRC(14) from five
parameters** `{hbonf, hMoment, hB, hsmall, hpartA}`. -/
theorem lrc14_from_momentfloor_nodes'
    (nuShape measGP : Shape → ℝ)
    (hbonf : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)
    (hMoment : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (momentBar (clusterSize s) : ℝ) ≤ nuShape s)
    (hB : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (capRat (clusterSize s) : ℝ) ≤ measGP s)
    (hsmall : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 7 → (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement :=
  lrc14_from_momentfloor_nodes nuShape measGP hbonf hMoment hB hsmall
    (fun v _ => clusterSize_shapeOf_le v) hpartA

end LRC14
end LonelyRunner
