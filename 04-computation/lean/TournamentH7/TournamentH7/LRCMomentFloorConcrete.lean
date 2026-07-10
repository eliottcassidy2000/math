/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-09-S190)
-/
import Mathlib.Tactic
import TournamentH7.LRCWitnessMomentFloor
import TournamentH7.LRCBonferroniMeasure

/-!
# Discharging the moment-floor legs against the CONCRETE `witnessG2` (opus-S190)

death-star-S4 de-opaqued the skeleton: `witnessG2 s = (slowμ (goodSet s.2 ∩ safeSet s.1)).toReal`,
`shapeOf` concrete. This unblocks the two OPACITY-blocked legs of `lrc14_from_momentfloor_nodes`
(opus-S186), which are now provable proof terms:

* **`hbonf`** — the Bonferroni bound `μ(GOOD) + μ(G_P) − 1 ≤ μ(GOOD ∩ G_P)`, from
  `BonferroniMeasure.toReal_bonferroni` (kps-S30) on the probability measure `slowμ` with the measurable
  small-part event `safeSet`;
* **`hsize`** — `clusterSize (shapeOf v) ≤ 13`, since the cluster is a filter of the 13 absolute speeds.

`lrc14_from_momentfloor_concrete` instantiates the S186 node with the concrete measures
`nuShape = (slowμ (goodSet)).toReal`, `measGP = (slowμ (safeSet)).toReal` and DISCHARGES `hbonf`, `hsize`.
The open surface shrinks to the genuine analytic content: `hMoment` (the density floor `μ(GOOD) ≥
momentBar`, = THM-661; the pure-cluster/diam≤75 case is death-star's `GoodSetBridge`), `hB` (Lemma B,
`μ(G_P) ≥ cap`), `hsmall` (`k ≤ 7` pigeonhole), and `hpartA` (the reach).

Kernel-pure: no `sorry`, no `native_decide`. Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LonelyRunner.LRC14.MomentFloor

open LonelyRunner.LRC14 MeasureTheory

/-- The concrete good-set density `nuShape s = (slowμ (goodSet s.2)).toReal`. -/
noncomputable def nuShapeConcrete (s : Shape) : ℝ :=
  (DenseCovers.slowμ (TournamentH7.GoodSet.goodSet s.2)).toReal

/-- The concrete small-part safe density `measGP s = (slowμ (safeSet s.1)).toReal`. -/
noncomputable def measGPConcrete (s : Shape) : ℝ :=
  (DenseCovers.slowμ (DenseCovers.safeSet s.1)).toReal

/-- **`hsize` DISCHARGED.** The cluster co-offset list is a filter of the 13 absolute speeds, so it has
length ≤ 13. -/
theorem clusterSize_shapeOf_le (v : Fin 13 → ℤ) : clusterSize (shapeOf v) ≤ 13 := by
  simp only [clusterSize, shapeOf, List.length_map]
  exact le_trans (List.length_filter_le _ _) (by simp)

/-- **`hbonf` DISCHARGED.** The Bonferroni bound against the concrete `witnessG2`: on the probability
measure `slowμ`, `μ(GOOD).toReal + μ(G_P).toReal − 1 ≤ μ(GOOD ∩ G_P).toReal = witnessG2`. -/
theorem bonferroni_concrete (s : Shape) :
    nuShapeConcrete s + measGPConcrete s - 1 ≤ witnessG2 s :=
  LonelyRunner.BonferroniMeasure.toReal_bonferroni DenseCovers.slowμ
    (TournamentH7.GoodSet.goodSet s.2) (DenseCovers.safeSet s.1)
    (DenseCovers.measurableSet_safeSet s.1)

/-- **The moment-floor route with the concrete legs discharged (opus-S190).**
LRC(14) from the four remaining ANALYTIC obligations — `hMoment` (density floor `μ(GOOD) ≥ momentBar`,
THM-661), `hB` (Lemma B `μ(G_P) ≥ cap`), `hsmall` (`k ≤ 7` pigeonhole), `hpartA` (reach) — stated on the
CONCRETE measures. `hbonf` (Bonferroni) and `hsize` (cluster length ≤ 13) are supplied as proof terms;
`hR0` was already discharged in `lrc14_from_momentfloor_nodes`. -/
theorem lrc14_from_momentfloor_concrete
    (hMoment : ∀ s : Shape, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (momentBar (clusterSize s) : ℝ) ≤ nuShapeConcrete s)
    (hB : ∀ s : Shape, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (capRat (clusterSize s) : ℝ) ≤ measGPConcrete s)
    (hsmall : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 7 → (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement :=
  lrc14_from_momentfloor_nodes nuShapeConcrete measGPConcrete
    bonferroni_concrete hMoment hB hsmall
    (fun v _ => clusterSize_shapeOf_le v) hpartA

end LonelyRunner.LRC14.MomentFloor
