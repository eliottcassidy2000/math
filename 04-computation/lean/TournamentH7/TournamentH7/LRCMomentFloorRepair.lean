/-
  TournamentH7.LRCMomentFloorRepair — the hB soundness repair + the safe-only certificate
  consumer (mac-mini-2026-07-09-S65 cont.22).

  SOUNDNESS FINDING (same genus as the cont.17 hfloor repair): the `hB` node of
  `lrc14_from_momentfloor_concrete` (opus-S190) is UNSATISFIABLE as stated — it quantifies
  over ALL `s : Shape`, and `s = ([0], List.replicate 13 0)` has `clusterSize s = 13`,
  `capRat 13 = 1`, but `measGPConcrete s = (slowμ (safeSet [0])).toReal = 0` (the speed `0`
  has `fract 0 = 0 ∉ [1/14, 13/14]`, so `safeSet [0] = ∅`).  On REACHABLE shapes
  (`shapeOf v`, nonzero speeds) the node is fine: `|P| = 13 − k ≤ 5` there, and the capRat
  ladder is EXACTLY the per-|P| safe-measure minima (engine-verified, all six rows equal:
  min over |S|=5 is 2243/5880 at {1,5,7,8,9}, |S|=4: 1979/4004 at {1,11,12,13},
  |S|=3: 55/91 at {1,12,13}, |S|=2: 66/91 at {1,13}, |S|=1: 6/7 uniform, |S|=0: 1).

  This file: (1) `lrc14_from_momentfloor_concrete_shapes` — the repaired consumer with `hB`
  narrowed to `shapeOf v` (hMoment left shape-quantified: its E-uniform scope is THM-661's
  claim; no falsity found there — but formalizers should be aware the same narrowing is
  available if needed).  (2) `safeSet_anti` — safeSet is antitone in list membership (the
  dedup/superset bridge for the certificate table).  (3) `measGP_ge_of_sorted_bands` — the
  safe-only sorted-interval floor: engine component lists ⟹ `Σ lengths ≤ measGPConcrete`,
  consuming brick (iii); this is the consumer for the per-family capRat certificates
  (LRCSafeCertSize*.lean).
-/
import Mathlib
import TournamentH7.LRCMomentFloorConcrete
import TournamentH7.LRCIntervalBridge
import TournamentH7.LRCUnionVolume

namespace LonelyRunner
namespace LRC14
namespace MomentFloorRepair

open DenseCovers MeasureTheory MomentFloor

/-- **The repaired moment-floor consumer.**  Same route as
`lrc14_from_momentfloor_concrete`, with `hB` quantified over REACHABLE shapes only
(`shapeOf v`, nonzero speeds) — the original ∀-Shape `hB` is unsatisfiable
(`s = ([0], [0,…,0])`).  `hbonf`/`hsize` remain internally discharged. -/
theorem lrc14_from_momentfloor_concrete_shapes
    (hMoment : ∀ s : Shape, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
      (momentBar (clusterSize s) : ℝ) ≤ nuShapeConcrete s)
    (hB : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      8 ≤ clusterSize (shapeOf v) → clusterSize (shapeOf v) ≤ 13 →
      (capRat (clusterSize (shapeOf v)) : ℝ) ≤ measGPConcrete (shapeOf v))
    (hsmall : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) →
      clusterSize (shapeOf v) ≤ 7 → (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement := by
  apply lrc14_from_witness_floor
  · apply witness_floor_from_cluster_cases hsmall
    · -- the large node, with hB applied at the reachable shape only
      intro v hv h8 h13
      set k := clusterSize (shapeOf v) with hk
      have hrat : momentBar k + capRat k - 1 = witnessMP := momentBar_add_capRat k
      have hcast : (momentBar k : ℝ) + (capRat k : ℝ) - 1 = (witnessMP : ℝ) := by
        have := congrArg (fun q : ℚ => (q : ℝ)) hrat
        push_cast at this ⊢
        linarith
      have hMk : (momentBar k : ℝ) ≤ nuShapeConcrete (shapeOf v) := hMoment (shapeOf v) h8 h13
      have hBk : (capRat k : ℝ) ≤ measGPConcrete (shapeOf v) := hB v hv h8 h13
      have hbk : nuShapeConcrete (shapeOf v) + measGPConcrete (shapeOf v) - 1
          ≤ witnessG2 (shapeOf v) := bonferroni_concrete (shapeOf v)
      linarith
    · exact fun v _ => clusterSize_shapeOf_le v
  · exact hpartA

/-- `safeSet` is ANTITONE in list membership: more speeds, smaller safe set.  The
dedup/superset bridge: a shape's P-list certificate reduces to any membership-superset's. -/
theorem safeSet_anti {P P' : List ℤ} (h : ∀ p ∈ P, p ∈ P') :
    safeSet P' ⊆ safeSet P := by
  intro x hx p hp
  exact hx p (h p hp)

/-- Measure form of the antitone bridge. -/
theorem measGP_anti {P P' : List ℤ} (h : ∀ p ∈ P, p ∈ P') :
    (slowμ (safeSet P')).toReal ≤ (slowμ (safeSet P)).toReal := by
  apply ENNReal.toReal_mono (measure_ne_top _ _)
  exact measure_mono (safeSet_anti h)

/-- **The safe-only certificate consumer**: a sorted-disjoint engine list whose intervals
pass the per-speed band checks floors `(slowμ (safeSet P)).toReal` by `Σ lengths`.
This is `measGPConcrete` at reachable shapes — the exact form the capRat certificate
tables (LRCSafeCertSize*.lean) discharge. -/
theorem measGP_ge_of_sorted_bands (P : List ℤ) (l : List (ℝ × ℝ))
    (hposP : ∀ p ∈ P, (0 : ℤ) < p)
    (hcert : ∀ q ∈ l, ∀ p ∈ P, ∃ j : ℤ, (j : ℝ) + 1 / 14 ≤ (p : ℝ) * q.1 ∧
      (p : ℝ) * q.2 ≤ (j : ℝ) + 13 / 14)
    (hin : ∀ q ∈ l, 0 ≤ q.1 ∧ q.1 ≤ q.2 ∧ q.2 ≤ 1)
    (hsorted : l.Pairwise (fun q r => q.2 ≤ r.1)) :
    (l.map (fun q => q.2 - q.1)).sum ≤ (slowμ (safeSet P)).toReal := by
  have hsub : ∀ q ∈ l, Set.Ico q.1 q.2 ⊆ safeSet P := fun q hq =>
    IntervalBridge.Ico_subset_safeSet_of_bounds hposP (hcert q hq)
  have hfloor := UnionVolume.slowμ_ge_sum_of_sorted_Ico_subset l hsub hin hsorted
  have hnn : 0 ≤ (l.map (fun q => q.2 - q.1)).sum := by
    apply List.sum_nonneg
    intro y hy
    obtain ⟨q, hq, rfl⟩ := List.mem_map.mp hy
    have h := hin q hq
    linarith [h.2.1]
  calc (l.map (fun q => q.2 - q.1)).sum
      = (ENNReal.ofReal ((l.map (fun q => q.2 - q.1)).sum)).toReal :=
        (ENNReal.toReal_ofReal hnn).symm
    _ ≤ (slowμ (safeSet P)).toReal := ENNReal.toReal_mono (measure_ne_top _ _) hfloor

end MomentFloorRepair
end LRC14
end LonelyRunner
