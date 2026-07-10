/-
  TournamentH7.LRCWitnessFloorRepair — the repaired witness-floor assembly
  (mac-mini-2026-07-09-S65 cont.17; implements the cont.16 soundness-flag repair).

  The engine (lrc14_interval_measure_engine_macmini_S65cont16) proved the skeleton's original
  `hfloor`/`hsmall` UNSATISFIABLE as stated: `v = (1,…,13)` has `clusterSize (shapeOf v) = 0`
  and `witnessG2 = slowμ(safeSet {1..13}) = 0 < witnessMP` (the AP's lonely set is one point).
  Exact failure boundary: the m_P floor fails precisely for `clusterSize ≤ 2` — THM-530's
  admissibility (`k ≥ 3`), dropped by the skeleton statement.

  THE REPAIR (this file, additive — no existing statement edited):
  * `k = 0` ⟹ every speed has `|v i| ≤ 13` ⟹ 14 divides no speed ⟹ the q = 14 sieve fires
    (`LonelyRunner.sieve_one_div`) — no floor needed at all.
  * `k ∈ {1, 2}` ⟹ only POSITIVITY of `witnessG2` is needed (to feed `hpartA`); the engine's
    exact floors are `7/858` and `313/9702`.
  * `3 ≤ k ≤ 7` and `8 ≤ k ≤ 13` keep the `witnessMP = m_P` floor (canon-true; engine-certified).
  The assembly `lrc14_from_repaired_nodes` derives `LRC14Statement` from the four repaired
  legs + `hpartA` — every hypothesis now SATISFIABLE, with the engine's table as certificate data.
-/
import Mathlib
import TournamentH7.LRCMomentFloorConcrete

namespace LonelyRunner
namespace LRC14
namespace Repair

open MomentFloor

/-- `clusterSize (shapeOf v) = 0` means the `> 13` filter of the speed list is empty, so every
speed satisfies `|v i| ≤ 13`. -/
theorem speeds_le_of_clusterSize_zero (v : Fin 13 → ℤ)
    (h0 : clusterSize (shapeOf v) = 0) (i : Fin 13) : |v i| ≤ 13 := by
  by_contra hgt
  push_neg at hgt
  have hmem : |v i| ∈ List.ofFn (fun j => |v j|) := List.mem_ofFn.mpr ⟨i, rfl⟩
  have hmemf : |v i| ∈ (List.ofFn fun j => |v j|).filter (fun a => decide (13 < a)) :=
    List.mem_filter.mpr ⟨hmem, by simpa using hgt⟩
  have hlen : (((List.ofFn fun j => |v j|).filter (fun a => decide (13 < a))).map
      (fun a => (List.ofFn fun j => |v j|).foldr max 0 - a)).length = 0 := by
    simpa [clusterSize, shapeOf] using h0
  rw [List.length_map] at hlen
  exact List.ne_nil_of_mem hmemf (List.eq_nil_of_length_eq_zero hlen)

/-- **The k = 0 leg: the q = 14 sieve.**  All speeds `≤ 13` in absolute value, hence 14 divides
none of them, hence `t = 1/14` is a lonely time.  (This is exactly why the AP `{1..13}` — the
shape that FALSIFIED the original `hfloor` — is nonetheless lonely: it is non-covering at 14.) -/
theorem lonely_of_clusterSize_zero (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (h0 : clusterSize (shapeOf v) = 0) : ∃ t : ℝ, Lonely 14 v t := by
  refine ⟨1 / 14, sieve_one_div 14 14 v le_rfl (by norm_num) ?_⟩
  intro i hdvd
  have habs : (14 : ℤ) ∣ |v i| := (dvd_abs _ _).mpr hdvd
  have hle : (14 : ℤ) ≤ |v i| := Int.le_of_dvd (abs_pos.mpr (hv i)) habs
  have := speeds_le_of_clusterSize_zero v h0 i
  omega

/-- **The repaired witness-floor assembly.**  LRC(14) from SATISFIABLE legs:
positivity at `k ∈ {1,2}` (engine floors `7/858`, `313/9702`), the m_P floor on `3 ≤ k ≤ 7`
and `8 ≤ k ≤ 13` (canon/THM-530 admissible range), the reach node `hpartA`, and the internal
`k = 0` sieve leg.  Replaces the unsatisfiable original `hfloor` quantification. -/
theorem lrc14_from_repaired_nodes
    (hk12 : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → 1 ≤ clusterSize (shapeOf v) →
      clusterSize (shapeOf v) ≤ 2 → 0 < witnessG2 (shapeOf v))
    (hsmall3 : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → 3 ≤ clusterSize (shapeOf v) →
      clusterSize (shapeOf v) ≤ 7 → (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hlarge : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → 8 ≤ clusterSize (shapeOf v) →
      clusterSize (shapeOf v) ≤ 13 → (witnessMP : ℝ) ≤ witnessG2 (shapeOf v))
    (hpartA : ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) →
      (1 : ℝ) / 14 ≤ Mreach v) :
    LRC14Statement := by
  intro v hv
  rcases Nat.eq_zero_or_pos (clusterSize (shapeOf v)) with h0 | h1
  · exact lonely_of_clusterSize_zero v hv h0
  · have hMP : (0 : ℝ) < (witnessMP : ℝ) := by
      have : (0 : ℚ) < witnessMP := by norm_num [witnessMP]
      exact_mod_cast this
    have hpos : 0 < witnessG2 (shapeOf v) := by
      rcases Nat.lt_or_ge (clusterSize (shapeOf v)) 3 with h2 | h3
      · exact hk12 v hv h1 (Nat.lt_succ_iff.mp h2)
      · rcases Nat.lt_or_ge (clusterSize (shapeOf v)) 8 with h7 | h8
        · exact lt_of_lt_of_le hMP (hsmall3 v hv h3 (Nat.lt_succ_iff.mp h7))
        · exact lt_of_lt_of_le hMP (hlarge v hv h8 (clusterSize_shapeOf_le v))
    exact lonely_of_Mreach_ge v hv (hpartA v hpos)

end Repair
end LRC14
end LonelyRunner
