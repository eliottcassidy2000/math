/-
  TournamentH7.LRCGoodSetSmall — brick (i) of the hk12 discharge chain
  (mac-mini-2026-07-09-S65 cont.19).

  `goodSet E = univ` for clusters of length ≤ 2: with one tooth the difference is 0 and
  `fract 0 = 0 ∉ (0, 1/7]`; with two teeth, if `fract(d·x) ∈ (0, 1/7]` then
  `fract(−d·x) = 1 − fract(d·x) ≥ 6/7 ∉ (0, 1/7]` — the two teeth cannot BOTH have occupied
  following-arcs, so one of them witnesses the union.  Composes with
  `LRCIntervalBridge.witnessG2_pos_of_anchor`: for `clusterSize ∈ {1,2}` the goodSet condition
  is free, and the engine's anchor table (lrc14_anchor_table_macmini_S65cont19) supplies the
  safeSet bounds — the hk12 leg of `lrc14_from_repaired_nodes` becomes anchors-in, proof-out.
-/
import Mathlib
import TournamentH7.LRCGoodSet

namespace LonelyRunner
namespace LRC14
namespace GoodSetSmall

open TournamentH7.GoodSet

/-- The pointwise sign flip: if `fract(d·x)` lies in `(0, 1/7]` then `fract(−d·x)` does not
(it equals `1 − fract(d·x) ≥ 6/7`). -/
theorem fract_neg_not_mem {d : ℤ} {x : ℝ}
    (h : Int.fract ((d : ℝ) * x) ∈ Set.Ioc (0 : ℝ) (1 / 7)) :
    Int.fract (((-d : ℤ) : ℝ) * x) ∉ Set.Ioc (0 : ℝ) (1 / 7) := by
  intro hmem
  have hne0 : Int.fract ((d : ℝ) * x) ≠ 0 := ne_of_gt h.1
  have hcast : (((-d : ℤ) : ℝ) * x) = -(((d : ℤ) : ℝ) * x) := by push_cast; ring
  rw [hcast, Int.fract_neg hne0] at hmem
  have h1 := h.2
  have h2 := hmem.2
  linarith

/-- The zero difference is never in the arc. -/
theorem fract_zero_not_mem (x : ℝ) :
    Int.fract (((0 : ℤ) : ℝ) * x) ∉ Set.Ioc (0 : ℝ) (1 / 7) := by
  simp

/-- **Brick (i): clusters of length ≤ 2 have `goodSet = univ`.** -/
theorem goodSet_univ_of_length_le_two (E : List ℤ) (hne : E ≠ []) (hlen : E.length ≤ 2) :
    goodSet E = Set.univ := by
  apply Set.eq_univ_of_forall
  intro x
  rcases E with _ | ⟨a, E'⟩
  · exact absurd rfl hne
  rcases E' with _ | ⟨b, E''⟩
  · -- E = [a]
    refine Set.mem_iUnion₂.mpr ⟨a, by simp, Set.mem_iInter₂.mpr ?_⟩
    intro c hc
    have hca : c = a := by simpa using hc
    subst hca
    simpa using fract_zero_not_mem x
  rcases E'' with _ | ⟨c, t⟩
  · -- E = [a, b]
    by_cases hcase : Int.fract (((b - a : ℤ) : ℝ) * x) ∈ Set.Ioc (0 : ℝ) (1 / 7)
    · -- witness b: the (a − b) difference flips out of the arc
      refine Set.mem_iUnion₂.mpr ⟨b, by simp, Set.mem_iInter₂.mpr ?_⟩
      intro d hd
      have hd' : d = a ∨ d = b := by simpa using hd
      rcases hd' with rfl | rfl
      · have : (d - b : ℤ) = -(b - d) := by ring
        rw [this]
        exact fract_neg_not_mem hcase
      · simpa using fract_zero_not_mem x
    · -- witness a
      refine Set.mem_iUnion₂.mpr ⟨a, by simp, Set.mem_iInter₂.mpr ?_⟩
      intro d hd
      have hd' : d = a ∨ d = b := by simpa using hd
      rcases hd' with rfl | rfl
      · simpa using fract_zero_not_mem x
      · exact hcase
  · -- length ≥ 3 contradicts hlen
    simp at hlen

end GoodSetSmall
end LRC14
end LonelyRunner
