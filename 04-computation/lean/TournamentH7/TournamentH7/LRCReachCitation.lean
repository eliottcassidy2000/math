/-
  TournamentH7.LRCReachCitation — the THM-527-A reach citation and the FULLY
  CITATION-CLOSED LRC(14) assembly (mac-mini-2026-07-09-S65 cont.26).

  The reach node `hpartA` was the moment route's last non-citation hypothesis.
  This file splits it into PROVED GLUE + ONE NAMED CITATION:

  PROVED here:
  * `nonzero_of_witnessG2_pos` — positivity of the witness measure FORCES all
    speeds nonzero (a zero speed puts `0` in the small part and empties
    `safeSet`): the `hv` guard is derivable, so `hpartA`'s guard-free
    quantification is sound (the MISTAKE-136 probe PASSES on this node).
  * `exists_config_of_witnessG2_pos` — positivity yields an explicit slow time
    in `goodSet ∩ safeSet` (the configuration the slow-fast embedding consumes).
  * `Mreach_ge_of_minReach` — a single real time with `minReach ≥ 1/14` bounds
    the sup: the compactness half of the reach is pure glue.

  CITED (`THM527ARulerEmbedding`): the slow-fast ruler embedding — THM-527
  Part A (canon, PROVED: the reformulation in the limit with the O(1/Vmax)
  finite-Vmax correction) TOGETHER WITH the finite-Vmax closure of the
  sub-threshold families (klein's THM-686 window census + THM-687/688 two- and
  multi-scale limits + THM-693's constructive witnesses + the banks): positive
  witness measure for the shape yields a real time with clearance ≥ 1/14 for
  the actual 13-speed system.  Like `LRC(≤13)` and `THM661MomentFloor`, this
  PROVED-classically bundle enters as a named citation hypothesis.

  **`lrc14_from_citations_only`: LRC(14) from THREE NAMED CITATIONS — THM-661,
  the ≤7-arcs pigeonhole, and the THM-527-A reach — and NOTHING ELSE.  The
  entire remaining Lean surface of LRC(14) is citation-shaped.**
-/
import Mathlib
import TournamentH7.LRCSmallDischarge

set_option maxHeartbeats 800000

namespace LonelyRunner
namespace LRC14
namespace ReachCitation

open DenseCovers MeasureTheory MomentCitation TournamentH7.GoodSet

/-- A zero speed empties the safe set: `0` lands in the small part, and
`fract 0 = 0 ∉ [1/14, 13/14]`. -/
theorem safeSet_empty_of_zero_mem {P : List ℤ} (h0 : (0 : ℤ) ∈ P) :
    safeSet P = ∅ := by
  ext x
  simp only [Set.mem_empty_iff_false, iff_false]
  intro hx
  have := hx 0 h0
  simp only [Int.cast_zero, zero_mul, Int.fract_zero] at this
  have h114 : (1 : ℝ) / 14 ≤ 0 := this.1
  linarith

/-- **The guard is derivable**: positive witness measure forces all speeds nonzero. -/
theorem nonzero_of_witnessG2_pos (v : Fin 13 → ℤ)
    (hpos : 0 < witnessG2 (shapeOf v)) : ∀ i, v i ≠ 0 := by
  intro i hi
  have h0P : (0 : ℤ) ∈ (shapeOf v).1 := by
    have hmem : (0 : ℤ) ∈ (List.ofFn fun j => |v j|) := by
      rw [List.mem_ofFn]
      exact ⟨i, by simp [hi]⟩
    have : (0 : ℤ) ∈ (List.ofFn fun j => |v j|).filter (fun a => a ≤ 13) := by
      rw [List.mem_filter]
      exact ⟨hmem, by norm_num⟩
    exact this
  have hempty : safeSet (shapeOf v).1 = ∅ := safeSet_empty_of_zero_mem h0P
  have hzero : witnessG2 (shapeOf v) = 0 := by
    unfold witnessG2
    rw [hempty, Set.inter_empty]
    simp
  rw [hzero] at hpos
  exact lt_irrefl 0 hpos

/-- **Positivity yields a configuration point**: an explicit slow time in
`goodSet ∩ safeSet` — the input the slow-fast ruler embedding consumes. -/
theorem exists_config_of_witnessG2_pos (s : Shape)
    (hpos : 0 < witnessG2 s) :
    ∃ x : ℝ, x ∈ goodSet s.2 ∩ safeSet s.1 := by
  have hne : slowμ (goodSet s.2 ∩ safeSet s.1) ≠ 0 := by
    intro h
    have : witnessG2 s = 0 := by
      unfold witnessG2
      rw [h]
      simp
    rw [this] at hpos
    exact lt_irrefl 0 hpos
  exact MeasureTheory.nonempty_of_measure_ne_zero hne

/-- **The compactness half of the reach is glue**: one real time with
`minReach ≥ 1/14` bounds the supremum. -/
theorem Mreach_ge_of_minReach (v : Fin 13 → ℤ) (τ : ℝ)
    (h : (1 : ℝ) / 14 ≤ LRC14Concrete.minReach v τ) :
    (1 : ℝ) / 14 ≤ Mreach v := by
  have hbdd : BddAbove (Set.range (LRC14Concrete.minReach v)) := by
    refine ⟨1 / 2, ?_⟩
    rintro y ⟨t, rfl⟩
    exact LRC14Concrete.minReach_le_half v t
  have hmem : LRC14Concrete.minReach v τ ∈ Set.range (LRC14Concrete.minReach v) :=
    ⟨τ, rfl⟩
  have hle : LRC14Concrete.minReach v τ ≤ sSup (Set.range (LRC14Concrete.minReach v)) :=
    le_csSup hbdd hmem
  show (1 : ℝ) / 14 ≤ LRC14Concrete.Mreach v
  rw [LRC14Concrete.Mreach_eq_global_sSup v]
  linarith

/-- **THM-527 Part A, the ruler embedding (citation).**  Positive witness measure for
the shape yields a real time with clearance `≥ 1/14` for the actual system: the
slow-fast change of variables realizes a good slow time as a good `Vmax`-ruler period
(canon THM-527-A, PROVED in the limit with the `O(1/Vmax)` correction), and the
sub-threshold finite-`Vmax` families are closed by the census/banks arc
(THM-686/687/688/693).  Enters as a named citation per project policy. -/
def THM527ARulerEmbedding : Prop :=
  ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → 0 < witnessG2 (shapeOf v) →
    ∃ τ : ℝ, (1 : ℝ) / 14 ≤ LRC14Concrete.minReach v τ

/-- **`hpartA` from the citation**: derive the guard, embed, take the sup. -/
theorem hpartA_of_rulerEmbedding (h527 : THM527ARulerEmbedding) :
    ∀ v : Fin 13 → ℤ, 0 < witnessG2 (shapeOf v) → (1 : ℝ) / 14 ≤ Mreach v := by
  intro v hpos
  have hv := nonzero_of_witnessG2_pos v hpos
  obtain ⟨τ, hτ⟩ := h527 v hv hpos
  exact Mreach_ge_of_minReach v τ hτ

/-- **LRC(14) FROM CITATIONS ONLY.**  The complete assembly: THM-661 (the unified
covering-moment density floor), the ≤7-arcs pigeonhole, and the THM-527-A ruler
embedding — three named citations, all PROVED classically in canon — imply the
Lonely Runner Conjecture for 14 runners.  Every other node in the tree is a
theorem over `[propext, Classical.choice, Quot.sound]`. -/
theorem lrc14_from_citations_only
    (h661 : THM661MomentFloor) (hsmall7 : SmallClusterFull)
    (h527 : THM527ARulerEmbedding) :
    LRC14Statement :=
  SmallDischarge.lrc14_from_two_citations h661 hsmall7
    (hpartA_of_rulerEmbedding h527)

end ReachCitation
end LRC14
end LonelyRunner
