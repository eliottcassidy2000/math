/-
  TournamentH7.LRCCoverBoxes -- assignment-box skeleton for the LRC(14)
  cover/decorrelation route (codex-2026-06-22-S87c).

  `coverSet E` says every one of the six inner sectors is hit by some speed in
  `E`.  The first formal step toward a Vitali/union-bound treatment is to
  choose one witnessing speed per sector.  This file packages that step as a
  union over finite sector-to-speed assignment boxes:

      coverSet E ⊆ ⋃ σ, sectorBox E σ.

  The final theorem is the corresponding outer-measure union bound for the
  slow-time measure.  It is deliberately loose: it does not prove independence,
  decorrelation, or the quasi-independence refinement needed for the tight
  margin.
-/

import TournamentH7.LRCDenseCovers

namespace LonelyRunner
namespace CoverBoxes

open MeasureTheory

/-- A choice of one speed from `E` for each of the six inner sectors. -/
abbrev Assignment (E : List ℤ) := Fin 6 → {e : ℤ // e ∈ E.toFinset}

/-- The box where assignment `σ` puts its `j`th chosen speed in inner sector
`[(j+1)/7,(j+2)/7)`, for `j : Fin 6`. -/
def sectorBox (E : List ℤ) (σ : Assignment E) : Set ℝ :=
  {x | ∀ j : Fin 6,
    ((((j : ℕ) + 1 : ℕ) : ℝ) / 7 ≤ Int.fract (((σ j).1 : ℝ) * x)) ∧
    Int.fract (((σ j).1 : ℝ) * x) <
      (((((j : ℕ) + 1 : ℕ) : ℝ) + 1) / 7)}

/-- `coverSet` is covered by sector-assignment boxes.  This is the formal
choice-of-witnesses skeleton behind the later union/Vitali cover. -/
theorem coverSet_subset_iUnion_sectorBox (E : List ℤ) :
    DenseCovers.coverSet E ⊆ ⋃ σ : Assignment E, sectorBox E σ := by
  classical
  intro x hx
  have hw : ∀ j : Fin 6, ∃ e : {e : ℤ // e ∈ E.toFinset},
      ((((j : ℕ) + 1 : ℕ) : ℝ) / 7 ≤ Int.fract (((e : ℤ) : ℝ) * x)) ∧
      Int.fract (((e : ℤ) : ℝ) * x) <
        (((((j : ℕ) + 1 : ℕ) : ℝ) + 1) / 7) := by
    intro j
    obtain ⟨e, heE, hle, hlt⟩ :=
      hx (((j : ℕ) + 1 : ℕ)) (by omega) (by have := j.isLt; omega)
    exact ⟨⟨e, List.mem_toFinset.mpr heE⟩, hle, hlt⟩
  choose σ hσ using hw
  exact Set.mem_iUnion.mpr ⟨σ, fun j => hσ j⟩

/-- The loose countable-union measure bound for the assignment-box cover.  The
tight HYP-2840 step is to replace this lossy union ledger by a Vitali/overlap
or decorrelation refinement. -/
theorem slowμ_coverSet_le_tsum_sectorBox (E : List ℤ) :
    DenseCovers.slowμ (DenseCovers.coverSet E) ≤
      ∑' σ : Assignment E, DenseCovers.slowμ (sectorBox E σ) := by
  calc
    DenseCovers.slowμ (DenseCovers.coverSet E)
        ≤ DenseCovers.slowμ (⋃ σ : Assignment E, sectorBox E σ) :=
      measure_mono (coverSet_subset_iUnion_sectorBox E)
    _ ≤ ∑' σ : Assignment E, DenseCovers.slowμ (sectorBox E σ) :=
      measure_iUnion_le _

/-! ## Axiom audit -/

#print axioms coverSet_subset_iUnion_sectorBox
#print axioms slowμ_coverSet_le_tsum_sectorBox

end CoverBoxes
end LonelyRunner
