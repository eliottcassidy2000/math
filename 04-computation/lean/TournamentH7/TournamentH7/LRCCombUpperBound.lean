/- LRCCombUpperBound.lean -- opus-2026-07-17-S351 (HYP-7340).
   THE SINGLE-COMB UPPER BOUND, sharp: `vol (badArcs w lam ∩ W) ≤ 2*lam`
   for a unit window `W`.

   The corpus's `fragmentation` bound is deliberately lossy at the edges
   (`(w*L + 2*lam + 1) * (2*lam/w)`, i.e. `2*lam + O(1/w)` at `L = 1`) — the
   `+1` counts a possible straddling arc at each end.  That extra `O(1/w)`
   is fatal here: seven combs at `lam = 1/14` must sum to `≤ 1` EXACTLY, with
   no slack to spare.

   THE FIX — a SHIFTED window.  Put the endpoints at half-cells,
   `W = Ico (1/(2w)) (1 + 1/(2w))`.  Arcs have radius `lam/w ≤ 1/(2w)`, so no
   arc straddles either endpoint: precisely the `w` arcs `j = 1..w` can meet
   `W`, each of measure `2*lam/w`, and the bound is `w · 2*lam/w = 2*lam` with
   a UNIFORM per-arc estimate and no boundary case analysis.

   Feeds `LRCSevenWallExistence.block_lonely_point_of_sum_le` (hypothesis
   `∑ μ (A i) ≤ 1`): seven combs at `lam = 1/14` give `7 · (1/7) = 1`.
   Kernel-pure: no sorry, no native_decide. -/
import Mathlib
import TournamentH7.FragmentationLemma

open MeasureTheory

namespace LonelyRunner.LRC14.Hunter

/-- Only the arcs `j = 1 … w` can meet the half-cell-shifted unit window. -/
theorem badArcs_window_subset {w : ℕ} {lam : ℝ} (hw : 0 < w) (hlam : 0 < lam)
    (hhalf : 2 * lam ≤ 1) :
    LRC14.badArcs w lam ∩ Set.Ico (1 / (2 * (w : ℝ))) (1 + 1 / (2 * (w : ℝ)))
      ⊆ ⋃ k ∈ Finset.range w,
          Set.Ioo ((((k : ℤ) + 1 : ℤ) : ℝ) / w - lam / w)
                  ((((k : ℤ) + 1 : ℤ) : ℝ) / w + lam / w) := by
  have hwR : (0 : ℝ) < w := by exact_mod_cast hw
  rintro t ⟨htA, htW⟩
  rw [LRC14.badArcs, Set.mem_iUnion] at htA
  obtain ⟨j, hj⟩ := htA
  rw [Set.mem_Ioo] at hj
  rw [Set.mem_Ico] at htW
  -- clear denominators: the window/arc inequalities times w
  have h2w : (0 : ℝ) < 2 * w := by linarith
  have e1 : t * w < (j : ℝ) + lam := by
    have h := hj.2
    have hrw : (j : ℝ) / w + lam / w = ((j : ℝ) + lam) / w := by ring
    rw [hrw, lt_div_iff₀ hwR] at h
    linarith
  have e2 : (1 : ℝ) / 2 ≤ t * w := by
    have h := htW.1
    rw [div_le_iff₀ h2w] at h
    linarith
  have e3 : (j : ℝ) - lam < t * w := by
    have h := hj.1
    have hrw : (j : ℝ) / w - lam / w = ((j : ℝ) - lam) / w := by ring
    rw [hrw, div_lt_iff₀ hwR] at h
    linarith
  have e4 : t * w < (w : ℝ) + 1 / 2 := by
    have h := htW.2
    have hrw : (1 : ℝ) + 1 / (2 * w) = ((w : ℝ) + 1 / 2) / w := by
      field_simp
    rw [hrw, lt_div_iff₀ hwR] at h
    linarith
  -- the arc index j is forced into [1, w]
  have hjlow : (0 : ℝ) < (j : ℝ) := by linarith
  have hjhigh : (j : ℝ) < (w : ℝ) + 1 := by linarith
  have hj1 : 1 ≤ j := by exact_mod_cast hjlow
  have hjw : j ≤ (w : ℤ) := by exact_mod_cast Int.lt_add_one_iff.mp (by exact_mod_cast hjhigh)
  -- reindex j as k + 1 with k < w
  have hidx : (j - 1).toNat ∈ Finset.range w := by
    rw [Finset.mem_range]
    omega
  refine Set.mem_biUnion hidx ?_
  · have hcast : (((((j - 1).toNat : ℤ) + 1 : ℤ)) : ℝ) = (j : ℝ) := by
      have : ((j - 1).toNat : ℤ) = j - 1 := Int.toNat_of_nonneg (by omega)
      rw [this]; push_cast; ring
    rw [Set.mem_Ioo, hcast]
    exact ⟨hj.1, hj.2⟩

/-- **THE SINGLE-COMB UPPER BOUND (sharp).** -/
theorem single_comb_le {w : ℕ} {lam : ℝ} (hw : 0 < w) (hlam : 0 < lam)
    (hhalf : 2 * lam ≤ 1) :
    volume (LRC14.badArcs w lam
        ∩ Set.Ico (1 / (2 * (w : ℝ))) (1 + 1 / (2 * (w : ℝ))))
      ≤ ENNReal.ofReal (2 * lam) := by
  have hwR : (0 : ℝ) < w := by exact_mod_cast hw
  have harc : ∀ k : ℕ,
      volume (Set.Ioo ((((k : ℤ) + 1 : ℤ) : ℝ) / w - lam / w)
                      ((((k : ℤ) + 1 : ℤ) : ℝ) / w + lam / w))
        = ENNReal.ofReal (2 * (lam / w)) := by
    intro k
    rw [Real.volume_Ioo]
    congr 1
    ring
  calc volume (LRC14.badArcs w lam
          ∩ Set.Ico (1 / (2 * (w : ℝ))) (1 + 1 / (2 * (w : ℝ))))
      ≤ volume (⋃ k ∈ Finset.range w,
          Set.Ioo ((((k : ℤ) + 1 : ℤ) : ℝ) / w - lam / w)
                  ((((k : ℤ) + 1 : ℤ) : ℝ) / w + lam / w)) :=
        measure_mono (badArcs_window_subset hw hlam hhalf)
    _ ≤ ∑ k ∈ Finset.range w,
          volume (Set.Ioo ((((k : ℤ) + 1 : ℤ) : ℝ) / w - lam / w)
                          ((((k : ℤ) + 1 : ℤ) : ℝ) / w + lam / w)) :=
        measure_biUnion_finset_le _ _
    _ = ∑ _k ∈ Finset.range w, ENNReal.ofReal (2 * (lam / w)) :=
        Finset.sum_congr rfl fun k _ => harc k
    _ = (w : ENNReal) * ENNReal.ofReal (2 * (lam / w)) := by
        rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul]
    _ = ENNReal.ofReal (2 * lam) := by
        rw [show ((w : ENNReal)) = ENNReal.ofReal (w : ℝ) by
          simp [ENNReal.ofReal_natCast],
          ← ENNReal.ofReal_mul (le_of_lt hwR)]
        congr 1
        field_simp

/-! ## Axiom audit -/
#print axioms badArcs_window_subset
#print axioms single_comb_le

end LonelyRunner.LRC14.Hunter
