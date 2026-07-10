/-
  TournamentH7.LRCGoodSetBand — the goodSet band lemma (mac-mini-2026-07-09-S65 cont.21).

  The last plumbing lemma of the witness-floor discharge pipeline: a checkable
  per-difference band condition forces `Ico α β ⊆ goodSet E` — the `hgood` input of
  cont.18's `witnessG2_pos_of_anchor`, and (composed with brick (iii)) the full
  m_P floors of `hsmall3`/`hlarge` become engine-lists-in, proof-out
  (`witnessG2_ge_of_sorted_bands` below is the leg-shaped composition).

  Band shapes (e := the positive difference; endpoints rational in practice):
  * up (a < b, e = b − a):  ∃ j, j + 1/7 < e·α  ∧  e·β ≤ j + 1
      — on [α,β): e·x ∈ (j + 1/7, j + 1), so ⌊e·x⌋ = j and fract > 1/7.
  * down (b < a, e = a − b): ∃ j, j ≤ e·α  ∧  e·β ≤ j + 6/7
      — e·x ∈ [j, j + 6/7), fract(e·x) < 6/7, and the reflection
        fract(−y) ∈ {0} ∪ {1 − fract y} lands outside (0, 1/7].
  * b = a: the zero difference, free (`fract 0 = 0`).
-/
import Mathlib
import TournamentH7.LRCFourteenSkeleton
import TournamentH7.LRCGoodSetSmall
import TournamentH7.LRCIntervalBridge
import TournamentH7.LRCUnionVolume

namespace LonelyRunner
namespace LRC14
namespace GoodSetBand

open DenseCovers MeasureTheory TournamentH7.GoodSet

/-- The per-element band certificate for tooth `a`, element `b`, interval `[α, β)`. -/
def bandCert (a b : ℤ) (α β : ℝ) : Prop :=
  b = a ∨
    (a < b ∧ ∃ j : ℤ, (j : ℝ) + 1 / 7 < ((b - a : ℤ) : ℝ) * α ∧
      ((b - a : ℤ) : ℝ) * β ≤ (j : ℝ) + 1) ∨
    (b < a ∧ ∃ j : ℤ, (j : ℝ) ≤ ((a - b : ℤ) : ℝ) * α ∧
      ((a - b : ℤ) : ℝ) * β ≤ (j : ℝ) + 6 / 7)

/-- **The goodSet band lemma.** A tooth `a ∈ E` whose every difference passes the band
check puts the whole interval inside `goodSet E`. -/
theorem Ico_subset_goodSet_of_bounds {E : List ℤ} {a : ℤ} (ha : a ∈ E) {α β : ℝ}
    (hband : ∀ b ∈ E, bandCert a b α β) :
    Set.Ico α β ⊆ goodSet E := by
  intro x hx
  refine Set.mem_iUnion₂.mpr ⟨a, List.mem_toFinset.mpr ha, Set.mem_iInter₂.mpr ?_⟩
  intro b hb
  have hbE : b ∈ E := List.mem_toFinset.mp hb
  simp only [Set.mem_setOf_eq]
  have h := hband b hbE
  simp only [bandCert] at h
  rcases h with rfl | ⟨hab', j, hj1, hj2⟩ | ⟨hba, j, hj1, hj2⟩
  · -- b = a: the zero difference
    rw [sub_self]
    exact GoodSetSmall.fract_zero_not_mem x
  · -- up branch: e = b − a > 0, e·x ∈ (j + 1/7, j + 1)
    have heR : (0 : ℝ) < ((b - a : ℤ) : ℝ) := by exact_mod_cast sub_pos.mpr hab'
    have hex1 : (j : ℝ) + 1 / 7 < ((b - a : ℤ) : ℝ) * x :=
      lt_of_lt_of_le hj1 (mul_le_mul_of_nonneg_left hx.1 heR.le)
    have hex2 : ((b - a : ℤ) : ℝ) * x < (j : ℝ) + 1 :=
      lt_of_lt_of_le (mul_lt_mul_of_pos_left hx.2 heR) hj2
    have hfloor : ⌊((b - a : ℤ) : ℝ) * x⌋ = j := by
      apply Int.floor_eq_iff.mpr
      constructor <;> [linarith; linarith]
    have hfract : Int.fract (((b - a : ℤ) : ℝ) * x) = ((b - a : ℤ) : ℝ) * x - (j : ℝ) := by
      rw [Int.fract, hfloor]
    intro hmem
    rw [hfract] at hmem
    have := hmem.2
    linarith
  · -- down branch: e = a − b > 0, e·x ∈ [j, j + 6/7), then reflect
    have hcast : ((b - a : ℤ) : ℝ) * x = -(((a - b : ℤ) : ℝ) * x) := by push_cast; ring
    have heR : (0 : ℝ) < ((a - b : ℤ) : ℝ) := by exact_mod_cast sub_pos.mpr hba
    have hex1 : (j : ℝ) ≤ ((a - b : ℤ) : ℝ) * x :=
      le_trans hj1 (mul_le_mul_of_nonneg_left hx.1 heR.le)
    have hex2 : ((a - b : ℤ) : ℝ) * x < (j : ℝ) + 6 / 7 :=
      lt_of_lt_of_le (mul_lt_mul_of_pos_left hx.2 heR) hj2
    have hfloor : ⌊((a - b : ℤ) : ℝ) * x⌋ = j := by
      apply Int.floor_eq_iff.mpr
      constructor <;> [linarith; linarith]
    have hfract : Int.fract (((a - b : ℤ) : ℝ) * x) = ((a - b : ℤ) : ℝ) * x - (j : ℝ) := by
      rw [Int.fract, hfloor]
    rw [hcast]
    by_cases h0 : Int.fract (((a - b : ℤ) : ℝ) * x) = 0
    · -- e·x is exactly the integer j: the reflection has fract 0
      have hexj : ((a - b : ℤ) : ℝ) * x = (j : ℝ) := by
        rw [hfract] at h0; linarith
      rw [hexj]
      have hjj : -(j : ℝ) = ((-j : ℤ) : ℝ) := by push_cast; ring
      rw [hjj, Int.fract_intCast]
      simp
    · -- fract(e·x) ∈ (0, 6/7): the reflection has fract = 1 − fract(e·x) > 1/7
      rw [Int.fract_neg h0]
      intro hmem
      have h3 : Int.fract (((a - b : ℤ) : ℝ) * x) < 6 / 7 := by
        rw [hfract]; linarith
      linarith [hmem.2]

/-- Both band checks at once: the interval sits inside the intersection `witnessG2` measures. -/
theorem Ico_subset_good_inter_safe {E P : List ℤ} {a : ℤ} (ha : a ∈ E) {α β : ℝ}
    (hgood : ∀ b ∈ E, bandCert a b α β)
    (hposP : ∀ p ∈ P, (0 : ℤ) < p)
    (hsafe : ∀ p ∈ P, ∃ j : ℤ, (j : ℝ) + 1 / 14 ≤ (p : ℝ) * α ∧
      (p : ℝ) * β ≤ (j : ℝ) + 13 / 14) :
    Set.Ico α β ⊆ goodSet E ∩ safeSet P :=
  Set.subset_inter (Ico_subset_goodSet_of_bounds ha hgood)
    (IntervalBridge.Ico_subset_safeSet_of_bounds hposP hsafe)

/-- **The leg-shaped floor**: a sorted-disjoint engine list, each interval carrying a tooth
and both band certificates, forces `witnessG2 s ≥ Σ lengths` — the full m_P floor shape of
`hsmall3`/`hlarge`, consuming brick (iii). -/
theorem witnessG2_ge_of_sorted_bands (s : Shape) (l : List (ℝ × ℝ))
    (hposP : ∀ p ∈ s.1, (0 : ℤ) < p)
    (hcert : ∀ q ∈ l, ∃ a ∈ s.2, (∀ b ∈ s.2, bandCert a b q.1 q.2) ∧
      (∀ p ∈ s.1, ∃ j : ℤ, (j : ℝ) + 1 / 14 ≤ (p : ℝ) * q.1 ∧
        (p : ℝ) * q.2 ≤ (j : ℝ) + 13 / 14))
    (hin : ∀ q ∈ l, 0 ≤ q.1 ∧ q.1 ≤ q.2 ∧ q.2 ≤ 1)
    (hsorted : l.Pairwise (fun q r => q.2 ≤ r.1)) :
    (l.map (fun q => q.2 - q.1)).sum ≤ witnessG2 s := by
  have hsub : ∀ q ∈ l, Set.Ico q.1 q.2 ⊆ goodSet s.2 ∩ safeSet s.1 := by
    intro q hq
    obtain ⟨a, haE, hgood, hsafe⟩ := hcert q hq
    exact Ico_subset_good_inter_safe haE hgood hposP hsafe
  have hfloor := UnionVolume.slowμ_ge_sum_of_sorted_Ico_subset l hsub hin hsorted
  have hnn : 0 ≤ (l.map (fun q => q.2 - q.1)).sum := by
    apply List.sum_nonneg
    intro y hy
    obtain ⟨q, hq, rfl⟩ := List.mem_map.mp hy
    have h := hin q hq
    linarith [h.2.1]
  have hne : slowμ (goodSet s.2 ∩ safeSet s.1) ≠ ⊤ := measure_ne_top _ _
  unfold witnessG2
  calc (l.map (fun q => q.2 - q.1)).sum
      = (ENNReal.ofReal ((l.map (fun q => q.2 - q.1)).sum)).toReal :=
        (ENNReal.toReal_ofReal hnn).symm
    _ ≤ (slowμ (goodSet s.2 ∩ safeSet s.1)).toReal := ENNReal.toReal_mono hne hfloor

end GoodSetBand
end LRC14
end LonelyRunner
