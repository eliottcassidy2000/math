/-
  TournamentH7.LRCMomentDual  (mac-mini-2026-06-22-S29)

  Standalone module (not yet root-imported), so it does not affect the team build.
  Target-build verified by codex-2026-06-26 with:
      lake build TournamentH7.LRCMomentDual

  The THM-534 moment-LP dual reduction for hp0cap, formalized:
      p0(E) = slowμ(coverSet E)  ≤  L_y(E) := ∫ g(missCount) dslowμ
  whenever the dual polynomial `g` is FEASIBLE: `g(t) ≥ 1[t=0]` for `t ∈ {0,…,6}`.

  This machine-checks the reduction `hp0cap ⟸ (L_y(E) ≤ cap_k)`, turning the deep
  cover bound into the scalar moment-extremality "consec maximizes L_y".  The proof
  is POINTWISE: `1_{coverSet}(x) = 1[missCount=0] ≤ g(missCount x)` (dual feasibility
  + `missCount ≤ 6`), then integral monotonicity.  `missCount E x` = #missed inner
  sectors, built on kps's `measurableSet_sector_hit`.
-/
import Mathlib
import TournamentH7.LRCDenseCovers

open MeasureTheory Set
open LonelyRunner.DenseCovers
open scoped Classical

namespace TournamentH7.MomentDual

/-- The set where sector `j` is hit by some speed in `E` (kps's coverSet building block). -/
def sectorHitSet (E : List ℤ) (j : ℕ) : Set ℝ :=
  ⋃ e ∈ E.toFinset, {x : ℝ | (j : ℝ) / 7 ≤ Int.fract ((e : ℝ) * x) ∧
    Int.fract ((e : ℝ) * x) < ((j : ℝ) + 1) / 7}

theorem measurableSet_sectorHitSet (E : List ℤ) (j : ℕ) :
    MeasurableSet (sectorHitSet E j) :=
  (E.toFinset.finite_toSet).measurableSet_biUnion
    (fun e _ => measurableSet_sector_hit e j)

/-- The number of inner sectors `{1,…,6}` NOT hit at `x`. -/
noncomputable def missCount (E : List ℤ) (x : ℝ) : ℕ :=
  ((Finset.Icc 1 6).filter (fun j => x ∉ sectorHitSet E j)).card

theorem missCount_le_six (E : List ℤ) (x : ℝ) : missCount E x ≤ 6 := by
  have h := Finset.card_filter_le (Finset.Icc 1 6) (fun j => x ∉ sectorHitSet E j)
  simpa [Nat.card_Icc] using h

theorem missCount_eq_sum (E : List ℤ) (x : ℝ) :
    missCount E x = ∑ j ∈ Finset.Icc 1 6, (if x ∉ sectorHitSet E j then 1 else 0) := by
  rw [missCount, Finset.card_filter]

theorem measurable_missCount (E : List ℤ) : Measurable (fun x => missCount E x) := by
  simp only [missCount_eq_sum]
  refine Finset.measurable_sum _ (fun j _ => ?_)
  have hset : (fun x => if x ∉ sectorHitSet E j then (1 : ℕ) else 0)
      = Set.indicator (sectorHitSet E j)ᶜ (fun _ => 1) := by
    funext x; by_cases hx : x ∈ sectorHitSet E j <;> simp [Set.indicator, hx]
  rw [hset]
  exact (measurable_const).indicator (measurableSet_sectorHitSet E j).compl

theorem measurable_g_missCount (E : List ℤ) (g : ℕ → ℝ) :
    Measurable (fun x => g (missCount E x)) :=
  (measurable_from_top (f := g)).comp (measurable_missCount E)

/-- `coverSet E = {x | missCount E x = 0}` (all inner sectors hit ⟺ none missed). -/
theorem coverSet_eq_missCount_zero (E : List ℤ) :
    coverSet E = {x | missCount E x = 0} := by
  ext x
  simp only [coverSet, missCount, Set.mem_setOf_eq, Finset.card_eq_zero,
    Finset.filter_eq_empty_iff, Finset.mem_Icc, not_not, sectorHitSet,
    Set.mem_iUnion, List.mem_toFinset, exists_prop]
  constructor
  · intro h j hj; exact h j hj.1 hj.2
  · intro h j hj1 hj6; exact h ⟨hj1, hj6⟩

/-- The witness floor functional `L_y(E) := ∫ g(missCount) dslowμ`. -/
noncomputable def Ly (E : List ℤ) (g : ℕ → ℝ) : ℝ :=
  ∫ x, g (missCount E x) ∂slowμ

theorem integrable_g_missCount (E : List ℤ) (g : ℕ → ℝ) :
    Integrable (fun x => g (missCount E x)) slowμ := by
  refine Integrable.mono' (integrable_const (∑ n ∈ Finset.range 7, |g n|))
    (measurable_g_missCount E g).aestronglyMeasurable
    (Filter.Eventually.of_forall (fun x => ?_))
  rw [Real.norm_eq_abs]
  exact Finset.single_le_sum (fun i _ => abs_nonneg (g i))
    (Finset.mem_range.mpr (Nat.lt_succ_of_le (missCount_le_six E x)))

/-- **THM-534 moment-LP dual, formalized.**  If `g` is dual-feasible
(`g(t) ≥ 1[t=0]` for `t ≤ 6`), then `p0(E) ≤ L_y(E)`. -/
theorem p0_le_Ly (E : List ℤ) (g : ℕ → ℝ)
    (hfeas : ∀ n : ℕ, n ≤ 6 → (if n = 0 then (1 : ℝ) else 0) ≤ g n) :
    (slowμ (coverSet E)).toReal ≤ Ly E g := by
  have hcov : (slowμ (coverSet E)).toReal
      = ∫ x, Set.indicator (coverSet E) (fun _ => (1 : ℝ)) x ∂slowμ := by
    rw [integral_indicator_const (1 : ℝ) (measurableSet_coverSet E)]
    simp [Measure.real]
  rw [Ly, hcov]
  refine integral_mono ((integrable_const (1 : ℝ)).indicator (measurableSet_coverSet E))
    (integrable_g_missCount E g) (fun x => ?_)
  by_cases hx : x ∈ coverSet E
  · have h0 : missCount E x = 0 := by rw [coverSet_eq_missCount_zero] at hx; exact hx
    rw [Set.indicator_of_mem hx, h0]
    simpa using hfeas 0 (by norm_num)
  · rw [Set.indicator_of_notMem hx]
    have hne : missCount E x ≠ 0 := by
      rw [coverSet_eq_missCount_zero] at hx; simpa using hx
    have h := hfeas (missCount E x) (missCount_le_six E x)
    rwa [if_neg hne] at h

#print axioms p0_le_Ly

end TournamentH7.MomentDual
