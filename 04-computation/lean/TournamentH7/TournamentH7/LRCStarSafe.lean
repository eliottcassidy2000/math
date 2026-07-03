/-
  TournamentH7.LRCStarSafe — THE STAR-SAFE POSITIVITY CAPSTONE
  (kind-pasteur-2026-07-03-S26).

  The err-free, measure-theoretic close of the seven-wall for a covering `c ≤ 7`
  block, on the REAL danger sets (Lebesgue measure on the circle `ℝ/ℤ`).

  Ingredients, all already in the corpus:
    * `LRCCommensuration.volume_danger`  — each runner's danger has measure EXACTLY
      `2·(1/14) = 1/7` (the TIGHT singles bound the region-length route lacked — kps-S25);
    * `LRCCommensuration.seven_commensuration'` — a 7-divisible center and a non-7 leaf
      have danger-overlap EXACTLY `(2·1/14)² = 1/49`, for every pair of phases;
    * `Hunter.star_union_le`  — the star-Bonferroni bound
      `μ(⋃) + (c−1)/49 ≤ c/7`.

  Combining: `μ(⋃ danger) ≤ (6c+1)/49 ≤ 43/49 < 1` for `c ≤ 7`, so the SAFE set
  (complement of the union) has measure `≥ (48−6c)/49 > 0` — hence is nonempty, and
  the whole `c`-block is simultaneously lonely at some phase.  This supplies BOTH the
  tight singles and the pair floor at once, exactly as the S25 analysis identified.

  Kernel-pure (`propext, Classical.choice, Quot.sound`); no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCHunterLedger
import TournamentH7.LRCCommensuration

open MeasureTheory Finset Set
open scoped ENNReal

namespace LonelyRunner
namespace StarSafe

open LRCCommensuration
open LonelyRunner.LRC14.Hunter

/-- **THE STAR-SAFE MEASURE BOUND**: for a `c ≤ 7` block whose center speed `v 0` is
`7`-divisible (nonzero) and whose leaves `v 1 … v (c−1)` are not `7`-divisible (all
nonzero), the union of the real danger sets has measure at most `(6c+1)/49`, so the
safe set has measure at least `(48−6c)/49 > 0`. -/
theorem star_safe_pos (c : ℕ) (hc1 : 1 ≤ c) (hc7 : c ≤ 7)
    (v : ℕ → ℤ) (φ : ℕ → UnitAddCircle)
    (hv0 : v 0 ≠ 0) (hv0_7 : (7 : ℤ) ∣ v 0)
    (hvne : ∀ i ∈ Finset.range c, v i ≠ 0)
    (hvi_n7 : ∀ i ∈ Finset.Ico 1 c, ¬ (7 : ℤ) ∣ v i) :
    0 < volume ((⋃ i ∈ Finset.range c, danger (v i) (φ i) (1 / 14))ᶜ) := by
  set A : ℕ → Set UnitAddCircle := fun i => danger (v i) (φ i) (1 / 14) with hA
  have hAmeas : ∀ i, MeasurableSet (A i) := fun i => measurableSet_danger _ _ _
  -- singles: each danger has measure 1/7
  have hs : ∀ i ∈ Finset.range c, volume (A i) = ENNReal.ofReal (1 / 7) := by
    intro i hi
    rw [hA]
    rw [volume_danger (hvne i hi) (φ i) (by norm_num) (by norm_num)]
    norm_num
  -- pairs: center ∩ leaf has measure 1/49
  have hp : ∀ i ∈ Finset.Ico 1 c, volume (A 0 ∩ A i) = ENNReal.ofReal (1 / 49) := by
    intro i hi
    rw [hA]
    exact seven_commensuration' hv0 hv0_7 (hvi_n7 i hi) (φ 0) (φ i)
  -- the star-Bonferroni bound
  have hunter := star_union_le volume A hAmeas c hs hp
  -- evaluate the two constant sums
  have hcard_p : (Finset.Ico 1 c).card = c - 1 := Nat.card_Ico 1 c
  have hcard_s : (Finset.range c).card = c := Finset.card_range c
  have hc1R : (1 : ℝ) ≤ (c : ℝ) := by exact_mod_cast hc1
  have hcast1 : ((c - 1 : ℕ) : ℝ) = (c : ℝ) - 1 := by
    rw [Nat.cast_sub hc1]; norm_num
  have hsum_p : ∑ _i ∈ Finset.Ico 1 c, ENNReal.ofReal (1 / 49)
      = ENNReal.ofReal (((c : ℝ) - 1) / 49) := by
    rw [Finset.sum_const, hcard_p, nsmul_eq_mul, ← ENNReal.ofReal_natCast,
      ← ENNReal.ofReal_mul (by positivity), hcast1]
    congr 1; ring
  have hsum_s : ∑ _i ∈ Finset.range c, ENNReal.ofReal (1 / 7)
      = ENNReal.ofReal ((c : ℝ) / 7) := by
    rw [Finset.sum_const, hcard_s, nsmul_eq_mul, ← ENNReal.ofReal_natCast,
      ← ENNReal.ofReal_mul (by positivity)]
    congr 1; ring
  rw [hsum_p, hsum_s] at hunter
  -- U := μ(⋃)
  set U : ℝ≥0∞ := volume (⋃ i ∈ Finset.range c, A i) with hU
  -- split the singles sum: c/7 = (6c+1)/49 + (c-1)/49
  have hsplit : ENNReal.ofReal ((c : ℝ) / 7)
      = ENNReal.ofReal ((6 * (c : ℝ) + 1) / 49) + ENNReal.ofReal (((c : ℝ) - 1) / 49) := by
    rw [← ENNReal.ofReal_add (by positivity) (div_nonneg (by linarith [hc1R]) (by norm_num))]
    congr 1; ring
  rw [hsplit] at hunter
  -- cancel the finite pair term
  have hUle : U ≤ ENNReal.ofReal ((6 * (c : ℝ) + 1) / 49) := by
    have hfin : ENNReal.ofReal (((c : ℝ) - 1) / 49) ≠ ⊤ := ENNReal.ofReal_ne_top
    exact (ENNReal.add_le_add_iff_right hfin).mp hunter
  -- (6c+1)/49 ≤ 43/49 < 1
  have hcR : (c : ℝ) ≤ 7 := by exact_mod_cast hc7
  have hbound : ENNReal.ofReal ((6 * (c : ℝ) + 1) / 49) ≤ ENNReal.ofReal (43 / 49) := by
    apply ENNReal.ofReal_le_ofReal
    rw [div_le_div_iff_of_pos_right (by norm_num)]
    linarith
  have hUlt : U < 1 := by
    calc U ≤ ENNReal.ofReal ((6 * (c : ℝ) + 1) / 49) := hUle
      _ ≤ ENNReal.ofReal (43 / 49) := hbound
      _ < ENNReal.ofReal 1 := by
          apply ENNReal.ofReal_lt_ofReal_iff_of_nonneg (by norm_num) |>.mpr
          norm_num
      _ = 1 := ENNReal.ofReal_one
  -- safe measure = 1 - U > 0
  have hUmeas : MeasurableSet (⋃ i ∈ Finset.range c, A i) :=
    Finset.measurableSet_biUnion _ (fun i _ => hAmeas i)
  have hUne : U ≠ ⊤ := by
    rw [hU]; exact ne_top_of_lt hUlt
  have huniv : volume (Set.univ : Set UnitAddCircle) = 1 := by
    have := AddCircle.measure_univ (T := 1)
    simpa using this
  rw [show (⋃ i ∈ Finset.range c, danger (v i) (φ i) (1 / 14))
      = (⋃ i ∈ Finset.range c, A i) from rfl]
  rw [MeasureTheory.measure_compl hUmeas hUne, huniv]
  rw [← hU]
  exact tsub_pos_of_lt hUlt

/-- **THE STAR LONELY POINT**: a covering `c ≤ 7` block with a `7`-divisible center is
simultaneously lonely at some phase `x` on the circle. -/
theorem exists_star_lonely (c : ℕ) (hc1 : 1 ≤ c) (hc7 : c ≤ 7)
    (v : ℕ → ℤ) (φ : ℕ → UnitAddCircle)
    (hv0 : v 0 ≠ 0) (hv0_7 : (7 : ℤ) ∣ v 0)
    (hvne : ∀ i ∈ Finset.range c, v i ≠ 0)
    (hvi_n7 : ∀ i ∈ Finset.Ico 1 c, ¬ (7 : ℤ) ∣ v i) :
    ∃ x : UnitAddCircle, ∀ i ∈ Finset.range c, x ∉ danger (v i) (φ i) (1 / 14) := by
  have hpos := star_safe_pos c hc1 hc7 v φ hv0 hv0_7 hvne hvi_n7
  have hne : ((⋃ i ∈ Finset.range c, danger (v i) (φ i) (1 / 14))ᶜ).Nonempty := by
    by_contra h
    rw [Set.not_nonempty_iff_eq_empty] at h
    rw [h] at hpos
    simp at hpos
  obtain ⟨x, hx⟩ := hne
  refine ⟨x, ?_⟩
  intro i hi
  rw [Set.mem_compl_iff, Set.mem_iUnion₂] at hx
  push_neg at hx
  exact hx i hi

end StarSafe
end LonelyRunner
