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

/-- **THE MODULAR STAR-SAFE CORE**: for `c ≤ 7` measurable sets each of measure
exactly `1/7`, whose center (`A 0`) meets every other with measure AT LEAST `1/49`,
the union misses a positive-measure safe set.  The pair floor is a LOWER bound, so
ANY source of a `≥ 1/49` center-overlap (seven-commensuration, or a drift floor)
plugs in.  This is the `star_hunter_add_le` Bonferroni + the `(48−6c)/49 > 0` budget,
applied to the REAL Lebesgue measure. -/
theorem star_safe_measure_pos_of_lb (c : ℕ) (hc1 : 1 ≤ c) (hc7 : c ≤ 7)
    (A : ℕ → Set UnitAddCircle) (hA : ∀ i, MeasurableSet (A i))
    (hs : ∀ i ∈ Finset.range c, volume (A i) = ENNReal.ofReal (1 / 7))
    (hp : ∀ i ∈ Finset.Ico 1 c, ENNReal.ofReal (1 / 49) ≤ volume (A 0 ∩ A i)) :
    0 < volume ((⋃ i ∈ Finset.range c, A i)ᶜ) := by
  have hunter := star_hunter_add_le volume A hA c
  have hc1R : (1 : ℝ) ≤ (c : ℝ) := by exact_mod_cast hc1
  have hcR : (c : ℝ) ≤ 7 := by exact_mod_cast hc7
  have hcast1 : ((c - 1 : ℕ) : ℝ) = (c : ℝ) - 1 := by rw [Nat.cast_sub hc1]; norm_num
  have hsum_p : ∑ _i ∈ Finset.Ico 1 c, ENNReal.ofReal (1 / 49)
      = ENNReal.ofReal (((c : ℝ) - 1) / 49) := by
    rw [Finset.sum_const, Nat.card_Ico, nsmul_eq_mul, ← ENNReal.ofReal_natCast,
      ← ENNReal.ofReal_mul (by positivity), hcast1]
    congr 1; ring
  have hsum_s : ∑ i ∈ Finset.range c, volume (A i) = ENNReal.ofReal ((c : ℝ) / 7) := by
    rw [Finset.sum_congr rfl hs, Finset.sum_const, Finset.card_range, nsmul_eq_mul,
      ← ENNReal.ofReal_natCast, ← ENNReal.ofReal_mul (by positivity)]
    congr 1; ring
  have hpsum : ENNReal.ofReal (((c : ℝ) - 1) / 49)
      ≤ ∑ i ∈ Finset.Ico 1 c, volume (A 0 ∩ A i) := by
    rw [← hsum_p]; exact Finset.sum_le_sum hp
  set U : ℝ≥0∞ := volume (⋃ i ∈ Finset.range c, A i) with hU
  have hkey : U + ENNReal.ofReal (((c : ℝ) - 1) / 49) ≤ ENNReal.ofReal ((c : ℝ) / 7) := by
    calc U + ENNReal.ofReal (((c : ℝ) - 1) / 49)
        ≤ U + ∑ i ∈ Finset.Ico 1 c, volume (A 0 ∩ A i) := by gcongr
      _ ≤ ∑ i ∈ Finset.range c, volume (A i) := hunter
      _ = ENNReal.ofReal ((c : ℝ) / 7) := hsum_s
  have hsplit : ENNReal.ofReal ((c : ℝ) / 7)
      = ENNReal.ofReal ((6 * (c : ℝ) + 1) / 49) + ENNReal.ofReal (((c : ℝ) - 1) / 49) := by
    rw [← ENNReal.ofReal_add (by positivity) (div_nonneg (by linarith [hc1R]) (by norm_num))]
    congr 1; ring
  rw [hsplit] at hkey
  have hUle : U ≤ ENNReal.ofReal ((6 * (c : ℝ) + 1) / 49) :=
    (ENNReal.add_le_add_iff_right ENNReal.ofReal_ne_top).mp hkey
  have hUlt : U < 1 := by
    calc U ≤ ENNReal.ofReal ((6 * (c : ℝ) + 1) / 49) := hUle
      _ ≤ ENNReal.ofReal (43 / 49) := by
          apply ENNReal.ofReal_le_ofReal
          rw [div_le_div_iff_of_pos_right (by norm_num)]; linarith
      _ < ENNReal.ofReal 1 := by
          apply (ENNReal.ofReal_lt_ofReal_iff_of_nonneg (by norm_num)).mpr; norm_num
      _ = 1 := ENNReal.ofReal_one
  have hUmeas : MeasurableSet (⋃ i ∈ Finset.range c, A i) :=
    Finset.measurableSet_biUnion _ (fun i _ => hA i)
  have hUne : U ≠ ⊤ := by rw [hU]; exact ne_top_of_lt hUlt
  have huniv : volume (Set.univ : Set UnitAddCircle) = 1 := by
    have := AddCircle.measure_univ (T := 1); simpa using this
  rw [MeasureTheory.measure_compl hUmeas hUne, huniv, ← hU]
  exact tsub_pos_of_lt hUlt

/-- **THE STAR-SAFE MEASURE BOUND** (commensurate instance): for a `c ≤ 7` block whose
center speed `v 0` is `7`-divisible (nonzero) and whose leaves `v 1 … v (c−1)` are not
`7`-divisible (all nonzero), the safe set (complement of the union of real danger sets)
has positive measure.  The pair floor comes from `seven_commensuration'` (`= 1/49`). -/
theorem star_safe_pos (c : ℕ) (hc1 : 1 ≤ c) (hc7 : c ≤ 7)
    (v : ℕ → ℤ) (φ : ℕ → UnitAddCircle)
    (hv0 : v 0 ≠ 0) (hv0_7 : (7 : ℤ) ∣ v 0)
    (hvne : ∀ i ∈ Finset.range c, v i ≠ 0)
    (hvi_n7 : ∀ i ∈ Finset.Ico 1 c, ¬ (7 : ℤ) ∣ v i) :
    0 < volume ((⋃ i ∈ Finset.range c, danger (v i) (φ i) (1 / 14))ᶜ) := by
  refine star_safe_measure_pos_of_lb c hc1 hc7 (fun i => danger (v i) (φ i) (1 / 14))
    (fun i => measurableSet_danger _ _ _) ?_ ?_
  · intro i hi
    rw [volume_danger (hvne i hi) (φ i) (by norm_num) (by norm_num)]; norm_num
  · intro i hi
    rw [seven_commensuration' hv0 hv0_7 (hvi_n7 i hi) (φ 0) (φ i)]

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

/-- **THE STAR LONELY TIME** (real-time bridge): a covering `c ≤ 7` block with a
`7`-divisible center is `1/14`-lonely at some REAL time `t` — every runner keeps
integer distance `≥ 1/14`.  This is the err-free measure-theoretic loneliness at the
critical band, with NO citation, NO window, NO singles-bound loss. -/
theorem exists_star_lonely_real (c : ℕ) (hc1 : 1 ≤ c) (hc7 : c ≤ 7)
    (v : ℕ → ℤ)
    (hv0 : v 0 ≠ 0) (hv0_7 : (7 : ℤ) ∣ v 0)
    (hvne : ∀ i ∈ Finset.range c, v i ≠ 0)
    (hvi_n7 : ∀ i ∈ Finset.Ico 1 c, ¬ (7 : ℤ) ∣ v i) :
    ∃ t : ℝ, ∀ i ∈ Finset.range c, ∀ m : ℤ, (1 : ℝ) / 14 ≤ |(v i : ℝ) * t - m| := by
  obtain ⟨x, hx⟩ := exists_star_lonely c hc1 hc7 v (fun _ => 0) hv0 hv0_7 hvne hvi_n7
  obtain ⟨t, rfl⟩ := QuotientAddGroup.mk_surjective x
  refine ⟨t, ?_⟩
  intro i hi m
  have hxi := hx i hi
  rw [mem_danger] at hxi
  have hconv : (v i • ((t : ℝ) : UnitAddCircle)) + (0 : UnitAddCircle)
      = (((v i : ℝ) * t : ℝ) : UnitAddCircle) := by
    rw [add_zero]
    rw [show ((v i : ℝ) * t : ℝ) = (v i • t : ℝ) by rw [zsmul_eq_mul]]
    exact (QuotientAddGroup.mk_zsmul _ t (v i)).symm
  rw [hconv, Metric.mem_ball, dist_zero_right] at hxi
  push_neg at hxi
  rw [UnitAddCircle.norm_eq] at hxi
  calc (1 : ℝ) / 14 ≤ |(v i : ℝ) * t - round ((v i : ℝ) * t)| := hxi
    _ ≤ |(v i : ℝ) * t - m| := round_le _ m

/-! ## The measure route at its limit — mutual independence closes ANY n -/

/-- The unit circle carries a probability measure. -/
instance : IsProbabilityMeasure (volume : Measure UnitAddCircle) :=
  ⟨by have := AddCircle.measure_univ (T := 1); simpa using this⟩

open scoped ProbabilityTheory in
/-- **THE MUTUAL-INDEPENDENCE CLOSER**: if the real danger sets of a family are
MUTUALLY INDEPENDENT as events on the circle, the whole family is lonely at the
critical band `1/14` — for ANY number of runners, including 14.  The safe set has
measure `∏ (1 − 1/7) = (6/7)ⁿ > 0`.  This is the measure route at its limit: no cap
at 7, because full independence supplies ALL the inclusion–exclusion credits, not
just a spanning tree's `n−1`.  The honest obstruction is the hypothesis itself —
near-equal integer speeds are strongly correlated, NOT independent, and that
correlation is precisely the difficulty of LRC. -/
theorem exists_iIndep_lonely {ι : Type*} [Fintype ι] (v : ι → ℤ) (φ : ι → UnitAddCircle)
    (hvne : ∀ i, v i ≠ 0)
    (hindep : ProbabilityTheory.iIndepSet
      (fun i => danger (v i) (φ i) (1 / 14)) volume) :
    ∃ x : UnitAddCircle, ∀ i, x ∉ danger (v i) (φ i) (1 / 14) := by
  set A : ι → Set UnitAddCircle := fun i => danger (v i) (φ i) (1 / 14) with hA
  have hAmeas : ∀ i, MeasurableSet (A i) := fun i => measurableSet_danger _ _ _
  have hindep' := (ProbabilityTheory.iIndepSet_iff_iIndep A volume).mp hindep
  -- complements are measurable in the generated σ-algebras
  have hcompl : ∀ i, MeasurableSet[MeasurableSpace.generateFrom {A i}] ((A i)ᶜ) :=
    fun i => (MeasurableSpace.measurableSet_generateFrom (Set.mem_singleton _)).compl
  have hprod := hindep'.meas_iInter hcompl
  -- each complement has measure 6/7
  have hcompl_meas : ∀ i, volume ((A i)ᶜ) = ENNReal.ofReal (6 / 7) := by
    intro i
    rw [MeasureTheory.measure_compl (hAmeas i) (measure_ne_top _ _)]
    rw [hA]
    rw [volume_danger (hvne i) (φ i) (by norm_num) (by norm_num)]
    have huniv : volume (Set.univ : Set UnitAddCircle) = 1 := by
      have := AddCircle.measure_univ (T := 1); simpa using this
    rw [huniv]
    rw [show (2 : ℝ) * (1 / 14) = 1 / 7 by norm_num]
    rw [show (1 : ℝ≥0∞) = ENNReal.ofReal 1 from ENNReal.ofReal_one.symm]
    rw [← ENNReal.ofReal_sub _ (by norm_num)]
    norm_num
  -- the safe set has positive measure
  have hpos : 0 < volume (⋂ i, (A i)ᶜ) := by
    rw [hprod]
    rw [Finset.prod_congr rfl (fun i _ => hcompl_meas i)]
    rw [Finset.prod_const]
    apply ENNReal.pow_pos
    simp [ENNReal.ofReal_pos]
  -- positive measure => nonempty
  have hne : (⋂ i, (A i)ᶜ).Nonempty := by
    by_contra h
    rw [Set.not_nonempty_iff_eq_empty] at h
    rw [h] at hpos
    simp at hpos
  obtain ⟨x, hx⟩ := hne
  refine ⟨x, fun i => ?_⟩
  have := Set.mem_iInter.mp hx i
  rwa [Set.mem_compl_iff, hA] at this

open scoped ProbabilityTheory in
/-- Real-time form of the independence closer: mutually independent danger sets give a
real time `t` at which every runner keeps integer-distance `≥ 1/14`. -/
theorem exists_iIndep_lonely_real {ι : Type*} [Fintype ι] (v : ι → ℤ)
    (hvne : ∀ i, v i ≠ 0)
    (hindep : ProbabilityTheory.iIndepSet
      (fun i => danger (v i) (0 : UnitAddCircle) (1 / 14)) volume) :
    ∃ t : ℝ, ∀ i, ∀ m : ℤ, (1 : ℝ) / 14 ≤ |(v i : ℝ) * t - m| := by
  obtain ⟨x, hx⟩ := exists_iIndep_lonely v (fun _ => 0) hvne hindep
  obtain ⟨t, rfl⟩ := QuotientAddGroup.mk_surjective x
  refine ⟨t, ?_⟩
  intro i m
  have hxi := hx i
  rw [mem_danger] at hxi
  have hconv : (v i • ((t : ℝ) : UnitAddCircle)) + (0 : UnitAddCircle)
      = (((v i : ℝ) * t : ℝ) : UnitAddCircle) := by
    rw [add_zero, show ((v i : ℝ) * t : ℝ) = (v i • t : ℝ) by rw [zsmul_eq_mul]]
    exact (QuotientAddGroup.mk_zsmul _ t (v i)).symm
  rw [hconv, Metric.mem_ball, dist_zero_right] at hxi
  push_neg at hxi
  rw [UnitAddCircle.norm_eq] at hxi
  calc (1 : ℝ) / 14 ≤ |(v i : ℝ) * t - round ((v i : ℝ) * t)| := hxi
    _ ≤ |(v i : ℝ) * t - m| := round_le _ m

end StarSafe
end LonelyRunner
