/-
  TournamentH7.LRCAP44Floor — the AP₄₄ density-floor certificate via the Farey roof,
  using TWO nodes (`0/1, 1/1` endpoints + the `1/2` node).  boxeph-2026-07-07-S2.

  Demonstrates MULTI-NODE roof aggregation — the pattern that scales to the tight AP₇₆.
  Four roof-superlevel intervals (Farey-44 cells adjacent to `0/1`, `1/2`, `1/1`):

    (0, 6/301)   ∪   (141/287, 1/2)   ∪   (1/2, 146/287)   ∪   (295/301, 1)

  each ⊆ Good(1/7)(AP₄₄) by one `good_of_roof_gt` call, pairwise disjoint (sorted),
  total length `6/301 + 5/574 + 5/574 + 6/301 = 101/1763 ≈ 0.05729 ≥ m_P`.

  So the systematic roof route reaches diameter ≤ 43 (vs ≤ 29 for the two endpoint
  intervals alone / ≤ 19 for death-star's hand-built AP₂₀).  Adding the `q=3,4,5,6`
  nodes (24 intervals total) reaches the tight AP₇₆ — same pattern, more cells.

  Kernel-pure: no sorry, no native_decide.
-/
import Mathlib
import TournamentH7.LRCFareyRoofBridge
import TournamentH7.LRCTailDiameter

namespace LonelyRunner
namespace AP44Floor

open TailDiameter FareyRoofBridge MeasureTheory
open scoped ENNReal

/-- Interval 1: near `x = 0`, Farey-44 cell `(0/1, 1/44)`. -/
theorem I1_good :
    Set.Ioo (0 : ℝ) (6 / 301) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 44) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hx0, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (0 : ℤ)) (q := (1 : ℤ)) (p' := (1 : ℤ)) (q' := (44 : ℤ))
    (k := (44 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hx0]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxu]),
    le_of_lt hx0, by nlinarith [hxu]⟩

/-- Interval 2: left of `x = 1/2`, Farey-44 cell `(21/43, 1/2)`. -/
theorem I2_good :
    Set.Ioo (141 / 287 : ℝ) (1 / 2) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 44) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (21 : ℤ)) (q := (43 : ℤ)) (p' := (1 : ℤ)) (q' := (2 : ℤ))
    (k := (44 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxl]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 3: right of `x = 1/2`, Farey-44 cell `(1/2, 22/43)`. -/
theorem I3_good :
    Set.Ioo (1 / 2 : ℝ) (146 / 287) ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 44) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hxu⟩ := hx
  refine ⟨good_of_roof_gt (p := (1 : ℤ)) (q := (2 : ℤ)) (p' := (22 : ℤ)) (q' := (43 : ℤ))
    (k := (44 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hxu]) (by push_cast; nlinarith [hxu]),
    by nlinarith [hxl], by nlinarith [hxu]⟩

/-- Interval 4: near `x = 1`, Farey-44 cell `(43/44, 1/1)`. -/
theorem I4_good :
    Set.Ioo (295 / 301 : ℝ) 1 ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 44) ∩ Set.Icc (0 : ℝ) 1 := by
  intro x hx; obtain ⟨hxl, hx1⟩ := hx
  refine ⟨good_of_roof_gt (p := (43 : ℤ)) (q := (44 : ℤ)) (p' := (1 : ℤ)) (q' := (1 : ℤ))
    (k := (44 : ℤ)) (x := x) (θ := (1 / 7 : ℝ))
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by push_cast; nlinarith [hxl]) (by push_cast; nlinarith [hx1]) (by push_cast; nlinarith [hxl]),
    by nlinarith [hxl], le_of_lt hx1⟩

/-- disjointness of two `Ioo` when the first ends at or before the second starts. -/
private theorem ioo_disj {a b c e : ℝ} (h : b ≤ c) :
    Disjoint (Set.Ioo a b) (Set.Ioo c e) := by
  apply Set.disjoint_left.mpr
  intro x hx1 hx2
  exact absurd (lt_of_lt_of_le hx1.2 (le_trans h (le_of_lt hx2.1))) (lt_irrefl x)

/-- **THE AP₄₄ CERTIFICATE (unconditional).**  `μ_{1/7}(AP₄₄) ≥ m_P`, from four
roof-derived intervals across the nodes `0, 1/2, 1`. -/
theorem ap44_certificate :
    ENNReal.ofReal ((14249 : ℝ) / 252252) ≤ muGood (1 / 7) (Finset.Icc (1 : ℤ) 44) := by
  set I1 := Set.Ioo (0 : ℝ) (6 / 301) with hI1
  set I2 := Set.Ioo (141 / 287 : ℝ) (1 / 2) with hI2
  set I3 := Set.Ioo (1 / 2 : ℝ) (146 / 287) with hI3
  set I4 := Set.Ioo (295 / 301 : ℝ) 1 with hI4
  -- pairwise disjointness (sorted: 6/301 ≤ 141/287 ≤ 1/2 ≤ 146/287 ≤ 295/301)
  have d12 : Disjoint I1 I2 := ioo_disj (by norm_num)
  have d13 : Disjoint I1 I3 := ioo_disj (by norm_num)
  have d14 : Disjoint I1 I4 := ioo_disj (by norm_num)
  have d23 : Disjoint I2 I3 := ioo_disj (le_refl _)
  have d24 : Disjoint I2 I4 := ioo_disj (by norm_num)
  have d34 : Disjoint I3 I4 := ioo_disj (by norm_num)
  -- union ⊆ Good ∩ [0,1]
  have hu : I1 ∪ I2 ∪ I3 ∪ I4 ⊆ Good (1 / 7) (Finset.Icc (1 : ℤ) 44) ∩ Set.Icc (0 : ℝ) 1 :=
    Set.union_subset (Set.union_subset (Set.union_subset I1_good I2_good) I3_good) I4_good
  -- measure of the disjoint union = sum of lengths
  have hvol : volume (I1 ∪ I2 ∪ I3 ∪ I4)
      = volume I1 + volume I2 + volume I3 + volume I4 := by
    rw [measure_union ((d14.union_left d24).union_left d34) measurableSet_Ioo,
        measure_union (d13.union_left d23) measurableSet_Ioo,
        measure_union d12 measurableSet_Ioo]
  have hbound : volume I1 + volume I2 + volume I3 + volume I4
      ≤ muGood (1 / 7) (Finset.Icc (1 : ℤ) 44) := by
    rw [← hvol]; exact measure_mono hu
  -- arithmetic: m_P ≤ 6/301 + 5/574 + 5/574 + 6/301 = 101/1763
  refine le_trans ?_ hbound
  simp only [hI1, hI2, hI3, hI4, Real.volume_Ioo]
  rw [← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num),
      ← ENNReal.ofReal_add (by norm_num) (by norm_num)]
  exact ENNReal.ofReal_le_ofReal (by norm_num)

/-- On `{0,…,43}` via rotation invariance, ready for `TailDiameter.muGood_ge_APD`
(diameter `≤ 43`). -/
theorem ap44_certificate_icc0 :
    ENNReal.ofReal ((14249 : ℝ) / 252252) ≤ muGood (1 / 7) (Finset.Icc (0 : ℤ) 43) := by
  have hE : (Finset.Icc (1 : ℤ) 44).image (fun e => e - 1) = Finset.Icc (0 : ℤ) 43 := by
    ext n; simp only [Finset.mem_image, Finset.mem_Icc]
    constructor
    · rintro ⟨a, ⟨h1, h2⟩, rfl⟩; omega
    · intro h; exact ⟨n + 1, ⟨by omega, by omega⟩, by omega⟩
  have htr : muGood (1 / 7) (Finset.Icc (0 : ℤ) 43) = muGood (1 / 7) (Finset.Icc (1 : ℤ) 44) := by
    have := muGood_translate (1 / 7) (Finset.Icc (1 : ℤ) 44) 1
    rw [hE] at this; exact this
  rw [htr]; exact ap44_certificate

end AP44Floor
end LonelyRunner
