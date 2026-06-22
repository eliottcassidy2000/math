/-
  TournamentH7.LRCMarginalUniform -- the marginal phase-uniformity ATOM of the
  resonance bound (kind-pasteur-2026-06-22-S31, HYP-2840).

  Each runner's phase `frac(w·x)` is marginally uniform: the slow-time measure of
  `{x : frac(w·x) ∈ I}` equals `|I|`.  For the decorrelation/Vitali UPPER bounds we
  only need the `≤` direction, which is the cheap half — the forward `w`-fold-cover
  inclusion `{frac(w·x) ∈ [a,b)} ∩ [0,1) ⊆ ⋃_{i<w} [(a+i)/w, (b+i)/w)` plus
  measure subadditivity (no disjointness, no reverse inclusion needed):

      slowμ {x : frac(w·x) ∈ [a,b)} ≤ b - a.

  This is the atom every "resonance box" measure is bounded by; the decorrelation
  product/IE structure (HYP-2840) is built on it.  Sorry-free.
-/

import TournamentH7.LRCDenseCovers

namespace LonelyRunner
namespace MarginalUniform

open MeasureTheory

/-- **Marginal uniformity (upper bound).**  For a positive integer speed `w` and a
sub-interval `[a,b) ⊆ [0,1)`, the slow-time measure of `{x : frac(w·x) ∈ [a,b)}` is
at most `b - a`.  Proof: the set (in one period) is contained in the `w`-fold cover
`⋃_{i<w} [(a+i)/w, (b+i)/w)`, so its measure is `≤ Σ_{i<w} (b-a)/w = b-a`. -/
theorem slowμ_fract_Ico_le (w : ℕ) (hw : 1 ≤ w) {a b : ℝ}
    (ha : 0 ≤ a) (hab : a ≤ b) (hb : b ≤ 1) :
    DenseCovers.slowμ {x | Int.fract ((w : ℝ) * x) ∈ Set.Ico a b} ≤
      ENNReal.ofReal (b - a) := by
  have hwR : (0 : ℝ) < w := by exact_mod_cast hw
  -- forward w-fold-cover inclusion
  have hsub : {x | Int.fract ((w : ℝ) * x) ∈ Set.Ico a b} ∩ Set.Ico (0 : ℝ) 1 ⊆
      ⋃ i ∈ Finset.Ico (0 : ℤ) (w : ℤ),
        Set.Ico ((a + (i : ℝ)) / w) ((b + (i : ℝ)) / w) := by
    rintro x ⟨⟨hflo, hfhi⟩, hx0, hx1⟩
    -- hflo : a ≤ frac(w x),  hfhi : frac(w x) < b
    -- `frac(w x) = w x - ⌊w x⌋`
    have hfeq : Int.fract ((w : ℝ) * x) = (w : ℝ) * x - (⌊(w : ℝ) * x⌋ : ℝ) := by
      have h := Int.floor_add_fract ((w : ℝ) * x); linarith
    rw [hfeq] at hflo hfhi
    -- index range: 0 ≤ ⌊w x⌋ < w
    have hwx0 : 0 ≤ (w : ℝ) * x := mul_nonneg (le_of_lt hwR) hx0
    have hi0 : (0 : ℤ) ≤ ⌊(w : ℝ) * x⌋ := Int.floor_nonneg.mpr hwx0
    have hiw : ⌊(w : ℝ) * x⌋ < (w : ℤ) := by
      rw [Int.floor_lt]; push_cast; nlinarith
    refine Set.mem_iUnion₂.mpr ⟨⌊(w : ℝ) * x⌋, Finset.mem_Ico.mpr ⟨hi0, hiw⟩, ?_⟩
    rw [Set.mem_Ico]
    constructor
    · rw [div_le_iff₀ hwR]; nlinarith
    · rw [lt_div_iff₀ hwR]; nlinarith
  -- measure bound via restrict + subadditivity
  have hrestr : DenseCovers.slowμ {x | Int.fract ((w : ℝ) * x) ∈ Set.Ico a b}
      = volume ({x | Int.fract ((w : ℝ) * x) ∈ Set.Ico a b} ∩ Set.Ico (0 : ℝ) 1) := by
    rw [DenseCovers.slowμ, Measure.restrict_apply' measurableSet_Ico]
  rw [hrestr]
  calc volume ({x | Int.fract ((w : ℝ) * x) ∈ Set.Ico a b} ∩ Set.Ico (0 : ℝ) 1)
      ≤ volume (⋃ i ∈ Finset.Ico (0 : ℤ) (w : ℤ),
          Set.Ico ((a + (i : ℝ)) / w) ((b + (i : ℝ)) / w)) := measure_mono hsub
    _ ≤ ∑ i ∈ Finset.Ico (0 : ℤ) (w : ℤ),
          volume (Set.Ico ((a + (i : ℝ)) / w) ((b + (i : ℝ)) / w)) :=
        measure_biUnion_finset_le _ _
    _ = ∑ _i ∈ Finset.Ico (0 : ℤ) (w : ℤ), ENNReal.ofReal ((b - a) / w) := by
        apply Finset.sum_congr rfl
        intro i _
        rw [Real.volume_Ico, div_sub_div_same]
        congr 1; ring
    _ = (Finset.Ico (0 : ℤ) (w : ℤ)).card • ENNReal.ofReal ((b - a) / w) := by
        rw [Finset.sum_const]
    _ = w • ENNReal.ofReal ((b - a) / w) := by rw [Int.card_Ico]; simp
    _ = ENNReal.ofReal (b - a) := by
        have hwne : (w : ℝ) ≠ 0 := ne_of_gt hwR
        rw [nsmul_eq_mul, ← ENNReal.ofReal_natCast w,
          ← ENNReal.ofReal_mul (Nat.cast_nonneg w)]
        congr 1
        field_simp

/-- **Single-speed sector atom.**  Each speed `w` lands in a single inner sector
`[j/7, (j+1)/7)` (`j ≤ 6`) with slow-time measure `≤ 1/7`.  This is the atom every
cover-bound box measure is bounded by (the `coverSet` sectors are exactly these). -/
theorem slowμ_fract_sector_le (w : ℕ) (hw : 1 ≤ w) (j : ℕ) (hj : j ≤ 6) :
    DenseCovers.slowμ {x | Int.fract ((w : ℝ) * x) ∈
        Set.Ico ((j : ℝ) / 7) (((j : ℝ) + 1) / 7)} ≤ ENNReal.ofReal (1 / 7) := by
  have hjR : (j : ℝ) ≤ 6 := by exact_mod_cast hj
  have hwidth : ((j : ℝ) + 1) / 7 - (j : ℝ) / 7 = 1 / 7 := by ring
  have h := slowμ_fract_Ico_le w hw (a := (j : ℝ) / 7) (b := ((j : ℝ) + 1) / 7)
    (by positivity) (by linarith) (by linarith)
  rwa [hwidth] at h

/-! ## Axiom audit -/

#print axioms slowμ_fract_Ico_le
#print axioms slowμ_fract_sector_le

end MarginalUniform
end LonelyRunner
