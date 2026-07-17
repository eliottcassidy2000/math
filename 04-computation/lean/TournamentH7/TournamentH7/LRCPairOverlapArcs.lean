/-
  TournamentH7.LRCPairOverlapArcs — THE ARC-MEASURE RENDERING OF THE PAIR
  OVERLAP (boxeph-2026-07-17-S74; LEM-042(A)'s Lean face, lower-bound form).

  Danger set on ℝ: `dangerR v = {t | ∃ m : ℤ, |v·t − m| < 1/14}` (1-periodic).
  For the consecutive pair (k, k+1) the Bezout witness is trivial — the j-th
  overlap component carries the SAME integer witness `m = −j` on both runners —
  so every component is one explicit rational interval

      comp k j = Ioo (max ((−j−ρ)/k) ((−j−ρ)/(k+1)))
                     (min ((−j+ρ)/k) ((−j+ρ)/(k+1))),   ρ = 1/14.

  On the window (−1/2, 1/2) the components (|j| ≤ (2k+1)/14) sit symmetrically
  around 0 with no wrap split.  Kernel-pure results:

    · `comp_subset_danger`  — every component lies in dangerR k ∩ dangerR (k+1);
    · `volume_comp_zero`    — the 0-component has length exactly 2ρ/(k+1)
                              (the capped `w_min` term of LEM-042's trapezoid);
    · `volume_comp_tent`    — the j ≠ 0 components have length ≥ the tent term
                              ((2k+1)/14 − |j|)/(k(k+1));
    · `comp_disjoint_of_lt` — the components are pairwise disjoint (no range
                              condition needed);
    · `consecutive_overlap_credit` — the assembled lower bound: for 14J ≤ 2k+1,
        Σ_{j ∈ Icc (−J) J} ofReal(tentAbs j) ≤ volume (D_k ∩ D_{k+1} ∩ window),
      the CONCRETE CREDIT for `good_pos_of_path_credits` (LRCTreeHunter).
      Identifying the left sum with `muNum k (k+1)/(14k(k+1))` is the named
      cast bridge (finite rational bookkeeping).

  Kernel-pure: no `native_decide`, no `sorry`.
-/
import Mathlib

open MeasureTheory Set

namespace LonelyRunner.LRC14.Arcs

noncomputable section

/-- The 1/14-danger set of runner `v` on the line. -/
def dangerR (v : ℕ) : Set ℝ := {t | ∃ m : ℤ, |(v : ℝ) * t - m| < 1/14}

/-- The `j`-th overlap component of the consecutive pair `(k, k+1)`. -/
def comp (k : ℕ) (j : ℤ) : Set ℝ :=
  Ioo (max (((-j : ℝ) - 1/14) / k) (((-j : ℝ) - 1/14) / (k + 1)))
      (min (((-j : ℝ) + 1/14) / k) (((-j : ℝ) + 1/14) / (k + 1)))

/-- The trapezoid term the `j`-th component is worth. -/
def tentAbs (k : ℕ) (j : ℤ) : ℝ :=
  if j = 0 then 2 * (1/14) / ((k : ℝ) + 1)
  else ((2 * k + 1 : ℝ) / 14 - |(j : ℝ)|) / ((k : ℝ) * (k + 1))

theorem comp_subset_danger (k : ℕ) (hk : 1 ≤ k) (j : ℤ) :
    comp k j ⊆ dangerR k ∩ dangerR (k + 1) := by
  have hkR : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk
  have hk1R : (0 : ℝ) < (k : ℝ) + 1 := by linarith
  rintro t ⟨hlo, hhi⟩
  rw [max_lt_iff] at hlo
  rw [lt_min_iff] at hhi
  refine ⟨⟨-j, ?_⟩, ⟨-j, ?_⟩⟩
  · rw [abs_lt]
    have h1 := (div_lt_iff₀ hkR).mp hlo.1
    have h2 := (lt_div_iff₀ hkR).mp hhi.1
    constructor
    · push_cast
      nlinarith
    · push_cast
      nlinarith
  · rw [abs_lt]
    have hcast : ((k + 1 : ℕ) : ℝ) = (k : ℝ) + 1 := by push_cast; ring
    rw [hcast]
    have h1 := (div_lt_iff₀ hk1R).mp hlo.2
    have h2 := (lt_div_iff₀ hk1R).mp hhi.2
    constructor
    · push_cast
      nlinarith
    · push_cast
      nlinarith

/-- The 0-component's exact length: `2ρ/(k+1)` — the capped trapezoid term. -/
theorem volume_comp_zero (k : ℕ) (hk : 1 ≤ k) :
    volume (comp k 0) = ENNReal.ofReal (2 * (1/14) / ((k : ℝ) + 1)) := by
  have hkR : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk
  have hk1R : (0 : ℝ) < (k : ℝ) + 1 := by linarith
  have hmax : max (((-(0:ℤ) : ℝ) - 1/14) / k) (((-(0:ℤ) : ℝ) - 1/14) / (k + 1))
      = (-(1/14) : ℝ) / ((k : ℝ) + 1) := by
    rw [max_eq_right]
    · norm_num
    · rw [div_le_div_iff₀ hkR hk1R]
      push_cast
      nlinarith
  have hmin : min (((-(0:ℤ) : ℝ) + 1/14) / k) (((-(0:ℤ) : ℝ) + 1/14) / (k + 1))
      = (1/14 : ℝ) / ((k : ℝ) + 1) := by
    rw [min_eq_right]
    · norm_num
    · rw [div_le_div_iff₀ hk1R hkR]
      push_cast
      nlinarith
  rw [comp, hmax, hmin, Real.volume_Ioo]
  congr 1
  field_simp
  ring

/-- The `j ≠ 0` components carry at least the tent term. -/
theorem volume_comp_tent (k : ℕ) (hk : 1 ≤ k) (j : ℤ) (hj : j ≠ 0)
    (hj14 : 14 * |(j : ℝ)| ≤ 2 * k + 1) :
    ENNReal.ofReal (((2 * k + 1 : ℝ) / 14 - |(j : ℝ)|) / ((k : ℝ) * (k + 1)))
      ≤ volume (comp k j) := by
  have hkR : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk
  have hk1R : (0 : ℝ) < (k : ℝ) + 1 := by linarith
  rw [comp, Real.volume_Ioo]
  apply ENNReal.ofReal_le_ofReal
  rcases lt_or_gt_of_ne hj with hneg | hpos
  · -- j ≤ −1: the component sits at t ≈ |j|/k > 0
    have hj1 : (1 : ℝ) ≤ -(j : ℝ) := by
      have : (j : ℤ) ≤ -1 := by omega
      have : (j : ℝ) ≤ -1 := by exact_mod_cast this
      linarith
    have habs : |(j : ℝ)| = -(j : ℝ) := abs_of_neg (by linarith)
    -- min ≥ (−j+ρ)/(k+1), max ≤ (−j−ρ)/k
    have hmin : ((-j : ℝ) + 1/14) / (k + 1)
        ≤ min (((-j : ℝ) + 1/14) / k) (((-j : ℝ) + 1/14) / (k + 1)) := by
      apply le_min _ (le_refl _)
      rw [div_le_div_iff₀ hk1R hkR]
      nlinarith
    have hmax : max (((-j : ℝ) - 1/14) / k) (((-j : ℝ) - 1/14) / (k + 1))
        ≤ ((-j : ℝ) - 1/14) / k := by
      apply max_le (le_refl _)
      rw [div_le_div_iff₀ hk1R hkR]
      nlinarith
    have hval : ((2 * k + 1 : ℝ) / 14 - |(j : ℝ)|) / ((k : ℝ) * (k + 1))
        = ((-j : ℝ) + 1/14) / (k + 1) - ((-j : ℝ) - 1/14) / k := by
      rw [habs]
      field_simp
      ring
    linarith
  · -- j ≥ 1: the component sits at t ≈ −j/k < 0
    have hj1 : (1 : ℝ) ≤ (j : ℝ) := by exact_mod_cast hpos
    have habs : |(j : ℝ)| = (j : ℝ) := abs_of_pos (by linarith)
    -- min ≥ (−j+ρ)/k, max ≤ (−j−ρ)/(k+1)
    have hmin : ((-j : ℝ) + 1/14) / k
        ≤ min (((-j : ℝ) + 1/14) / k) (((-j : ℝ) + 1/14) / (k + 1)) := by
      apply le_min (le_refl _)
      rw [div_le_div_iff₀ hkR hk1R]
      nlinarith
    have hmax : max (((-j : ℝ) - 1/14) / k) (((-j : ℝ) - 1/14) / (k + 1))
        ≤ ((-j : ℝ) - 1/14) / (k + 1) := by
      apply max_le _ (le_refl _)
      rw [div_le_div_iff₀ hkR hk1R]
      nlinarith
    have hval : ((2 * k + 1 : ℝ) / 14 - |(j : ℝ)|) / ((k : ℝ) * (k + 1))
        = ((-j : ℝ) + 1/14) / k - ((-j : ℝ) - 1/14) / (k + 1) := by
      rw [habs]
      field_simp
      ring
    linarith

/-- Components are pairwise disjoint: for `i < j`, `comp k j` lies entirely to the
left of `comp k i` (no range condition needed). -/
theorem comp_disjoint_of_lt (k : ℕ) (hk : 1 ≤ k) {i j : ℤ} (hij : i < j) :
    Disjoint (comp k i) (comp k j) := by
  have hkR : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk
  have hk1R : (0 : ℝ) < (k : ℝ) + 1 := by linarith
  have hijR : (i : ℝ) + 1 ≤ (j : ℝ) := by exact_mod_cast hij
  -- same-denominator endpoints: upper end of comp j ≤ (−j+ρ)/(k+1) ≤ (−i−ρ)/(k+1)
  -- ≤ lower end of comp i, since i + 2/14 ≤ j
  have hkey : ((-j : ℝ) + 1/14) / (k + 1) ≤ ((-i : ℝ) - 1/14) / (k + 1) := by
    rw [div_le_div_iff₀ hk1R hk1R]
    nlinarith
  rw [Set.disjoint_left]
  rintro t ⟨hlo_i, _⟩ ⟨_, hhi_j⟩
  rw [max_lt_iff] at hlo_i
  rw [lt_min_iff] at hhi_j
  have h1 : t < ((-j : ℝ) + 1/14) / (k + 1) := hhi_j.2
  have h2 : ((-i : ℝ) - 1/14) / (k + 1) < t := hlo_i.2
  linarith

/-- Components stay inside the unit window `(−1/2, 1/2)` whenever `14|j| ≤ 2k+1`. -/
theorem comp_subset_window (k : ℕ) (hk : 1 ≤ k) (j : ℤ)
    (hj14 : 14 * |(j : ℝ)| ≤ 2 * k + 1) :
    comp k j ⊆ Ioo (-(1:ℝ)/2) (1/2) := by
  have hkR : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk
  have hk1 : (1 : ℝ) ≤ (k : ℝ) := by exact_mod_cast hk
  have hk1R : (0 : ℝ) < (k : ℝ) + 1 := by linarith
  have hjabs : -(|(j : ℝ)|) ≤ (j : ℝ) ∧ (j : ℝ) ≤ |(j : ℝ)| := ⟨neg_abs_le _, le_abs_self _⟩
  rintro t ⟨hlo, hhi⟩
  rw [max_lt_iff] at hlo
  rw [lt_min_iff] at hhi
  constructor
  · -- −1/2 < t : from the k-tooth lower end
    have h := hlo.1
    rw [div_lt_iff₀ hkR] at h
    nlinarith [hjabs.2, abs_nonneg ((j : ℝ)), hk1]
  · -- t < 1/2 : from the k-tooth upper end
    have h := hhi.1
    rw [lt_div_iff₀ hkR] at h
    nlinarith [hjabs.1, abs_nonneg ((j : ℝ)), hk1]

/-- **The assembled credit**: for every `J` with `14J ≤ 2k+1`, the trapezoid sum is a
lower bound for the pair-overlap measure on the unit window. -/
theorem consecutive_overlap_credit (k J : ℕ) (hk : 1 ≤ k)
    (hJ : 14 * (J : ℝ) ≤ 2 * k + 1) :
    ∑ j ∈ Finset.Icc (-(J : ℤ)) (J : ℤ), ENNReal.ofReal (tentAbs k j)
      ≤ volume ((dangerR k ∩ dangerR (k + 1)) ∩ Ioo (-(1:ℝ)/2) (1/2)) := by
  have habsJ : ∀ j ∈ Finset.Icc (-(J : ℤ)) (J : ℤ), 14 * |(j : ℝ)| ≤ 2 * k + 1 := by
    intro j hj
    rw [Finset.mem_Icc] at hj
    have h1 : |(j : ℤ)| ≤ (J : ℤ) := abs_le.mpr hj
    have h2 : |(j : ℝ)| ≤ (J : ℝ) := by exact_mod_cast (by exact_mod_cast h1 : ((|j| : ℤ) : ℝ) ≤ (J : ℝ))
    linarith
  calc ∑ j ∈ Finset.Icc (-(J : ℤ)) (J : ℤ), ENNReal.ofReal (tentAbs k j)
      ≤ ∑ j ∈ Finset.Icc (-(J : ℤ)) (J : ℤ), volume (comp k j) := by
        apply Finset.sum_le_sum
        intro j hj
        by_cases h0 : j = 0
        · subst h0
          rw [tentAbs, if_pos rfl, volume_comp_zero k hk]
        · rw [tentAbs, if_neg h0]
          exact volume_comp_tent k hk j h0 (habsJ j hj)
    _ = volume (⋃ j ∈ Finset.Icc (-(J : ℤ)) (J : ℤ), comp k j) := by
        rw [measure_biUnion_finset]
        · intro i hi j hj hij
          rcases lt_or_gt_of_ne hij with h | h
          · exact comp_disjoint_of_lt k hk h
          · exact (comp_disjoint_of_lt k hk h).symm
        · intro j _
          exact measurableSet_Ioo
    _ ≤ volume ((dangerR k ∩ dangerR (k + 1)) ∩ Ioo (-(1:ℝ)/2) (1/2)) := by
        apply measure_mono
        refine Set.iUnion₂_subset fun j hj => ?_
        intro t ht
        exact ⟨comp_subset_danger k hk j ht, comp_subset_window k hk j (habsJ j hj) ht⟩

#print axioms consecutive_overlap_credit

end

end LonelyRunner.LRC14.Arcs
