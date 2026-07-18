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

/-! ## The cast bridge (boxeph-S75): the credit in closed form

The symmetric trapezoid sum evaluates to `1/49 + r(6−r)/(49k(k+1))`, `r = k mod 7`
— the LEM-044 closed form, now as a volume lower bound.  The six-cycle of nonzero
residues mod 7 organizes the bookkeeping: `tentAbs` pairs antipodal residues
(`tentAbs_neg`), and the excess factor `r(6−r)` is the product of antipodes. -/

theorem tentAbs_neg (k : ℕ) (j : ℤ) : tentAbs k (-j) = tentAbs k j := by
  rcases eq_or_ne j 0 with rfl | hj
  · simp
  · rw [tentAbs, tentAbs, if_neg (neg_ne_zero.mpr hj), if_neg hj]
    have habs : |((-j : ℤ) : ℝ)| = |(j : ℝ)| := by push_cast; exact abs_neg _
    rw [habs]

theorem tentAbs_nonneg (k : ℕ) (hk : 1 ≤ k) (j : ℤ)
    (hj14 : 14 * |(j : ℝ)| ≤ 2 * k + 1) : 0 ≤ tentAbs k j := by
  have hkR : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk
  rw [tentAbs]
  split_ifs with h
  · positivity
  · apply div_nonneg _ (by positivity)
    linarith [abs_nonneg ((j : ℝ))]

theorem sum_tent_symm (k : ℕ) : ∀ J : ℕ,
    ∑ j ∈ Finset.Icc (-(J : ℤ)) (J : ℤ), tentAbs k j
      = tentAbs k 0 + 2 * ∑ j ∈ Finset.Icc (1 : ℤ) (J : ℤ), tentAbs k j := by
  intro J
  induction J with
  | zero => simp
  | succ n ih =>
    have hins1 : Finset.Icc (-((n + 1 : ℕ) : ℤ)) ((n + 1 : ℕ) : ℤ)
        = insert (-((n + 1 : ℕ) : ℤ))
            (insert ((n + 1 : ℕ) : ℤ) (Finset.Icc (-(n : ℤ)) (n : ℤ))) := by
      ext x
      simp only [Finset.mem_Icc, Finset.mem_insert]
      push_cast
      omega
    have hins2 : Finset.Icc (1 : ℤ) ((n + 1 : ℕ) : ℤ)
        = insert ((n + 1 : ℕ) : ℤ) (Finset.Icc (1 : ℤ) (n : ℤ)) := by
      ext x
      simp only [Finset.mem_Icc, Finset.mem_insert]
      push_cast
      omega
    rw [hins1]
    rw [Finset.sum_insert (by
      simp only [Finset.mem_insert, Finset.mem_Icc]
      push_cast
      omega)]
    rw [Finset.sum_insert (by
      simp only [Finset.mem_Icc]
      push_cast
      omega)]
    rw [hins2]
    rw [Finset.sum_insert (by
      simp only [Finset.mem_Icc]
      push_cast
      omega)]
    rw [ih, tentAbs_neg]
    ring

theorem sum_lin_Icc (c : ℝ) : ∀ q : ℕ,
    ∑ j ∈ Finset.Icc (1 : ℤ) (q : ℤ), (c - (j : ℝ))
      = (q : ℝ) * c - (q : ℝ) * ((q : ℝ) + 1) / 2 := by
  intro q
  induction q with
  | zero => simp
  | succ n ih =>
    have hins : Finset.Icc (1 : ℤ) ((n + 1 : ℕ) : ℤ)
        = insert ((n + 1 : ℕ) : ℤ) (Finset.Icc (1 : ℤ) (n : ℤ)) := by
      ext x
      simp only [Finset.mem_Icc, Finset.mem_insert]
      push_cast
      omega
    rw [hins]
    rw [Finset.sum_insert (by
      simp only [Finset.mem_Icc]
      push_cast
      omega)]
    rw [ih]
    push_cast
    ring

/-- **The evaluated credit**: the symmetric trapezoid sum in closed form. -/
theorem sum_tent_eval (k : ℕ) (hk : 1 ≤ k) :
    ∑ j ∈ Finset.Icc (-((k / 7 : ℕ) : ℤ)) ((k / 7 : ℕ) : ℤ), tentAbs k j
      = 1 / 49 + ((k % 7 : ℕ) : ℝ) * (6 - ((k % 7 : ℕ) : ℝ))
          / (49 * (k : ℝ) * ((k : ℝ) + 1)) := by
  have hkR : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk
  have hk0 : ((k : ℝ)) ≠ 0 := ne_of_gt hkR
  have hk10 : ((k : ℝ) + 1) ≠ 0 := by positivity
  have hkqr : (k : ℝ) = 7 * ((k / 7 : ℕ) : ℝ) + ((k % 7 : ℕ) : ℝ) := by
    exact_mod_cast (Nat.div_add_mod k 7).symm
  rw [sum_tent_symm]
  have hcong : ∀ j ∈ Finset.Icc (1 : ℤ) ((k / 7 : ℕ) : ℤ), tentAbs k j
      = ((2 * (k : ℝ) + 1) / 14 - (j : ℝ)) / ((k : ℝ) * ((k : ℝ) + 1)) := by
    intro j hj
    rw [Finset.mem_Icc] at hj
    have hj0 : j ≠ 0 := by omega
    have hjpos : (0 : ℝ) < (j : ℝ) := by exact_mod_cast (by omega : (0 : ℤ) < j)
    rw [tentAbs, if_neg hj0, abs_of_pos hjpos]
  rw [Finset.sum_congr rfl hcong, ← Finset.sum_div, sum_lin_Icc]
  rw [tentAbs, if_pos rfl]
  have h70 : (7 * ((k / 7 : ℕ) : ℝ) + ((k % 7 : ℕ) : ℝ)) ≠ 0 := by
    rw [← hkqr]; exact hk0
  have h71 : (7 * ((k / 7 : ℕ) : ℝ) + ((k % 7 : ℕ) : ℝ) + 1) ≠ 0 := by
    rw [← hkqr]; exact hk10
  rw [hkqr]
  field_simp [h70, h71]
  ring

/-- **The closed-form credit theorem**: the LEM-044 consecutive value is a volume
lower bound for the pair overlap on the unit window — the concrete credit for
`good_pos_of_path_credits`, end to end. -/
theorem consecutive_credit_closed (k : ℕ) (hk : 1 ≤ k) :
    ENNReal.ofReal (1 / 49 + ((k % 7 : ℕ) : ℝ) * (6 - ((k % 7 : ℕ) : ℝ))
        / (49 * (k : ℝ) * ((k : ℝ) + 1)))
      ≤ volume ((dangerR k ∩ dangerR (k + 1)) ∩ Ioo (-(1:ℝ)/2) (1/2)) := by
  have hdiv : ((k / 7 : ℕ) : ℝ) * 7 ≤ (k : ℝ) := by
    exact_mod_cast Nat.div_mul_le_self k 7
  have hJ : 14 * ((k / 7 : ℕ) : ℝ) ≤ 2 * k + 1 := by linarith
  have hnn : ∀ j ∈ Finset.Icc (-((k / 7 : ℕ) : ℤ)) ((k / 7 : ℕ) : ℤ),
      0 ≤ tentAbs k j := by
    intro j hj
    rw [Finset.mem_Icc] at hj
    apply tentAbs_nonneg k hk
    have hcl : ((j : ℝ)) ≤ ((k / 7 : ℕ) : ℝ) := by
      have h : ((j : ℤ) : ℝ) ≤ ((((k / 7 : ℕ) : ℤ)) : ℝ) := Int.cast_le.mpr hj.2
      push_cast at h
      exact h
    have hcg : (-((k / 7 : ℕ) : ℝ)) ≤ (j : ℝ) := by
      have h : (((-((k / 7 : ℕ) : ℤ)) : ℤ) : ℝ) ≤ ((j : ℤ) : ℝ) := Int.cast_le.mpr hj.1
      push_cast at h
      exact h
    have h2 : |(j : ℝ)| ≤ ((k / 7 : ℕ) : ℝ) := abs_le.mpr ⟨hcg, hcl⟩
    linarith
  calc ENNReal.ofReal (1 / 49 + ((k % 7 : ℕ) : ℝ) * (6 - ((k % 7 : ℕ) : ℝ))
          / (49 * (k : ℝ) * ((k : ℝ) + 1)))
      = ENNReal.ofReal (∑ j ∈ Finset.Icc (-((k / 7 : ℕ) : ℤ)) ((k / 7 : ℕ) : ℤ),
          tentAbs k j) := by rw [sum_tent_eval k hk]
    _ = ∑ j ∈ Finset.Icc (-((k / 7 : ℕ) : ℤ)) ((k / 7 : ℕ) : ℤ),
          ENNReal.ofReal (tentAbs k j) := ENNReal.ofReal_sum_of_nonneg hnn
    _ ≤ volume ((dangerR k ∩ dangerR (k + 1)) ∩ Ioo (-(1:ℝ)/2) (1/2)) :=
        consecutive_overlap_credit k (k / 7) hk hJ

/-! ## The per-runner mass bound (boxeph-S76): the skeleton's other hypothesis -/

/-- One tooth of runner `v`. -/
def tooth (v : ℕ) (m : ℤ) : Set ℝ :=
  Ioo (((m : ℝ) - 1/14) / v) (((m : ℝ) + 1/14) / v)

theorem mem_tooth_of_danger {v : ℕ} (hv : 1 ≤ v) {t : ℝ} {m : ℤ}
    (h : |(v : ℝ) * t - m| < 1/14) : t ∈ tooth v m := by
  have hvR : (0 : ℝ) < (v : ℝ) := by exact_mod_cast hv
  rw [abs_lt] at h
  constructor
  · rw [div_lt_iff₀ hvR]
    linarith [h.1]
  · rw [lt_div_iff₀ hvR]
    linarith [h.2]

theorem danger_window_cover (v : ℕ) (hv : 1 ≤ v) :
    dangerR v ∩ Ioo (-(1:ℝ)/2) (1/2)
      ⊆ ⋃ m ∈ Finset.Icc (-(((v + 1) / 2 : ℕ) : ℤ)) (((v + 1) / 2 : ℕ) : ℤ),
          (tooth v m ∩ Ioo (-(1:ℝ)/2) (1/2)) := by
  rintro t ⟨⟨m, hm⟩, htW⟩
  have hvR : (0 : ℝ) < (v : ℝ) := by exact_mod_cast hv
  have h2M : v ≤ 2 * ((v + 1) / 2) := by omega
  have h2MR : (v : ℝ) ≤ 2 * (((v + 1) / 2 : ℕ) : ℝ) := by exact_mod_cast h2M
  obtain ⟨hW1, hW2⟩ := htW
  have habs := abs_lt.mp hm
  have hvt1 : (v : ℝ) * (-(1:ℝ)/2) < (v : ℝ) * t := by
    apply mul_lt_mul_of_pos_left _ hvR
    linarith
  have hvt2 : (v : ℝ) * t < (v : ℝ) * (1/2) := by
    apply mul_lt_mul_of_pos_left _ hvR
    linarith
  simp only [Set.mem_iUnion, Finset.mem_Icc]
  refine ⟨m, ⟨?_, ?_⟩, mem_tooth_of_danger hv hm, hW1, hW2⟩
  · by_contra hcon
    push_neg at hcon
    have hmZ : (m : ℤ) ≤ -(((v + 1) / 2 : ℕ) : ℤ) - 1 := by omega
    have hmR : (m : ℝ) ≤ -(((v + 1) / 2 : ℕ) : ℝ) - 1 := by exact_mod_cast hmZ
    nlinarith [habs.1, habs.2]
  · by_contra hcon
    push_neg at hcon
    have hmZ : (((v + 1) / 2 : ℕ) : ℤ) + 1 ≤ m := by omega
    have hmR : (((v + 1) / 2 : ℕ) : ℝ) + 1 ≤ (m : ℝ) := by
      have h0 : (((((v + 1) / 2 : ℕ) : ℤ) + 1 : ℤ) : ℝ) ≤ ((m : ℤ) : ℝ) :=
        Int.cast_le.mpr hmZ
      push_cast at h0
      exact h0
    nlinarith [habs.1, habs.2]

theorem tooth_measure_le (v : ℕ) (hv : 1 ≤ v) (m : ℤ) :
    volume (tooth v m ∩ Ioo (-(1:ℝ)/2) (1/2)) ≤ ENNReal.ofReal (1 / (7 * v)) := by
  have hvR : (0 : ℝ) < (v : ℝ) := by exact_mod_cast hv
  calc volume (tooth v m ∩ Ioo (-(1:ℝ)/2) (1/2))
      ≤ volume (tooth v m) := measure_mono Set.inter_subset_left
    _ = ENNReal.ofReal (((m : ℝ) + 1/14) / v - ((m : ℝ) - 1/14) / v) := Real.volume_Ioo
    _ = ENNReal.ofReal (1 / (7 * v)) := by
        congr 1
        field_simp
        ring

theorem tooth_top_le (v : ℕ) (_hv : 1 ≤ v) (M : ℤ) :
    volume (tooth v M ∩ Ioo (-(1:ℝ)/2) (1/2))
      ≤ ENNReal.ofReal (1/2 - ((M : ℝ) - 1/14) / v) := by
  have hsub : tooth v M ∩ Ioo (-(1:ℝ)/2) (1/2)
      ⊆ Ioo (((M : ℝ) - 1/14) / v) (1/2) := by
    rintro t ⟨⟨h1, _⟩, ⟨_, h4⟩⟩
    exact ⟨h1, h4⟩
  calc volume (tooth v M ∩ Ioo (-(1:ℝ)/2) (1/2))
      ≤ volume (Ioo (((M : ℝ) - 1/14) / v) (1/2)) := measure_mono hsub
    _ = ENNReal.ofReal (1/2 - ((M : ℝ) - 1/14) / v) := Real.volume_Ioo

theorem tooth_bot_le (v : ℕ) (_hv : 1 ≤ v) (M : ℤ) :
    volume (tooth v (-M) ∩ Ioo (-(1:ℝ)/2) (1/2))
      ≤ ENNReal.ofReal (1/2 - ((M : ℝ) - 1/14) / v) := by
  have hsub : tooth v (-M) ∩ Ioo (-(1:ℝ)/2) (1/2)
      ⊆ Ioo (-(1:ℝ)/2) ((((-M : ℤ) : ℝ) + 1/14) / v) := by
    rintro t ⟨⟨_, h2⟩, ⟨h3, _⟩⟩
    exact ⟨h3, h2⟩
  calc volume (tooth v (-M) ∩ Ioo (-(1:ℝ)/2) (1/2))
      ≤ volume (Ioo (-(1:ℝ)/2) ((((-M : ℤ) : ℝ) + 1/14) / v)) := measure_mono hsub
    _ = ENNReal.ofReal ((((-M : ℤ) : ℝ) + 1/14) / v - (-(1:ℝ)/2)) := Real.volume_Ioo
    _ = ENNReal.ofReal (1/2 - ((M : ℝ) - 1/14) / v) := by
        congr 1
        push_cast
        ring

/-- **The per-runner mass bound**: each danger set carries at most `1/7` of the
unit window — the second hypothesis of `good_pos_of_path_credits`. -/
theorem danger_measure_le (v : ℕ) (hv : 1 ≤ v) :
    volume (dangerR v ∩ Ioo (-(1:ℝ)/2) (1/2)) ≤ ENNReal.ofReal (1/7) := by
  have hvR : (0 : ℝ) < (v : ℝ) := by exact_mod_cast hv
  set M : ℕ := (v + 1) / 2 with hM
  have hM1 : 1 ≤ M := by omega
  have hins : Finset.Icc (-(M : ℤ)) (M : ℤ)
      = insert (-(M : ℤ)) (insert (M : ℤ)
          (Finset.Icc (-(M : ℤ) + 1) ((M : ℤ) - 1))) := by
    ext x
    simp only [Finset.mem_Icc, Finset.mem_insert]
    omega
  have hcard : (Finset.Icc (-(M : ℤ) + 1) ((M : ℤ) - 1)).card = 2 * M - 1 := by
    rw [Int.card_Icc]
    omega
  have hsum : ∑ m ∈ Finset.Icc (-(M : ℤ)) (M : ℤ),
        volume (tooth v m ∩ Ioo (-(1:ℝ)/2) (1/2))
      ≤ ENNReal.ofReal (((2 * M - 1 : ℕ) : ℝ) * (1 / (7 * v)))
        + 2 * ENNReal.ofReal (1/2 - ((M : ℝ) - 1/14) / v) := by
    rw [hins]
    rw [Finset.sum_insert (by
      simp only [Finset.mem_insert, Finset.mem_Icc]
      omega)]
    rw [Finset.sum_insert (by
      simp only [Finset.mem_Icc]
      omega)]
    have hbody : ∑ m ∈ Finset.Icc (-(M : ℤ) + 1) ((M : ℤ) - 1),
          volume (tooth v m ∩ Ioo (-(1:ℝ)/2) (1/2))
        ≤ ENNReal.ofReal (((2 * M - 1 : ℕ) : ℝ) * (1 / (7 * v))) := by
      calc ∑ m ∈ Finset.Icc (-(M : ℤ) + 1) ((M : ℤ) - 1),
            volume (tooth v m ∩ Ioo (-(1:ℝ)/2) (1/2))
          ≤ ∑ _m ∈ Finset.Icc (-(M : ℤ) + 1) ((M : ℤ) - 1),
              ENNReal.ofReal (1 / (7 * v)) :=
            Finset.sum_le_sum fun m _ => tooth_measure_le v hv m
        _ = ENNReal.ofReal (((2 * M - 1 : ℕ) : ℝ) * (1 / (7 * v))) := by
            rw [Finset.sum_const, hcard, nsmul_eq_mul,
              ← ENNReal.ofReal_natCast (2 * M - 1),
              ← ENNReal.ofReal_mul (by positivity)]
    have htop := tooth_top_le v hv (M : ℤ)
    have hbot := tooth_bot_le v hv (M : ℤ)
    calc volume (tooth v (-(M : ℤ)) ∩ Ioo (-(1:ℝ)/2) (1/2))
          + (volume (tooth v (M : ℤ) ∩ Ioo (-(1:ℝ)/2) (1/2))
            + ∑ m ∈ Finset.Icc (-(M : ℤ) + 1) ((M : ℤ) - 1),
                volume (tooth v m ∩ Ioo (-(1:ℝ)/2) (1/2)))
        ≤ ENNReal.ofReal (1/2 - ((M : ℝ) - 1/14) / v)
          + (ENNReal.ofReal (1/2 - ((M : ℝ) - 1/14) / v)
            + ENNReal.ofReal (((2 * M - 1 : ℕ) : ℝ) * (1 / (7 * v)))) := by
          exact add_le_add hbot (add_le_add htop hbody)
      _ = ENNReal.ofReal (((2 * M - 1 : ℕ) : ℝ) * (1 / (7 * v)))
          + 2 * ENNReal.ofReal (1/2 - ((M : ℝ) - 1/14) / v) := by
          ring
  have hchain := (measure_mono (danger_window_cover v hv)).trans
    ((measure_biUnion_finset_le _ _).trans hsum)
  refine hchain.trans ?_
  -- final arithmetic, split by parity
  rcases Nat.even_or_odd v with ⟨w, hw⟩ | ⟨w, hw⟩
  · -- even v = w + w : M = w, boundary teeth carry exactly 1/(14v) each
    have hMw : M = w := by omega
    have hw1 : 1 ≤ w := by omega
    have hw0 : (0:ℝ) < (w : ℝ) := by exact_mod_cast hw1
    have hwR : ((w : ℝ)) * 2 = (v : ℝ) := by
      exact_mod_cast (by omega : w * 2 = v)
    have hB : (1:ℝ)/2 - ((M : ℝ) - 1/14) / v = 1 / (14 * v) := by
      rw [hMw, ← hwR]
      field_simp
      ring
    have hA : ((2 * M - 1 : ℕ) : ℝ) = (v : ℝ) - 1 := by
      have h1 : 2 * M - 1 = v - 1 := by omega
      rw [h1, Nat.cast_sub hv, Nat.cast_one]
    have hv1R : (1:ℝ) ≤ (v : ℝ) := by exact_mod_cast hv
    have hfin : ((v:ℝ) - 1) * (1 / (7 * (v:ℝ))) + (1 / (14 * (v:ℝ)) + 1 / (14 * (v:ℝ)))
        = 1/7 := by
      field_simp
      ring
    rw [hB, hA, two_mul,
      ← ENNReal.ofReal_add (by positivity) (by positivity),
      ← ENNReal.ofReal_add (mul_nonneg (by linarith) (by positivity)) (by positivity),
      hfin]
  · -- odd v = 2w + 1 : M = w + 1, boundary teeth are empty (nonpositive length)
    have hMw : M = w + 1 := by omega
    have hvw : (v : ℝ) = 2 * (w : ℝ) + 1 := by exact_mod_cast hw
    have hB : (1:ℝ)/2 - ((M : ℝ) - 1/14) / v ≤ 0 := by
      rw [hMw, sub_nonpos, le_div_iff₀ hvR, hvw]
      push_cast
      linarith
    have hA : ((2 * M - 1 : ℕ) : ℝ) = (v : ℝ) := by
      exact_mod_cast (by omega : 2 * M - 1 = v)
    have hfin : (v:ℝ) * (1 / (7 * (v:ℝ))) = 1/7 := by
      field_simp
    rw [ENNReal.ofReal_of_nonpos hB, mul_zero, add_zero, hA, hfin]

/-! ## The Lipschitz margin converter (boxeph-S79): M-floors become μ₀-floors

The missing link of the certificate chain: a margin-`M` instant spawns a good
interval of radius `(M − 1/14)/vmax`, so ANY quantitative bound `M > 1/14`
yields `μ₀ ≥ 2(M − 1/14)/vmax` — which the weight-1 sampling bridge
(LRCGridSampling + THM-984) converts into a census certificate at an explicit
modulus.  M-floor → μ₀-floor → liveCount floor → certificate, all mechanical. -/

theorem good_interval_of_margin (n : ℕ) (v : Fin n → ℕ) (vmax : ℕ)
    (hmax : ∀ i, v i ≤ vmax) (hvmax : 1 ≤ vmax) (t₀ M : ℝ)
    (hM : ∀ (i : Fin n) (m : ℤ), M ≤ |(v i : ℝ) * t₀ - m|) (hM14 : 1/14 < M) :
    Set.Ioo (t₀ - (M - 1/14) / vmax) (t₀ + (M - 1/14) / vmax)
      ⊆ (⋃ i, dangerR (v i))ᶜ := by
  have hvmaxR : (0 : ℝ) < (vmax : ℝ) := by exact_mod_cast hvmax
  rintro t ⟨ht1, ht2⟩
  simp only [Set.mem_compl_iff, Set.mem_iUnion, not_exists]
  intro i hdanger
  obtain ⟨m, hm⟩ := hdanger
  have hviR : ((v i : ℕ) : ℝ) ≤ (vmax : ℝ) := by exact_mod_cast hmax i
  have habs : |t - t₀| < (M - 1/14) / vmax := by
    rw [abs_lt]
    constructor <;> linarith
  have hkey : |(v i : ℝ) * t₀ - m|
      ≤ |(v i : ℝ) * t - m| + (v i : ℝ) * |t - t₀| := by
    have hdec : (v i : ℝ) * t₀ - m
        = ((v i : ℝ) * t - m) - (v i : ℝ) * (t - t₀) := by ring
    calc |(v i : ℝ) * t₀ - m|
        = |((v i : ℝ) * t - m) - (v i : ℝ) * (t - t₀)| := by rw [hdec]
      _ ≤ |(v i : ℝ) * t - m| + |(v i : ℝ) * (t - t₀)| := abs_sub _ _
      _ = |(v i : ℝ) * t - m| + (v i : ℝ) * |t - t₀| := by
          rw [abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ ((v i : ℕ) : ℝ))]
  have hsmall : (v i : ℝ) * |t - t₀| < M - 1/14 := by
    calc (v i : ℝ) * |t - t₀| ≤ (vmax : ℝ) * |t - t₀| :=
          mul_le_mul_of_nonneg_right hviR (abs_nonneg _)
      _ < (vmax : ℝ) * ((M - 1/14) / vmax) :=
          mul_lt_mul_of_pos_left habs hvmaxR
      _ = M - 1/14 := by field_simp
  have := hM i m
  linarith

/-- The measure corollary: `μ₀ ≥ 2(M − 1/14)/vmax`. -/
theorem good_measure_of_margin (n : ℕ) (v : Fin n → ℕ) (vmax : ℕ)
    (hmax : ∀ i, v i ≤ vmax) (hvmax : 1 ≤ vmax) (t₀ M : ℝ)
    (hM : ∀ (i : Fin n) (m : ℤ), M ≤ |(v i : ℝ) * t₀ - m|) (hM14 : 1/14 < M) :
    ENNReal.ofReal (2 * (M - 1/14) / vmax) ≤ volume ((⋃ i, dangerR (v i))ᶜ) := by
  have hIoo : volume (Set.Ioo (t₀ - (M - 1/14) / vmax) (t₀ + (M - 1/14) / vmax))
      = ENNReal.ofReal (2 * (M - 1/14) / vmax) := by
    rw [Real.volume_Ioo]
    congr 1
    ring
  rw [← hIoo]
  exact measure_mono (good_interval_of_margin n v vmax hmax hvmax t₀ M hM hM14)

#print axioms consecutive_overlap_credit
#print axioms consecutive_credit_closed
#print axioms danger_measure_le
#print axioms good_measure_of_margin

end

end LonelyRunner.LRC14.Arcs
