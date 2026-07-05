/-
  TournamentH7.LRCMergeGridAttainment — THE MERGE-GRID ATTAINMENT THEOREM
  (klein-2026-07-05-S137, HYP-4114; THM-592's formalization, kps-S3's suggested-next).

  The max-min profile  M(V) = max_t min_i dist(vᵢt, ℤ)  of a finite positive family is
  ATTAINED, and attained at a MERGE-GRID point  t* = m/(vᵢ + vⱼ)  (i = j giving the
  half-integer peaks).  This is the theorem beneath every exact-M sweep in the project:
  with it, scanning the finite grid `{m/(vᵢ+vⱼ)}` is a KERNEL-CHECKABLE evaluation of M.

  Stage A (this section): the distance-to-ℤ toolkit (`distZ`, nearest-point property,
  1-Lipschitz continuity, integer periodicity), the profile as a finite inf, its
  continuity/periodicity/positivity, and ATTAINMENT via compactness + period transport.
-/
import Mathlib
import TournamentH7.LonelyRunner

namespace LonelyRunner
namespace MergeGrid

/-- Distance to the nearest integer, through `round`. -/
noncomputable def distZ (x : ℝ) : ℝ := |x - round x|

theorem distZ_nonneg (x : ℝ) : 0 ≤ distZ x := abs_nonneg _

theorem distZ_le_half (x : ℝ) : distZ x ≤ 1/2 := by
  simpa [distZ] using abs_sub_round x

/-- `round` is the nearest integer: `distZ x ≤ |x − m|` for every integer `m`. -/
theorem distZ_le (x : ℝ) (m : ℤ) : distZ x ≤ |x - m| := by
  rcases eq_or_ne m (round x) with rfl | hne
  · simp [distZ]
  · have h1 : (1 : ℝ) ≤ |(m : ℝ) - round x| := by
      have hne' : m - round x ≠ 0 := sub_ne_zero.mpr hne
      have hz : (1 : ℤ) ≤ |m - round x| := Int.one_le_abs hne'
      have hz' : ((1 : ℤ) : ℝ) ≤ ((|m - round x| : ℤ) : ℝ) := by exact_mod_cast hz
      rw [Int.cast_abs] at hz'
      push_cast at hz'
      linarith
    have h2 : |x - round x| ≤ 1/2 := abs_sub_round x
    have htri : |(m : ℝ) - round x| ≤ |(m : ℝ) - x| + |x - round x| := abs_sub_le _ _ _
    have hsw : |(m : ℝ) - x| = |x - (m : ℝ)| := abs_sub_comm _ _
    unfold distZ
    linarith
  
/-- `distZ` is 1-Lipschitz. -/
theorem distZ_lipschitz (x y : ℝ) : distZ x ≤ |x - y| + distZ y := by
  have h1 : distZ x ≤ |x - round y| := distZ_le x (round y)
  have h2 : |x - (round y : ℝ)| ≤ |x - y| + |y - round y| := abs_sub_le _ _ _
  unfold distZ at *
  linarith

theorem distZ_continuous : Continuous distZ := by
  have hlip : LipschitzWith 1 distZ := by
    apply LipschitzWith.of_dist_le_mul
    intro x y
    rw [Real.dist_eq, Real.dist_eq, NNReal.coe_one, one_mul]
    rw [abs_le]
    constructor
    · have := distZ_lipschitz y x
      have habs : |y - x| = |x - y| := abs_sub_comm _ _
      linarith [habs ▸ this]
    · linarith [distZ_lipschitz x y]
  exact hlip.continuous

/-- Integer shift invariance. -/
theorem distZ_add_int (x : ℝ) (k : ℤ) : distZ (x + k) = distZ x := by
  unfold distZ
  rw [round_add_intCast]
  push_cast
  ring_nf

/-- The min-distance profile of a positive family (nonempty index). -/
noncomputable def profile {n : ℕ} (v : Fin (n+1) → ℤ) (t : ℝ) : ℝ :=
  Finset.univ.inf' Finset.univ_nonempty fun i => distZ ((v i : ℝ) * t)

theorem profile_continuous {n : ℕ} (v : Fin (n+1) → ℤ) : Continuous (profile v) := by
  apply Continuous.finset_inf'_apply
  intro i _
  exact distZ_continuous.comp (continuous_const.mul continuous_id)

theorem profile_periodic {n : ℕ} (v : Fin (n+1) → ℤ) :
    Function.Periodic (profile v) 1 := by
  intro t
  unfold profile
  congr 1
  funext i
  have hexp : (v i : ℝ) * (t + 1) = (v i : ℝ) * t + ((v i : ℤ) : ℝ) := by
    push_cast
    ring
  rw [hexp, distZ_add_int]

/-- The profile is somewhere positive: at `t₀ = 1/(2W+1)` every runner sits in
`(0, 1/2)`. -/
theorem profile_pos {n : ℕ} (v : Fin (n+1) → ℤ) (hv : ∀ i, 0 < v i) (W : ℤ)
    (hW : ∀ i, v i ≤ W) :
    (1 : ℝ) / (2 * (W : ℝ) + 1) ≤ profile v (1 / (2 * (W : ℝ) + 1)) := by
  have hWpos : (0 : ℝ) < (W : ℝ) := by
    have := hv 0
    have := hW 0
    exact_mod_cast lt_of_lt_of_le (hv 0) (hW 0)
  have hden : (0 : ℝ) < 2 * (W : ℝ) + 1 := by linarith
  unfold profile
  rw [Finset.le_inf'_iff]
  intro i _
  have hviR : (0 : ℝ) < (v i : ℝ) := by exact_mod_cast hv i
  have hviW : (v i : ℝ) ≤ (W : ℝ) := by exact_mod_cast hW i
  set x : ℝ := (v i : ℝ) * (1 / (2 * (W : ℝ) + 1)) with hx
  have hxpos : 0 < x := by
    rw [hx]
    positivity
  have hxlt : x < 1/2 := by
    rw [hx]
    rw [mul_one_div, div_lt_iff₀ hden]
    linarith
  have hround : round x = 0 := by
    apply round_eq_zero_iff.mpr
    exact Set.mem_Ico.mpr ⟨by linarith, by linarith⟩
  unfold distZ
  rw [hround]
  push_cast
  rw [sub_zero, abs_of_pos hxpos, hx]
  have h1 : (1 : ℝ) ≤ (v i : ℝ) := by exact_mod_cast hv i
  calc (1 : ℝ) / (2 * (W : ℝ) + 1) = 1 * (1 / (2 * (W : ℝ) + 1)) := by ring
    _ ≤ (v i : ℝ) * (1 / (2 * (W : ℝ) + 1)) :=
        mul_le_mul_of_nonneg_right h1 (le_of_lt (div_pos one_pos hden))

/-- **ATTAINMENT**: the profile attains its global maximum. -/
theorem profile_attains_max {n : ℕ} (v : Fin (n+1) → ℤ) :
    ∃ tstar : ℝ, ∀ t : ℝ, profile v t ≤ profile v tstar := by
  obtain ⟨tstar, -, hmax⟩ :=
    (isCompact_Icc (a := (0:ℝ)) (b := 1)).exists_isMaxOn
      (Set.nonempty_Icc.mpr (by norm_num))
      (profile_continuous v).continuousOn
  refine ⟨tstar, fun t => ?_⟩
  have hper : profile v (t - (⌊t⌋ : ℝ) * 1) = profile v t :=
    (profile_periodic v).sub_int_mul_eq ⌊t⌋
  rw [mul_one] at hper
  rw [← hper]
  have hf : t - (⌊t⌋ : ℝ) = Int.fract t := Int.self_sub_floor t
  refine hmax (Set.mem_Icc.mpr ⟨?_, ?_⟩)
  · rw [hf]
    exact Int.fract_nonneg t
  · rw [hf]
    exact le_of_lt (Int.fract_lt_one t)

/-! ## Stage B/C: the merge-grid classification at the max -/

/-- An integer within `1/2` of `x` is `round x`. -/
theorem round_eq_of_abs_lt_half {x : ℝ} {k : ℤ} (h : |x - k| < 1/2) : round x = k := by
  have h2 : |x - round x| ≤ 1/2 := abs_sub_round x
  have h3 : |(k : ℝ) - round x| < 1 := by
    calc |(k : ℝ) - round x| ≤ |(k : ℝ) - x| + |x - round x| := abs_sub_le _ _ _
      _ < 1/2 + 1/2 := by
          rw [abs_sub_comm (k : ℝ) x]
          linarith
      _ = 1 := by norm_num
  have h4 : |k - round x| < 1 := by
    have : ((|k - round x| : ℤ) : ℝ) < 1 := by
      rw [Int.cast_abs]
      push_cast
      exact h3
    exact_mod_cast this
  rw [abs_lt] at h4
  omega

/-- **THE MERGE-GRID ATTAINMENT THEOREM** (THM-592, Lean face): the profile maximum of
a positive family is attained at a merge-grid point `t* = m/(vᵢ + vⱼ)` — with `i = j`
delivering the half-integer peaks.  Consequently every exact-M sweep needs only the
finite grid `{m/(vᵢ+vⱼ)}`. -/
theorem merge_grid_attainment {n : ℕ} (v : Fin (n+1) → ℤ) (hv : ∀ i, 0 < v i) :
    ∃ (tstar : ℝ) (i j : Fin (n+1)) (m : ℤ),
      ((v i + v j : ℤ) : ℝ) * tstar = (m : ℝ) ∧
      ∀ t : ℝ, profile v t ≤ profile v tstar := by
  obtain ⟨tstar, hmax⟩ := profile_attains_max v
  set M : ℝ := profile v tstar with hM
  -- positivity of M
  obtain ⟨imax, -, himax⟩ :=
    Finset.exists_max_image (Finset.univ : Finset (Fin (n+1))) v ⟨0, Finset.mem_univ 0⟩
  have hW : ∀ i, v i ≤ v imax := fun i => himax i (Finset.mem_univ i)
  have hWposR : (0 : ℝ) < (v imax : ℝ) := by exact_mod_cast hv imax
  have hMpos : 0 < M := by
    have h1 := profile_pos v hv (v imax) hW
    have h2 := hmax (1 / (2 * ((v imax : ℤ) : ℝ) + 1))
    have h3 : (0 : ℝ) < 1 / (2 * ((v imax : ℤ) : ℝ) + 1) := by positivity
    linarith
  -- M ≤ 1/2
  have hMhalf : M ≤ 1/2 := by
    have h1 : M ≤ distZ ((v 0 : ℝ) * tstar) := by
      rw [hM]
      unfold profile
      exact Finset.inf'_le _ (Finset.mem_univ 0)
    linarith [distZ_le_half ((v 0 : ℝ) * tstar)]
  rcases eq_or_lt_of_le hMhalf with hhalf | hlt
  · -- PEAK CASE: some active runner sits at distance exactly 1/2: t* = odd/(2vᵢ)
    obtain ⟨i, -, hieq⟩ := Finset.exists_mem_eq_inf'
      (Finset.univ_nonempty) (fun i => distZ ((v i : ℝ) * tstar))
    have hdi : distZ ((v i : ℝ) * tstar) = 1/2 := by
      rw [← hhalf, hM]
      unfold profile
      rw [hieq]
    set x : ℝ := (v i : ℝ) * tstar with hx
    have habs : |x - round x| = 1/2 := hdi
    rcases abs_eq (by norm_num : (0:ℝ) ≤ 1/2) |>.mp habs with hpm | hpm
    · refine ⟨tstar, i, i, 2 * round x + 1, ?_, hmax⟩
      push_cast
      rw [show ((v i : ℝ) + (v i : ℝ)) * tstar = 2 * x by rw [hx]; ring]
      linarith
    · refine ⟨tstar, i, i, 2 * round x - 1, ?_, hmax⟩
      push_cast
      rw [show ((v i : ℝ) + (v i : ℝ)) * tstar = 2 * x by rw [hx]; ring]
      linarith
  · -- INTERIOR CASE: M < 1/2; extract a right-descender and a left-ascender
    set ε₀ : ℝ := min ((1/2 - M) / (2 * (v imax : ℝ))) (M / (2 * (v imax : ℝ))) with hε₀
    have hε₀pos : 0 < ε₀ := by
      rw [hε₀]
      apply lt_min <;> positivity
    -- the per-k witness on the RIGHT
    have hright : ∀ k : ℕ, ∃ i : Fin (n+1),
        distZ ((v i : ℝ) * (tstar + ε₀ / (k+1))) ≤ M := by
      intro k
      have h1 : profile v (tstar + ε₀ / (k+1)) ≤ M := hmax _
      unfold profile at h1
      obtain ⟨i, -, hieq⟩ := Finset.exists_mem_eq_inf'
        (Finset.univ_nonempty) (fun i => distZ ((v i : ℝ) * (tstar + ε₀ / (k+1))))
      rw [hieq] at h1
      exact ⟨i, h1⟩
    -- the per-k witness on the LEFT
    have hleft : ∀ k : ℕ, ∃ j : Fin (n+1),
        distZ ((v j : ℝ) * (tstar - ε₀ / (k+1))) ≤ M := by
      intro k
      have h1 : profile v (tstar - ε₀ / (k+1)) ≤ M := hmax _
      unfold profile at h1
      obtain ⟨i, -, hieq⟩ := Finset.exists_mem_eq_inf'
        (Finset.univ_nonempty) (fun i => distZ ((v i : ℝ) * (tstar - ε₀ / (k+1))))
      rw [hieq] at h1
      exact ⟨i, h1⟩
    choose FR hFR using hright
    choose FL hFL using hleft
    obtain ⟨iR, hiR'⟩ := Finite.exists_infinite_fiber FR
    obtain ⟨iL, hiL'⟩ := Finite.exists_infinite_fiber FL
    have hiR : (FR ⁻¹' {iR}).Infinite := Set.infinite_coe_iff.mp hiR'
    have hiL : (FL ⁻¹' {iL}).Infinite := Set.infinite_coe_iff.mp hiL'
    have hVbound : ∀ i : Fin (n+1),
        (v i : ℝ) * ε₀ ≤ (1/2 - M)/2 ∧ (v i : ℝ) * ε₀ ≤ M/2 := by
      intro i
      have hvile : (v i : ℝ) ≤ (v imax : ℝ) := by exact_mod_cast hW i
      have hmain : ∀ C : ℝ, 0 < C → ε₀ ≤ C / (2 * (v imax : ℝ)) →
          (v i : ℝ) * ε₀ ≤ C/2 := by
        intro C hC hle
        calc (v i : ℝ) * ε₀ ≤ (v imax : ℝ) * ε₀ :=
              mul_le_mul_of_nonneg_right hvile (le_of_lt hε₀pos)
          _ ≤ (v imax : ℝ) * (C / (2 * (v imax : ℝ))) :=
              mul_le_mul_of_nonneg_left hle (le_of_lt hWposR)
          _ = C/2 := by field_simp
      exact ⟨hmain _ (by linarith) (by rw [hε₀]; exact min_le_left _ _),
             hmain _ (by linarith [hMpos]) (by rw [hε₀]; exact min_le_right _ _)⟩
    have hεk : ∀ (i : Fin (n+1)) (k : ℕ),
        0 < (v i : ℝ) * (ε₀/(k+1)) ∧ (v i : ℝ) * (ε₀/(k+1)) ≤ (v i : ℝ) * ε₀ := by
      intro i k
      have hvpos : (0:ℝ) < (v i : ℝ) := by exact_mod_cast hv i
      constructor
      · positivity
      · apply mul_le_mul_of_nonneg_left ?_ (le_of_lt hvpos)
        apply div_le_self (le_of_lt hε₀pos)
        have : (0:ℝ) ≤ (k : ℝ) := Nat.cast_nonneg k
        linarith
    -- THE RIGHT-DESCENDER: round(x_R) − x_R = M
    have hRfact : (round ((v iR : ℝ) * tstar) : ℝ) - (v iR : ℝ) * tstar = M := by
      set x : ℝ := (v iR : ℝ) * tstar with hxdef
      have hvpos : (0:ℝ) < (v iR : ℝ) := by exact_mod_cast hv iR
      have hd0 : M ≤ |x - round x| := by
        have h1 : M ≤ distZ x := by
          rw [hM]
          unfold profile
          exact Finset.inf'_le _ (Finset.mem_univ iR)
        simpa [distZ] using h1
      have hcore : ∀ k : ℕ, FR k = iR →
          |x + (v iR : ℝ) * (ε₀/(k+1)) - round x| ≤ M := by
        intro k hk
        have h1 := hFR k
        rw [hk] at h1
        have hexp : (v iR : ℝ) * (tstar + ε₀/(k+1)) = x + (v iR : ℝ) * (ε₀/(k+1)) := by
          rw [hxdef]
          ring
        rw [hexp] at h1
        unfold distZ at h1
        have hεp := (hεk iR k).1
        have hεle := (hεk iR k).2
        have hxm : |x - (round (x + (v iR : ℝ) * (ε₀/(k+1))) : ℤ)| < 1/2 := by
          have htri : |x - (round (x + (v iR : ℝ) * (ε₀/(k+1))) : ℝ)|
              ≤ |x - (x + (v iR : ℝ) * (ε₀/(k+1)))|
                + |x + (v iR : ℝ) * (ε₀/(k+1))
                    - round (x + (v iR : ℝ) * (ε₀/(k+1)))| := abs_sub_le _ _ _
          have hxy : |x - (x + (v iR : ℝ) * (ε₀/(k+1)))| = (v iR : ℝ) * (ε₀/(k+1)) := by
            rw [show x - (x + (v iR : ℝ) * (ε₀/(k+1))) = -((v iR : ℝ) * (ε₀/(k+1))) by ring,
              abs_neg, abs_of_pos hεp]
          have hb := (hVbound iR).1
          calc |x - (round (x + (v iR : ℝ) * (ε₀/(k+1))) : ℝ)|
              ≤ _ + _ := htri
            _ ≤ (v iR : ℝ) * (ε₀/(k+1)) + M := by rw [hxy]; linarith
            _ ≤ (v iR : ℝ) * ε₀ + M := by linarith
            _ ≤ (1/2 - M)/2 + M := by linarith
            _ < 1/2 := by linarith
        have hpin := round_eq_of_abs_lt_half hxm
        rw [hpin]
        exact h1
      obtain ⟨k₀, hk₀mem⟩ := hiR.nonempty
      have hk₀ : FR k₀ = iR := by simpa using hk₀mem
      have hsign : x - (round x : ℝ) ≤ -M := by
        by_contra hcon
        push Not at hcon
        have hge : M ≤ x - round x := by
          rcases abs_cases (x - (round x : ℝ)) with ⟨he, -⟩ | ⟨he, -⟩
          · linarith [hd0]
          · linarith [hd0]
        have hc := hcore k₀ hk₀
        have hεp := (hεk iR k₀).1
        have habs : |x + (v iR : ℝ) * (ε₀/(k₀+1)) - round x|
            = x + (v iR : ℝ) * (ε₀/(k₀+1)) - round x := abs_of_pos (by linarith)
        linarith [hc, habs]
      have hupper : ∀ k : ℕ, FR k = iR →
          (round x : ℝ) - x ≤ M + (v iR : ℝ) * (ε₀/(k+1)) := by
        intro k hk
        have hc := hcore k hk
        have h1 := (abs_le.mp hc).1
        linarith
      have hlim : (round x : ℝ) - x ≤ M := by
        apply le_of_forall_pos_le_add
        intro ε' hε'
        have hvv : (0:ℝ) < (v iR : ℝ) * ε₀ := by positivity
        obtain ⟨K, hK⟩ := exists_nat_gt ((v iR : ℝ) * ε₀ / ε')
        obtain ⟨k, hkmem, hkgt⟩ := hiR.exists_gt K
        have hkfib : FR k = iR := by simpa using hkmem
        have hup := hupper k hkfib
        have hKpos : (0:ℝ) < (K : ℝ) := by
          calc (0:ℝ) < (v iR : ℝ) * ε₀ / ε' := div_pos hvv hε'
            _ < K := hK
        have hkR : (K : ℝ) < (k : ℝ) := by exact_mod_cast hkgt
        have hfrac : (v iR : ℝ) * (ε₀/(k+1)) < ε' := by
          rw [mul_div_assoc']
          rw [div_lt_iff₀ (by positivity : (0:ℝ) < (k:ℝ)+1)]
          have h2 : (v iR : ℝ) * ε₀ < ε' * K := by
            rw [div_lt_iff₀ hε'] at hK
            linarith
          nlinarith
        linarith
      have hfinal : (round x : ℝ) - x = M := by
        have := hsign
        linarith [le_antisymm hlim (by linarith : M ≤ (round x : ℝ) - x)]
      exact hfinal
    -- THE LEFT-ASCENDER: x_L − round(x_L) = M
    have hLfact : (v iL : ℝ) * tstar - (round ((v iL : ℝ) * tstar) : ℝ) = M := by
      set x : ℝ := (v iL : ℝ) * tstar with hxdef
      have hvpos : (0:ℝ) < (v iL : ℝ) := by exact_mod_cast hv iL
      have hd0 : M ≤ |x - round x| := by
        have h1 : M ≤ distZ x := by
          rw [hM]
          unfold profile
          exact Finset.inf'_le _ (Finset.mem_univ iL)
        simpa [distZ] using h1
      have hcore : ∀ k : ℕ, FL k = iL →
          |x - (v iL : ℝ) * (ε₀/(k+1)) - round x| ≤ M := by
        intro k hk
        have h1 := hFL k
        rw [hk] at h1
        have hexp : (v iL : ℝ) * (tstar - ε₀/(k+1)) = x - (v iL : ℝ) * (ε₀/(k+1)) := by
          rw [hxdef]
          ring
        rw [hexp] at h1
        unfold distZ at h1
        have hεp := (hεk iL k).1
        have hεle := (hεk iL k).2
        have hxm : |x - (round (x - (v iL : ℝ) * (ε₀/(k+1))) : ℤ)| < 1/2 := by
          have htri : |x - (round (x - (v iL : ℝ) * (ε₀/(k+1))) : ℝ)|
              ≤ |x - (x - (v iL : ℝ) * (ε₀/(k+1)))|
                + |x - (v iL : ℝ) * (ε₀/(k+1))
                    - round (x - (v iL : ℝ) * (ε₀/(k+1)))| := abs_sub_le _ _ _
          have hxy : |x - (x - (v iL : ℝ) * (ε₀/(k+1)))| = (v iL : ℝ) * (ε₀/(k+1)) := by
            rw [show x - (x - (v iL : ℝ) * (ε₀/(k+1))) = (v iL : ℝ) * (ε₀/(k+1)) by ring,
              abs_of_pos hεp]
          have hb := (hVbound iL).1
          calc |x - (round (x - (v iL : ℝ) * (ε₀/(k+1))) : ℝ)|
              ≤ _ + _ := htri
            _ ≤ (v iL : ℝ) * (ε₀/(k+1)) + M := by rw [hxy]; linarith
            _ ≤ (v iL : ℝ) * ε₀ + M := by linarith
            _ ≤ (1/2 - M)/2 + M := by linarith
            _ < 1/2 := by linarith
        have hpin := round_eq_of_abs_lt_half hxm
        rw [hpin]
        exact h1
      obtain ⟨k₀, hk₀mem⟩ := hiL.nonempty
      have hk₀ : FL k₀ = iL := by simpa using hk₀mem
      have hsign : M ≤ x - (round x : ℝ) := by
        by_contra hcon
        push Not at hcon
        have hge : M ≤ (round x : ℝ) - x := by
          rcases abs_cases (x - (round x : ℝ)) with ⟨he, -⟩ | ⟨he, -⟩
          · linarith [hd0]
          · linarith [hd0]
        have hc := hcore k₀ hk₀
        have hεp := (hεk iL k₀).1
        have habs : |x - (v iL : ℝ) * (ε₀/(k₀+1)) - round x|
            = -(x - (v iL : ℝ) * (ε₀/(k₀+1)) - round x) := by
          apply abs_of_neg
          linarith
        linarith [hc, habs]
      have hupper : ∀ k : ℕ, FL k = iL →
          x - (round x : ℝ) ≤ M + (v iL : ℝ) * (ε₀/(k+1)) := by
        intro k hk
        have hc := hcore k hk
        have h1 := (abs_le.mp hc).2
        linarith
      have hlim : x - (round x : ℝ) ≤ M := by
        apply le_of_forall_pos_le_add
        intro ε' hε'
        have hvv : (0:ℝ) < (v iL : ℝ) * ε₀ := by positivity
        obtain ⟨K, hK⟩ := exists_nat_gt ((v iL : ℝ) * ε₀ / ε')
        obtain ⟨k, hkmem, hkgt⟩ := hiL.exists_gt K
        have hkfib : FL k = iL := by simpa using hkmem
        have hup := hupper k hkfib
        have hKpos : (0:ℝ) < (K : ℝ) := by
          calc (0:ℝ) < (v iL : ℝ) * ε₀ / ε' := div_pos hvv hε'
            _ < K := hK
        have hkR : (K : ℝ) < (k : ℝ) := by exact_mod_cast hkgt
        have hfrac : (v iL : ℝ) * (ε₀/(k+1)) < ε' := by
          rw [mul_div_assoc']
          rw [div_lt_iff₀ (by positivity : (0:ℝ) < (k:ℝ)+1)]
          have h2 : (v iL : ℝ) * ε₀ < ε' * K := by
            rw [div_lt_iff₀ hε'] at hK
            linarith
          nlinarith
        linarith
      exact le_antisymm hlim hsign
    refine ⟨tstar, iR, iL,
      round ((v iR : ℝ) * tstar) + round ((v iL : ℝ) * tstar), ?_, hmax⟩
    push_cast
    have hexp : ((v iR : ℝ) + (v iL : ℝ)) * tstar
        = (v iR : ℝ) * tstar + (v iL : ℝ) * tstar := by ring
    rw [hexp]
    linarith [hRfact, hLfact]

#print axioms profile_attains_max
#print axioms merge_grid_attainment

end MergeGrid
end LonelyRunner
