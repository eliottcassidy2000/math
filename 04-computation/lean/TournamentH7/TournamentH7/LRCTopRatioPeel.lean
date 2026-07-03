/-
  TournamentH7.LRCTopRatioPeel — THE 91-RATIO PEEL (kind-pasteur-2026-07-02-S18, HYP-3975).

  The intermediate band's uniform theorem.  The co-star (12,1) data route was sized and
  rejected this session (646,646 base subsets × finite windows ≈ 2·10⁶ instance rows —
  see lrc14_costar_r1_generator_kps_S18.py and its probe): per-base data cannot scale.
  Instead ONE theorem with NO census data closes every TOP-HEAVY family:

      |v i₀| ≥ 91 · max_{i ≠ i₀} |v i|   ⟹   v is 14-lonely,

  for ANY nonzero 13-family — no covering, no distinctness, no boundedness hypotheses.

  MECHANISM (three layers):
   * the LRC(≤13) citation node on the other TWELVE runners gives a point t₀ that is
     1/13-lonely for them;
   * the 1/13 − 1/14 margin transports along speeds ≤ B to the whole interval
     [t₀ − δ, t₀ + δ], δ = 1/(182·B) (Lipschitz) — consumed EXACTLY:
     1/13 − 1/182 = 1/14;
   * the far runner's dial covers a unit interval over the window (w·2δ ≥ 1), so a
     ruler point τ = (j + 1/2)/w lies inside — where the far runner sits at distance
     exactly 1/2 from ℤ.

  CONSEQUENCE: the families remaining after `lrc14_of_peel20` and this peel are the
  COMPRESSED ones — top ratio below 91 all the way down while exceeding the band.
-/
import Mathlib
import TournamentH7.LRC13Citation
import TournamentH7.LRC14ConcreteSurface

namespace LonelyRunner
namespace LRC14

/-- A half-integer is at distance at least `1/2` from every integer. -/
theorem half_int_dist (k : ℤ) : (1 : ℝ) / 2 ≤ |(k : ℝ) + 1 / 2| := by
  rcases le_or_gt 0 k with hk | hk
  · have : (0 : ℝ) ≤ (k : ℝ) := by exact_mod_cast hk
    rw [abs_of_nonneg (by linarith)]
    linarith
  · have hk1 : k ≤ -1 := by omega
    have : (k : ℝ) ≤ -1 := by exact_mod_cast hk1
    rw [abs_of_nonpos (by linarith)]
    linarith

/-- **THE 91-RATIO PEEL.**  If one runner is at least `91×` faster than every other,
the family is lonely — from the LRC(≤13) citation node alone.  The margin arithmetic
is exactly tight: `1/13 − 1/182 = 1/14`. -/
theorem top_ratio_lonely (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (i0 : Fin 13) {B : ℤ} (hB : 0 < B) (hrest : ∀ i, i ≠ i0 → |v i| ≤ B)
    (hfar : 91 * B ≤ |v i0|) : ∃ t : ℝ, Lonely 14 v t := by
  -- the twelve other runners, cited at gap 1/13
  set w : Fin 12 → ℤ := fun j => v (i0.succAbove j) with hw
  have hwne : ∀ j, w j ≠ 0 := fun j => hv _
  obtain ⟨t₀, h13⟩ := cite 12 (by omega) w hwne
  have hBR : (0 : ℝ) < (B : ℝ) := by exact_mod_cast hB
  set δ : ℝ := 1 / (182 * (B : ℝ)) with hδ
  have hδpos : 0 < δ := by positivity
  set W : ℝ := |(v i0 : ℝ)| with hW
  have hWpos : (0 : ℝ) < W := by
    rw [hW]
    exact abs_pos.mpr (by exact_mod_cast hv i0)
  have hWfar : 91 * (B : ℝ) ≤ W := by
    rw [hW, ← Int.cast_abs]
    exact_mod_cast hfar
  have hlen : 1 ≤ W * (2 * δ) := by
    rw [hδ]
    have h1 : W * (2 * (1 / (182 * (B : ℝ)))) = W / (91 * (B : ℝ)) := by
      field_simp
      ring
    rw [h1, le_div_iff₀ (by positivity)]
    linarith
  -- the ruler index
  obtain ⟨j, hj1, hj2⟩ : ∃ j : ℤ, W * (t₀ - δ) - 1 / 2 ≤ (j : ℝ) ∧
      (j : ℝ) ≤ W * (t₀ + δ) - 1 / 2 := by
    refine ⟨⌈W * (t₀ - δ) - 1 / 2⌉, Int.le_ceil _, ?_⟩
    have h1 := Int.ceil_lt_add_one (W * (t₀ - δ) - 1 / 2)
    have hspread : W * (t₀ + δ) - W * (t₀ - δ) = W * (2 * δ) := by ring
    linarith
  set τ : ℝ := ((j : ℝ) + 1 / 2) / W with hτ
  have hτlo : t₀ - δ ≤ τ := by
    rw [hτ, le_div_iff₀ hWpos]
    linarith
  have hτhi : τ ≤ t₀ + δ := by
    rw [hτ, div_le_iff₀ hWpos]
    linarith
  refine ⟨τ, ?_⟩
  intro i m
  rcases eq_or_ne i i0 with rfl | hne
  · -- the far runner sits at a half-integer
    have hWne : W ≠ 0 := hWpos.ne'
    rcases abs_cases ((v i : ℝ)) with ⟨habs, _⟩ | ⟨habs, _⟩
    · -- v i = W
      have hval : (v i : ℝ) * τ = (j : ℝ) + 1 / 2 := by
        rw [hτ, ← hW] at *
        rw [show (v i : ℝ) = W by rw [hW, ← habs]]
        field_simp
      rw [hval]
      have := half_int_dist (j - m)
      calc (1 : ℝ) / (14 : ℕ) ≤ 1 / 2 := by norm_num
        _ ≤ |((j - m : ℤ) : ℝ) + 1 / 2| := half_int_dist (j - m)
        _ = |(j : ℝ) + 1 / 2 - m| := by
            congr 1
            push_cast
            ring
    · -- v i = −W
      have hval : (v i : ℝ) * τ = -((j : ℝ) + 1 / 2) := by
        rw [show (v i : ℝ) = -W by rw [hW]; linarith [habs]]
        rw [hτ]
        field_simp
      rw [hval]
      calc (1 : ℝ) / (14 : ℕ) ≤ 1 / 2 := by norm_num
        _ ≤ |((-j - 1 - m : ℤ) : ℝ) + 1 / 2| := half_int_dist _
        _ = |-((j : ℝ) + 1 / 2) - m| := by
            congr 1
            push_cast
            ring
  · -- one of the twelve: Lipschitz transport of the 1/13 margin
    obtain ⟨j', hj'⟩ := Fin.exists_succAbove_eq hne
    have h0 : (1 : ℝ) / 13 ≤ |(v i : ℝ) * t₀ - m| := by
      have := h13 j' m
      rw [hw] at this
      simpa [hj'] using this
    have hvB : |(v i : ℝ)| ≤ (B : ℝ) := by
      rw [← Int.cast_abs]
      exact_mod_cast hrest i hne
    have htri : |(v i : ℝ) * t₀ - m| ≤ |(v i : ℝ) * τ - m| + |(v i : ℝ)| * |t₀ - τ| := by
      calc |(v i : ℝ) * t₀ - m|
          = |((v i : ℝ) * τ - m) + (v i : ℝ) * (t₀ - τ)| := by congr 1; ring
        _ ≤ |(v i : ℝ) * τ - m| + |(v i : ℝ) * (t₀ - τ)| := abs_add_le _ _
        _ = |(v i : ℝ) * τ - m| + |(v i : ℝ)| * |t₀ - τ| := by rw [abs_mul]
    have hwin : |t₀ - τ| ≤ δ := by
      rw [abs_le]
      constructor <;> linarith
    have hlip : |(v i : ℝ)| * |t₀ - τ| ≤ (B : ℝ) * δ := by
      apply mul_le_mul hvB hwin (abs_nonneg _) hBR.le
    have hBδ : (B : ℝ) * δ = 1 / 182 := by
      rw [hδ]
      field_simp
    calc (1 : ℝ) / (14 : ℕ) = 1 / 13 - 1 / 182 := by norm_num
      _ ≤ |(v i : ℝ) * τ - m| := by
          rw [hBδ] at hlip
          linarith

/-- Convenience form: peelability with the max-entry witness (the bound is the
floor-division `|v i₀| / 91`). -/
theorem top_ratio_lonely' (cite : LRCUpTo13) (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (i0 : Fin 13) (hpeel : ∀ i, i ≠ i0 → 91 * |v i| ≤ |v i0|) :
    ∃ t : ℝ, Lonely 14 v t := by
  set B : ℤ := |v i0| / 91 with hB
  have hBpos : 0 < B := by
    -- some other index exists; its entry is nonzero, so |v i0| ≥ 91
    have h1 : (0 : Fin 13) ≠ i0 ∨ (1 : Fin 13) ≠ i0 := by
      by_cases h : (0 : Fin 13) = i0
      · right
        rw [← h]
        decide
      · left
        exact h
    have hex : ∃ i, i ≠ i0 := h1.elim (fun h => ⟨0, h⟩) (fun h => ⟨1, h⟩)
    obtain ⟨i, hi⟩ := hex
    have h91 := hpeel i hi
    have habs : 1 ≤ |v i| := Int.one_le_abs (hv i)
    rw [hB]
    omega
  apply top_ratio_lonely cite v hv i0 hBpos
  · intro i hi
    have := hpeel i hi
    rw [hB]
    omega
  · rw [hB]
    omega

/-- **THE REMAINING SURFACE, exactly**: LRC(14) from the citation node, the `W = 22`
band data, and loneliness of the COMPRESSED covering families — covering, some entry
beyond the band, and no entry dominating all others by the factor 91.  Everything
else is machine-checked. -/
theorem lrc14_of_compressed_lonely (cite : LRCUpTo13)
    (hwindow : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → (∀ i, |v i| ≤ 22) →
      ∃ t : ℝ, Lonely 14 v t)
    (hcompressed : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → CoveringFamily v →
      ¬ (∀ i, |v i| ≤ 22) →
      (∀ i0 : Fin 13, ∃ i, i ≠ i0 ∧ |v i0| < 91 * |v i|) →
      ∃ t : ℝ, Lonely 14 v t) :
    LRC14Statement := by
  apply lrc14_of_covering_lonely
  intro v hv hcov
  by_cases hbnd : ∀ i, |v i| ≤ 22
  · exact hwindow v hv hbnd
  · by_cases hpeel : ∃ i0 : Fin 13, ∀ i, i ≠ i0 → 91 * |v i| ≤ |v i0|
    · obtain ⟨i0, hp⟩ := hpeel
      exact top_ratio_lonely' cite v hv i0 hp
    · push Not at hpeel
      exact hcompressed v hv hcov hbnd (fun i0 => by
        obtain ⟨i, hi, hlt⟩ := hpeel i0
        exact ⟨i, hi, by omega⟩)

end LRC14
end LonelyRunner
