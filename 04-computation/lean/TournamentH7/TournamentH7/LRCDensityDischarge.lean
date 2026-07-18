/-
  TournamentH7.LRCDensityDischarge — the density-route discharge for a separated
  far element (the geometric core of boxeph-S96–S100).  boxeph-2026-07-18-S107.

  The density route closes a far-element family by keeping the good set NONEMPTY:
  a frame that is 1/13-lonely at `t₀` is 1/14-lonely on the whole interval
  `[t₀−δ, t₀+δ]`, `δ = 1/(182·V)` (`V` bounds the frame speeds), by the reverse
  triangle inequality (each frame runner moves at rate ≤ V, losing ≤ V·δ = 1/182,
  and 1/13 − 1/182 = 1/14).  If the far element `d = v vstar` is separated
  (`d ≥ 91·V`, i.e. the good interval has length `2δ ≥ 1/d`), that interval contains
  a half-integer point `t = (k+1/2)/d` where the far runner sits at distance 1/2 ≥ 1/14
  from every integer.  So the whole family is 1/14-lonely at `t`.

  `density_far_extension` (PROVED, kernel-pure): frame 1/13-lonely + far element
  `≥ 91·V` ⟹ `∃ t, Lonely 14 v t`.  This is the density route's Φ>0 mechanism
  (interval-completion), a proof DISTINCT from the descent floor's round+kick.
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRC13Citation

namespace LonelyRunner

variable {ι : Type*}

/-- A half-integer sits at distance ≥ 1/2 from every integer. -/
theorem half_integer_far (k m : ℤ) : (1 : ℝ) / 2 ≤ |(k : ℝ) + 1 / 2 - m| := by
  have hj : (k : ℝ) + 1 / 2 - m = ((k - m : ℤ) : ℝ) + 1 / 2 := by push_cast; ring
  rw [hj]
  rcases le_or_gt 0 (k - m) with h | h
  · have : (0 : ℝ) ≤ ((k - m : ℤ) : ℝ) := by exact_mod_cast h
    rw [abs_of_nonneg (by linarith)]; linarith
  · have : ((k - m : ℤ) : ℝ) ≤ -1 := by exact_mod_cast (by omega : k - m ≤ -1)
    rw [abs_of_neg (by linarith)]; linarith

/-- **The density-route discharge (separated far element) — PROVED.**
If the non-`vstar` speeds are `1/13`-lonely at `t₀` and bounded by `V`, and the far
element `v vstar ≥ 91·V`, then the whole family is `1/14`-lonely — at a point where the
far runner sits at distance `1/2` inside the frame's good interval `[t₀−δ, t₀+δ]`. -/
theorem density_far_extension
    (v : ι → ℤ) (vstar : ι) (t0 V : ℝ) (hV : 0 < V)
    (hframe : ∀ i, i ≠ vstar → ∀ m : ℤ, (1 : ℝ) / 13 ≤ |(v i : ℝ) * t0 - m|)
    (hbound : ∀ i, i ≠ vstar → |(v i : ℝ)| ≤ V)
    (hfar : (91 : ℝ) * V ≤ (v vstar : ℝ)) :
    ∃ t : ℝ, Lonely 14 v t := by
  set d : ℝ := (v vstar : ℝ) with hd
  have hdpos : (0 : ℝ) < d := lt_of_lt_of_le (by positivity) hfar
  set δ : ℝ := 1 / (182 * V) with hδ
  have hδpos : (0 : ℝ) < δ := by rw [hδ]; positivity
  -- an integer `k` with `d(t0−δ) − 1/2 ≤ k ≤ d(t0+δ) − 1/2`
  set k : ℤ := ⌈d * (t0 - δ) - 1 / 2⌉ with hk
  have hklo : d * (t0 - δ) - 1 / 2 ≤ (k : ℝ) := Int.le_ceil _
  have hkhi : (k : ℝ) ≤ d * (t0 + δ) - 1 / 2 := by
    have hlt : (k : ℝ) < d * (t0 - δ) - 1 / 2 + 1 := Int.ceil_lt_add_one _
    -- 2·d·δ ≥ 1 from d ≥ 91V and δ = 1/(182V)
    have h2δ : (1 : ℝ) ≤ d * (2 * δ) := by
      have heq : d * (2 * δ) = d / (91 * V) := by rw [hδ]; field_simp; ring
      rw [heq, le_div_iff₀ (by positivity)]
      nlinarith [hfar]
    nlinarith [hlt, h2δ]
  set t : ℝ := ((k : ℝ) + 1 / 2) / d with ht
  have hdt : d * t = (k : ℝ) + 1 / 2 := by rw [ht]; field_simp
  -- `t ∈ [t0−δ, t0+δ]`
  have hta : t0 - δ ≤ t := by
    rw [ht, le_div_iff₀ hdpos]; nlinarith [hklo]
  have htb : t ≤ t0 + δ := by
    rw [ht, div_le_iff₀ hdpos]; nlinarith [hkhi]
  have hdist : |t - t0| ≤ δ := by rw [abs_le]; constructor <;> linarith
  refine ⟨t, ?_⟩
  intro i m
  by_cases hi : i = vstar
  · -- the far runner sits at distance 1/2
    rw [hi]
    have : (v vstar : ℝ) * t - m = ((k : ℝ) + 1 / 2 - m) := by rw [← hd, hdt]
    rw [this]
    calc (1 : ℝ) / 14 ≤ 1 / 2 := by norm_num
      _ ≤ _ := half_integer_far k m
  · -- a frame runner stays lonely across the good interval
    have hsplit : (v i : ℝ) * t - m = ((v i : ℝ) * t0 - m) + (v i : ℝ) * (t - t0) := by ring
    have hkick : |(v i : ℝ) * (t - t0)| ≤ 1 / 182 := by
      rw [abs_mul]
      have hb := hbound i hi
      have : |(v i : ℝ)| * |t - t0| ≤ V * δ :=
        mul_le_mul hb hdist (abs_nonneg _) (le_of_lt hV)
      have hVδ : V * δ = 1 / 182 := by rw [hδ]; field_simp
      linarith [this, hVδ.le]
    have hrt := abs_sub_abs_le_abs_sub ((v i : ℝ) * t0 - m) (-((v i : ℝ) * (t - t0)))
    rw [abs_neg, sub_neg_eq_add, ← hsplit] at hrt
    have hf := hframe i hi m
    linarith [hf, hrt, hkick]

/-- **The density discharge as a complete LRC(14) rung.**  A 13-family whose 12 non-max
speeds are bounded by `V` and whose far element `v vstar ≥ 91·V` is `1/14`-lonely: the
LRC(≤13) citation makes the 12 speeds `1/13`-lonely at some `t₀`, and
`density_far_extension` completes the far runner inside their good interval.  This is the
density-route discharge for a separated far element, self-contained on the citation. -/
theorem density_far_bridge (cite : LRCUpTo13) (v : Fin 13 → ℤ)
    (hpos : ∀ i, 0 < v i) (vstar : Fin 13) (V : ℝ) (hV : 0 < V)
    (hbound : ∀ i, i ≠ vstar → |(v i : ℝ)| ≤ V)
    (hfar : (91 : ℝ) * V ≤ (v vstar : ℝ)) :
    ∃ t : ℝ, Lonely 14 v t := by
  set w : Fin 12 → ℤ := fun j => v (vstar.succAbove j) with hw
  have hwne : ∀ j, w j ≠ 0 := fun j => (hpos _).ne'
  obtain ⟨t0, ht0⟩ := cite 12 (le_refl 12) w hwne
  have hframe : ∀ i, i ≠ vstar → ∀ m : ℤ, (1 : ℝ) / 13 ≤ |(v i : ℝ) * t0 - m| := by
    intro i hi m
    obtain ⟨j, hj⟩ := Fin.exists_succAbove_eq hi
    have h := ht0 j m
    simp only [hw, hj] at h
    have hcast : (1 : ℝ) / ((12 + 1 : ℕ) : ℝ) = (1 : ℝ) / 13 := by norm_num
    rwa [hcast] at h
  exact density_far_extension v vstar t0 V hV hframe hbound hfar

end LonelyRunner

#print axioms LonelyRunner.half_integer_far
#print axioms LonelyRunner.density_far_extension
#print axioms LonelyRunner.density_far_bridge
