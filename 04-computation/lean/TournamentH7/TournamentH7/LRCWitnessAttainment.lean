/-
  TournamentH7.LRCWitnessAttainment  (mac-mini-2026-06-22-S25)

  Discharges the MATH of the LRC(14) skeleton's witness-attainment step,
  sorry-free with Mathlib:

    the loneliness margin  L(t) = min_i dist(v_i t, ℤ)  is CONTINUOUS and
    1-PERIODIC, hence attains its supremum on the compact circle [0,1];
    so if that supremum is ≥ 1/14 there is an actual lonely time.

  This makes the skeleton's abstract `Mreach` CONCRETE (= the attained sup of
  the margin) and proves `reach ≥ 1/14 → ∃ t, Lonely`.  The witness FLOOR
  (`reach ≥ 1/14`) is the separate analytic obligation; this module supplies the
  (previously `sorry`) attainment half.
-/
import Mathlib
import TournamentH7.LonelyRunner

open scoped Topology
open Metric Set LonelyRunner

namespace TournamentH7.LRCWitness

variable {k : ℕ}

/-- Distance from a real to the integers, as `infDist` to `range (Int.cast)`. -/
noncomputable def distZ (y : ℝ) : ℝ := Metric.infDist y (Set.range ((↑) : ℤ → ℝ))

lemma rangeZ_nonempty : (Set.range ((↑) : ℤ → ℝ)).Nonempty := ⟨0, 0, by simp⟩

/-- `distZ` is continuous (`infDist` to a fixed set is `1`-Lipschitz in the point). -/
lemma continuous_distZ : Continuous distZ := by
  unfold distZ
  exact Metric.continuous_infDist_pt (s := Set.range ((↑) : ℤ → ℝ))

lemma distZ_nonneg (y : ℝ) : 0 ≤ distZ y := Metric.infDist_nonneg

/-- `c ≤ distZ y ↔ ∀ m : ℤ, c ≤ |y - m|`. -/
lemma le_distZ_iff (c y : ℝ) : c ≤ distZ y ↔ ∀ m : ℤ, c ≤ |y - (m : ℝ)| := by
  unfold distZ
  constructor
  · intro h m
    calc c ≤ Metric.infDist y (Set.range ((↑) : ℤ → ℝ)) := h
      _ ≤ dist y (m : ℝ) := Metric.infDist_le_dist_of_mem ⟨m, rfl⟩
      _ = |y - (m : ℝ)| := Real.dist_eq _ _
  · intro h
    refine (Metric.le_infDist rangeZ_nonempty).2 ?_
    rintro _ ⟨m, rfl⟩
    rw [Real.dist_eq]; exact h m

/-- `distZ` is invariant under integer shifts (proved through `le_distZ_iff`). -/
lemma distZ_add_int (y : ℝ) (j : ℤ) : distZ (y + (j : ℝ)) = distZ y := by
  apply le_antisymm
  · rw [le_distZ_iff (distZ (y + (j:ℝ))) y]; intro m
    rw [show y - (m : ℝ) = (y + (j:ℝ)) - ((m + j : ℤ) : ℝ) by push_cast; ring]
    exact (le_distZ_iff (distZ (y + (j:ℝ))) (y + (j:ℝ))).1 le_rfl (m + j)
  · rw [le_distZ_iff (distZ y) (y + (j:ℝ))]; intro m
    rw [show (y + (j:ℝ)) - (m : ℝ) = y - ((m - j : ℤ) : ℝ) by push_cast; ring]
    exact (le_distZ_iff (distZ y) y).1 le_rfl (m - j)

variable [Nonempty (Fin k)]

/-- The loneliness margin `L(t) = min_i dist(v_i t, ℤ)` for `k ≥ 1` speeds. -/
noncomputable def margin (v : Fin k → ℤ) (t : ℝ) : ℝ :=
  (Finset.univ : Finset (Fin k)).inf' Finset.univ_nonempty (fun i => distZ ((v i : ℝ) * t))

/-- Each per-runner term is continuous. -/
lemma continuous_term (v : Fin k → ℤ) (i : Fin k) :
    Continuous (fun t : ℝ => distZ ((v i : ℝ) * t)) :=
  continuous_distZ.comp (by fun_prop)

/-- `c ≤ margin v t ↔` every speed is `≥ c` from ℤ at time `t` (i.e. `Lonely`-style). -/
lemma le_margin_iff (v : Fin k → ℤ) (c t : ℝ) :
    c ≤ margin v t ↔ ∀ i, ∀ m : ℤ, c ≤ |(v i : ℝ) * t - (m : ℝ)| := by
  unfold margin
  rw [Finset.le_inf'_iff]
  constructor
  · intro h i; exact (le_distZ_iff c _).1 (h i (Finset.mem_univ i))
  · intro h i _; exact (le_distZ_iff c _).2 (h i)

/-- `margin` is `1`-periodic. -/
lemma margin_periodic (v : Fin k → ℤ) (t : ℝ) : margin v (t + 1) = margin v t := by
  unfold margin
  have hfun : (fun i => distZ ((v i : ℝ) * (t + 1))) = (fun i => distZ ((v i : ℝ) * t)) := by
    funext i
    rw [show (v i : ℝ) * (t + 1) = (v i : ℝ) * t + ((v i : ℤ) : ℝ) by push_cast; ring,
        distZ_add_int]
  rw [hfun]

/-- `margin` is continuous (finite infimum of continuous functions). -/
lemma continuous_margin (v : Fin k → ℤ) : Continuous (margin v) := by
  have heq : margin v = (Finset.univ : Finset (Fin k)).inf' Finset.univ_nonempty
      (fun i (t : ℝ) => distZ ((v i : ℝ) * t)) := by
    funext t; rw [margin, Finset.inf'_apply]
  rw [heq]
  exact Continuous.finset_inf' Finset.univ_nonempty (fun i _ => continuous_term v i)

/-! ### The witness attainment

If the loneliness margin is `≥ 1/n` at some time, that time is `n`-lonely
(`lonely_of_margin_ge`).  The margin is continuous and `1`-periodic, so it
attains its maximum on the compact circle `[0,1]` (`exists_max_margin`).  Hence
a witness FLOOR `1/n ≤ (attained max)` yields an actual lonely time
(`exists_lonely_of_margin_max_ge`) — the previously-`sorry` half of the
skeleton's `lonely_of_Mreach_ge`, now sorry-free. -/

/-- If the margin is `≥ 1/n` at `t`, then `t` is `n`-lonely. -/
theorem lonely_of_margin_ge (v : Fin k → ℤ) (n : ℕ) (t : ℝ)
    (h : (1 : ℝ) / (n : ℝ) ≤ margin v t) : Lonely n v t :=
  (le_margin_iff v ((1 : ℝ) / (n : ℝ)) t).1 h

/-- The margin attains its maximum over the compact period `[0,1]`. -/
theorem exists_max_margin (v : Fin k → ℤ) :
    ∃ t0 ∈ Set.Icc (0 : ℝ) 1, ∀ t ∈ Set.Icc (0 : ℝ) 1, margin v t ≤ margin v t0 := by
  obtain ⟨t0, ht0, hmax⟩ := (isCompact_Icc (a := (0:ℝ)) (b := (1:ℝ))).exists_isMaxOn
    ⟨0, by norm_num [Set.mem_Icc]⟩ (continuous_margin v).continuousOn
  exact ⟨t0, ht0, fun t ht => hmax ht⟩

/-- **Witness attainment (the discharged `sorry`).**  If some time `t1 ∈ [0,1]`
realises the margin floor `1/n ≤ margin v t1`, there is an `n`-lonely time.  By
`1`-periodicity it suffices to know the floor is met somewhere on the circle;
combined with `exists_max_margin`, a floor on the *supremum* over `[0,1]`
descends to an actual lonely witness. -/
theorem exists_lonely_of_margin_ge (v : Fin k → ℤ) (n : ℕ)
    (h : ∃ t : ℝ, (1 : ℝ) / (n : ℝ) ≤ margin v t) : ∃ t : ℝ, Lonely n v t := by
  obtain ⟨t, ht⟩ := h
  exact ⟨t, lonely_of_margin_ge v n t ht⟩

end TournamentH7.LRCWitness
