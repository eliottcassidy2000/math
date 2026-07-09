/-
  TournamentH7.LRCAPTight — the tight (≤) half, completing `M(AP) = 1/14` exactly
  (kind-pasteur-2026-07-09-S111).

  kps-S110 proved `Mreach(AP) ≥ 1/14` (the loneliness the LRC asserts).  This file proves the matching
  `Mreach(AP) ≤ 1/14` via **Dirichlet's approximation theorem** with `n = 13`: for every time `τ` there
  is `k ∈ {1,…,13}` with `‖k·τ‖ ≤ 1/14` (14 points `{0,τ,…,13τ}` force two within `1/14`), so the AP
  runner `v = k` is within `1/14` of the origin — `minReach ≤ 1/14` at every `τ`, hence `Mreach ≤ 1/14`.
  Together with kps-S110: **`M(AP) = 1/14` exactly** — the AP is the LRC(14) extremal, achieving the
  bound with equality.  Uses Mathlib's `Real.exists_nat_abs_mul_sub_round_le`.
-/
import Mathlib
import TournamentH7.LRCAPExtremal
import TournamentH7.LRCWitnessAttainmentBridge

namespace LonelyRunner
namespace LRC14Concrete

/-- **`nearInt x ≤ |x − m|` for every integer `m`.**  The nearest-integer distance is `≤` the distance
to any particular integer (`nearInt = distZ = infDist` to `ℤ`). -/
theorem nearInt_le_abs_sub (x : ℝ) (m : ℤ) : nearInt x ≤ |x - (m : ℝ)| := by
  rw [← TournamentH7.LRCWitness.distZ_eq_nearInt, ← Real.dist_eq]
  exact Metric.infDist_le_dist_of_mem ⟨m, rfl⟩

/-- `minReach v t ≤ nearInt (v i · t)` for each runner `i` (the finite infimum is `≤` any term). -/
theorem minReach_le (v : Fin 13 → ℤ) (t : ℝ) (i : Fin 13) :
    minReach v t ≤ nearInt ((v i : ℝ) * t) := by
  unfold minReach
  exact ciInf_le ⟨0, by rintro _ ⟨j, rfl⟩; exact nearInt_nonneg _⟩ i

/-- The AP min-reach is `≤ 1/14` at every time (Dirichlet: some `k∈{1..13}` has `‖kτ‖ ≤ 1/14`). -/
theorem minReach_AP_le (τ : ℝ) : minReach (fun i : Fin 13 => (i.val : ℤ) + 1) τ ≤ 1 / 14 := by
  obtain ⟨k, hk0, hk13, hround⟩ := Real.exists_nat_abs_mul_sub_round_le τ (n := 13) (by norm_num)
  set i : Fin 13 := ⟨k - 1, by omega⟩ with hi
  have hvi : ((((fun j : Fin 13 => (j.val : ℤ) + 1) i) : ℤ) : ℝ) = (k : ℝ) := by
    have hk1 : 1 ≤ k := hk0
    simp only [hi]
    push_cast [Nat.cast_sub hk1]
    ring
  have hround' : |(k : ℝ) * τ - round ((k : ℝ) * τ)| ≤ 1 / 14 := by
    have h14 : (1 : ℝ) / ((13 : ℕ) + 1) = 1 / 14 := by norm_num
    rwa [h14] at hround
  calc minReach (fun j : Fin 13 => (j.val : ℤ) + 1) τ
      ≤ nearInt ((((fun j : Fin 13 => (j.val : ℤ) + 1) i : ℤ) : ℝ) * τ) := minReach_le _ _ i
    _ = nearInt ((k : ℝ) * τ) := by rw [hvi]
    _ ≤ |(k : ℝ) * τ - (round ((k : ℝ) * τ) : ℝ)| := nearInt_le_abs_sub _ _
    _ ≤ 1 / 14 := hround'

/-- **`Mreach(AP) ≤ 1/14`** — the tight half. -/
theorem mreach_AP_le : Mreach (fun i : Fin 13 => (i.val : ℤ) + 1) ≤ 1 / 14 := by
  rw [Mreach_eq_global_sSup]
  refine csSup_le ⟨minReach _ 0, 0, rfl⟩ ?_
  rintro _ ⟨τ, rfl⟩
  exact minReach_AP_le τ

/-- **`M(AP) = 1/14` exactly** — the AP is the LRC(14) extremal, achieving the bound with equality
(`≥` = kps-S110 loneliness, `≤` = Dirichlet tightness). -/
theorem mreach_AP_eq : Mreach (fun i : Fin 13 => (i.val : ℤ) + 1) = 1 / 14 :=
  le_antisymm mreach_AP_le mreach_AP_ge

end LRC14Concrete
end LonelyRunner
