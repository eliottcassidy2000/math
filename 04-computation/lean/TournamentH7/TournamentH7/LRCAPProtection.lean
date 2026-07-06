/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S103)
-/
import Mathlib

/-!
# The AP-11 protection lemma: the inductive core of (C)/hdich (HYP-4356)

The crux (C)/hdich reduces (S100) to the AP-completion step: adding a 12th runner to the
dilated AP `{1,…,11}` gives either the tight completion (`v = 12`, `M = 1/13`) or a loose
family (`M ≥ 2/25`) — never an in-window value.  The clean mechanism (verified
`ap_completion_opus_S103.out`):

* At `t = 1/12`, every runner of `{1,…,11}` sits at circle-distance `≥ 1/12` (the AP is
  self-protecting), and **any** integer speed `v` with `12 ∤ v` ALSO sits at `≥ 1/12`
  (because `dist(v/12, ℤ) = |v − 12m|/12 ≥ 1/12` when `12 ∤ v`).
* So the family `{1,…,11, v}` is lonely at `1/12 > 2/25` whenever `12 ∤ v` — LOOSE.
* Only `12 ∣ v` can break the protection, i.e. `v = 12w`, which is exactly the l = 1 lift
  stratum (`LRCLiftRigidityRows`: tight at `w = 1`, `M ≥ 2/25` otherwise).

Hence the AP-completion step is: `{1,…,11, v}` is tight (`v = 12`) or loose — the induction's
"add-a-runner" core.
-/

namespace LonelyRunner
namespace APProtection

/-- **The protection core.**  If `12 ∤ s`, the speed `s` is at circle-distance `≥ 1/12`
from `0` at time `1/12`: `dist(s/12, ℤ) ≥ 1/12`. -/
theorem int_far_of_not_dvd (s : ℤ) (hs : ¬ (12 : ℤ) ∣ s) :
    ∀ m : ℤ, (1 : ℝ) / 12 ≤ |(s : ℝ) / 12 - m| := by
  intro m
  have hne : s - 12 * m ≠ 0 := by
    intro h
    exact hs ⟨m, by linarith [sub_eq_zero.mp h]⟩
  have h1 : (1 : ℤ) ≤ |s - 12 * m| := Int.one_le_abs (by exact_mod_cast hne)
  have h1R : (1 : ℝ) ≤ |((s - 12 * m : ℤ) : ℝ)| := by
    rw [← Int.cast_abs]; exact_mod_cast h1
  have hval : (s : ℝ) / 12 - m = ((s - 12 * m : ℤ) : ℝ) / 12 := by push_cast; ring
  rw [hval, abs_div, abs_of_pos (by norm_num : (0:ℝ) < 12),
    le_div_iff₀ (by norm_num : (0:ℝ) < 12)]
  linarith

/-- **AP-11 protection**: the family `{1,…,11, v}` with `12 ∤ v` is lonely at `t = 1/12`
with margin `1/12` — hence LOOSE (`M ≥ 1/12 > 2/25`), outside the gap window.  Every runner
avoids `0` because none is divisible by `12` (`1..11` trivially, `v` by hypothesis). -/
theorem ap11_loose_of_not_dvd (v : ℤ) (hv : ¬ (12 : ℤ) ∣ v) :
    ∃ t : ℝ, (∀ s : ℤ, s ∈ ({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, v} : Finset ℤ) →
      ∀ m : ℤ, (1 : ℝ) / 12 ≤ |(s : ℝ) * t - m|) := by
  refine ⟨1 / 12, fun s hs m => ?_⟩
  -- reduce `|s·(1/12) − m|` to `|s/12 − m|`
  have hval : (s : ℝ) * (1 / 12) - m = (s : ℝ) / 12 - m := by ring
  rw [hval]
  -- every listed speed is not divisible by 12
  have hnd : ¬ (12 : ℤ) ∣ s := by
    simp only [Finset.mem_insert, Finset.mem_singleton] at hs
    rcases hs with rfl|rfl|rfl|rfl|rfl|rfl|rfl|rfl|rfl|rfl|rfl|rfl
    all_goals first | decide | exact hv
  exact int_far_of_not_dvd s hnd m

#print axioms int_far_of_not_dvd
#print axioms ap11_loose_of_not_dvd

end APProtection
end LonelyRunner
