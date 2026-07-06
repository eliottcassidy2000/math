/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S104)
-/
import Mathlib

/-!
# The divisor-protection lemma: the mechanism behind the Farey ladder (HYP-4366)

Generalizing the AP-11 protection (S103, `LRCAPProtection`) to every modulus `k`: a family
of integer speeds NONE of which is a multiple of `k` is lonely at `t = 1/k` with margin
`1/k`.  Since `1/k > 2/(2k+1)` (the mediant, the second Farey rung), such a family is LOOSE
— its `M` sits above the whole gap window `(1/(k+1), 2/(2k+1))`.

This is the structural mechanism behind the universal Farey ladder `j/(kj+1)` (S100,
HYP-4306; kps HYP-4357): it cleanly separates the spectrum.  A family with NO `k`-multiple
is loose.  The tight family `{1,…,k}` and EVERY window-candidate therefore MUST contain a
`k`-multiple (the "big" runner) — `{1,…,k}` contains `k`; the second-best `{1,…,k−1, 2k}`
contains `2k`; the deep well contains `(k−1)k`.  Loneliness at the tight value `1/(k+1)` is
exactly the price of admitting one `k`-multiple to sit at `0` when `t = 1/k`.
-/

namespace LonelyRunner
namespace DivisorProtection

/-- **Divisor distance.**  If `k ∤ s` (and `k ≥ 1`), the speed `s` is at circle-distance
`≥ 1/k` from `0` at time `1/k`: `dist(s/k, ℤ) ≥ 1/k`, because `|s − k·m| ≥ 1`. -/
theorem int_far_of_not_dvd_k (k : ℕ) (hk : 0 < k) (s : ℤ) (hs : ¬ (k : ℤ) ∣ s) :
    ∀ m : ℤ, (1 : ℝ) / k ≤ |(s : ℝ) / k - m| := by
  intro m
  have hkR : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk
  have hne : s - (k : ℤ) * m ≠ 0 := by
    intro h
    apply hs
    rw [sub_eq_zero.mp h]
    exact dvd_mul_right _ _
  have h1 : (1 : ℤ) ≤ |s - (k : ℤ) * m| := Int.one_le_abs hne
  -- s/k − m = (s − k·m)/k  (needs k ≠ 0)
  have hid : (s : ℝ) / k - m = ((s : ℝ) - (k : ℝ) * m) / k := by
    rw [sub_div, mul_div_cancel_left₀ _ hkR.ne']
  have h1R : (1 : ℝ) ≤ |(s : ℝ) - (k : ℝ) * (m : ℝ)| := by
    have e : (s : ℝ) - (k : ℝ) * (m : ℝ) = ((s - (k : ℤ) * m : ℤ) : ℝ) := by push_cast; ring
    rw [e, ← Int.cast_abs]; exact_mod_cast h1
  rw [hid, abs_div, abs_of_pos hkR, le_div_iff₀ hkR]
  have hkk : (1 : ℝ) / k * k = 1 := by field_simp
  linarith

/-- **Divisor protection.**  A family `V` of integer speeds with NO multiple of `k` is
lonely at `t = 1/k` with margin `1/k`: every `s ∈ V` stays `≥ 1/k` from every integer. -/
theorem lonely_of_all_not_dvd (k : ℕ) (hk : 0 < k) (V : Finset ℤ)
    (hV : ∀ s ∈ V, ¬ (k : ℤ) ∣ s) :
    ∀ s ∈ V, ∀ m : ℤ, (1 : ℝ) / k ≤ |(s : ℝ) * (1 / k) - m| := by
  intro s hs m
  have hval : (s : ℝ) * (1 / k) - m = (s : ℝ) / k - m := by ring
  rw [hval]
  exact int_far_of_not_dvd_k k hk s (hV s hs) m

/-- **The ladder separation** (mediant form).  A `k`-multiple-free family is loose above the
mediant: its margin `1/k` strictly exceeds `2/(2k+1)`.  Hence every family inside the gap
window `(1/(k+1), 2/(2k+1))` — and the tight family at `1/(k+1)` — contains a `k`-multiple. -/
theorem margin_gt_mediant (k : ℕ) (hk : 0 < k) :
    (2 : ℝ) / (2 * k + 1) < 1 / k := by
  have hkR : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk
  rw [div_lt_div_iff₀ (by positivity) hkR]
  push_cast; nlinarith

#print axioms int_far_of_not_dvd_k
#print axioms lonely_of_all_not_dvd
#print axioms margin_gt_mediant

end DivisorProtection
end LonelyRunner
