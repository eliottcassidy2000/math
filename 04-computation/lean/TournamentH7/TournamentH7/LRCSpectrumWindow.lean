/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S116)
-/
import Mathlib

/-!
# The gap window in the Lonely Runner spectrum framework (HYP-4486)

Our gap-emptiness obligation (G) is exactly the **first-gap case of the Lonely Runner Spectrum
Conjecture**.  Kravitz conjectured that the achievable maximum-loneliness `ML(v)` below `1/n` is
always a *rung* `s/(ns+1)`; Fan–Sun disproved this (counterexamples at `n = 4, 6, 7`, all
generalized arithmetic progressions) and **amended** it to: `ML(v) = s/(ns+k)` with `k ≤ n`, or
`ML(v) ≥ 1/n`.  For our `n`-speed problem the extremal rung is `1/(n+1)` (`s=1`) and the next is
`2/(2n+1)` (`s=2`); the first gap window is `(1/(n+1), 2/(2n+1))`.

This file pins down which amended-spectrum values `s/(ns+k)` can land in that window:

  **`s/(ns+k) ∈ (1/(n+1), 2/(2n+1))  ⟺  k < s < 2k`.**

Consequences: the Kravitz rungs (`k = 1`) are *never* strictly inside (no integer `s` with
`1 < s < 2`), so a gap member must have order `k ≥ 2`; and the minimal admissible form is
`k = 2, s = 3`, the mediant `3/(3n+2)` (matching the Farey clearance `q ≥ 3n+2` of
`LRCFareyGap`).  So (G) sharpens to: *no `n`-speed family attains an `ML` of the form `s/(ns+k)`
with `2 ≤ k` and `k < s < 2k`* — and by the Fan–Sun evidence the achievable order `k` is small,
so this is a finite family of exceptional generalized-AP shapes to rule out.
-/

namespace LonelyRunner
namespace SpectrumWindow

/-- **Which amended-spectrum forms hit the first-gap window.**  For positive `n, s, k`, the value
`s/(ns+k)` lies strictly in `(1/(n+1), 2/(2n+1))` iff `k < s < 2k`. -/
theorem form_in_window_iff (n s k : ℤ) (hn : 0 < n) (hs : 0 < s) (hk : 0 < k) :
    ((1 : ℚ) / ((n : ℚ) + 1) < (s : ℚ) / ((n : ℚ) * s + k) ∧
      (s : ℚ) / ((n : ℚ) * s + k) < 2 / (2 * (n : ℚ) + 1)) ↔ (k < s ∧ s < 2 * k) := by
  have hnq : (0 : ℚ) < (n : ℚ) := by exact_mod_cast hn
  have hsq : (0 : ℚ) < (s : ℚ) := by exact_mod_cast hs
  have hkq : (0 : ℚ) < (k : ℚ) := by exact_mod_cast hk
  have hd : (0 : ℚ) < (n : ℚ) * s + k := by positivity
  have hn1 : (0 : ℚ) < (n : ℚ) + 1 := by linarith
  have hn2 : (0 : ℚ) < 2 * (n : ℚ) + 1 := by linarith
  rw [div_lt_div_iff₀ hn1 hd, div_lt_div_iff₀ hd hn2]
  constructor
  · rintro ⟨h1, h2⟩
    refine ⟨?_, ?_⟩
    · have : (k : ℚ) < s := by nlinarith
      exact_mod_cast this
    · have : (s : ℚ) < 2 * k := by nlinarith
      exact_mod_cast this
  · rintro ⟨h1, h2⟩
    have h1q : (k : ℚ) < s := by exact_mod_cast h1
    have h2q : (s : ℚ) < 2 * k := by exact_mod_cast h2
    exact ⟨by nlinarith, by nlinarith⟩

/-- **Kravitz rungs are never strictly inside the window.**  No rung `s/(ns+1)` (`k = 1`) lies in
`(1/(n+1), 2/(2n+1))`, since that would need `1 < s < 2`.  A first-gap member must have order
`k ≥ 2`. -/
theorem rung_not_in_window (n s : ℤ) (hn : 0 < n) (hs : 0 < s) :
    ¬ ((1 : ℚ) / ((n : ℚ) + 1) < (s : ℚ) / ((n : ℚ) * s + 1) ∧
      (s : ℚ) / ((n : ℚ) * s + 1) < 2 / (2 * (n : ℚ) + 1)) := by
  intro h
  have h' := (form_in_window_iff n s 1 hn hs one_pos).mp (by push_cast at h ⊢; exact h)
  omega

/-- **Numerator–denominator locking (the complexity is one parameter).**  For an in-window
spectrum value `M = s/(Ns+k)` (so `k < s`, from `form_in_window_iff`), the denominator `q = Ns+k`
is locked to the numerator by `Ns < q < (N+1)s` — i.e. `N < q/s < N+1`.  So the numerator `s`
(mac-mini's `c`), the denominator `q` (the witness denominator, opus-S109 `q ≤ 2·max`), and the
order `k` all grow together: bounding any one bounds the rest, and via the far-element structure
(mac-mini-S25) bounds the height.  Order `k`, numerator `s`, and Stern-Brocot depth are one
complexity parameter, and the height upper bound is equivalent to bounding it. -/
theorem window_num_denom_locked (N s k : ℤ) (hk : 0 < k) (hks : k < s) :
    N * s < N * s + k ∧ N * s + k < (N + 1) * s :=
  ⟨by linarith, by nlinarith⟩

#print axioms form_in_window_iff
#print axioms rung_not_in_window
#print axioms window_num_denom_locked

end SpectrumWindow
end LonelyRunner
