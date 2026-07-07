/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S119)
-/
import Mathlib

/-!
# Binder infeasibility and the parity kill (HYP-4516)

The mediant-attainer trichotomy (mac-mini-S28, opus-S118): for the canonical family
`F(N) = {1,…,N-2, N, 3(N-1)}`, the far element `3(N-1)` binds the smallest feasible small speed
`b ∈ {2,3,5}` at the competing denominator `Q = 3(N-1) + b ∈ {3N-1, 3N, 3N+2}`, and `M(F(N)) =
3/Q`.  A speed `b` can bind (at residue `3`, i.e. distance `3/Q`) only if the congruence
`b·x ≡ 3 (mod Q)` is solvable; when it is not, that branch dies and a tighter one takes over.
Since `3/(3N-1) > 1/N > 3/(3N+2)`, the mediant `3/(3N+2)` (an in-gap value) wins exactly when the
looser branches `b = 2` and `b = 3` are both killed.

This file proves the arithmetic core:

* `no_solution_of_gcd_not_dvd`: `b·x ≡ r (mod Q)` is unsolvable whenever `gcd(b,Q) ∤ r`.
* `parity_kill`: the speed-2 branch `2·x ≡ 3 (mod Q)` is unsolvable for **even** `Q` (`2 ∤ 3`).

Consequence for LRC(14): our case is `N = 12`, so the speed-2 competitor sits at
`Q = 3N-1 = 35` — **odd** — so speed-2 *is* feasible and `M(F(12)) = 3/35 > 2/25` is loose; the
canonical mediant `3/38` is missed.  The kill fires instead at `N` odd (`Q = 3N-1` even), which is
why the mediant is attained at `N ≡ 1 (mod 6)` (and, keeping the `b=5` branch alive,
`N ≢ 1 (mod 5)`).  So the gap being missed at `N = 12` is a *parity* fact about `N`, not the
compositeness of `38 = 2·19`.
-/

namespace LonelyRunner
namespace BinderInfeasible

/-- **Binder infeasibility.**  If `gcd(b,Q) ∤ r`, the congruence `b·x ≡ r (mod Q)` has no
solution: no integer `x` gives `Q ∣ (b·x − r)`.  (If it did, `gcd(b,Q)` would divide both `b·x`
and `b·x − r`, hence `r`.) -/
theorem no_solution_of_gcd_not_dvd (Q b r : ℤ) (h : ¬ (Int.gcd b Q : ℤ) ∣ r) (x : ℤ) :
    ¬ (Q ∣ (b * x - r)) := by
  intro hd
  refine h ?_
  have hgb : (Int.gcd b Q : ℤ) ∣ b := Int.gcd_dvd_left b Q
  have hgQ : (Int.gcd b Q : ℤ) ∣ Q := Int.gcd_dvd_right b Q
  have hgbx : (Int.gcd b Q : ℤ) ∣ b * x := hgb.mul_right x
  have hgd : (Int.gcd b Q : ℤ) ∣ (b * x - r) := hgQ.trans hd
  have hr : (Int.gcd b Q : ℤ) ∣ (b * x - (b * x - r)) := dvd_sub hgbx hgd
  simpa using hr

/-- **The parity kill.**  For even `Q`, `2·x ≡ 3 (mod Q)` has no solution: `Q ∣ (2x − 3)` would
force `2x − 3` even, but it is odd.  This is the mechanism that kills the speed-2 branch of the
mediant trichotomy when the competing denominator `Q = 3N−1` is even (i.e. `N` odd). -/
theorem parity_kill (Q : ℤ) (hQ : Even Q) (x : ℤ) : ¬ (Q ∣ (2 * x - 3)) := by
  rintro ⟨k, hk⟩
  have h2 : Even (2 * x - 3) := by rw [hk]; exact hQ.mul_right k
  obtain ⟨m, hm⟩ := h2
  omega

#print axioms no_solution_of_gcd_not_dvd
#print axioms parity_kill

end BinderInfeasible
end LonelyRunner
