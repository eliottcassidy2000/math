/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: kind-pasteur (LRC multi-agent project, 2026-07-10-S127)
-/
import Mathlib
import TournamentH7.LRCTwoDetunedClearing

/-!
# The mod-2g lift for the `(2,2)` half-harmonic residual — the terminating base case (kps-S127)

The last detuned residual (opus's THM-678 / my `ExceptionalDetunedDispatch`): a family
`v = g·H ∪ {δ₁, δ₂}` with `q₁ = q₂ = 2`, i.e. `δ₁ ≡ δ₂ ≡ g/2 (mod g)` — two **congruent half-harmonics**.
The branch count at `g` saturates (two `q=2` coordinates each cover `g/2` branches, filling `[0,g)`), so no
single branch clears both.

**monad's THM-678 mod-2g lift is a 2-adic DESCENT, not a closure.** At modulus `2g` the two detuned become
`q = 4` (`gcd(δᵢ, 2g) = g/2`), where the count would fire — BUT every ODD harmonic multiplier `m` becomes a
NEW half-harmonic of `2g` (`g·m ≡ g mod 2g`, i.e. `q = 2` at `2g`). So the obstruction descends one level:
lifting requires handling the odd-`m` harmonics, the same problem again. Full closure lands at the open core
(mac-mini's THM-676 descent-burden; klein's pair-sum / dead-zone route). **It is NOT independently
dischargeable from the LRC(≤13) citation** — I do not claim it is.

**What the lift DOES close — the terminating base case (verified `lrc14_two_detuned_lift_kps_S127`).** When
there are NO odd harmonic multipliers — every non-detuned coordinate divisible by `2g`, not merely `g` — the
descent halts in one step: at `2g` the ONLY non-multiples are the two detuned, both `q = 4`, and opus's
generic `d = 2` dispatch fires (`½·(¼+¼) < 1`). This is `lonely14_of_two_detuned_lift2g` below, kernel-pure,
reusing opus's `lonely14_of_two_detuned'` verbatim at the doubled modulus.  The proof needs only `¬ g ∣ δᵢ`
(a `q = 2` coordinate is automatically not `2g`-divisible and has `q(2g) ≠ 2`) plus the `2g`-divisibility of
the harmonic part — no gcd valuation arithmetic.

So the `(2,2)` residual splits cleanly: **base case (no odd multipliers) — CLOSED here; descent case (some
odd multiplier) — the open-core 2-adic tower, cited.**

Kernel-pure: no `sorry`, no `native_decide`.  Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LonelyRunner
namespace DetunedD2

open LonelyRunner

/-- **The mod-2g lift (terminating base case).**  A `d = 2` detuned family whose entire harmonic part is
divisible by `2g` (no odd multipliers) is lonely — apply opus's generic `d = 2` dispatch at the doubled
modulus `2g`, where both detuned speeds sit at `q = 4`.

This closes the `(2,2)` congruent-half-harmonic residual exactly when the descent terminates immediately
(no odd harmonic multiplier).  The proof only uses that each detuned speed is not divisible by `g`: then it
is not divisible by `2g` either, and `2g / gcd(δ, 2g) = 2` would force `g ∣ δ` — impossible.  So opus's
`(q₁,q₂) ≠ (2,2)` hypothesis holds at `2g` for free. -/
theorem lonely14_of_two_detuned_lift2g (cite : LRCUpTo13)
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0) (g : ℤ) (hg : 2 ≤ g)
    (i₁ i₂ : Fin 13) (hne : i₁ ≠ i₂)
    (hlift : ∀ j, j ≠ i₁ → j ≠ i₂ → 2 * g ∣ v j)
    (hδ1 : ¬ g ∣ v i₁) (hδ2 : ¬ g ∣ v i₂) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hg2 : 2 ≤ 2 * g := by omega
  have hgdvd : g ∣ 2 * g := ⟨2, by ring⟩
  have hδ1' : ¬ 2 * g ∣ v i₁ := fun h => hδ1 (dvd_trans hgdvd h)
  have hδ2' : ¬ 2 * g ∣ v i₂ := fun h => hδ2 (dvd_trans hgdvd h)
  -- a `q(2g) = 2` coordinate would have `gcd(δ, 2g) = g`, forcing `g ∣ δ`
  have hq_ne : ∀ δ : ℤ, ¬ g ∣ δ → 2 * g / (Int.gcd δ (2 * g) : ℤ) ≠ 2 := by
    intro δ hδ h1
    apply hδ
    have hdvd2g : ((Int.gcd δ (2 * g) : ℤ)) ∣ 2 * g := Int.gcd_dvd_right _ _
    have hcancel := Int.ediv_mul_cancel hdvd2g
    rw [h1] at hcancel                          -- 2 * ↑gcd = 2 * g
    have hgcd_eq : (Int.gcd δ (2 * g) : ℤ) = g := by omega
    have hdvdδ : ((Int.gcd δ (2 * g) : ℤ)) ∣ δ := Int.gcd_dvd_left _ _
    rwa [hgcd_eq] at hdvdδ
  have hq : ¬ (2 * g / (Int.gcd (v i₁) (2 * g) : ℤ) = 2 ∧ 2 * g / (Int.gcd (v i₂) (2 * g) : ℤ) = 2) :=
    fun h => hq_ne (v i₁) hδ1 h.1
  exact lonely14_of_two_detuned' cite v hv (2 * g) hg2 i₁ i₂ hne hlift hδ1' hδ2' hq

end DetunedD2
end LonelyRunner
