/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-05-S84)
-/
import TournamentH7.LRCLiftRowsL7

/-!
# The tower lift: kernel rows migrate up the 13-adic tower for free (HYP-4136)

The tower's one-way traffic (HYP-4126, verified numerically S83) at the GATE level:
`speedOK13 s num den → speedOK13 s (13·num) (13·den)`.  A strict witness `num/den` is the
same rational as `13·num/(13·den)`, and the strict-gate inequality survives the rescaling
exactly — so every kernel row at denominator `13^ℓ` is automatically a row at `13^{ℓ+1}`:
covers project DOWN the tower, witnesses lift UP.  Consequently the S77/S82 row corpus at
169 feeds any future level-3 (2197) table without recomputation, and the tower-limit
dichotomy (HYP-4126: only 13-adic dilations cover at every level) is the right asymptotic
form of the shadow question.

`strictLonely13_of_kernelWitness` consumes lifted rows unchanged (it is
denominator-generic), so no new consumer is needed.
-/

namespace LonelyRunner
namespace TowerLift

open KernelGate13

/-- **The gate-level tower lift**: the strict `1/13`-gate at witness `num/den` transports
verbatim to `13·num/(13·den)` — rows migrate up the 13-adic tower for free. -/
theorem speedOK13_lift {s num : ℤ} {den : ℕ} (h : speedOK13 s num den) :
    speedOK13 s (13 * num) (13 * den) := by
  unfold speedOK13 at h ⊢
  have hmod : (s * (13 * num)) % ((13 * den : ℕ) : ℤ) = 13 * ((s * num) % den) := by
    push_cast
    rw [show s * (13 * num) = 13 * (s * num) by ring]
    exact Int.mul_emod_mul_of_pos _ _ (by norm_num)
  rw [hmod]
  push_cast
  omega

/-- Demo: the S82 pattern-A row at 169 is a row at 2197 with no new kernel check. -/
theorem rowA_check_2197 : ∀ i, speedOK13 (LiftRowsL7.rowA i) (13 * 6) (13 * 169) :=
  fun i => speedOK13_lift (LiftRowsL7.rowA_check i)

#print axioms speedOK13_lift

end TowerLift
end LonelyRunner
