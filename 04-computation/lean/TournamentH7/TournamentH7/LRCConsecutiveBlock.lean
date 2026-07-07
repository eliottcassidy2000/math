/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S128)
-/
import TournamentH7.LRCWitnessAttainment

/-!
# Consecutive blocks are loose: the uniform-lift escape of the crux (HYP-4616)

The `(C)`-crux endgame (opus-S127) isolates the families that evade the finite mod-`q` covering to
the `lcm`-lifts of the AP `{i + L·k_i}`.  The **uniform** case (all `k_i` equal) is a *translate*
of the AP — a block of 12 consecutive integers `{m,…,m+11}` — and this file proves those are loose:

> for `m ≥ 2`, `M({m,…,m+11}) ≥ 2/25`, witnessed at `t = 1/(2m+11)`.

Indeed `M({m,…,m+11}) = m/(2m+11)` exactly (verified), with maximizer `t = 1/(2m+11)`: at that
time each speed `m+j` sits at `(m+j)/(2m+11)`, whose distance to `ℤ` is `min(m+j, m+11−j)/(2m+11) ≥
m/(2m+11)`, and `m/(2m+11) ≥ 2/25 ⟺ m ≥ 2`.  Only `m = 1` (the AP) dips below `2/25` (to `1/13`).
So the translate escape contributes only the AP; every other consecutive block is loose.
-/

namespace LonelyRunner
namespace ConsecutiveBlock

open TournamentH7.LRCWitness

/-- **Lattice-distance bound.**  For `0 ≤ j ≤ 11` and `m ≥ 2`, the point `m+j` is at distance `≥ m`
from every multiple of `q = 2m+11` (its nearest multiples `0` and `q` are `m+j ≥ m` and
`m+11−j ≥ m` away). -/
theorem block_dist_ge (m : ℤ) (hm : 2 ≤ m) (j : ℤ) (hj0 : 0 ≤ j) (hj : j ≤ 11) (n : ℤ) :
    m ≤ |m + j - n * (2 * m + 11)| := by
  by_cases hn : n ≤ 0
  · -- n ≤ 0 : m + j - n·q ≥ m + j ≥ m  (subtracting a nonpositive multiple only grows it)
    have hnn : 0 ≤ (-n) * (2 * m + 11) := mul_nonneg (by omega) (by omega)
    rw [abs_of_nonneg (by nlinarith)]
    nlinarith
  · -- n ≥ 1 : n·q − (m + j) ≥ q − (m+11) = m
    rw [not_le] at hn
    have hq : (2 * m + 11) ≤ n * (2 * m + 11) :=
      by nlinarith [mul_nonneg (show (0:ℤ) ≤ n - 1 by omega) (show (0:ℤ) ≤ 2 * m + 11 by omega)]
    rw [abs_of_nonpos (by nlinarith)]
    nlinarith

variable {k : ℕ}

/-- **The consecutive block `{m,…,m+11}` is loose for `m ≥ 2`.**  At `t = 1/(2m+11)` every speed is
`≥ 2/25` from `ℤ`, so the loneliness margin — and hence `M` — is `≥ 2/25`.  Thus no translate of the
AP other than the AP itself lies in the gap `(1/13, 2/25)`. -/
theorem block_margin_ge (m : ℤ) (hm : 2 ≤ m) :
    (2 : ℝ) / 25 ≤ margin (fun j : Fin 12 => m + (j : ℤ)) (1 / (2 * (m : ℝ) + 11)) := by
  have hqpos : (0 : ℝ) < 2 * (m : ℝ) + 11 := by positivity
  rw [le_margin_iff]
  intro i n
  -- |(m+i)/(2m+11) − n| = |m+i − n(2m+11)| / (2m+11) ≥ m/(2m+11) ≥ 2/25
  have hj0 : (0 : ℤ) ≤ (i : ℤ) := Int.natCast_nonneg i
  have hj : (i : ℤ) ≤ 11 := by have := i.2; omega
  have hlat : (m : ℤ) ≤ |m + (i : ℤ) - n * (2 * m + 11)| := block_dist_ge m hm i hj0 hj n
  -- cast the integer lattice bound to ℝ
  have hlatR : (m : ℝ) ≤ |(m : ℝ) + (i : ℝ) - (n : ℝ) * (2 * (m : ℝ) + 11)| := by
    exact_mod_cast hlat
  -- the target absolute value, cleared of the denominator
  have hcomb : (↑(m + (i : ℤ)) : ℝ) * (1 / (2 * (m : ℝ) + 11)) - (n : ℝ)
      = ((m : ℝ) + (i : ℝ) - (n : ℝ) * (2 * (m : ℝ) + 11)) / (2 * (m : ℝ) + 11) := by
    push_cast
    field_simp
  rw [hcomb, abs_div, abs_of_pos hqpos, le_div_iff₀ hqpos]
  -- 2/25 · (2m+11) ≤ |…|,  and |…| ≥ m ≥ 2/25·(2m+11) since 25m ≥ 4m+22 (m ≥ 2)
  have hmR : (2 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hm
  nlinarith [hlatR, hmR]

#print axioms block_dist_ge
#print axioms block_margin_ge

end ConsecutiveBlock
end LonelyRunner
