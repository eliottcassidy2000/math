/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-18)
-/
import TournamentH7.FragmentationLemma

/-!
# Multi-deletion harmonic-crown consumers

This module formalizes two downstream arithmetic steps from THM-1153.

The first specializes `LRC14.killer_budget` to the actual lonely-runner radius
`1 / 14`.  If at most six danger combs cover an interval of length `L`, then
any certified lower bound `ell0 ≤ L` pays the reciprocal endpoint budget

`ell0 * (7 - r) ≤ ∑ i, 1 / ws i`.

The second starts from the five cleared harmonic-crown inequalities for the
top seven normalized speeds.  The compact seed `x13 = 1`, `x12 < 13` then
implies the exact ceiling

`1350 * x7 < 613466231`.

The analytic production of the interval lower bound and of the five harmonic
inequalities remains outside this module; both are explicit hypotheses here.
-/

open scoped BigOperators

namespace LRC14
namespace MultiDeletionCrown

/-- The `1 / 14` multi-killer budget in the normalization used by THM-1153.

The hypothesis `r ≤ 6` is used only to make the coefficient `7 - r`
nonnegative, so that an interval lower bound `ell0 ≤ L` can be substituted
into the fragmentation inequality. -/
theorem harmonic_budget_fourteenth
    (r : ℕ) (hr : r ≤ 6) (ws : Fin r → ℕ) (hws : ∀ i, 1 ≤ ws i)
    (ell0 L x : ℝ) (hell0 : ell0 ≤ L) (hL : 0 ≤ L)
    (hcover : Set.Icc x (x + L) ⊆ ⋃ i, badArcs (ws i) (1 / 14)) :
    ell0 * (7 - (r : ℝ)) ≤ ∑ i, (1 : ℝ) / ws i := by
  have hr_real : (r : ℝ) ≤ 6 := by
    exact_mod_cast hr
  have hcoefficient : 0 ≤ (7 : ℝ) - r := by
    linarith
  have hinterval : ell0 * (7 - (r : ℝ)) ≤ L * (7 - (r : ℝ)) :=
    mul_le_mul_of_nonneg_right hell0 hcoefficient
  have hkb :=
    killer_budget r ws hws (1 / 14) L x (by norm_num) hL hcover
  have hscaled :=
    mul_le_mul_of_nonneg_left hkb (show (0 : ℝ) ≤ 7 by norm_num)
  have hbudget : L * (7 - (r : ℝ)) ≤ ∑ i, (1 : ℝ) / ws i := by
    convert hscaled using 1 <;> ring
  exact hinterval.trans hbudget

/-- The exact cleared triangular recurrence for the compact top-seven crown.

The five middle hypotheses are the denominator-cleared forms of

* `(5 / 42) * x11 ≤ x12 + x13`,
* `(12 / 77) * x10 ≤ x11 + x12 + x13`,
* `(6 / 35) * x9 ≤ x10 + x11 + x12 + x13`,
* `(10 / 63) * x8 ≤ x9 + x10 + x11 + x12 + x13`, and
* `(3 / 28) * x7 ≤ x8 + x9 + x10 + x11 + x12 + x13`.

No positivity assumptions are needed for this linear-arithmetic consumer. -/
theorem compact_top_seven_ceiling
    (x7 x8 x9 x10 x11 x12 x13 : ℝ)
    (hx13 : x13 = 1) (hx12 : x12 < 13)
    (h11 : 5 * x11 ≤ 42 * (x12 + x13))
    (h10 : 12 * x10 ≤ 77 * (x11 + x12 + x13))
    (h9 : 6 * x9 ≤ 35 * (x10 + x11 + x12 + x13))
    (h8 : 10 * x8 ≤ 63 * (x9 + x10 + x11 + x12 + x13))
    (h7 : 3 * x7 ≤ 28 * (x8 + x9 + x10 + x11 + x12 + x13)) :
    1350 * x7 < 613466231 := by
  linarith

end MultiDeletionCrown
end LRC14

#print axioms LRC14.MultiDeletionCrown.harmonic_budget_fourteenth
#print axioms LRC14.MultiDeletionCrown.compact_top_seven_ceiling
