/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Arithmetic consumers for THM-1237 positioned seven-wall Hunter transfer

THM-1166 supplies the periodic one-interval discrepancies and THM-1221
supplies the global `15/154` spanning-tree credit.  This module checks the
forest-Hunter count and the exact arithmetic downstream of those providers.
-/

namespace LRC14
namespace PositionedSevenWallHunter

/-- A nonempty induced subgraph of a forest has enough vertex surplus for the
pointwise Hunter inequality. -/
theorem forest_active_set_hunter {vertices edges : ℕ}
    (hnonempty : 1 ≤ vertices) (hforest : edges ≤ vertices - 1) :
    (1 : ℤ) ≤ vertices - edges := by
  omega

/-- The exact pair-positioning error is maximized at pair mass `1/14` on the
proved range `0 ≤ rho ≤ 1/14`. -/
theorem pair_position_error_bound {rho : ℝ}
    (h0 : 0 ≤ rho) (h14 : rho ≤ 1 / 14) :
    rho * (1 - rho) ≤ 13 / 196 := by
  have hsecond : 0 ≤ 13 / 14 - rho := by linarith
  have hproduct := mul_nonneg (sub_nonneg.mpr h14) hsecond
  nlinarith

/-- Abstract algebraic endpoint of the local forest-Hunter proof. -/
theorem positioned_forest_transfer
    {forestWeight harmonic safeMass : ℝ}
    (hlocal : forestWeight - (6 / 49) * harmonic ≤ safeMass) :
    forestWeight - (6 / 49) * harmonic ≤ safeMass := by
  exact hlocal

/-- A covered protected needle must pay the weighted harmonic/gcd debt
`24mH+13mG ≥ 30/11`.  The hypotheses are the denominator-cleared needle
length and the THM-1221 positioned-tree cover inequality. -/
theorem protected_needle_debt
    {L m harmonic gcdDebt : ℝ}
    (hm : 0 < m) (hneedle : 1 ≤ 7 * m * L)
    (hcovered : (15 / 154) * L ≤
      (6 / 49) * harmonic + (13 / 196) * gcdDebt) :
    30 / 11 ≤ 24 * m * harmonic + 13 * m * gcdDebt := by
  have hmul := mul_le_mul_of_nonneg_left hcovered (le_of_lt hm)
  nlinarith

/-- Subtracting a periodic forbidden event from a positioned safe-mass lower
bound incurs only its local upper mass. -/
theorem periodic_forbidden_transfer
    {forestWeight harmonic safeMass forbiddenMass beta L periodError
      outsideMass : ℝ}
    (hsafe : forestWeight - (6 / 49) * harmonic ≤ safeMass)
    (hforbidden : forbiddenMass ≤ beta * L + periodError)
    (houtside : safeMass - forbiddenMass ≤ outsideMass) :
    forestWeight - (6 / 49) * harmonic - beta * L - periodError ≤
      outsideMass := by
  linarith

/-- Exact strict-spectrum/BAD incidence margins and periodic discrepancy
constants. -/
theorem bad_transfer_constants :
    (15 / 154 : ℚ) - 2 / 21 = 1 / 462 ∧
      (2 / 21 : ℚ) * (1 - 2 / 21) = 38 / 441 ∧
      (15 / 154 : ℚ) - 60 / 637 = 45 / 14014 ∧
      (60 / 637 : ℚ) * (1 - 60 / 637) = 34620 / 405769 := by
  norm_num

/-- A rank-`r` threshold forest with positioned profit at least `lambda` per
edge forces a safe point unless its credit is absorbed by harmonic debt. -/
theorem threshold_rank_covered
    {rank lambda harmonic safeMass : ℝ}
    (hsafe : rank * lambda - (6 / 49) * harmonic ≤ safeMass)
    (hcovered : safeMass ≤ 0) :
    rank * lambda ≤ (6 / 49) * harmonic := by
  linarith

/-- Constant ledger for the protected-needle rescaling. -/
theorem protected_needle_constants :
    (1 / 14 : ℚ) * (1 - 1 / 14) = 13 / 196 ∧
      (15 / 154 : ℚ) * (1 / 7) * 196 = 30 / 11 ∧
      (6 / 49 : ℚ) * 196 = 24 := by
  norm_num

#print axioms forest_active_set_hunter
#print axioms pair_position_error_bound
#print axioms positioned_forest_transfer
#print axioms protected_needle_debt
#print axioms periodic_forbidden_transfer
#print axioms bad_transfer_constants
#print axioms threshold_rank_covered
#print axioms protected_needle_constants

end PositionedSevenWallHunter
end LRC14
