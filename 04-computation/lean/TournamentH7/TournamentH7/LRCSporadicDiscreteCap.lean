/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-17)
-/
import Mathlib

/-!
# A denominator-quantized cap for the twelve-speed sporadic branch

This file isolates the elementary arithmetic at the end of two canonical
reductions:

* THM-668, Part 2 (pair-sum ruler): for a finite positive core `P` with
  `b = max P`, its maximum is attained at `p/q`, where
  `q = u+v <= 2*b`; putting
  `m = min_{x in P} dist(x*p mod q,{0,q})` gives `mu=M(P)=m/q`.
  The displayed `q` is the pair-sum ruler itself and need not be the reduced
  denominator of `m/q`.
* THM-759 (tight-instance ratio bound): if adjoining a top speed `w` makes an
  `n`-speed tight family, then
  `w <= b/((n+1)*mu-1)`.  On the sporadic branch `mu>1/n`, the denominator is
  positive, so this is equivalently
  `w*((n+1)*mu-1) <= b`.

Thus, if an eleven-speed core has maximum `b` and rational maximum
`mu = m/q > 1/n`, THM-668 supplies the inclusive bound `q <= 2*b`.
Integrality then
forces the quantitative gap

`mu >= 1/n + 1/(2*n*b)`.

Combining that gap with the tight-completion ratio inequality

`w * ((n+1)*mu - 1) <= b`

gives

`w <= 2*n*b^2/(2*b+n+1)`.

At `n=12`, this is the inclusive rational bound
`w <= 24*b^2/(2*b+13)` (hence the corresponding integer floor), a strict
improvement over the coarse `w < 12*b`.  Strictness enters only in
`mu>1/n`, equivalently `q<n*m`; both the ruler bound and final cap are
inclusive.  The theorem is deliberately conditional on the rational ruler
data and the completion inequality: the THM-668/759 bridges themselves are
not formalized here, and this module does not assert the still-open global
emptiness of the sporadic branch.

The module is standalone and not added to the umbrella import.
-/

namespace LonelyRunner
namespace SporadicDiscreteCap

/-- A strict rational inequality above `1/n` is separated from it by at least
`1/(2*n*b)` once the positive denominator is at most `2*b`. -/
theorem rational_gap_of_denominator_le_two_mul
    (n b m q : ℕ) (hn : 0 < n) (hb : 0 < b) (hq : 0 < q)
    (hqle : q ≤ 2 * b) (hstrict : q < n * m) :
    (1 : ℚ) / n + 1 / (2 * n * b) ≤ (m : ℚ) / q := by
  have hstep : q + 1 ≤ n * m := by omega
  have hnq : (0 : ℚ) < (n : ℚ) * q := mul_pos (by positivity) (by positivity)
  have hnb : (0 : ℚ) < 2 * (n : ℚ) * b := by positivity
  have hqleQ : (q : ℚ) ≤ 2 * b := by exact_mod_cast hqle
  have hstepQ : (q : ℚ) + 1 ≤ n * m := by exact_mod_cast hstep
  have hfirst : (1 : ℚ) / n + 1 / (n * q) ≤ (m : ℚ) / q := by
    have hnum : (0 : ℚ) ≤ n * m - q - 1 := by nlinarith
    rw [← sub_nonneg]
    have hid : (m : ℚ) / q - (1 / n + 1 / (n * q)) =
        (n * m - q - 1) / (n * q) := by
      field_simp
      ring
    rw [hid]
    exact div_nonneg hnum (le_of_lt hnq)
  have hsecond : (1 : ℚ) / (2 * n * b) ≤ 1 / (n * q) := by
    apply one_div_le_one_div_of_le hnq
    nlinarith
  linarith

/-- Abstract completion step: the denominator gap and the tight-completion
ratio inequality imply the improved universal cap. -/
theorem completion_cap_of_gap
    (n b w mu : ℚ) (hn : 0 < n) (hb : 0 < b) (hw : 0 ≤ w)
    (hgap : 1 / n + 1 / (2 * n * b) ≤ mu)
    (hcompletion : w * ((n + 1) * mu - 1) ≤ b) :
    w ≤ 2 * n * b ^ 2 / (2 * b + n + 1) := by
  have hden : 0 < 2 * b + n + 1 := by positivity
  have hcoef : (2 * b + n + 1) / (2 * n * b) ≤ (n + 1) * mu - 1 := by
    have hn' : n ≠ 0 := ne_of_gt hn
    have hb' : b ≠ 0 := ne_of_gt hb
    calc
      (2 * b + n + 1) / (2 * n * b)
          = (n + 1) * (1 / n + 1 / (2 * n * b)) - 1 := by
              field_simp
              ring
      _ ≤ (n + 1) * mu - 1 := by nlinarith
  have hmul : w * ((2 * b + n + 1) / (2 * n * b)) ≤ b := by
    nlinarith [mul_le_mul_of_nonneg_left hcoef hw]
  apply (le_div_iff₀ hden).2
  have hn' : n ≠ 0 := ne_of_gt hn
  have hb' : b ≠ 0 := ne_of_gt hb
  field_simp at hmul ⊢
  nlinarith

/-- The complete arithmetic brick, from integral ruler data to the improved
cap. -/
theorem completion_cap_of_ruler
    (n b w m q : ℕ) (hn : 0 < n) (hb : 0 < b) (hq : 0 < q)
    (hqle : q ≤ 2 * b) (hstrict : q < n * m)
    (hcompletion :
      (w : ℚ) * (((n + 1 : ℕ) : ℚ) * ((m : ℚ) / q) - 1) ≤ b) :
    (w : ℚ) ≤ 2 * n * b ^ 2 / (2 * b + n + 1) := by
  apply completion_cap_of_gap (n := n) (b := b) (w := w)
      (mu := (m : ℚ) / q) (by exact_mod_cast hn) (by exact_mod_cast hb)
      (by positivity)
  · simpa using rational_gap_of_denominator_le_two_mul n b m q hn hb hq hqle hstrict
  · simpa using hcompletion

/-- Twelve-speed specialization: a nonextremal eleven-speed ruler core forces
`w <= 24*b^2/(2*b+13)`. -/
theorem twelve_speed_sporadic_cap
    (b w m q : ℕ) (hb : 0 < b) (hq : 0 < q)
    (hqle : q ≤ 2 * b) (hstrict : q < 12 * m)
    (hcompletion :
      (w : ℚ) * (13 * ((m : ℚ) / q) - 1) ≤ b) :
    (w : ℚ) ≤ 24 * b ^ 2 / (2 * b + 13) := by
  have h := completion_cap_of_ruler 12 b w m q (by norm_num) hb hq hqle hstrict
    (by norm_num at hcompletion ⊢; exact hcompletion)
  norm_num at h
  convert h using 1 <;> ring

#print axioms rational_gap_of_denominator_le_two_mul
#print axioms completion_cap_of_gap
#print axioms completion_cap_of_ruler
#print axioms twelve_speed_sporadic_cap

end SporadicDiscreteCap
end LonelyRunner
