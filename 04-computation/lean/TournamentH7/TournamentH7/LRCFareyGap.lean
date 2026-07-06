/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S113)
-/
import Mathlib

/-!
# The loneliness gap is a Farey gap: witness denominators are large (HYP-4456)

The loneliness minimax gap window `(1/(k+1), 2/(2k+1))` has endpoints that are **Farey
neighbors**: `2·(k+1) − 1·(2k+1) = 1`.  Consequently no fraction of small denominator lies
strictly inside — the coarsest is the mediant `3/(3k+2)`.  So any minimax value `M = p/q`
landing in the gap has `q ≥ 3k+2`.

Combined with the witness-denominator lever (opus-S109, `LRCWitnessDenominator`: `M = c/q` in
lowest terms ⟹ `q ≤ 2·max|vᵢ|`), this is the **clearance-depth wall** (mac-mini S17): a gap
member has `max|vᵢ| ≥ (3k+2)/2`, a height *lower* bound.  The upper bound (single-cluster
height) is the remaining analytic piece; together they pin a gap member to a finite window —
which shrinks with `k`, the arithmetic face of the gap's n-specific emptiness.

This file proves the general Farey-gap denominator bound and specializes it to the window.
-/

namespace LonelyRunner
namespace FareyGap

/-- **Farey-gap denominator bound.**  If `a/b < p/q < c/d` with positive denominators and the
endpoints are Farey neighbors (`c·b − a·d = 1`), then `b + d ≤ q`: nothing of denominator below
the mediant denominator `b + d` lies strictly between two Farey neighbors.

Proof is the one-line integer identity `(pb−aq)·d + (cq−pd)·b = q·(cb−ad)`, with each bracket
`≥ 1` (strict integer inequalities) and `cb−ad = 1`. -/
theorem denom_ge_of_between {a p c b d q : ℤ} (hb : 0 < b) (hd : 0 < d)
    (h1 : a * q < p * b) (h2 : p * d < c * q) (hdet : c * b - a * d = 1) :
    b + d ≤ q := by
  have e1 : 1 ≤ p * b - a * q := by omega
  have e2 : 1 ≤ c * q - p * d := by omega
  have key : (p * b - a * q) * d + (c * q - p * d) * b = q * (c * b - a * d) := by ring
  rw [hdet, mul_one] at key
  have hle1 : d ≤ (p * b - a * q) * d := le_mul_of_one_le_left hd.le e1
  have hle2 : b ≤ (c * q - p * d) * b := le_mul_of_one_le_left hb.le e2
  linarith

/-- **The LRC gap window forces a large witness denominator.**  For the loneliness-gap window
`(1/(k+1), 2/(2k+1))` (Farey neighbors), any `M = p/q` strictly inside has `q ≥ 3k+2`.  With the
witness lever `q ≤ 2·max|vᵢ|` this yields the clearance-depth height lower bound
`max|vᵢ| ≥ (3k+2)/2`. -/
theorem gap_witness_denom_ge (k : ℕ) {p q : ℤ}
    (h1 : (1 : ℤ) * q < p * ((k : ℤ) + 1))
    (h2 : p * (2 * (k : ℤ) + 1) < 2 * q) :
    (3 * (k : ℤ) + 2) ≤ q := by
  have h := denom_ge_of_between (a := 1) (b := (k : ℤ) + 1) (c := 2) (d := 2 * (k : ℤ) + 1)
    (p := p) (q := q) (by positivity) (by positivity) h1 h2 (by ring)
  linarith

#print axioms denom_ge_of_between
#print axioms gap_witness_denom_ge

end FareyGap
end LonelyRunner
