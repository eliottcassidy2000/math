/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-02-S41)
-/
import TournamentH7.RatIntervals
import TournamentH7.RatIntervalsWrap

/-!
# THM-604(a) over `Region`: the origin-nest attainment on the full simplex census

THM-604 (origin-nest max law): for a primitive pattern `m₁ ≤ … ≤ m_d` with `Σ mᵢ ≤ 7`, the
phased d-fold comb overlap has `max_θ ov = 2r/m_d`, attained at the aligned (origin) phases.
This file machine-checks the ATTAINMENT half on the complete census — all 30 primitive
patterns of depth 2–5 — as one `native_decide` over exact `Region` arithmetic:

    length (⋂ᵢ wrap (comb mᵢ r 0)) = 2r / m_d    at r = 1/14.

The shift-uniform upper bound (`ov(θ⃗) ≤ 2r/m_d` for ALL phases — the coprime-pair count
argument of THM-604) is the remaining [PAPER→LEAN] half; with it, the THM-602(C) deviation
charges `max − mean = 2r/m_d − (2r)^d` become closed-form at every depth, no stored tables.
-/

namespace LonelyRunner
namespace RatIntervals

/-- The origin-aligned nest of a pattern: the wrapped intersection of all its combs at
phase `0`. -/
def originNestQ (m : List ℕ) (r : ℚ) : ℚ :=
  match m with
  | [] => 0
  | v :: rest =>
      length (rest.foldl (fun L u => inter L (wrap (comb u r 0))) (wrap (comb v r 0)))

/-- The complete THM-604 simplex census: primitive nondecreasing patterns, depth 2–5,
`Σ mᵢ ≤ 7` (30 patterns). -/
def censusPatterns : List (List ℕ) :=
  [[1,1],[1,2],[1,3],[1,4],[1,5],[1,6],[2,3],[2,5],[3,4],
   [1,1,1],[1,1,2],[1,1,3],[1,1,4],[1,1,5],[1,2,2],[1,2,3],[1,2,4],[1,3,3],[2,2,3],
   [1,1,1,1],[1,1,1,2],[1,1,1,3],[1,1,1,4],[1,1,2,2],[1,1,2,3],[1,2,2,2],
   [1,1,1,1,1],[1,1,1,1,2],[1,1,1,1,3],[1,1,1,2,2]]

/-- **THM-604 attainment, machine-checked on the full census**: the origin nest of every
census pattern has length exactly `2r/m_d` at `r = 1/14`. -/
theorem originNest_attained_census :
    censusPatterns.all
      (fun m => decide (originNestQ m (1/14) = 2 * (1/14) / ((m.getLast?.getD 1 : ℕ) : ℚ)))
      = true := by
  native_decide

/-- Quantified corollary: extraction from the census check. -/
theorem originNest_attained {m : List ℕ} (hm : m ∈ censusPatterns) :
    originNestQ m (1/14) = 2 * (1/14) / ((m.getLast?.getD 1 : ℕ) : ℚ) := by
  have h := List.all_eq_true.mp originNest_attained_census m hm
  exact of_decide_eq_true h

end RatIntervals
end LonelyRunner
