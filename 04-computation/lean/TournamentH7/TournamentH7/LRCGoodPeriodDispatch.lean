/-
  TournamentH7.LRCGoodPeriodDispatch — assembling the good-period leg (kind-pasteur-2026-07-09-S99).

  The good-period existence for the covering case of LRC(14) is a DICHOTOMY on the longest AP in the
  co-offset set `E` (LEM-012 ∪ LEM-013):

    * near-AP branch  `longest-AP ≥ k−6`  — klein's LEM-012 (Dirichlet clustering).
    * dissociated     `longest-AP ≤ k−7`  — kps-S94 LEM-013 (the maxgap margin; klein-S200's
                                            `worst7Struct_hasGoodPeriod` is the finite-check template).

  These two ranges TILE all `L` with no gap (`k−7` and `k−6` are consecutive), so a good period exists
  for every cluster.  This file formalizes (i) that dispatch over klein's concrete decidable
  `HasGoodPeriod` predicate, and (ii) the geometric BRIDGE core `IsGoodPeriod ⟹ clearance > 1/14`
  that klein's file asserts ("THM-527: this forces `M(S) ≥ 1/14`") but leaves to cite: a circular
  gap of normalized length `> 1/7` has a midpoint at distance `> 1/14` from every phase, so the
  observer runner `Vmax` sits `> 1/14` from the whole cluster.  Builds on `LRCGoodPeriodMaxgap`.
-/
import Mathlib
import TournamentH7.LRCGoodPeriodMaxgap

namespace LRC14

/-- **The good-period dichotomy dispatch.**  If a cluster's longest AP has length `L`, and BOTH
branches are supplied — near-AP (`k−6 ≤ L`, LEM-012) and dissociated (`L ≤ k−7`, LEM-013) each
yielding a good period — then a good period exists unconditionally, because the two ranges tile all
`L` (for `k ≥ 7`, every `L` satisfies `k−6 ≤ L` or `L ≤ k−7`).  This is the exact surface the two
lemmas discharge. -/
theorem hasGoodPeriod_of_dichotomy (E : List ℕ) (Vmax k L : ℕ) (hk : 7 ≤ k)
    (hNearAP : k - 6 ≤ L → HasGoodPeriod E Vmax)
    (hDissoc : L ≤ k - 7 → HasGoodPeriod E Vmax) :
    HasGoodPeriod E Vmax := by
  by_cases h : k - 6 ≤ L
  · exact hNearAP h
  · exact hDissoc (by omega)

/-- **The bridge core (clearance).**  A good period (`IsGoodPeriod`, i.e. `7 · maxCircGap > Vmax`)
gives a circular gap whose HALF — the observer's clearance from the nearest cluster phase — exceeds
`1/14` in normalized units: `(maxCircGap)/(2·Vmax) > 1/14`.  This is the arithmetic heart of
"good period ⟹ `M(S) ≥ 1/14`". -/
theorem isGoodPeriod_clearance (E : List ℕ) (Vmax j : ℕ) (hV : 0 < Vmax)
    (h : IsGoodPeriod E Vmax j) :
    (1 : ℚ) / 14 < (maxCircGap E Vmax j : ℚ) / (2 * (Vmax : ℚ)) := by
  have hVQ : (0 : ℚ) < (Vmax : ℚ) := by exact_mod_cast hV
  have hq : (Vmax : ℚ) < 7 * (maxCircGap E Vmax j : ℚ) := by
    have : Vmax < 7 * maxCircGap E Vmax j := h
    exact_mod_cast this
  rw [lt_div_iff₀ (by positivity)]
  linarith

/-- **The geometric bridge (abstract).**  The midpoint of an empty circular gap `(lo, hi)` of length
`> 1/7` clears BOTH gap endpoints by more than `1/14`; every other phase lies beyond an endpoint, so
the midpoint clears the whole cluster by `> 1/14` — the observer's `1/14`-loneliness witness.  Stated
on the two endpoints (the nearest phases); the outer phases are farther by construction. -/
theorem gap_midpoint_clears (lo hi : ℚ) (hgap : (1 : ℚ) / 7 < hi - lo) :
    (1 : ℚ) / 14 < (lo + hi) / 2 - lo ∧ (1 : ℚ) / 14 < hi - (lo + hi) / 2 :=
  ⟨by linarith, by linarith⟩

/-- **Sanity:** the two `native_decide` clusters of `LRCGoodPeriodMaxgap` (the hard 7-structured
sets where the arc-count route fails) feed the dissociated branch of the dispatch — a good period
exists, hence `M(S) ≥ 1/14` via the clearance bridge. -/
example : HasGoodPeriod [0,7,14,21,26,29,37,44,51,58,67,75,82] 91 :=
  worst7Struct_hasGoodPeriod

end LRC14
