/-
  TournamentH7.LRCGoodPeriodNonStrict — the NON-STRICT good-period certificate (mac-mini-2026-07-09-S64).

  LRC(14) loneliness is `M(S) ≥ 1/14` — NON-strict: equality `M = 1/14` SATISFIES the conjecture
  (the observer reaches distance `1/n`, it need not exceed it).  The strict good-period predicate
  `IsGoodPeriod` (`7·maxCircGap > Vmax`, klein-S200/kps-S99) proves the STRICT `M > 1/14` and so
  MISSES the **knife-edge** clusters where `7·maxCircGap = Vmax` exactly (`M = 1/14`).

  Those clusters are real covering-case inputs — the wraparound-boundary sets `spread = 6·Vmax/7`
  (mac-mini-S64).  Concretely `{0,7,10,14,18,20,21,26,28,35,36,37,42}` at `Vmax = 49 = 7²` has
  `spread = 42 = 6·49/7`, so at `j = 1` the phases fill `[0, 6/7]` and the wraparound gap is `1/7`
  exactly: `7·maxCircGap = 49 = Vmax`.  It has **no** strict good period yet is lonely at `1/14`.
  `hasGoodPeriod_of_dichotomy` (kps-S99) would drop it.

  This file adds the non-strict predicate, the `M ≥ 1/14` clearance bridge, the strict⟹non-strict
  upgrade (so every existing `native_decide` cert still feeds it), the non-strict dispatch, and a
  `native_decide` witness that the knife-edge is caught by non-strict and missed by strict.  Builds
  on `LRCGoodPeriodMaxgap`; the `j = 1` mechanism is `good_period_j1_wraparound_nonstrict`
  (`7·spread ≤ 6·Vmax ⟹ gapLen ≥ 1/7`) in `LRCGoodPeriodJ1`.
-/
import Mathlib
import TournamentH7.LRCGoodPeriodMaxgap

namespace LRC14

/-- **Non-strict good period.** `j` leaves a circular gap `≥ Vmax/7` (`7·maxCircGap ≥ Vmax`) — the
non-strict form matching LRC's `M(S) ≥ 1/14`.  Includes the knife-edge `7·maxCircGap = Vmax`
(`M = 1/14` exactly) that the strict `IsGoodPeriod` excludes. -/
def IsGoodPeriodNonStrict (E : List ℕ) (Vmax j : ℕ) : Prop := Vmax ≤ 7 * maxCircGap E Vmax j

instance (E : List ℕ) (Vmax j : ℕ) : Decidable (IsGoodPeriodNonStrict E Vmax j) :=
  inferInstanceAs (Decidable (Vmax ≤ 7 * maxCircGap E Vmax j))

/-- A cluster **has a non-strict good period** if some `j ∈ {1,…,Vmax−1}` gives a gap `≥ Vmax/7`. -/
def HasGoodPeriodNonStrict (E : List ℕ) (Vmax : ℕ) : Prop :=
  ∃ j ∈ Finset.Ioo 0 Vmax, IsGoodPeriodNonStrict E Vmax j

instance (E : List ℕ) (Vmax : ℕ) : Decidable (HasGoodPeriodNonStrict E Vmax) :=
  inferInstanceAs (Decidable (∃ j ∈ Finset.Ioo 0 Vmax, IsGoodPeriodNonStrict E Vmax j))

/-- **The non-strict clearance bridge — `IsGoodPeriodNonStrict ⟹ M ≥ 1/14`.**  A circular gap
`≥ Vmax/7` has a midpoint at distance `≥ 1/14` (normalized) from the nearest phase, so the observer
runner `Vmax` sits `≥ 1/14` from the whole cluster: `maxCircGap/(2·Vmax) ≥ 1/14`.  Equality is the
knife-edge.  The equality-tolerant twin of `isGoodPeriod_clearance`. -/
theorem isGoodPeriodNonStrict_clearance (E : List ℕ) (Vmax j : ℕ) (hV : 0 < Vmax)
    (h : IsGoodPeriodNonStrict E Vmax j) :
    (1 : ℚ) / 14 ≤ (maxCircGap E Vmax j : ℚ) / (2 * (Vmax : ℚ)) := by
  have hVQ : (0 : ℚ) < (Vmax : ℚ) := by exact_mod_cast hV
  have hq : (Vmax : ℚ) ≤ 7 * (maxCircGap E Vmax j : ℚ) := by
    have : Vmax ≤ 7 * maxCircGap E Vmax j := h
    exact_mod_cast this
  rw [le_div_iff₀ (by positivity)]
  linarith

/-- **Strict upgrades to non-strict** (pointwise). -/
theorem isGoodPeriod_imp_nonStrict (E : List ℕ) (Vmax j : ℕ)
    (h : IsGoodPeriod E Vmax j) : IsGoodPeriodNonStrict E Vmax j :=
  le_of_lt h

/-- **Every strict good-period certificate feeds the non-strict dispatch.**  So the existing
`native_decide` strict certs (`worst7Struct_hasGoodPeriod`, klein-S200) upgrade for free. -/
theorem hasGoodPeriod_imp_nonStrict (E : List ℕ) (Vmax : ℕ)
    (h : HasGoodPeriod E Vmax) : HasGoodPeriodNonStrict E Vmax := by
  obtain ⟨j, hj, hgp⟩ := h
  exact ⟨j, hj, isGoodPeriod_imp_nonStrict E Vmax j hgp⟩

/-- **The non-strict dispatch** — the dichotomy of `hasGoodPeriod_of_dichotomy` on the non-strict
predicate, so it covers the knife-edge.  Near-AP (`k−6 ≤ L`, LEM-012) and dissociated (`L ≤ k−7`,
LEM-013) tile all `L` for `k ≥ 7`; each branch need only supply a NON-strict good period (e.g. the
`j = 1` wraparound `good_period_j1_wraparound_nonstrict` at `spread ≤ 6·Vmax/7`). -/
theorem hasGoodPeriodNonStrict_of_dichotomy (E : List ℕ) (Vmax k L : ℕ) (hk : 7 ≤ k)
    (hNearAP : k - 6 ≤ L → HasGoodPeriodNonStrict E Vmax)
    (hDissoc : L ≤ k - 7 → HasGoodPeriodNonStrict E Vmax) :
    HasGoodPeriodNonStrict E Vmax := by
  by_cases h : k - 6 ≤ L
  · exact hNearAP h
  · exact hDissoc (by omega)

/-- **The knife-edge the strict predicate MISSES (mac-mini-S64).**  The wraparound-boundary cluster
`{0,7,10,14,18,20,21,26,28,35,36,37,42}` at `Vmax = 49 = 7²` (`spread = 42 = 6·49/7`) has best
`7·maxCircGap = 49 = Vmax` — a NON-strict good period (`M = 1/14`, at `j = 1`), but NO strict one. -/
theorem knifeEdge_hasGoodPeriodNonStrict :
    HasGoodPeriodNonStrict [0,7,10,14,18,20,21,26,28,35,36,37,42] 49 := by
  native_decide

/-- …and it has NO strict good period — so `hasGoodPeriod_of_dichotomy` (strict) would drop it,
while the non-strict dispatch catches it.  This is exactly why the non-strict layer is needed. -/
theorem knifeEdge_not_hasGoodPeriod :
    ¬ HasGoodPeriod [0,7,10,14,18,20,21,26,28,35,36,37,42] 49 := by
  native_decide

end LRC14
