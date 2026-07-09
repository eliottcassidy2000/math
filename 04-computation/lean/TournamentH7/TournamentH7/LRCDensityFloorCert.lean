/-
  TournamentH7.LRCDensityFloorCert — the a-priori certification of `∫W > 0`, done honestly
  (kind-pasteur-2026-07-09-S113).

  Goal: certify that the smooth uncovered-measure surrogate `W` is positive somewhere at the continuum,
  a-priori, so that the bridge `LRCSmoothBridge.exists_pos_of_integral_pos` fires.

  **The route is the density floor, NOT the grid average** (klein-S201).  The grid-average identity
  `E_grid[W] = E_x[W] + R_grid` is *not* a valid a-priori certificate of `E_x[W] = ∫W > 0`, because the
  discrepancy bound `|R_grid| ≤ C/V²` FAILS at resonant small rulers: at the tight AP `{0,…,12}` and its
  own ruler `V = 13`, the grid `j/13` lands on `maxgap`'s equidistribution *nulls*, so `E_grid = 1/13 <
  1/7` while `E_x = 0.211 > 1/7` — a discrepancy of `0.134`, nowhere near `O(1/V²)` (klein-S201,
  MISTAKE-127/128/129 lineage: "existence is a max, never a mean").  The grid route is sound only for
  `V ≥ Q+1` (non-resonant), which for the good-period leg is automatic (`spread ≥ 858 > Q`).

  What IS a-priori and grid-free is the **continuum density floor**: `ρ* = vol{x : maxgap > 1/7} =
  vol{W > 0} ≥ m_P = 14249/252252 > 0` (THM-530 Bonferroni, holds even at the AP where `μ_good = 0.44`).
  And a **positive-measure set is nonempty** — so the density floor hands a good point directly, with no
  grid and no integral-positivity machinery (`exists_pos_of_measure_support_pos`).  This file states that
  certification and the honest arithmetic skeleton of the grid↔continuum transfer.  Self-contained apart
  from `LRCReachWitness`.
-/
import Mathlib
import TournamentH7.LRCReachWitness

namespace LonelyRunner
namespace LRC14Concrete

open MeasureTheory

/-- **Positive good-measure ⟹ good point (the density-floor certification).**  If the good set
`{x : 0 < W x}` has positive Lebesgue measure — which the Bonferroni density floor `vol{W>0} ≥ m_P > 0`
provides a-priori, with no grid — then `W` is strictly positive at some point.  A positive-measure set is
nonempty; that is the whole content.  This is the honest replacement for the grid-average route, which
klein-S201 showed fails at resonant rulers. -/
theorem exists_pos_of_measure_support_pos {W : ℝ → ℝ}
    (h : 0 < volume {x | 0 < W x}) : ∃ x, 0 < W x :=
  nonempty_of_measure_ne_zero h.ne'

/-- **The density floor wires to the reach (continuum, grid-free).**  Given the density floor
`0 < vol{W > 0}` (Bonferroni `≥ m_P`, a-priori, holds at every cluster including the AP) and the
reformulation `good x ⟹ ∃ lonely instant`, the runner set is lonely `Mreach ≥ 1/14`.  No grid, no
`R_grid`, no resonant-ruler pathology. -/
theorem mreach_ge_of_good_measure_pos (v : Fin 13 → ℤ) (W : ℝ → ℝ)
    (hpos : 0 < volume {x | 0 < W x})
    (hrefl : (∃ x, 0 < W x) → ∃ τ : ℝ, ∀ i, (1 : ℝ) / 14 ≤ nearInt ((v i : ℝ) * τ)) :
    (1 : ℝ) / 14 ≤ Mreach v :=
  Mreach_ge_of_lonely_instant v (hrefl (exists_pos_of_measure_support_pos hpos))

/-- **Continuum positivity from a grid computation (honest skeleton).**  If the grid average `Eg` and
the continuum integral `Ex` differ by at most `B`, and the *computed* grid value beats the error
(`B < Eg`), then the continuum integral is positive.  Pure arithmetic — always valid; the content is
entirely in whether the input `B` is a correct discrepancy bound.  The analytic bound is
monad-explorer's THM-665 (`|R_grid| ≤ TV(W')/(12·V²)`, Poisson aliasing + bounded variation).
**Caveat (klein-S201):** that bound holds only for non-resonant `V ≥ Q+1`; at resonant rulers (e.g. the
tight AP at `V = 13`) the true discrepancy is `O(1)`, so this skeleton must be fed a regime-correct `B`.
The density floor above is the grid-free alternative that has no such restriction. -/
theorem continuum_pos_of_grid (Eg Ex B : ℝ) (hdisc : |Eg - Ex| ≤ B) (hgrid : B < Eg) : 0 < Ex := by
  have h : Eg - Ex ≤ B := (abs_le.mp hdisc).2
  linarith

/-- **Grid positivity from the continuum floor (the good-period-EXISTS transfer, large-`Vmax`).**  Dual
of the above: if the continuum floor beats the error (`B < Ex`, e.g. `Ex ≥ m_P` and `B ≤ C/V² < m_P`),
the grid average is positive — a good period exists on the ruler.  Sound for `V ≥ Q+1`; the density
floor supplies `Ex`, klein-S205's `V ≳ 1.41·spread` embedding turns the good period into loneliness. -/
theorem grid_pos_of_continuum (Eg Ex B : ℝ) (hdisc : |Eg - Ex| ≤ B) (hcont : B < Ex) : 0 < Eg := by
  have h : Ex - Eg ≤ B := by rw [abs_le] at hdisc; linarith [hdisc.1]
  linarith

end LRC14Concrete
end LonelyRunner
