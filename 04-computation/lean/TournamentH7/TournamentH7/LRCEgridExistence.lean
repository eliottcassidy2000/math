/-
  TournamentH7.LRCEgridExistence — the grid-mean existence route for the dissociated good-period
  branch (kind-pasteur-2026-07-09-S96).

  A good period exists ⟺ some dilation `j` leaves uncovered measure `W j > 0` ⟺ the grid mean
  `E_grid[W] = (1/Vmax)·Σ_{j<Vmax} W j > 0`.  By LEM-011 (klein-S194), `E_grid[W] = (6/7)^k + R`
  with `R = Σ_{Vmax | n·e} 𝒲̂(n)` the resonance residual.  Hence a good period exists as soon as

      |R| < (6/7)^k.

  This SIDESTEPS the arc-count `#arcs` bound (opus-S169's one open item): it needs only that the
  resonance residual is smaller than the main term `(6/7)^k`, an inequality on the NEAR-RESONANCE
  COUNT (kps-S93) — few low-height `n` with `Vmax | n·e` for dissociated (Sidon-like) `E`, which is
  Mertens-SAFE (a count, not a cancellation).  Verified (kps-S96): `|R|/(6/7)^k ≤ 0.69 < 1` over
  dissociated clusters at the critical `Vmax`, INCLUDING mac-mini's hard `7`-structured case (diffs
  `≡ 0 mod 7`) that broke the `c<D3` certificate.

  This file formalizes the ELEMENTARY chain `|R| < mean ⟹ E_grid > 0 ⟹ ∃ good j`.  The two inputs —
  the LEM-011 identity `E_grid = mean + R` and the residual bound `|R| < mean` — enter as hypotheses.
  Self-contained (imports only Mathlib).
-/
import Mathlib

namespace LRC14Egrid

/-- **The grid-mean is positive when the residual is small.** If the grid sum decomposes as
`Σ_{j<N} W j = N·(mean + R)` (LEM-011) and `|R| < mean`, then the sum is strictly positive. -/
theorem gridsum_pos_of_residual_small (W : ℕ → ℚ) (N : ℕ) (mean R : ℚ)
    (hN : 0 < N) (hid : ∑ j ∈ Finset.range N, W j = (N : ℚ) * (mean + R))
    (hlt : |R| < mean) :
    0 < ∑ j ∈ Finset.range N, W j := by
  have hR : -mean < R := (abs_lt.mp hlt).1
  have hpos : (0 : ℚ) < mean + R := by linarith
  have hNQ : (0 : ℚ) < (N : ℚ) := by exact_mod_cast hN
  rw [hid]
  exact mul_pos hNQ hpos

/-- **Grid-mean existence.** If the uncovered measures `W j` are nonnegative and their grid sum is
positive, then some dilation `j < N` is a good period (`W j > 0`).  A sum of nonnegatives is positive
exactly when a summand is. -/
theorem exists_good_of_gridsum_pos (W : ℕ → ℚ) (N : ℕ)
    (_hnn : ∀ j, 0 ≤ W j) (hsum : 0 < ∑ j ∈ Finset.range N, W j) :
    ∃ j ∈ Finset.range N, 0 < W j := by
  by_contra h
  push_neg at h
  have hle : ∑ j ∈ Finset.range N, W j ≤ 0 := Finset.sum_nonpos (fun j hj => h j hj)
  linarith

/-- **The full E_grid route (kps-S96): `|R| < (6/7)^k ⟹ a good period exists`.** Combining the two
steps: with the LEM-011 identity `Σ W = N·((6/7)^k + R)`, nonnegative `W`, and the residual bound
`|R| < (6/7)^k`, some dilation `j < N` has `W j > 0` — a good period.  No `#arcs`. -/
theorem exists_good_of_residual_small (W : ℕ → ℚ) (N k : ℕ) (R : ℚ)
    (hN : 0 < N) (hnn : ∀ j, 0 ≤ W j)
    (hid : ∑ j ∈ Finset.range N, W j = (N : ℚ) * ((6 / 7 : ℚ) ^ k + R))
    (hlt : |R| < (6 / 7 : ℚ) ^ k) :
    ∃ j ∈ Finset.range N, 0 < W j :=
  exists_good_of_gridsum_pos W N hnn
    (gridsum_pos_of_residual_small W N ((6 / 7 : ℚ) ^ k) R hN hid hlt)

-- Sanity: k = 13, N = 93, residual R = 0.09 (≈ the measured 0.69·(6/7)^13) is `< (6/7)^13 ≈ 0.1348`,
-- so a good period exists whenever the grid sum has this shape.
example (W : ℕ → ℚ) (hnn : ∀ j, 0 ≤ W j)
    (hid : ∑ j ∈ Finset.range 93, W j = (93 : ℚ) * ((6 / 7 : ℚ) ^ 13 + (9/100 : ℚ))) :
    ∃ j ∈ Finset.range 93, 0 < W j :=
  exists_good_of_residual_small W 93 13 (9/100) (by norm_num) hnn hid (by norm_num)

end LRC14Egrid
