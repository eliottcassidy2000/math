/-
  TournamentH7.LRCEgridExistenceCorrected — the E_grid good-period route with the trivial-shift
  correction and the R₀-signed / R_grid-absolute split (klein-2026-07-09-S202).

  CORRECTS a blind spot in `LRCEgridExistence.exists_good_of_residual_small` (kps-S96): it concludes
  `∃ j ∈ Finset.range N, 0 < W j`, but `Finset.range N` INCLUDES `j = 0`, the trivial shift, where
  `W 0 = 6/7 > 0` ALWAYS.  So that conclusion is a tautology and does NOT certify a good period
  (which needs `j ∈ {1,…,N−1}`).  Likewise `E_grid[W] > 0` is trivially true (the `j=0` term alone
  contributes `6/7`).  This is exactly the tight-AP failure (klein-S201): `E = {0,…,12}` at `N = 13`
  has `E_grid[W] = (6/7)/13` — every nonzero shift gives `W = 0`, no good period — yet `E_grid[W] > 0`.

  THE CORRECTION: the good-period threshold is `E_grid[W] > (6/7)/N`, i.e. the grid sum must EXCEED
  the trivial head term `W 0 = 6/7`, and the conclusion is a NONTRIVIAL `j ≠ 0` with `W j > 0`.
  With the klein-S202 split `E_grid[W] = E[W] + R_grid` (`E[W] = (6/7)^k + R₀`, the `n·e=0` exact
  relations / density floor kept SIGNED; `R_grid` the pure grid resonances `n·e = mN`), and
  `|R_grid| ≤ R_grid_abs`, the threshold becomes the SOUND, verified condition

      R_grid_abs  <  E[W] − (6/7)/N   ⟹   a good period exists.

  Keeping `R₀` signed inside `E[W]` (bounded `≥ 0` by the density floor) sidesteps the total-absolute
  wall (kps-S98: `Σ_{N|n·e}|𝒲̂(n)|` is NOT uniformly `< (6/7)^k`).  Self-contained (imports Mathlib).
-/
import Mathlib

namespace LRC14EgridCorrected

/-- **Nontrivial grid existence (the trivial-shift correction).** If the `W j` are nonnegative, the
trivial shift has value `W 0 = h`, and the grid sum STRICTLY EXCEEDS `h`, then some NONTRIVIAL
dilation `j ∈ {1,…,N−1}` has `W j > 0` — a genuine good period, not the trivial `j = 0`. -/
theorem exists_nontrivial_good_of_gridsum_gt_head
    (W : ℕ → ℚ) (N : ℕ) (h : ℚ)
    (hN : 0 < N) (_hnn : ∀ j, 0 ≤ W j) (hhead : W 0 = h)
    (hgt : h < ∑ j ∈ Finset.range N, W j) :
    ∃ j ∈ Finset.range N, j ≠ 0 ∧ 0 < W j := by
  have h0 : (0 : ℕ) ∈ Finset.range N := Finset.mem_range.mpr hN
  have hsplit : W 0 + ∑ j ∈ (Finset.range N).erase 0, W j = ∑ j ∈ Finset.range N, W j :=
    Finset.add_sum_erase _ W h0
  have hpos : 0 < ∑ j ∈ (Finset.range N).erase 0, W j := by
    rw [hhead] at hsplit; linarith
  by_contra hc
  push_neg at hc
  have hle : ∑ j ∈ (Finset.range N).erase 0, W j ≤ 0 := by
    apply Finset.sum_nonpos
    intro j hj
    rw [Finset.mem_erase] at hj
    exact hc j hj.2 hj.1
  linarith

/-- **The threshold is `E_grid[W] > (6/7)/N`, via the R₀-signed / R_grid-absolute split.** Given the
klein-S202 identity `Σ_{j<N} W j = N·(EW + Rgrid)` (`EW = E[W]`, the continuum mean with the exact
`n·e=0` relations kept signed inside it), the trivial head `W 0 = 6/7`, a bound `|Rgrid| ≤ RgridAbs`
on the pure-grid residual, and the SOUND condition `RgridAbs < EW − (6/7)/N`, a nontrivial good period
exists. -/
theorem exists_good_of_grid_residual
    (W : ℕ → ℚ) (N : ℕ) (EW Rgrid RgridAbs : ℚ)
    (hN : 0 < N) (hnn : ∀ j, 0 ≤ W j) (hhead : W 0 = (6 / 7 : ℚ))
    (hid : ∑ j ∈ Finset.range N, W j = (N : ℚ) * (EW + Rgrid))
    (hbound : |Rgrid| ≤ RgridAbs)
    (hthr : RgridAbs < EW - (6 / 7 : ℚ) / N) :
    ∃ j ∈ Finset.range N, j ≠ 0 ∧ 0 < W j := by
  have hNQ : (0 : ℚ) < (N : ℚ) := by exact_mod_cast hN
  have hRlow : -RgridAbs ≤ Rgrid := (abs_le.mp hbound).1
  -- N·(EW + Rgrid) ≥ N·(EW − RgridAbs) > N·((6/7)/N) = 6/7
  have hstep : (6 / 7 : ℚ) < (N : ℚ) * (EW + Rgrid) := by
    have h1 : (N : ℚ) * (EW - RgridAbs) ≤ (N : ℚ) * (EW + Rgrid) := by
      apply mul_le_mul_of_nonneg_left _ (le_of_lt hNQ); linarith
    have h2 : (6 / 7 : ℚ) < (N : ℚ) * (EW - RgridAbs) := by
      have := (mul_lt_mul_of_pos_left hthr hNQ)
      rw [mul_sub, mul_div_cancel₀ _ (ne_of_gt hNQ)] at this
      linarith
    linarith
  exact exists_nontrivial_good_of_gridsum_gt_head W N (6 / 7 : ℚ) hN hnn hhead
    (by rw [hid]; exact hstep)

/-- Sanity — the TIGHT AP `{0,…,12}` at `N = 13` is correctly EXCLUDED.  There `E_grid[W] = (6/7)/13`
exactly (grid sum `= 6/7`, the head), so `hgt : 6/7 < 6/7` is FALSE — no nontrivial good period is
claimed (klein-S201).  The corrected lemma cannot be misapplied to it. -/
example : ¬ ((6 / 7 : ℚ) < (6 / 7 : ℚ)) := lt_irrefl _

/-- Sanity — a genuine case fires.  `k = 13`, `N = 91`, `E[W] = (6/7)^13 + R₀` with `R₀ = 0` (take the
density-floor value at its main term) and `RgridAbs = 1/100 < (6/7)^13 − (6/7)/91 ≈ 0.1348 − 0.0094`.
A nontrivial good period exists. -/
example (W : ℕ → ℚ) (hnn : ∀ j, 0 ≤ W j) (hhead : W 0 = (6 / 7 : ℚ))
    (hid : ∑ j ∈ Finset.range 91, W j
         = (91 : ℚ) * (((6 / 7 : ℚ) ^ 13) + (1 / 200 : ℚ))) :
    ∃ j ∈ Finset.range 91, j ≠ 0 ∧ 0 < W j :=
  exists_good_of_grid_residual W 91 ((6 / 7 : ℚ) ^ 13) (1 / 200) (1 / 100)
    (by norm_num) hnn hhead hid (by norm_num) (by norm_num)

end LRC14EgridCorrected
