/-
  TournamentH7.LRCArcCountExistence — the a-priori GLUE of the dissociated good-period branch
  (kind-pasteur-2026-07-09-S95).

  The dissociated branch of THM-527-A (LRC(14) covering leg) closes by mac-mini-S61/opus-S168's
  ARC-COUNT existence: the good periods number `#good ≥ ρ*·Vmax − #arcs`, so a good period EXISTS
  as soon as `#arcs < ρ*·Vmax`.  Route (c) supplies this a-priori via `ρ* ≥ D3(E)` (THM-661) and
  `#arcs ≤ c·Vmax` with `c := #arcs/spread < D3(E)`.

  This file formalizes the ELEMENTARY LOGIC of that closure — the implications that turn the two
  cited analytic inputs (`ρ* ≥ D3` and the arc-count reduction) into `0 < #good`.  The analytic
  inputs themselves (the covering-moment bound `ρ* ≥ D3`, mac-mini/opus's `#arcs = o(Vmax)`) enter
  as hypotheses, exactly as they are separately established.  Also formalizes the complementary
  AVERAGING existence route (kps-S95): if the mean largest-gap over dilations beats `1/7`, some
  dilation is good (`max ≥ mean`), and the second-moment (contraharmonic) lower bound on `maxgap`
  that the lemniscate `r⁴ = r²cos2θ` metaphor suggests.  Self-contained (imports only Mathlib).
-/
import Mathlib

namespace LRC14ArcCount

/-- **Arc-count existence (the reduction's payoff).** Given the arc-count reduction
`#good ≥ ρ*·Vmax − #arcs` (mac-mini-S58/S61) and `#arcs < ρ*·Vmax`, a good period exists
(`0 < #good`). -/
theorem arccount_existence (rho arcs good Vmax : ℚ)
    (hred : rho * Vmax - arcs ≤ good) (hlt : arcs < rho * Vmax) :
    0 < good := by linarith

/-- **`c < D3` forces the arc-count inequality.** With `Vmax > 0`, the covering-moment bound
`D3 ≤ ρ*` (THM-661), the arc bound `#arcs ≤ c·Vmax`, and Route-(c)'s `c < D3`, we get
`#arcs < ρ*·Vmax`. -/
theorem c_lt_D3_forces_arccount (rho D3 c arcs Vmax : ℚ)
    (hV : 0 < Vmax) (hrho : D3 ≤ rho) (harcs : arcs ≤ c * Vmax) (hc : c < D3) :
    arcs < rho * Vmax := by
  have h1 : c * Vmax < D3 * Vmax := mul_lt_mul_of_pos_right hc hV
  have h2 : D3 * Vmax ≤ rho * Vmax := mul_le_mul_of_nonneg_right hrho (le_of_lt hV)
  linarith

/-- **Route (c): `c < D3` ⟹ a good period exists.** The full a-priori glue of the dissociated
branch: `Vmax > 0`, `D3 ≤ ρ*`, `#arcs ≤ c·Vmax`, `c < D3`, and the reduction
`#good ≥ ρ*·Vmax − #arcs` together give `0 < #good`. -/
theorem c_lt_D3_existence (rho D3 c arcs good Vmax : ℚ)
    (hV : 0 < Vmax) (hrho : D3 ≤ rho) (harcs : arcs ≤ c * Vmax) (hc : c < D3)
    (hred : rho * Vmax - arcs ≤ good) :
    0 < good :=
  arccount_existence rho arcs good Vmax hred
    (c_lt_D3_forces_arccount rho D3 c arcs Vmax hV hrho harcs hc)

/-- **Averaging existence (kps-S95, complementary route).** If the MEAN largest-gap over the
dilations `j ∈ J` exceeds the threshold `thr = Vmax/7` — i.e. `|J|·thr < Σ_j maxgap(j)` — then some
dilation `j` has `maxgap(j) > thr`: a good period.  Pure `max ≥ mean`.  (kps-S95 verified the mean
ratio `≥ 1.047 > 1` adversarially for dissociated 13-clusters across all spreads.) -/
theorem averaging_existence (J : Finset ℕ) (maxgap : ℕ → ℚ) (thr : ℚ)
    (havg : (J.card : ℚ) * thr < ∑ j ∈ J, maxgap j) :
    ∃ j ∈ J, thr < maxgap j := by
  by_contra h
  push_neg at h
  have hsum : ∑ j ∈ J, maxgap j ≤ ∑ _j ∈ J, thr := Finset.sum_le_sum (fun j hj => h j hj)
  rw [Finset.sum_const, nsmul_eq_mul] at hsum
  linarith

/-- **The second-moment (contraharmonic) lower bound on `maxgap`** — the lemniscate `r⁴ = r²cos2θ`
second-moment cue.  For nonnegative gaps `g i` summing to `total > 0`, each `≤ mg`, the
contraharmonic mean `Σ g² / total ≤ mg`.  (kps-S95: this bound alone gives mean ratio `≈ 0.85 < 1`,
so it is necessary but not sufficient — the true `maxgap` is needed; recorded as the clean
elementary inequality behind the averaging route.) -/
theorem maxgap_ge_contraharmonic (G : Finset ℕ) (g : ℕ → ℚ) (mg total : ℚ)
    (hg : ∀ i ∈ G, 0 ≤ g i) (hmax : ∀ i ∈ G, g i ≤ mg)
    (htot : total = ∑ i ∈ G, g i) (hpos : 0 < total) :
    (∑ i ∈ G, (g i)^2) / total ≤ mg := by
  have hsq : ∑ i ∈ G, (g i)^2 ≤ mg * total := by
    rw [htot, Finset.mul_sum]
    apply Finset.sum_le_sum
    intro i hi
    rw [pow_two]
    exact mul_le_mul_of_nonneg_right (hmax i hi) (hg i hi)
  rw [div_le_iff₀ hpos]
  linarith

-- Sanity: c = 0.5 < D3 = 0.6 ≤ ρ* = 0.9, #arcs ≤ 0.5·Vmax, Vmax = 100 ⟹ a good period exists.
example : (0 : ℚ) < 40 :=
  c_lt_D3_existence 0.9 0.6 0.5 50 40 100 (by norm_num) (by norm_num) (by norm_num)
    (by norm_num) (by norm_num)

end LRC14ArcCount
