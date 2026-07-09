---
id: THM-665
title: The sharp aliasing bound for the smooth surrogate — E_grid[W](V) = Σ_m Ŵ(mV) exactly, and |E_grid[W] − ∫W| ≤ TV(W′)/(12 V²); with the exact TV(W′) ledger (TV ≈ 12.2·spread², so V₀ ≈ 2.8·spread — kps-S108's window CONFIRMED from the proved side) and the structural corollary that the a-priori existence route never fires on covering clusters (V/spread ≈ 1 there): the bounded window IS the covering case
status: PROVED (three steps, below; the identity and both constants are rigorous). Machine-verified: bound holds 12/12 rows on kps-S108's two clusters at V = 200..6400 (actual slack 8–418×); aliasing identity verified to 1e-5 (|m| ≤ 59 partial sums); exact TV(W′) ledger over 6 clusters.
source: monad-explorer-2026-07-09-S1 (HYP-5707) — the brick named in kps-S108/kps-S112 NEXT ("formalize E_grid[W] = E_x[W] + R_grid + the O(1/V²) bound / sharpen to shrink the window").
depends_on: []
related:
  - THM-661   # the moment-bound density floor whose E[W] ingredient is ∫W here
  - THM-663   # the covering-case closure this bound's window analysis feeds
  - HYP-5157  # the S11 cell engine computing the exact ledger; W = the S11 excess mass V_θ
external: Poisson summation / aliasing for absolutely convergent Fourier series; BV coefficient bounds. Standard analysis, self-contained proof below.
---

# THM-665 — the sharp aliasing bound for the smooth surrogate

## Setting

`E = {e_1 < … < e_k} ⊂ Z` a co-offset cluster, `θ = 1/7`,
`W(x) = Σ_i (g_i(x) − θ)₊` the smooth surrogate (opus-S170 = the S11 excess mass), where
`g_i(x)` are the circular gaps of `{frac(e_i x)}`. `W` is 1-periodic, **continuous**
(sorted gaps vary continuously through collisions), **piecewise linear** (breakpoints:
phase collisions `m/|e_i − e_j|`, wraps `m/e_i`, θ-crossings), and `W′` is **piecewise
constant with integer values** (off breakpoints, `W′ = Σ_{g_i > θ} (slope of g_i)`, each
slope a difference of two co-offsets). `TV(W′)` denotes the total variation of `W′` over
the circle (the sum of |slope jumps|, including the periodic seam).

The Vmax-ruler grid average is `E_grid[W](V) = (1/V) Σ_{j=0}^{V−1} W(j/V)`.

## Statement

**(i) (Exact aliasing identity.)**  `E_grid[W](V) = Σ_{m ∈ Z} Ŵ(mV)`, absolutely
convergent; in particular `R_grid := E_grid[W](V) − ∫₀¹W = Σ_{m ≠ 0} Ŵ(mV)`.

**(ii) (Coefficient bound.)**  `|Ŵ(r)| ≤ TV(W′)/(4π²r²)` for every integer `r ≠ 0`.

**(iii) (The sharp bound.)**
> `|E_grid[W](V) − ∫₀¹W| ≤ TV(W′) / (12 V²)`.

Hence `E_grid[W](V) > 0` (a good period EXISTS on the V-ruler) whenever
`V > V₀ := sqrt( TV(W′) / (12 ∫W) )`.

## Proof

**(i).** `W` is continuous and piecewise `C¹` with `W′` of bounded variation, so by (ii)
its Fourier coefficients are `O(1/m²)` and the Fourier series converges absolutely,
hence pointwise to `W` everywhere. Summing the series at the `V` grid points and using
`(1/V) Σ_{j<V} e(rj/V) = 1[V | r]` (finite geometric sum), the interchange justified by
absolute convergence:
`(1/V)Σ_j W(j/V) = Σ_r Ŵ(r)·1[V|r] = Σ_m Ŵ(mV)`.  The `m = 0` term is `∫W`. ∎

**(ii).** Two integrations by parts. First, `W` continuous, periodic, piecewise `C¹`:
`Ŵ(r) = Ŵ′(r)/(2πir)` (no boundary terms). Second, `W′` is a periodic BV function
(piecewise constant, jump measure `dW′ = Σ_c σ_c δ_{x_c}` with `Σ_c |σ_c| = TV(W′)`):
`Ŵ′(r) = (1/(2πir)) ∫ e(−rx) dW′(x)`, so `|Ŵ′(r)| ≤ TV(W′)/(2π|r|)`. Compose. ∎

**(iii).** `|R_grid| ≤ Σ_{m≠0} TV(W′)/(4π²m²V²) = (TV(W′)/(4π²V²))·2ζ(2)
= TV(W′)·(π²/3)/(4π²V²) = TV(W′)/(12V²)`. ∎

## The exact TV(W′) ledger (S11 cell engine; per-shape exact integers)

| cluster | spread | TV(W′) | TV/spread² | ∫W | V₀ (sharp) | V₀/spread |
|---|---|---|---|---|---|---|
| kps-S108 dissociated | 35 | 14904 | 12.17 | 0.13970 | 94.3 | 2.69 |
| kps-S108 7-structured (= mac-mini @91) | 82 | 81860 | 12.17 | 0.12930 | 229.7 | 2.80 |
| klein-S206 worst covering (Vmax 17) | 16 | 2628 | 10.27 | 0.12185 | 42.4 | 2.65 |
| tight AP (non-covering reference) | 12 | 1640 | 11.39 | 0.12662 | 32.9 | 2.74 |
| GW co-offsets | 23 | 2984 | 5.64 | 0.11488 | 46.5 | 2.02 |

**The TV law (empirical, exact values):** `TV(W′) ≈ 12.2·spread²` with remarkable
stability across structure classes (10.3–12.2; GW lower at 5.6). A crude provable
closed form is `TV(W′) ≤ 2 Σ_{i<j}(e_i−e_j)² + (θ-crossing term of the same order)`
(each collision of the pair `(i,j)` occurs `|e_i−e_j|` times per period and jumps `W′`
by at most `2|e_i−e_j|` — and only when a neighbor gap exceeds θ, which is why the true
constant is ~12, far below the crude `2·C(k,2)·avg ≈ k²`-scale bound).

## Consequences and honest reconciliation

1. **kps-S108's window is CONFIRMED from the proved side.** With the sharp constant
   `1/12` and the measured law `TV ≈ 12.2·spread²`, `V₀ = sqrt(TV/(12∫W)) ≈
   spread/sqrt(∫W) ≈ 2.7–2.8·spread` — exactly the "V ≳ 2.8·spread" kps-S108 reported.
   The `π²/3` envelope written in their prose is `4π² ≈ 39.5×` looser than the true
   constant (it would have put V₀ at ≈ 17·spread); their numeric window was right
   because the empirical fit, not the envelope, set it. No result is damaged; the
   *formalizable* statement is (iii) above with `1/12`.
2. **The a-priori existence route NEVER fires on covering clusters.** Covering velocity
   sets contain small speeds, so co-offsets reach `≈ Vmax` and `spread ≈ Vmax`, i.e.
   `V/spread ≈ 1 < 2.7`. The entire covering case (= the entire remaining LRC(14)
   content, klein-S206/THM-663) lives INSIDE the bounded window. The aliasing bound's
   role there is not existence-certification but (a) the formal `1/V²` convergence
   statement kps-S112 asked for, and (b) making the per-(E,V) finite checks certified
   (E_grid[W](V) is a finite sum of exact rationals — positivity is decidable).
3. **The named room: signed cancellation is square-root.** Measured `|R_grid|` runs
   8–418× BELOW the bound (iii): the corner phases `e(−mV·x_c)` cancel like
   `sqrt(#corners)` (≈ 570 corners for spread 35, `sqrt ≈ 24` — matching the observed
   generic slack). Proving that square-root cancellation is an equidistribution
   statement about the corner positions relative to the `1/V` grid — i.e. **the same
   Diophantine/Kronecker node** (kps-S112) the realization already needs. The two open
   quantitative questions are one question.
4. **What would shrink the window below the drift threshold `1.41·spread`:** a proven
   `|Σ_m Ŵ(mV)| ≤ TV/(c·V²)` with `c ≥ 47` (= 12·(2.8/1.41)²) — i.e. any proven
   factor-4 of the measured square-root cancellation. Then existence (this bound) +
   drift embed (klein-S205) would cover all `V > 1.41·spread` and the finite window
   would shrink to `(spread, 1.41·spread]`.

## Verification & files

`04-computation/lrc14_sharp_aliasing_tv_monad_S1.py` (+ `.out`): exact TV/∫W ledger
(6 clusters); bound validation 12/12 rows at V = 200..6400 with slack 8–418×; aliasing
identity partial sums agree with measured `R_grid` to 1e-5; the window table vs the
1.41·spread drift threshold.

---

> **ADDENDUM (klein-2026-07-09-S208, HYP-5691 — convention refinement of consequence 2;
> the theorem itself is untouched).** Consequence 2's "the a-priori existence route NEVER
> fires on covering clusters (spread ≈ Vmax)" is correct in the ALL-RUNNER co-offset
> convention (E = {Vmax − v : v ∈ S}), which the ledger rows use (e.g. the klein-S206
> worst covering row has spread 16 at Vmax 17 = all-runner). Under THM-527's actual
> architecture, the split S = P ∪ L sends the small speeds to G_P (measure side) and the
> cluster spread is Vmax − min L, which the covering constraint does NOT force to be
> proportional: covering duty for the q's P misses is dischargeable inside
> [(9/14)Vmax, Vmax] once Vmax ≥ ~40. CENSUS (exact, 60 instances, k ≥ 8):
> on the confined stratum (r > 2.8) the aliasing existence AND the Lean-faithful drift
> embed fire 20/20 — the a-priori chain is complete there; on the proportional stratum
> (r ≤ 1.41) they fire 0/13 and 1/34. So the corrected scope statement is: **the
> a-priori route fires exactly on the confined-L stratum of the covering family; the
> residual is the proportional stratum**, which THM-667's ladder further localizes to
> the mid-band (Vmax/14, 9·Vmax/14). Files: lrc14_split_ratio_census_klein_S208.py(+out).
