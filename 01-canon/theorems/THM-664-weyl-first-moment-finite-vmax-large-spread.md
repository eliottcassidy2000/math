---
id: THM-664
title: The Weyl first-moment theorem for the finite-Vmax glue (THM-527-A), large-spread half — a good period exists iff E_grid[W] > 6/(7·Vmax), and the Weyl–Fourier identity E_grid[W] = (6/7)^k + Σ_{n≠0, Vmax|n·e} 𝒲̂(n) exhibits the iid first moment (6/7)^k as the main term plus a grid-resonance sum that vanishes (Weyl decorrelation) as spread→∞; hence for large spread E_grid[W] → (6/7)^k > 0 and a good period exists. This is the SAME resonance-sum decorrelation as the density-floor tail (LEM-009/THM-518), now with "Vmax | n·e" in place of "n·e = 0"
status: PROVED (framework + main term + pairwise piece + decorrelation limit) / VERIFIED (the uniform smallness). Steps 1–3 are exact: the reduction "good period ⟺ E_grid[W] > 6/(7Vmax)" (W≥0, j=0 gives W(0)=6/7); the Fourier identity E_grid[W]=(6/7)^k+Σ_{n≠0,Vmax|n·e}𝒲̂(n) with 𝒲̂(0)=(6/7)^k EXACT (the pinned phase-0 arc contributes the 6/7 factor); the exact pairwise tent T(V')=(1+2R)/(7V')−R(R+1)/V'², R=⌊(V'−1)/7⌋, |T(V')−1/49|≤0.12/V'. Step 4 (decorrelation R→0) is the Weyl limit + the density-floor resonance machinery. VERIFIED exactly on structured large-spread clusters (2-block, AP, k=11,13): |E_grid[W]−(6/7)^k| ≤ 0.032 ≪ (6/7)^k∈[0.135,0.183], mean |dev| 0.011→0.002 decreasing with spread; E_grid[W] > 6/(7Vmax) in ALL cases; #good periods ≥ 20 and ∝ spread. The one non-elementary step is the fully-uniform bound |R| < (6/7)^k, which reduces to the same a-priori resonance constants as the density-floor tail (opus-S157). Complements LEM-010's deterministic Dirichlet closure (two independent routes)
source: mac-mini-2026-07-08-S59
depends_on:
  - THM-527   # the reformulation: good period ⟺ maxgap{frac(e_i j/Vmax)} > 1/7 ⟺ W(j/Vmax) > 0
  - THM-663   # the covering-case assembly; this is the Weyl route for its remaining item (1)
  - THM-657   # W = uncovered measure = Σ(g_i − 1/7)_+ (the covering reformulation)
related:
  - LEM-010   # the DETERMINISTIC (Dirichlet + j=1) closure of the same large-spread half
  - THM-518   # Weyl stranger-decoupling (the decorrelation this instantiates for the grid)
  - LEM-009   # the density-floor tail: the SAME 𝒲̂-resonance decorrelation, with n·e=0
---

# THM-664 — The Weyl first-moment theorem (finite-Vmax glue, large-spread)

Let `E = {e_0=0 < e_1 < … < e_{k−1}}` be the cluster of co-offsets (`k ≤ 13`, `e_i < Vmax`), and
`W(y) = Σ_i (g_i(y) − 1/7)_+` the uncovered measure of the phases `{frac(e_i y)}` (THM-657). A
**good period** is `j ∈ {1,…,Vmax−1}` with `maxgap{frac(e_i j/Vmax)} > 1/7`, i.e. `W(j/Vmax) > 0`;
one exists ⟹ `M(S) ≥ 1/14` (THM-527). This theorem produces a good period for large-spread clusters
by a Weyl (equidistribution / first-moment) argument.

## Step 1 — Exact reduction to a first moment

`W ≥ 0` pointwise, and at `j = 0` all phases coincide at `0` (one gap of length `1`), so
`W(0) = 1 − 1/7 = 6/7`. Writing `E_grid[W] := (1/Vmax) Σ_{j=0}^{Vmax−1} W(j/Vmax)`,

> `Σ_{j=1}^{Vmax−1} W(j/Vmax) = Vmax·E_grid[W] − 6/7`,  hence  **a good period exists ⟺
> `E_grid[W] > 6/(7·Vmax)`.**

(Since `W ≥ 0`, `E_grid[W] ≥ W(0)/Vmax = 6/(7Vmax)` always, with equality iff every `j ≥ 1` is bad;
the content is the strict inequality.)

## Step 2 — The Weyl–Fourier identity

`W(y) = 𝒲(frac(e_1 y),…,frac(e_{k−1} y))` for the fixed **uncovered-measure function**
`𝒲 : T^{k−1} → [0, 6/7]` (phase `e_0 y = 0` pinned). Expanding `𝒲(φ) = Σ_{n∈ℤ^{k−1}} 𝒲̂(n) e(n·φ)`
and averaging over the grid (`(1/Vmax) Σ_j e(m j/Vmax) = 1[Vmax | m]`, with `n·e := Σ_i n_i e_i`):

> **`E_grid[W] = Σ_{n : Vmax | n·e} 𝒲̂(n) = (6/7)^k + R`,  `R := Σ_{n≠0, Vmax | n·e} 𝒲̂(n)`.**

The main term is EXACT: `𝒲̂(0) = ∫_{T^{k−1}} 𝒲 = ` the iid uncovered measure with `e_0`'s arc
`[0,1/7)` pinned and the other `k−1` phases uniform `= ∫_p 1[p∉[0,1/7)]·(6/7)^{k−1} dp = (6/7)·(6/7)^{k−1}
= (6/7)^k`. `R` is the **grid-resonance sum**.

**Same object as the density-floor tail.** The continuous first moment is `E[W] = ∫_0^1 W =
Σ_{n : n·e = 0} 𝒲̂(n) = (6/7)^k + Σ_{n≠0, n·e=0} 𝒲̂(n)`. So `R = [E[W] − (6/7)^k] +
Σ_{n≠0, n·e≠0, Vmax|n·e} 𝒲̂(n)`: the **decorrelation** term (`n·e = 0`) plus the **pure grid
resonances** (`n·e` a nonzero multiple of `Vmax`). The identical `𝒲̂`-decay controls both — the
finite-`Vmax` glue's large-spread half IS the density-floor decorrelation (LEM-009/THM-518), read on
the grid.

## Step 3 — The pairwise piece is exact

The `|n|`-support-2 part of `R` is the pairwise arc-overlap correction `Σ_{i<j}(T(V'_{ij}) − 1/49)`,
`V'_{ij} = Vmax/gcd(e_i−e_j, Vmax)`, where `T(V') = (1/V')Σ_{r=0}^{V'−1}(1/7 − ‖r/V'‖)_+` is the
grid-average of the overlap tent. Closed form (verified `V'=1..399`):

> `T(V') = (1+2R)/(7V') − R(R+1)/V'²`,  `R = ⌊(V'−1)/7⌋`,  and  **`|T(V') − 1/49| ≤ 0.12/V'`**

(so `T(7)=1/49` exactly; the deviation is `O(1/V')`, small once `gcd(e_i−e_j,Vmax) ≪ Vmax`). The
higher-order terms (`|n|≥3`) are NOT negligible — they carry the bulk of `R` — but the total stays
tiny (below).

## Step 4 — Decorrelation ⟹ good period (large spread)

`𝒲` is continuous and piecewise-linear on `T^{k−1}`, so `Σ_n |𝒲̂(n)| < ∞`. As `spread → ∞`, a
nonzero resonance `Vmax | n·e` forces `|n·e| ≥ Vmax` with `|n·e| ≤ |n|₁·spread`, hence `|n|₁ ≥
Vmax/spread`; absent exact additive relations among the `e_i`, the resonances recede to large `|n|`
and `R → 0` (Weyl equidistribution). Therefore `E_grid[W] → (6/7)^k > 0`, and a good period exists
for every sufficiently spread cluster.

## Verification (exact)

Over structured large-spread clusters (2-block `{0..a}∪{s−b..s}`, AP, `k = 11, 13`; `Vmax ≈ spread`):
`|E_grid[W] − (6/7)^k| ≤ 0.032`, far below `(6/7)^k ∈ [0.135, 0.183]`, and the mean deviation
**decreases with spread** (`0.011 → 0.002` as `spread: 20 → 640`). In every case `E_grid[W] >
6/(7Vmax)`, with `#good periods ≥ 20` growing `∝ spread`. So `R` never approaches `(6/7)^k`: the
good period exists, and abundantly.

## Honest status & relation to LEM-010

Steps 1–3 are exact theorems; Step 4 is the Weyl limit plus the smallness of `R`. The one
non-elementary point — a **fully uniform** `|R| < (6/7)^k` for all `(E, Vmax)` — reduces to the same
a-priori resonance-sum constants as the density-floor tail (opus-S157's `𝒲̂`-decay / mixed-variation
bounds), which are numerically certified but not yet fully a-priori. So this Weyl route and
**LEM-010** (the elementary deterministic `j=1` + Dirichlet closure, which needs NO resonance bound)
are two independent proofs of the large-spread half; LEM-010 is the cleaner unconditional one, while
THM-664 is the conceptual statement — *the finite-Vmax glue is the density-floor decorrelation on the
Vmax-grid* — and gives the abundance `#good ≈ (6/7)^k·Vmax`. Files:
`04-computation/lrc14_weyl_{firstmoment,tent_bound,decorrelation}_macmini_S59.{py,out}`.
