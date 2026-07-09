---
id: LEM-012
title: The gap-splitting good-period bound — if a k-element co-offset set E contains an arithmetic progression of length L ≥ k−5, then a good period exists at j ≤ ⌈7(L−1)/(L−k+6)⌉ = O(k), by an ELEMENTARY argument (Dirichlet clustering of the sub-AP + a pigeonhole gap-split): cluster the L-term AP into a circular arc of span < (L−k+6)/7 so its complement is a single gap of length > (k−L+1)/7, which the remaining m=k−L ≤ 5 points can split into at most m+1 pieces, the largest > 1/7. This elementarily closes the near-AP branch (longest-AP ≥ k−5) of the general j*=O(k) capstone — including mac-mini's exact-AP case (L=k, bound = ⌈7(k−1)/6⌉) — leaving only the deeply-dissociated case (longest-AP ≤ k−6) to kps-S91; the two ranges tile all L with no gap
status: PROVED (elementary: Dirichlet's simultaneous-approximation theorem in ONE dimension + a pigeonhole on the split gap; no Weyl / equidistribution / resonance sum). VERIFIED: the Dirichlet-cluster j leaves a >1/7 gap in 100% of constructed long-AP clusters (L=k..k−5, m=0..5, k=11,12,13, Vmax=91/200/400), and j* ≤ Q in all sampled hard clusters with L ≥ k−5. Closes THM-527-A's finite check for the near-AP branch WITHOUT the 𝒲̂-resonance bound
source: klein-2026-07-09-S196
depends_on:
  - THM-527   # good period j ⟺ maxgap{frac(e_i j/Vmax)} > 1/7; the finite-Vmax glue THM-527-A
  - LEM-010   # the deterministic good-period lemma; this bounds its j* for the near-AP branch
related:
  - THM-664   # the Weyl route (the analytic alternative this sidesteps for L ≥ k−5)
  - LEM-011   # the 𝒲̂ route (needed only for the dissociated L ≤ k−6 complement, kps-S91)
external: Dirichlet's approximation theorem; the three-distance / pigeonhole.
---

# LEM-012 — The gap-splitting good-period bound (elementary, near-AP branch)

## Statement

Let `E = {e_0=0, e_1, …, e_{k−1}} ⊂ ℤ` (`k ≤ 13`) be a co-offset set, `V > max E` the ruler. A
**good period** is `j ∈ {1,…,V−1}` with `maxgap{ e·j mod V : e∈E } > V/7` (THM-527: `⟹ M(S) ≥ 1/14`).
Suppose `E` contains an **arithmetic progression of length `L`** (common difference `d`), i.e.
`{a, a+d, …, a+(L−1)d} ⊆ E`, with

> **`L ≥ k − 5`.**   Then a good period exists at some `1 ≤ j ≤ Q := ⌈7(L−1)/(L−k+6)⌉ = O(k)`.

(`Q ≤ 49` for all `k ≤ 13`; for `L = k`, `Q = ⌈7(k−1)/6⌉` — mac-mini-S59's exact-AP bound, recovered.)

## Proof (elementary)

Write `m := k − L ≤ 5` (the number of points of `E` off the AP) and `S := (L−k+6)/7 = 1 − (m+1)/7 > 0`
(the hypothesis `L ≥ k−5` is exactly `S > 0`).

**Step 1 — cluster the sub-AP (Dirichlet).** By Dirichlet's approximation theorem, for the integer
`Q = ⌈(L−1)/S⌉ = ⌈7(L−1)/(L−k+6)⌉` there is `j ∈ {1,…,Q}` with `‖j d / V‖ < 1/Q ≤ S/(L−1)` (this is
UNCONDITIONAL in `d, V`). Put `s := (jd \bmod V)` in `(−V/2, V/2]`, so `|s| < V·S/(L−1)`. At this
dilation the AP maps to `{(a+id)j \bmod V} = {aj + i·s \bmod V : i=0,…,L−1}` — an AP of step `s`
spanning `|(L−1)s| < V·S < V`, so it lies in a **single circular arc `A` of length `< S·V`** (no wrap).

**Step 2 — the complement is one big gap.** Since the `L` AP-points sit inside the arc `A` of length
`< S·V`, the complementary arc `Γ := (circle) ∖ A` is a single gap of length `> (1−S)·V = (m+1)V/7`,
containing NO AP-point.

**Step 3 — the `m` stray points cannot fill it (pigeonhole).** The remaining `m = k−L` points of `E`,
at this dilation, sit somewhere on the circle. Only those falling inside `Γ` subdivide it; there are at
most `m` of them, so they cut `Γ` into at most `m+1` sub-arcs whose lengths sum to `|Γ| > (m+1)V/7`.
Hence the **largest sub-arc `> |Γ|/(m+1) > V/7`** — a gap `> V/7`. So `maxgap{e·j mod V} > V/7`: `j` is a
good period, and `j ≤ Q`. ∎

(Coincidences only help: if stray points merge with each other or with AP-points, there are fewer
cuts and the surviving gap is only larger.)

## Verification

Machine-checked (`lrc14_gapsplit_mechanism_klein_S196`): over constructed clusters `E = AP_L ∪ {m
random points}` for every `L ∈ {k,…,k−5}` (`m = 0,…,5`), `k = 11,12,13`, `V = 91/200/400`, hard
(`spread ≥ 6V/7`, `j=1` fails) — the Dirichlet-cluster `j` leaves a `> 1/7` gap in **100%** of cases,
and `j* ≤ Q` always. Also (`lrc14_gapsplit_klein_S196`) over sampled hard clusters, `j* ≤ Q` holds
whenever `L ≥ k−5` (0 failures).

## Consequence — the near-AP branch of `j*=O(k)` is ELEMENTARY

The general `j*=O(k)` capstone (the last math gap of the covering case, THM-527-A / LEM-010) splits on
the **longest AP** `L(E)`:

- **`L ≥ k−5` (near-AP):** LEM-012 — `j* ≤ ⌈7(L−1)/(L−k+6)⌉`, **elementary** (Dirichlet + pigeonhole,
  no equidistribution). This is the *hard* branch — the near-AP is where `j*` is largest (`≈ k`), and
  where mac-mini had only the exact AP and opus needed the `r_N` (𝒲̂) verification. Now it is a
  one-page elementary proof.
- **`L ≤ k−6` (deeply dissociated):** kps-S91 — `j* ≤ 3` (verified; a-priori via the LEM-011
  𝒲̂-smallness, which holds because a short longest-AP means few resonances).

The two ranges **tile `{2,…,k}` with no gap** (`k−5` and `k−6` are consecutive). So `j*=O(k)` holds for
every cluster, closing THM-527-A's finite check to `V ≤ O(k)` (trivial) — **modulo only** the a-priori
`𝒲̂`-smallness of the dissociated branch (or an elementary argument for it). The near-AP branch — the
one everyone expected to need the analysis — is now elementary.

## Remaining (the dissociated branch)

`L ≤ k−6` (only reachable for `k ≥ 8`; `L ≤ 7` at `k=13` down to `L=2` Sidon at `k=8`). Empirically
`j* ≤ 3` (kps-S91). Two honest routes: (a) the 𝒲̂ partial-sum bound (LEM-011, the same object as the
density-floor tail) — `dissociated ⟹ few small-`n` resonances `n·E≡0` ⟹ the correction is small ⟹
`E_grid[W] ≈ (6/7)^k > 0` at small `N`); (b) an elementary gap-argument for short-longest-AP hard sets
(open). Either closes `j*=O(k)` and hence LRC(14)'s covering case outright.

## Route (c) — the dissociated branch closes by ARC-COUNT existence (mac-mini-S61)

The finite-`Vmax` glue needs only that a good period **EXISTS** (any `j ∈ {1,…,Vmax−1}`), not the
small-`j` bound `j*=O(k)`. For the dissociated branch this is cleaner than routes (a)/(b): the
**arc-count pigeonhole** (klein-S192, kps-S90)

> `#{good grid j} ≥ ρ*·Vmax − #arcs(Good_E)`,  so a good period exists whenever  `#arcs < ρ*·Vmax`.

Since `spread ≤ Vmax`, this holds iff `c := #arcs/spread < ρ*`. Verified over dissociated
spread-dense clusters (`L ≤ k−6`, `k=11,13`): **`c ≤ 0.19 (k=11) / 0.51 (k=13)`** while **`ρ* ≥
0.99 / 0.96`** — `c < ρ*` with a large margin (`≥ 0.45`). So the dissociated branch closes by
existence. Both inputs are a-priori-backed by established machinery, for LOW longest-AP:

- **`ρ* ≥ ρ_min` (bounded below):** `ρ* = μ_{1/7}(E) ≥ D3(E)` (THM-661), and `D3` is HIGH for low
  longest-AP — the decorrelation limits `D3_∞^{(L)} = 0.86/0.85/0.84/0.81/0.76/0.68/0.60` at
  `L=2..8` (opus-S158) DECREASE in `L`, so dissociated (`L ≤ k−6 ≤ 7`) gives `ρ* ≳ 0.6`.
- **`#arcs ≤ c(L)·spread` (bounded above), `c(L)` increasing in `L`, `Vmax`-INDEPENDENT** (the
  bounded-arc-count lemma, mac-mini-S58): `#arcs` depends only on the cluster's internal differences;
  for low longest-AP the coherent gap-events are few, `c(L) ≤ 0.5` for `L ≤ 7`.

The margin `c(k−6) ≤ 0.5 < ρ_min ≈ 0.6` closes the dissociated branch a-priori (non-tight), leaving
only the explicit `c(L)`/`ρ_min` constants (the same density-floor-tail + arc-count objects, both
already computed). Since the near-AP branch (`L ≥ k−5`) is LEM-012 (elementary) and the dissociated
branch (`L ≤ k−6`) is arc-count existence, **both branches of THM-527-A now have a-priori closure
routes** — LRC(14)'s covering case is closed modulo the two explicit constants + Lean. File:
`04-computation/lrc14_dissociated_arccount_macmini_S61.{py,out}`. (NB the arc-count pigeonhole was
VACUOUS for the near-AP branch — MISTAKE-127 — but is the RIGHT tool here: low longest-AP ⟹ `c` small
⟹ `c < ρ*`. The two branches use OPPOSITE tools, cf. the dichotomy below.)

> **⚠ CORRECTION (klein-S199, MISTAKE-128).** The a-priori form `#arcs/spread < D3(E)` is FALSE — and
> it fails at **ALL spreads by dilation**, not just in a finite-check window. `c=#arcs/spread`, `D3`,
> `μ` are all **dilation-invariant**. Counterexample: `E={0,7,14,21,26,29,37,44,51,58,67,75,82}` (hard,
> longest-AP=4, four co-offsets `≡0 mod 7`) has `#arcs=72` (grid-stable), spread 82 ⟹ `c=0.878`, but
> `D3(E)=0.629`, so `c/D3=1.40`. Dilating (`t·E`, `#arcs=72t`, spread `82t`) keeps `c/D3=1.40` at spread
> `≥200` — so the "spread ≥200 ⟹ c<D3" claim (mac-mini-S62) is contradicted; its "c/D3 decreasing in
> spread" was over random sets, missing this 7-structured dilated family. CAUSE: co-offset differences
> `≡0 mod 7` (step-7 sub-AP) resonate with the `1/7` threshold, spiking `#arcs`, while `D3(E)` (a moment
> LOWER bound) stays ≪ the true `μ`. **Route (c) still holds via the ACTUAL `ρ*=μ`** (`c=0.878 < μ=0.915`;
> over 4691 hard low-L sets `c<μ` always, min margin 0.081). So the a-priori closure needs `μ ≥ c` (a
> decorrelation lower bound `μ≳0.9`), NOT the D3 proxy — route (c) inherits a μ-lower-bound difficulty and
> is NOT the clean D3 closure it appeared. Existence is unaffected (LEM-010(ii) Dirichlet is the
> unconditional fallback; LEM-012 covers structured `L≥k−6`).

## Extension to L = k−6 (m=6): the ×7 collapse (klein-S197)

The gap-split needs `S = (L−k+6)/7 > 0`, i.e. `L ≥ k−5`. The boundary `L = k−6` (`m = 6`, `S = 0`) is
recovered by a companion **×7 argument**. Cluster the `L`-AP TIGHTLY: choose `Q = ⌈49(L−1)/3⌉` so the
Dirichlet `j ≤ Q` gives super-point width `δ = (L−1)·(jd mod V) < 3V/49`. Then:

- **If that `j` is good, done** (`j* ≤ Q`).
- **If it is bad** (`maxgap ≤ V/7`): the config is `1` super-point `+ 6` stray `= 7` clumps with all
  gaps `≤ V/7` summing to `V−δ`, so every gap `∈ [V/7−δ, V/7]` — the clumps are FORCED onto a `V/7`-grid,
  `q_c = p + c·V/7 + ε_c` with `|ε_c| ≤ δ`. **At dilation `7j`** the grid collapses:
  `7q_c = 7p + cV + 7ε_c ≡ 7p + 7ε_c \pmod V`, so the whole config lands in an arc of span `≤ 14δ <
  14·(3V/49) = 6V/7` — a gap `> V/7`. So `7j` is good: **`j* ≤ 7Q = 7⌈49(L−1)/3⌉ = O(k)`.**

Machine-confirmed (`lrc14_times7_klein_S197`): the perfect `V/7`-grid (`maxgap = 1/7` exactly, bad) maps
entirely to one point under `×7` (`maxgap → 1`); perturbed grids are already good directly. (For `m = 6`
the clustering-`j` is in fact good `~99%` of the time — the exact-grid bad case, which `×7` handles, is
rare.) The `×7` works ONLY at `m = 6` (`7` clumps force the grid); for `m ≥ 7` there are `≥ 8` clumps and
no grid is forced.

**Updated split (the clean dichotomy).** `L ≥ k−6` is now ELEMENTARY (`L ≥ k−5` gap-split + `L = k−6`
×7 collapse — the "concentrate" mechanism, Dirichlet only). `L ≤ k−7` (deeply dissociated / Sidon-like,
only `k ≥ 9`, longest-AP `≤ 6` at `k=13`) is the "spread" mechanism, closed by **Route (c)** (arc-count
existence, `c < ρ*` — mac-mini-S61) or **Route (a)** (𝒲̂ few-resonances, LEM-011). For `k = 8` there are
NO `L ≤ k−6` sets, so LEM-012 closes `k=8` outright. The two branches use OPPOSITE tools (concentrate a
long AP vs. exploit the spread of a short-AP set), matching opus-S166's two-mechanism framing.

### Route (c) sharpened — the clean per-cluster inequality `#arcs/spread < D3(E)` (mac-mini-S61)

Route (c) collapses to a SINGLE exact per-cluster inequality. Since `ρ* ≥ D3(E)` (THM-661, exact
degree-3 covering-moment bound) and `#arcs = c·spread ≤ c·Vmax` (`spread ≤ Vmax`), a good period
EXISTS whenever

> **`c := #arcs(E)/spread(E) < D3(E)`**   (then `#arcs ≤ c·Vmax < D3·Vmax ≤ ρ*·Vmax`, so
> `#{good grid} ≥ ρ*·Vmax − #arcs > 0`).

BOTH sides are exact and a-priori — `D3` via Farey-cell integration, `#arcs` via the arc structure
(`Vmax`-independent, S58). Verified over dissociated (`L ≤ k−6`) clusters: **`c < D3` holds for k=11
ALWAYS** (min margin `+0.58`), and for **k=13 except a narrow small-spread sliver** (`spread ≈ 80`:
`c=0.675` vs `D3=0.659`, margin `−0.016`). But small spread ⟹ the hard case has `Vmax ≤ 7·spread/6 ≈
93` ⟹ **inside the finite check** (kps-S30, `Vmax ≤ 1001`). For large spread `c` falls (toward its
asymptotic `≤ 0.5`) while `D3` rises (toward the decorrelation limit `≥ 0.6` for `L ≤ 7`, opus-S158),
so `c < D3` holds with margin. So the dissociated branch closes by **`[c < D3` for spread `≥ ~150] +
[small-spread finite check]`** — a clean, mostly-a-priori route with NO equidistribution or resonance
sum. Remaining: the explicit `c < D3` inequality over the large-spread dissociated range (both sides
exact — a finite/verifiable statement, not an analytic wall). File:
`04-computation/lrc14_dissociated_D3_vs_c_macmini_S61.{py,out}`.

### Route (c) — the concrete finite/verifiable closure (mac-mini-S62)

The inequality `c := #arcs/spread < D3(E)` is dilation-aware: `D3` is dilation-invariant, `c` shrinks
as spread grows (large spread ⟹ small `c`). Machine sweep over dissociated (`L ≤ k−6`) clusters shows
`c/D3` is **monotone decreasing in spread** and `< 1` throughout — `max c/D3` by spread `∈
{80,120,200,350,600,1000}`:

| spread | 80 | 120 | 200 | 350 | 600 | 1000 |
|---|---|---|---|---|---|---|
| k=11 (`L≤5`) | 0.38 | 0.30 | 0.16 | 0.15 | 0.12 | 0.10 |
| k=13 (`L≤7`) | 0.90 | 0.86 | 0.64 | 0.58 | 0.50 | 0.44 |

So the dissociated branch closes CONCRETELY as **`[spread ≥ 200: c < D3]` + `[spread < 200: finite
check]`**:
- **spread ≥ 200:** `c/D3 ≤ 0.64` (`≥ 35%` margin, robust; both sides EXACT — `D3` via Farey, `#arcs`
  via the arc structure) ⟹ `#arcs ≤ c·Vmax < D3·Vmax ≤ ρ*·Vmax` ⟹ a good period exists.
- **spread < 200:** the hard case (`j=1` fails) has `Vmax ≤ 7·spread/6 < 234` — a bounded FINITE
  CHECK, well inside kps-S30's exact `M(S) ≥ 1/14` sweep (`Vmax ≤ 1001`).

This is the finite/verifiable form the covering case reduces to: an exact per-cluster inequality
`c < D3` over the (dilation-normalized, hence effectively finite) large-spread dissociated shape
space — verified with `≥ 35%` margin and monotone in spread toward the decorrelation limits
(`c → c_∞ ≤ 0.44`, `D3 → D3_∞ ≥ 0.76` for `k=13`, `L ≤ 7`, opus-S158) — plus a `Vmax ≤ 234` finite
check. No equidistribution, no resonance sum, no analytic wall. File:
`04-computation/lrc14_dissociated_threshold_macmini_S62.{py,out}`.
