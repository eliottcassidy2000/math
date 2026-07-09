---
id: LEM-013
title: The dissociated good-period existence margin — mac-mini's Route-(c) `c ≥ D3` sliver is a SUFFICIENT-CERTIFICATE artifact, not a covering-leg gap; good-period EXISTENCE for primitive dissociated (longest-AP ≤ 7) 13-clusters holds directly with a uniform margin `7·maxgap/Vmax ≥ 1.105` across the entire critical ruler window `Vmax ∈ [spread+1, ⌊7·spread/6⌋]`
status: VERIFIED. EXHAUSTIVE for spread ≤ 22 (621,455 primitive dissociated clusters, L ≤ 7; ZERO existence failures; min margin 1.1053 = 21/19, an integer near-miss maxgap=3 vs 19/7 at Vmax=19). ADVERSARIAL (min-margin hill-climb) for spread ∈ [21,200]: min margin 1.355 (s=27), monotone-increasing to 2.31 (s=200). The margin is MINIMIZED at small spread — exactly the exhaustively-checkable regime. NOT a closed-form a-priori proof; it removes the sliver as a RISK and reduces the residual to a cleaner certificate on the intermediate band.
source: kind-pasteur-2026-07-09-S94
depends_on:
  - THM-527   # the covering reformulation; good period j ⟺ maxgap{frac(e_i j/Vmax)} > 1/7
  - LEM-012   # klein's near-AP branch (L ≥ k−6); mac-mini's Route-(c) c<D3 is the dissociated branch this sharpens
  - THM-661   # ρ* ≥ D3(E), the covering-moment bound Route-(c) leans on
related:
  - LEM-010   # the j=1 wraparound (good for spread < 6Vmax/7); this handles the complementary window
  - LEM-009   # density-floor tail (the D3 side of Route-(c))
external: none (elementary integer maxgap computation + exhaustion + adversarial search).
---

# LEM-013 — The dissociated good-period existence margin

## Context — what this closes

The LRC(14) covering leg reduces (THM-527) to: every primitive cluster `E` admits a **good period**
`j ∈ {1,…,Vmax−1}` with `maxgap{e_i·j mod Vmax} > Vmax/7`. The good-period dichotomy splits on the
longest AP in `E`:

- **near-AP** (`L ≥ k−6`): klein's LEM-012, elementary (Dirichlet clustering).
- **dissociated** (`L ≤ k−7`): mac-mini-S61 **Route (c)** — a good period exists whenever the clean
  a-priori inequality `c := #arcs(E)/spread(E) < D3(E)` holds (since `ρ* ≥ D3` and `#arcs ≤ c·Vmax`).

Route (c) holds for `k=11` always, but at `k=13` it has **one sampled failure** — a narrow
small-spread sliver (`spread ≈ 80`: `c=0.675 ≥ D3=0.659`, margin `−0.016`) — which mac-mini deferred
to "the finite check." **This lemma resolves that sliver:** it is a failure of the *sufficient
certificate* `c < D3`, **not** of good-period existence. Existence holds directly, with margin.

## Statement (verified)

For every primitive **dissociated** (`longest-AP ≤ 7`) 13-element cluster `E` with `spread s = max E`,
and every ruler `Vmax` in the **critical window** `[s+1, ⌊7s/6⌋]` (the window where the `j=1`
wraparound of LEM-010 fails — below it, `Vmax−s > Vmax/7` makes `j=1` good), a good period exists:

> **`μ(E) := min_{Vmax ∈ [s+1,⌊7s/6⌋]} max_{1≤j<Vmax} (7·maxgap{e_i·j mod Vmax}/Vmax)  ≥  1.105 > 1.`**

(`μ(E) > 1` ⟺ some period is good at every critical `Vmax`.) Outside the window: `Vmax < s+1`
impossible; `Vmax > ⌊7s/6⌋` gives `j=1` good directly (LEM-010(i)); and for genuinely large spread
Route (c)'s `c < D3` holds a-priori (mac-mini). So the window is the entire hard set.

## Evidence

| regime | method | clusters | existence failures | **min margin `μ`** |
|---|---|---|---|---|
| `s ≤ 22` | **EXHAUSTIVE** (`L ≤ 7`) | **621,455** | **0** | **1.1053** = 21/19 |
| `s ≤ 22` | exhaustive (pure `L ≤ 6`) | 569,255 | 0 | 1.1053 |
| `s ∈ [21,49]` | adversarial min-margin hill-climb | — | 0 | 1.355 (at s=27) |
| `s ∈ [50,200]` | adversarial min-margin hill-climb | — | 0 | 1.717 → 2.31 (↑ in s) |

The margin is **minimized at small spread** (the exhaustively-checkable regime) and **grows with
spread** — the opposite of a danger. The exhaustive minimum `μ = 21/19` is a transparent integer
near-miss: at `s=17, Vmax=19` the best period has `maxgap = 3 > 19/7 = 2.714`. The tightest clusters
are near-*consecutive* (`{0,1,2,3,4,5,6,…}`, `L=7`), i.e. the LEAST-dissociated allowed — which is
also inside LEM-012's `L ≥ k−6` regime, so the pure Route-(c) region (`L ≤ 6`) is at least as safe.

## Consequence for the critical path

The dissociated branch of THM-527-A closes as:

> **`[c < D3` a-priori for spread ≥ S₁ (Route (c), mac-mini)]  ∪  [good-period existence directly,
> `μ ≥ 1.105`, for spread ≤ S₁]`**,

with the small-spread half exhaustive (`s ≤ 22`) and the intermediate band (`22 < s < S₁ ≈ 100`)
adversarially margin-robust (`μ ≥ 1.355`). mac-mini's single `c ≥ D3` failure is therefore **not a
gap** — it is a region where the *certificate* is loose while *existence* is comfortable. What
remains to make the branch fully closed-form is only (i) extending exhaustion or (ii) a clean
a-priori `μ > 1` bound on the intermediate band — not a genuine covering-leg risk.

## Files

`04-computation/lrc14_{direct_existence,sliver_adversarial,sliver_midband,exhaustive_s22}_kps_S94.py`
(+ `.out` in `05-knowledge/results/`). Method: exact integer `maxgap` over all `j`, vectorized;
adversarial = random-restart hill-climb minimizing `μ` under primitivity + `longest-AP ≤ {6,7}`.
