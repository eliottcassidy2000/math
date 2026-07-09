---
id: THM-667
title: The clamped grid-port — the full density-floor measure μ(θ′) transfers to the Vmax-ruler grid at the 1/V² aliasing rate; grid-fraction{j : maxgap(j/V) ≥ θ″} ≥ μ(θ′) − TV(C′)/(12V²) for the clamp C of (maxgap − θ″)/(θ′ − θ″)
status: PROVED (reduction to THM-665's three steps applied to the clamp; below). Machine-verified: grid fraction ≥ bound on the zoo at V = 400/1600/6400 (see companion .out); TV(C′) exact per shape via the extended cell engine (maxgap = max of affines, subdivided at argmax and band crossings).
source: monad-explorer-2026-07-09-S2 (HYP-5717) — the upgrade THM-665's Markov-fraction weakness demanded: the (7/6)·∫W ≈ 0.15 grid fraction becomes the full closed-leg μ ≈ 0.3–0.94.
depends_on:
  - THM-665   # the aliasing identity + BV coefficient bound the proof reuses verbatim
related:
  - THM-661   # the closed density-floor legs whose μ this ports to the grid
  - LRCDriftEmbed (klein-S205)  # the consumer: teeth-gap at a grid index j ⟹ minReach ≥ 1/14
external: none beyond THM-665's.
---

# THM-667 — the clamped grid-port

## Statement

Let `E` be a cluster of co-offsets, `maxgap(x)` the largest circular gap of
`{frac(e·x) : e ∈ E}` (continuous, piecewise linear in `x`), and fix thresholds
`1/7 ≤ θ″ < θ′`. Define the **clamp**

> `C(x) = clamp( (maxgap(x) − θ″)/(θ′ − θ″), 0, 1 )`,

a continuous piecewise-linear 1-periodic function with

> `1[maxgap(x) > θ′] ≤ C(x) ≤ 1[maxgap(x) ≥ θ″]`.

Then for every integer `V ≥ 1`:

> **`(1/V)·#{ 0 ≤ j < V : maxgap(j/V) ≥ θ″ } ≥ μ(θ′) − TV(C′)/(12 V²)`,**

where `μ(θ′) = meas{x : maxgap(x) > θ′}` is the density-floor measure and `TV(C′)` is
the total variation of `C′` over the circle (exactly computable: `C′` is piecewise
constant with values `(slope of the active gap)/(θ′ − θ″)` inside the band and `0`
outside).

## Proof

By the sandwich, `#{j : maxgap(j/V) ≥ θ″}/V ≥ E_grid[C](V)` and `∫C ≥ μ(θ′)`.
`C` is continuous piecewise linear, so THM-665 (i)–(iii) apply verbatim with `C` in
place of `W`: `|E_grid[C](V) − ∫C| ≤ TV(C′)/(12V²)`. Chain the three inequalities. ∎

## Remarks

1. **What it buys.** THM-665 alone gives a grid fraction only via Markov
   (`E_grid[W]/Wmax ≈ 0.15`). The clamp ports the FULL measure — at the closed legs'
   values (μ ≥ 0.309 at k=13 up to 0.94 at near-AP clusters) — to the ruler grid.
   Combined with klein-S205's drift embed (which consumes a teeth-gap at a single grid
   index), this is the "existence half" of the P∪L composition at full strength.
2. **The cost dial.** `TV(C′) ≈ TV-of-maxgap-slope-in-band/(θ′ − θ″)`: the price of a
   narrow clamp band is `1/(θ′ − θ″)` in the aliasing term. The band must also cover
   the drift margin (`θ″ ≥ 1/7 + drift`), so the dial trades grid-error against the
   drift-embed's gap requirement — the optimization is per-shape and exactly computable.
3. **Same corner-cancellation room as THM-665:** the measured grid error again sits far
   below the bound (square-root over the band's corners); the unified Kronecker
   statement (HYP-5707 finding 4) would sharpen both simultaneously.
