---
id: THM-669
title: The reflection availability lemma — the φ-interval composition's low-anchor availability dominates the smooth surrogate at the lifted threshold, Av(E, r) ≥ ∫W_{1/7+r}(E)/(1−r²); with the layer cake and the parametric tent this gives unconditional explicit availability floors (|L| ≤ 6 needs no structure at all; k_L ≤ ~10 explicit via the tent; 11–13 reduce to the parametric-θ moment floors)
status: PROVED (three short steps, below). Machine-verified: exact Av (cell engine) vs the bound across the cluster zoo and an r-grid, 0 violations; the parametric tent floor re-derivation checked against THM-651 at θ = 1/7; end-to-end proven-availability vs the S2/S3 composition batteries (companion .out).
source: monad-explorer-2026-07-09-S4 (HYP-5737) — the last named analytic piece of the φ-interval composition (HYP-5717).
depends_on:
  - THM-667   # the clamped grid port that transfers the continuum availability to the ruler grid
related:
  - THM-651 / THM-661   # the tent/moment density floors this plugs into (parametric-θ versions)
  - HYP-5717 / THM-668  # the composition and the dispatch this completes
  - kps-S68 anchored-gap subset lemma  # a different anchoring (fixed anchor point vs position-weighted gaps) — complementary, cited for scope
---

# THM-669 — the reflection availability lemma

## Setting

Cluster `E` (co-offsets, `0 ∈ E`, `|E| = n`), teeth `{frac(e·x)}` at slow time `x`; gaps
`(a, a+g)` (the anchor `0` is always a tooth, so gaps do not wrap). The φ-interval
composition (HYP-5717) admits, per gap, valid fast phases of length
`(g(1−r) − 2ra − 1/7)₊/(1−r²)`, `r = S_L/Vmax`. Define the **availability functional**

> `Av(E, r) = (1−r²)^{-1} ∫₀¹ Σ_gaps ( g(1−r) − 2r·a − 1/7 )₊ dx`.

`Av(E, 0) = ∫W` (the smooth surrogate) exactly; for `r > 0` low-anchored gaps are worth
more — the composition's low-anchor structure.

## Statement

**(i) (Reflection reduction.)**  For every cluster `E` and `0 ≤ r < 1`:

> **`Av(E, r) ≥ (1−r²)^{-1} · ∫₀¹ Σ_gaps ( g − 1/7 − r )₊ dx = ∫W_{1/7+r}(E)/(1−r²)`.**

The availability dominates the smooth surrogate at the **lifted threshold** `1/7 + r`.

**(ii) (Layer cake.)**  `∫W_θ′(E) ≥ ∫ (maxgap − θ′)₊ dx = ∫_{θ′}^{1} μ_s(E) ds`, so any
per-threshold density floor integrates to an availability floor.

**(iii) (Parametric tent floor; THM-651's proof at general θ.)**  For `n` teeth and
`1/n < θ ≤ 2/n`:  `μ_θ(E) ≥ 1 − 2(n−1)(nθ − 1)/n`  (kink `β = 2/n − θ`; at `θ = 1/7`,
`n = 8` this is THM-651's `3/4`). Together with (i)–(ii) this yields **unconditional,
structure-free availability floors**, e.g.:

- `|E| ≤ 6`: pointwise `maxgap ≥ 1/|E|` gives `Av ≥ (1/|E| − 1/7 − r)₊/(1−r²)` — positive
  for every `r < 1/|E| − 1/7`, with no hypotheses whatsoever;
- `|E| = n ≤ 2/(1/7 + r)`: `Av ≥ (1−r²)^{-1} ∫_{1/7+r}^{2/n} (1 − 2(n−1)(ns−1)/n)₊ ds`,
  an explicit rational-in-`(n, r)` positive floor;
- `n ∈ {11, 12, 13}`: (i)–(ii) reduce the availability to the parametric-θ moment floors
  (THM-661's D3/B4 machinery re-run at `θ′ = 1/7 + r`) — flagged as the remaining input,
  machinery existing, constants per (n, θ′).

With THM-667 (the clamped grid port), every floor transfers to the Vmax-ruler grid at
cost `TV/(12V²)`, making the composition's mass criterion checkable from proved
quantities alone.

## Proof

**(i).** The map `x ↦ −x` is measure-preserving and reflects each configuration about
the anchor: `frac(e(−x)) = 1 − frac(ex)` (for nonzero phases), so the gap multiset of
`config(−x)` is `{(1−a−g, g) : (a, g) ∈ config(x)}` — same lengths, reflected positions.
Pair the integrand at `x` and `−x` gap-by-gap and use `(u)₊ + (v)₊ ≥ (u+v)₊`:

`(g(1−r) − 2ra − 1/7)₊ + (g(1−r) − 2r(1−a−g) − 1/7)₊ ≥ (2g(1−r) − 2r(1−g) − 2/7)₊
= 2·(g − 1/7 − r)₊`,

since `2g(1−r) − 2r + 2rg = 2g − 2r`. Integrate and divide by 2:
`Av(E,r) ≥ (1−r²)^{-1} ∫ Σ (g − 1/7 − r)₊`. ∎

**(ii).** `Σ_gaps (g − θ′)₊ ≥ (maxgap − θ′)₊` pointwise; `∫(maxgap − θ′)₊ dx =
∫_{θ′}^{∞} meas{maxgap > s} ds` (Fubini/layer cake). ∎

**(iii).** THM-651's three steps with `1/7` replaced by `θ` and tent kink `β`:
`E[F] = n(n−1)(θ−β)²/2` (pair equidistribution), `min_S F = 1 − nβ` on the safe event
(gap-sum budget), Markov. Optimize `β`: `2(1−nβ) = n(θ−β)` gives `β = 2/n − θ`
(admissible for `θ ≤ 2/n`), floor `1 − 2(n−1)(nθ−1)/n`. ∎

## Remarks

1. **What is genuinely new** is (i): the position-weight (the composition's low-anchor
   discount) is exchanged for a threshold lift by one reflection — no anchored-measure
   machinery needed. kps-S68's anchored-gap subset lemma anchors at a fixed point
   (gap@0); here the weight is linear in the gap's position and the symmetrization is
   exact. The factor `(1−r²)^{-1} ≥ 1` is kept (it only helps).
2. **The r-dial.** The composition wants small `r` (larger availability threshold room)
   but small `r` shrinks `L` (larger `|P|`, more G_P mass to beat). The assembled
   criterion — `Av`-floor(n, r) vs `|P|/7` + grid cost — is now an explicit finite
   optimization per shape class, with every ingredient proved except the `n ≥ 11`
   parametric moment floors (existing machinery, to be re-run at `θ′`).
3. **Honest gap to the observed 100%:** the proven floors are conservative (reflection
   + layer cake + tent each lose a constant); the companion computation quantifies
   proven-vs-empirical availability per battery cluster. The composition's empirical
   coverage is far above the proven floor — as expected, and safe.
