# The 11-runner even-part M-spectrum has a gap above 1/12 — the residual "near-AP" is discrete, not a continuum

*klein-2026-07-04-S126 (HYP-4080). Owner: work the single remaining sliver (the small-tightener ×
near-AP corner of the m=2 confinement). I probed the actual crux question — is the residual a
continuum accumulating at the tight value, or discrete? Full-enumeration answer: DISCRETE. This is
the M-value companion of kps's measure-quantization ([[the-lrc-inf-is-a-quantized-gap-from-the-tight-locus-and-the-extremizer-is-a-sporadic-perturbation-kps]]); an ingredient that sharpens the residual, not a closure.*

## The finding (exact, full-enumeration verified)

The confinement's last sliver is `M(2U ∪ {w₁,w₂}) ≥ 1/12` for **loose** 11-runner even parts `U`
(`M(U) > 1/12`) with small tighteners. The fear was a *continuum* of `U` with `M(U) → 1/12`, making
the good band `{g_E ≥ 1/12}` shrink to nothing (so any small tightener could hide it). **There is
no such continuum:**

> **`M(U) = 1/12` is isolated in the 11-runner spectrum.** Over ALL primitive 11-subsets of `{1..17}`
> (12 376 sets) plus dilated/large-speed near-AP families, **zero** have `M(U) ∈ (1/12, 2/23)`.
> The smallest loose value is `2/23`; the gap is **`δ₀ = 2/23 − 1/12 = 1/276 > 0`**.

The near-bottom spectrum is *exactly* the covering-min / Ostrowski ladder rungs (mac-mini S38)
`M = k/(11k+1)` — `1/12, 2/23, 3/34, …` accumulating at `1/11` — realized cleanly by
**`{1,…,10, 11k}`** (spread the top runner to `11k` gives rung `k`; `k=1→1/12` tight, `k=2→2/23`,
`k=3→3/34`). Nothing lives *between* the rungs at the bottom. This is the Lagrange-spectrum-is-
discrete-at-the-bottom phenomenon, for LRC(12).

## Why it matters for the sliver (honest: an ingredient)

Loose even parts have `M(U) ≥ 2/23`, a **definite margin** `M(U) − 1/12 ≥ 1/276` — not an
infinitesimal. Consequences for the residual:

- **The near-1/12 continuum is gone.** "Near-AP loose `U`" is not `M(U) = 1/12 + ε`; it is a
  discrete rung `≥ 2/23`. The residual `U` are ladder-type (spread-one-runner) configurations, the
  same shape kps closed by an exact residue-table formula (kps-S2, `{1..11,13,12k}`).
- **opus's Lemma 3 (large tightener, klein-S125 Lean) now has a uniform trigger:** it fires for
  `max(w₁,w₂) > u_max/(6·(M(U)−1/12)) ≤ u_max/(6·(1/276)) = 46·u_max`. So the residual is confined
  to `tighteners ≤ 46·u_max` — the margin `δ₀` makes the "large tightener" threshold an explicit
  multiple of `u_max`, no longer blowing up.

So after the gap, the residual is: *ladder-type loose `U` (`M(U) ≥ 2/23`) with tighteners
`≤ 46·u_max`, forming a tight covering family.* The near-AP continuum — the thing that made the
band vacuously small — is removed.

## What it does NOT do (the crux stays)

`u_max` is still unbounded, so `46·u_max` does not bound the tighteners absolutely; the residual is
sharpened, not closed. Bounding `u_max` (= kps's "extremizers have bounded lcm" = mac-mini's "bound
`v_max(U)`") is the LRC(14)-equivalent crux, unchanged. The gap is the **value-side** companion of
kps's **measure-side** quantization (`L ≥ 1/(14·lcm)`): both say the extremal data is quantized/
discrete off the tight locus, and both reduce to the same bounded-lcm/`u_max` wall.

## The path this opens

Because the residual `U` are the discrete ladder rungs `{1..10, 11k}` (spread-one-runner) rather
than a continuum, and kps closed exactly such a family (`{1..11,13,12k}`) by an *exact residue-table
formula* (no search, all `k`), the residual may be a **finite union of formula-closable ladder
families** — one closed form per rung shape, each valid for all `k` (all `u_max`), bypassing the
`u_max` bound. That is the concrete next move: enumerate the ladder-rung even-part shapes in the
residual and close each by kps's residue-table method, as an infinite family.

## Links

- Scripts: `04-computation/lrc14_residual_spectral_gap_klein_S126.py`,
  `lrc14_residual_gap_verify_klein_S126b.py` (+ `.out`). HYP-4080.
- Credits/threads: kps quantized-gap-from-tight-locus (measure side) + kps-S2 residue-liar formula;
  mac-mini S38 Ostrowski ladder `[0;11,k]` + S40 Chebyshev-equioscillation; opus THM-615 Lemma 3
  (large tightener) = klein LRCLargeTightener/LRCFolding/LRCMarginMeasure. The `u_max`/lcm bound is
  the unchanged open crux. See also [[fibonacci-is-the-covering-mins-foil-not-its-lever-the-anti-golden-eisenstein-sibling-klein-S124]] (the anti-golden ladder is this same ladder).
