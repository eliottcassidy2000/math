---
id: THM-531
title: The LRC(14) S3 residual reduces to ONE clean three-distance extremal statement — "consecutive minimizes μ_{1/7}" — with the exact 1/7 union-bound closure PROVED (m_k+μ_{1/7}(consec k)−1>0 for all k=7..13, exact rationals) and the AP-orbit invariance PROVED (μ_θ is translation- and scale-invariant, so every arithmetic progression has μ_θ=μ_θ(consecutive), killing the large-spread single-run tail)
status: MIXED. PROVED: the exact union-bound closure (finite rational computation, all k=7..13); the AP-orbit invariance (translation+scale, elementary). VERIFIED (exact/exhaustive): consecutive minimizes μ_{1/7} at k=8 (792 shapes), k=9,10,11 (bounded spread); all large-spread/multi-block/dissociated ≥ consecutive. OPEN (the isolated crux, now LRC-free): the extremal statement μ_{1/7}(E) ≥ μ_{1/7}(consecutive_k) for every non-AP integer set E. LRC(14) NOT proved.
source: mac-mini-2026-06-18-S5
depends_on:
  - THM-530   # the two-branch global-witness floor (k≤7 pigeonhole; k≥8 union bound)
  - THM-527   # the lonely-density reformulation (good period ⟺ x∈G_P and cluster gap)
  - THM-523   # covering-set reduction; meas(G_P) floor
related:
  - OPEN-Q-108
  - HYP-2592   # the 1/7 threshold slack ladder (codex)
  - HYP-2602   # consecutive minimizes μ_{1/7} (the reduced crux)
external: Lonely Runner Conjecture (first open case = 13 speeds). Steinhaus three-gap theorem; Weyl equidistribution.
---

# THM-531 — The 1/7 union closure, AP-orbit invariance, and the reduced crux

Building on THM-530 (the correct object is the **global-witness density**
`ρ*_{1/7}(P,E) = meas(G_P ∩ {x: maxgap{frac(e_i x)} > 1/7})`, with the threshold `1/7` —
a single safe **point** needs the cluster teeth, each of half-width `1/14`, to leave a gap
`> 2·(1/14) = 1/7` between adjacent centers — not the via-max `2/7`). THM-530's two-branch
floor: `k≤7` PROVED unconditional (pigeonhole: `7` gaps can't all be `≤1/7`, so
`μ_{1/7}=1`); `k≥8` the union bound `ρ*_{1/7}(P,E) ≥ meas(G_P) + μ_{1/7}(E) − 1`.

## A. The exact 1/7 union-bound closure (PROVED — finite rational computation)

The per-`k` minimum of the union bound, `m_k + μ_{1/7}(consec_k) − 1` where
`m_k := min_{|P|=13−k} meas(G_P)` and `μ_{1/7}(consec_k) = min_E μ_{1/7}(E)` (consecutive
is the verified minimiser, part C), is a **positive exact rational for every `k=7..13`**:

| k | μ_{1/7}(consec k) | m_k = min meas(G_P) | margin = m_k+μ−1 |
|---|---|---|---|
| 7 | 1 | 3029/10780 | **3029/10780 ≈ +0.281** |
| 8 | 691/735 | 2243/5880 | **1891/5880 ≈ +0.322** |
| 9 | 247/294 | 1979/4004 | **28117/84084 ≈ +0.334** |
| 10 | 38/49 | 55/91 | **242/637 ≈ +0.380** |
| 11 | 1381/2205 | 66/91 | **10078/28665 ≈ +0.352** |
| 12 | 13823/24255 | 6/7 | **10358/24255 ≈ +0.427** |
| 13 | 477/1078 | 1 | **477/1078 ≈ +0.443** |

So `ρ*_{1/7}(P,E) > 0` — hence `M(S) ≥ 1/14` for the reconstructed covering set `S` — for
**every** admissible `(P,E)`, *provided* `μ_{1/7}(E) ≥ μ_{1/7}(consec_k)`. The entire S3
residual is now reduced to that one inequality (the margin is the slack
`μ_{1/7}(consec_k) − thr_k`, `thr_k = 1 − m_k`, all `≥ 0.28`).

## B. AP-orbit invariance (PROVED — kills the large-spread single-run tail)

> **Proposition.** For any threshold `θ` and any integers `a, d≥1`, the arithmetic
> progression `E = {a, a+d, …, a+(k−1)d}` satisfies `μ_θ(E) = μ_θ({0,1,…,k−1})`.

**Proof.** `μ_θ(E) = meas{x∈[0,1): maxgap{frac((a+jd)x) : j=0..k−1} > θ}`.
*(Translation.)* `frac((a+jd)x) = frac(frac(ax) + frac(jd·x))`, i.e. the `k` points
`{frac(jd·x)}` all **rotated** by the common angle `frac(ax)`. A rigid rotation preserves
every circular gap, so `maxgap{frac((a+jd)x)} = maxgap{frac(jd·x)}` for all `x`; hence
`μ_θ(E) = meas{x: maxgap{frac(j·(dx))} > θ}`. *(Scale.)* The map `x ↦ frac(dx)` is
measure-preserving (`d`-to-`1`) on `[0,1)`, so for the measurable set
`F = {y: maxgap{frac(jy):j} > θ}`, `meas{x: frac(dx)∈F} = meas(F)`. Thus
`μ_θ(E) = meas(F) = μ_θ({0,…,k−1})`. ∎

**Consequence.** Every single arithmetic run — at *any* spread `d` — has exactly the
consecutive value. This rigorously disposes of the "S3 is infinite" obstruction
(THM-526/HYP-2581c: the dilated-AP family `{t,2t,…,12t,V}` with `t→∞`): such clusters are
APs, so their `μ_{1/7}` is pinned to `μ_{1/7}(consec_k)`, independent of the unbounded
spread. The non-compactness in `Vmax` is illusory at the level of `μ_{1/7}`.

## C. The reduced crux (VERIFIED, the isolated open statement)

> **CRUX (HYP-2602).** `μ_{1/7}(E) ≥ μ_{1/7}({0,1,…,k−1})` for every integer set `E`,
> `0∈E`, `|E|=k` — i.e. **the consecutive cluster (equivalently any AP, by part B)
> minimises the `1/7`-gap measure.**

VERIFIED: exhaustive over all shapes at `k=8` (792, spread ≤12; 0 violations), `k=9,10,11`
(bounded spread; 0 violations); every large-spread, multi-block, perforated, and
dissociated/Sidon shape tested has `μ_{1/7} ≥ μ_{1/7}(consec_k)` (often `=1`). This is a
**clean, LRC-free extremal statement about three-distance orbits**: among all `k`-point
integer dilation-orbits `{frac(e_i x)}`, the AP `{0,x,…,(k−1)x}` (the Steinhaus orbit) is
the most uniform and hence leaves the smallest `1/7`-gap measure. Contrast the `2/7` object,
where consecutive is NOT extremal (perforation wins, THM-527-F / HYP-2585) — the `1/7`
threshold is exactly where the AP becomes the clean minimiser.

## D. The remaining tail and why it is soft

Non-AP shapes split into bounded-spread (a finite exact check per `k`, done for `k≤11`) and
large-spread multi-block. As blocks separate, `frac(g·x)` decorrelates and `μ_{1/7}(E)`
tends to a block-average that the data keeps `≥ μ_{1/7}(consec_k)` (the separated blocks add
points, only shrinking gaps relative to a single block of the same size; the per-`k`
minimum is the maximally-coupled single run = AP). The only non-elementary ingredient is a
quantitative rate for this convergence — but the union slack (`≥ 0.28`, part A) means a
*crude* rate suffices, unlike the sharp floor the `2/7` object demanded.

## Honest status

PROVED: the exact `1/7` union closure (A) and the AP-orbit invariance (B). The S3 residual
is reduced to the single **LRC-free** extremal inequality (C), verified exhaustively at the
relevant `k`, with the only soft spot a crude multi-block equidistribution rate (D) that the
`0.28` slack makes far easier than the original floor. **LRC(14) is NOT proved**, but its
last inequality is now a clean three-distance statement — "the arithmetic progression
minimises the `1/7`-gap measure" — divorced from the runner problem.
