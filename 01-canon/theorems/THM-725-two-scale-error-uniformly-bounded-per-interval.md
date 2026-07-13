---
id: THM-725
title: The two-scale far-element correction is UNIFORMLY BOUNDED (not merely O(Σe'/w)) — via the per-interval reduction (the missed sector is constant on each maximal miss-exactly-one interval) and the MIN of the antiderivative bound and the trivial length bound; |Error| ≤ min(R_ct·6/(49w), 6/7), unifying THM-700's decay branch with an absolute cap
status: PROVED (elementary; the per-interval reduction is topological, the min-bound is two one-line estimates + Σ|I|≤1). SCOPE — HONEST: this proves BOUNDEDNESS and recovers THM-700's O(Σe'/w) as the large-w branch, but it does NOT provide a Σe'-free 1/w DECAY constant (the min caps at O(1) for w<Σe', it does not shrink). The Σe'-free decay (needed to close the density-row tails at moderate w) remains OPEN.
source: klein-2026-07-13-S274
depends_on:
  - THM-700   # the plateau decorrelation lemma (kps) — this adds the absolute cap + the per-interval reduction
  - THM-699   # the TV far-element contraction (mac-mini) — same object, weight side
related:
  - HYP-6285  # klein-S273 — the cross-sector cancellation lemma (this session works it)
  - HYP-6305  # klein-S274 — R_ct interval-count + the min-argument
---

# THM-725 — the two-scale far-element correction is uniformly bounded (per-interval reduction)

## Setup

For co-offsets `E'` and a far offset `w`, the cover-measure two-scale error (THM-700) is
`Error = Σ_{s=0}^6 ∫₀¹ f_s(x) g_s(wx) dx`, with `f_s = 1{E'` misses exactly sector `s}` and
`g_s(y) = 1{y∈[s/7,(s+1)/7)} − 1/7` (mean zero). THM-700 bounds this by `V(E')/(6w) ∝ Σe'/w`.

## The per-interval reduction

The `f_s` are disjoint (`E'` misses at most one sector "exactly"). Let `R = {x : E'` misses exactly
one sector`}`, a union of maximal intervals.

> **Lemma (missed sector constant per interval).** On each maximal interval `I ⊆ R`, the missed
> sector `s_I` is constant.

*Proof.* To change the missed sector from `s` to `s'` on a connected piece of `R`, sector `s` must
gain a phase (`0→1`) and `s'` must lose its last (`1→0`). Occupancy changes one boundary-crossing at
a time (generic offsets ⇒ no simultaneous crossings), so between the two events the miss-count is `0`
(if `s` fills first) or `2` (if `s'` empties first) — either way `x` leaves `R`, contradicting that `I`
is a single component of `R`. ∎

Hence `Error = Σ_{I⊆R} ∫_I g_{s_I}(wx)\,dx`, a **single fixed sector** per interval.

## The min-bound (boundedness)

Bound each interval term two ways:
- **antiderivative:** `|∫_I g_{s_I}(wx)dx| = |G_{s_I}(wb_I) − G_{s_I}(wa_I)|/w ≤ osc(G_s)/w = (6/49)/w`;
- **trivial:** `|∫_I g_{s_I}(wx)dx| ≤ |I|·‖g_s‖_∞ = (6/7)|I|`.

Take the min. With `R_ct = #{I}` and `Σ_I |I| = |R| ≤ 1`, and using `Σ_I min(a_I,b_I) ≤ min(Σa_I, Σb_I)`,

> **`|Error| ≤ min( R_ct · (6/49)/w , (6/7) )`.**

- **Large `w` branch** (`w ≳ R_ct`): `≤ R_ct·(6/49)/w`, i.e. THM-700's `O(Σe'/w)` decay (measured
  `R_ct ≈ 0.81·Σe'`, so this is `≈ 0.099·Σe'/w`, ~24× tighter than THM-700's crude `14Σe'/(6w)`).
- **Small `w` branch** (`w ≲ R_ct`): `≤ 6/7` — an **absolute cap**, independent of `E'` and `w`.

For the degree-3 majorant `Φ = E_x[ψ(N)]` the same argument, weighted by `Δψ(N)=ψ(N)−ψ(N−1)`
(`max_N|Δψ| = |Δψ(1)| = 2/3`), gives `|Error_Φ| ≤ (2/3)·min(R_ct·(6/49)/w, 6/7)`.

## What this does and does NOT do

- **DOES:** unifies THM-700 (decay) with an absolute cap `6/7`; the min shows THM-700's `Σe'/w` is
  only meaningful for `w ≳ Σe'` (beyond that the error is trivially `≤ 6/7`). The per-interval
  reduction is a clean structural sharpening.
- **Does NOT:** give a `Σe'`-free **decay** `|Error| ≤ C_abs/w`. The min-bound is `O(1)` for `w<Σe'`
  (it caps, it does not shrink). Closing the density-row tails at *moderate* separation (`w` a few
  times the cluster diameter, not `≫ Σe'`) needs the interval terms to **cancel** to `√R_ct` or
  better — which the min-argument does not supply. That decay constant is the open lemma (HYP-6285).

## Empirical (klein-S274, PRIME grid Ng≫w — the `Ng∝w` grid aliases and must not be used)

`err·w` (the true decay constant `C_abs`): **~0.2 at clean `w`, uniform in `Σe'`** (Σe'=21→420);
**~3 at resonant `w = lcm(E')`** (2-block clusters). The resonance is harmless for the row: `w=lcm ≥
diam`, so `w ≫ Σe'` and THM-700's `Σe'/w` branch already gives a small error there.

Files: `04-computation/lrc14_Rct_cross_sector_klein_S274.py`, `lrc14_prime_grid_klein_S274.py`,
`lrc14_Phi_min_constant_klein_S274.py` (+ outs).
