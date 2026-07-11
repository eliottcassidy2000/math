---
id: THM-700
title: "The plateau decorrelation lemma — the wide-spread wall's inductive step, PROVED: for E = E' ∪ {w} with integer offsets and w = max E, |meas(S7(E)) − Plat(E')| ≤ V(E')/(6w) where Plat(E') = p_0(E') + (1/7)p_1(E') and V(E') is the total variation of the 'E' misses exactly one sector' indicators; via the EXACT cover decomposition 1_cover(E) = 1_cover(E') + Σ_s f_s·1_{frac(wx)∈s} and Fourier decay Σ_ℓ |f̂_s(−ℓw)||ĝ_s(ℓ)| ≤ V(f_s)/(6w). This is the 1-D Weyl estimate HYP-2644 reduced the unbounded direction to; it converts the S_c reciprocal-sum wall from conditionally-convergent lattice geometry into an elementary bounded-variation Fourier bound"
status: PROVED (the decorrelation bound, rigorous — exact decomposition + BV/mean-zero Fourier estimate; verified O(1/w) with the bound holding on consec_8 core, w up to 1601). SCOPE: this is HYP-2644's inductive step (the single open analytic lemma of the wide-spread half); the FULL recursion closure (error accumulation across peels + the p_1 decorrelation + the sharp constant) remains, and the tight margin still lives in the finite check.
source: kind-pasteur-2026-07-11-S127 (cont.22) — attacking the S_c reciprocal-sum wall via its HYP-2644 Weyl form
depends_on:
  - HYP-2644   # the far-element plateau recursion (this proves its core decorrelation estimate)
  - THM-532    # meas(S7(E)) = M7(k) + corr(E), the seven-sector reduction
related:
  - THM-538    # the support-6 kernel / the lattice-sum form of the same wall
  - THM-699    # the support-6 kernel is a zero-mean permanent (the weight side; this is the S_c/analytic side)
  - HYP-2645   # Poisson dual = the finite x-cell integral (the convergent representation this bounds)
  - MISTAKE-078 # the divergent absolute envelope this sidesteps
external: Weyl equidistribution; Erdős–Turán; the Fourier coefficient decay of bounded-variation functions.
---

# THM-700 — the plateau decorrelation lemma

## Statement

Let `E = E' ∪ {w}` be a set of integer offsets with `0 ∈ E'` and `w = max E`. Write `p_0(F) = meas(S7(F))`
for the seven-sector cover measure, and let `Plat(E') := p_0(E') + (1/7)·p_1(E')`, where
`p_1(E') = meas{x : E' hits exactly six of the seven sectors}` (misses exactly one). Then

> **`|p_0(E) − Plat(E')| ≤ V(E') / (6w)`**,

where `V(E') = Σ_{s=0}^6 V(f_s)` is the total variation of the seven indicators
`f_s(x) = 1{E' misses exactly sector s}` (each `f_s` piecewise-constant, so `V(f_s) < ∞`).

In particular `p_0(E) → Plat(E')` at rate `O(1/w)` as the largest offset grows — the far element
decorrelates and acts as an independent uniform `1/7` fill of any singly-missed sector.

## Proof

### The exact cover decomposition (no error term)
`E = E' ∪ {w}` covers all seven sectors iff no sector is missed by *both* `E'` and the single point
`frac(wx)`. Splitting on the coverage of `E'`:
- `E'` covers all seven ⟹ `E` covers all seven (`w` irrelevant);
- `E'` misses ≥ 2 sectors ⟹ `E` misses ≥ 1 (the one point `frac(wx)` fills at most one) ⟹ not covered;
- `E'` misses *exactly one* sector `s` ⟹ `E` covers all seven iff `frac(wx) ∈ [s/7,(s+1)/7)`.

These three cases are disjoint and exhaustive, so **pointwise and exactly**
`1_cover(E,x) = 1_cover(E',x) + Σ_{s=0}^6 f_s(x)·1{frac(wx) ∈ [s/7,(s+1)/7)}`, `f_s = 1{E' misses exactly s}`.
Integrating, `p_0(E) = p_0(E') + Σ_s ∫_0^1 f_s(x)·1{frac(wx)∈s}\,dx`, and since `Σ_s f_s = 1{E' misses exactly 1}`
has integral `p_1(E')`, subtracting `Plat(E') = p_0(E') + (1/7)p_1(E')` gives the **centered error**
> `p_0(E) − Plat(E') = Σ_{s=0}^6 ∫_0^1 f_s(x)·g_s(wx)\,dx`, `g_s(y) := 1{y∈[s/7,(s+1)/7)} − 1/7`.

Each `g_s` is `1`-periodic with **mean zero**, and each `f_s` is bounded-variation.

### The Fourier-decay bound
Fix `s`. With `f_s(x) = Σ_m \hat f_s(m) e^{2πimx}` and `g_s(wx) = Σ_ℓ \hat g_s(ℓ) e^{2πiℓwx}` (`w∈ℤ`),
orthogonality collapses the integral to the diagonal `m + ℓw = 0`:
`∫_0^1 f_s(x)g_s(wx)dx = Σ_{ℓ} \hat f_s(−ℓw)\hat g_s(ℓ) = Σ_{ℓ≠0} \hat f_s(−ℓw)\hat g_s(ℓ)`
(the `ℓ=0` term drops since `\hat g_s(0)=0`, the mean-zero fact). Now:
- `|\hat g_s(ℓ)| = |\sin(πℓ/7)| / (π|ℓ|)` (the coefficient of a length-`1/7` interval), so `|\hat g_s(ℓ)| ≤ 1/(π|ℓ|)`;
- `|\hat f_s(m)| ≤ V(f_s) / (2π|m|)` for `m≠0` (the classical bound for bounded-variation functions).

Therefore, using `|ℓw| = |ℓ|w`,
`|∫ f_s g_s(wx)| ≤ Σ_{ℓ≠0} \frac{V(f_s)}{2π|ℓ|w}·\frac{1}{π|ℓ|} = \frac{V(f_s)}{2π^2 w} Σ_{ℓ≠0}\frac{1}{ℓ^2}
= \frac{V(f_s)}{2π^2 w}·\frac{π^2}{3} = \frac{V(f_s)}{6w}`.
Summing over the seven sectors, `|p_0(E) − Plat(E')| ≤ Σ_s V(f_s)/(6w) = V(E')/(6w)`. ∎

`V(E')` is finite and bounded by the core's total frequency: each `f_s` jumps only where some `frac(e'x)`
(`e'∈E'`) crosses the boundary of sector `s`, so `V(E') ≤ 14·Σ_{e'∈E'} e'` (a crude but explicit bound).

## Why it matters

- **It is the wall's inductive step, proved.** HYP-2644 reduced the entire LRC(14)-S3 unbounded direction to
  exactly this `|p_0 − Plat| ≤ C/w` estimate (with comfortable margin 0.13–0.18, an order looser than the
  tight `consec_k` finite check). THM-700 supplies the estimate with an explicit constant, via elementary
  bounded-variation Fourier decay — no lattice geometry, no conditional convergence.
- **It dissolves the reciprocal-sum wall's convergence problem.** The support-6 lattice sum `corr(E) =
  Σ_c D7(c) S_c(E)` is only conditionally convergent (MISTAKE-078); the honest convergent representation is
  the finite x-cell integral (HYP-2645). THM-700 works *directly* on that integral — the correlation
  `∫ f_s(x)g_s(wx)dx` IS the x-cell object — and the mean-zero `g_s` is precisely the centering that THM-699
  exhibited on the weight side (`Σ_c D7(c) = 0`). The two results meet: the weight has zero mean, the
  oscillation has zero mean, and the far element decorrelates them at rate `1/w`.

## Scope — what remains

THM-700 is the single-peel inductive step. To close the wide-spread half outright three things remain, all
downstream of this lemma:
1. **The `p_1` decorrelation** — the recursion tracks `Q(m) = max[p_0 + (1/7)p_1]`, so `p_1` needs the same
   far-element treatment (an identical BV/mean-zero correlation; mechanical from this proof).
2. **The accumulation** — peeling reduces a wide `k`-set to a bounded core through `k − k_0` steps; the
   per-step errors `V(E^{(i)})/(6 w_i)` must telescope below the margin. Since each `w_i` dominates its core
   and `V` grows only linearly in the core frequency, this is a summable-geometric bookkeeping, but it must
   be written.
3. **The sharp constant** — the crude `V(E') ≤ 14Σe'` overshoots the measured `≈0.2·(1/w)` by ~10²; the
   `|\sin(πℓ/7)|` numerator (which vanishes on `7∣ℓ`) and cancellation among the seven `f_s` recover it.

The tight margin (0.0014 at k=9) is **not** here — it lives entirely in the finite `consec_k` check. THM-700
governs only the loose wide-spread direction.

## Files

`04-computation/lrc14_plateau_decorrelation_kps_S127.py` (+ `.out`): the exact cover decomposition check,
`Plat(consec_8) = 0.36210` (matches HYP-2644), and the `O(1/w)` decorrelation with the `V/(6w)` bound holding
for `w` up to 1601.
