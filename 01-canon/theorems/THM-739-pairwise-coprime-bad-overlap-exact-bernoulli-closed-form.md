---
id: THM-739
title: The pairwise coprime bad-arc overlap in EXACT closed form — for coprime speeds c,c', the two 1/14-bad sets overlap with |bad_c ∩ bad_{c'}| = 1/49 + (1/cc')·[B₂({(c'−c)/14}) − B₂({(c'+c)/14})], a B₂-Bernoulli correction to independence; hence |bad_c ∩ bad_{c'}| ≤ 1/49 + 1/(4cc') and → 1/49 as cc'→∞. This is the pairwise coprime-overlap bound (klein-S292): the two bad sets are independent up to an explicit O(1/cc') Farey-scale term
status: PROVED (klein-2026-07-13-S293). Rigorous 5-line Fourier derivation: expand both bad-arc indicators, integrate, only frequencies nc+mc'=0 survive; coprimality forces n=c'k, m=−ck; the k-sum is a cos/k² Bernoulli series = 2π²B₂. Verified numerically (NG=2²², |direct−formula| ≤ 8e-7 across 12 coprime pairs incl. (2,3),(90,101),(100,171); every overlap ≤ 1/49+1/(4cc')). The |B₂|≤1/6 range gives the difference ∈[−1/4,1/4], hence the ±1/(4cc') envelope.
source: klein-2026-07-13-S293
depends_on:
  - THM-731   # the one-interval / disc route this feeds (the milder cancellation, klein-S292)
related:
  - THM-732   # kps's exact Bernoulli edge-pair disc form — SAME B₂ kernel, a sibling identity
  - THM-736   # mac-mini's Farey/three-gap deep-well closed form — same Farey-scale correction structure
  - HYP-6550  # klein-S292 (the one-interval bound this pairwise overlap sharpens: threshold 6/49 → 0.105)
  - HYP-6560  # klein-S293 (this result)
---

# THM-739 — the pairwise coprime bad-arc overlap is `1/49` up to a `B₂`-Bernoulli term

The exact pairwise version of the "bad sets are asymptotically independent" equidistribution that the
one-interval / milder-cancellation leg (klein-S292) needs.

## Setup

For a speed `c`, the `1/14`-bad set is `bad_c = {t∈[0,1) : ‖ct‖ < 1/14}` — `c` arcs of width `1/(7c)`,
total measure `1/7`. Its indicator is `1_{bad_c}(t) = 1_B(ct)` with `B=(−1/14,1/14)`,
`\hat{1_B}(0)=1/7`, `\hat{1_B}(n)=\sin(\pi n/7)/(\pi n)` for `n≠0`.

## The identity (PROVED)

For `gcd(c,c')=1`,
$$\boxed{\;|bad_c \cap bad_{c'}| \;=\; \tfrac1{49} \;+\; \tfrac1{cc'}\Big[\,B_2\!\big(\{\tfrac{c'-c}{14}\}\big) - B_2\!\big(\{\tfrac{c'+c}{14}\}\big)\Big]\;}$$
where `B_2(x)=x^2-x+\tfrac16` is the second Bernoulli polynomial and `{·}` the fractional part.

**Proof.** `|bad_c\cap bad_{c'}| = \int_0^1 1_B(ct)1_B(c't)\,dt = \sum_{n,m}\hat{1_B}(n)\hat{1_B}(m)\int_0^1 e((nc+mc')t)\,dt = \sum_{nc+mc'=0}\hat{1_B}(n)\hat{1_B}(m).`
With `gcd(c,c')=1`, `nc+mc'=0 \iff n=c'k,\ m=-ck` (`k\in\mathbb Z`). The `k=0` term is `\hat{1_B}(0)^2=1/49`.
For `k\ne0`, `\hat{1_B}(c'k)\hat{1_B}(ck)=\frac{\sin(\pi c'k/7)\sin(\pi ck/7)}{\pi^2 cc' k^2}`, and by
product-to-sum plus the classical `\sum_{k\ne0}\cos(2\pi\alpha k)/k^2 = 2\pi^2 B_2(\{\alpha\})` (at
`\alpha=(c'\mp c)/14`), the `k\ne0` sum collapses to `\frac1{cc'}[B_2(\{\frac{c'-c}{14}\})-B_2(\{\frac{c'+c}{14}\})]`. ∎

## The bound

`B_2` ranges over `[-\tfrac1{12},\tfrac16]` on `[0,1]`, so the bracket lies in `[-\tfrac14,\tfrac14]` and
$$\tfrac1{49}-\tfrac1{4cc'}\ \le\ |bad_c\cap bad_{c'}|\ \le\ \tfrac1{49}+\tfrac1{4cc'}.$$
So the two `1/14`-bad sets are **independent up to an explicit `O(1/cc')`** term: `|bad_c\cap bad_{c'}|\to
1/49=(1/7)^2` as `cc'\to\infty`. This is exactly the pairwise coprime-overlap bound `≤ 1/49` (klein-S292),
now exact. (Small coprime speeds carry a real correction — e.g. `(2,3)`: overlap `0.0476` vs `1/49=0.0204`;
`(90,101)`: `0.02041`.)

## Why it matters, and the two next steps

The one-interval leg (klein-S292/HYP-6550) bounds `|G(C)\cap[0,1/14)|` and gets stuck at threshold `6/49`
with a single speed because `\text{leading }6/7\approx1`. A **two-speed** inclusion–exclusion
`|G(\{c,c'\})| = 1-\tfrac27+|bad_c\cap bad_{c'}|` uses this identity to push the effective threshold down
toward `0.105`; the `k`-speed version (all `\binom{k}{2}` pairwise overlaps `\approx1/49`, higher-order
`\approx(1/7)^j`) is the route to the full multi-speed equidistribution — the milder cancellation.

**Two open extensions (same Fourier method):** (a) the **windowed** version
`|bad_c\cap bad_{c'}\cap[0,1/14)|` (multiply by `1_{[0,1/14)}` — a convolution of the Fourier coefficients,
still Bernoulli-shaped) is what the local one-interval bound literally needs; (b) `gcd(c,c')=g>1` replaces
`n=c'k,m=-ck` by `n=(c'/g)k, m=-(c/g)k`, giving a `g`-fold-denser correction. Both are mechanical.

## Provenance in the fleet's B₂ line

This is a sibling of kps's **THM-732** exact edge-pair Bernoulli disc form and mac-mini's **THM-736** Farey
three-gap closed form: the same `B_2({·/14})` kernel governs the disc, the deep-well far peel, and now the
pairwise bad-overlap. The recurring object across the covering endgame is `B_2` evaluated at Farey points
`k/14`.

## Addendum (klein-S294) — the microscopic resonance form, and why the WINDOWED overlap is not clean

The full-circle identity above has an exact *microscopic* companion. `bad_c ∩ bad_{c'}` is a union of arcs,
one per resonant pair `(j,k)` (`0≤j<c`, `k` = nearest, `|jc'−kc| < (c+c')/14`); by elementary
interval-overlap geometry each contributes length exactly
$$\ell_{j,k} = \frac1{cc'}\max\!\Big(0,\ \tfrac{c+c'}{14} - |jc'-kc|\Big).$$
Summing over all `c` residues `m_j = jc'-kc` (each hit once, coprimality) gives
`|bad_c∩bad_{c'}| = (1/cc')Σ_{|m|<(c+c')/14}(\tfrac{c+c'}{14}-|m|)`, whose leading term is
`(c+c')²/(196cc')` (`=1/49` at `c=c'`) with the `B_2` corrections above — a consistency check.

**The windowed overlap `W = |bad_c∩bad_{c'}∩[0,1/14)|` is the SAME sum restricted to the resonances whose
arc lies in `[0,1/14)` — a Farey/partial sum with NO one-line closed form**, unlike the full circle (the
window breaks the Fourier orthogonality; `\hat f` convolves against `\hat{1_{[0,1/14)}}`, so every
frequency contributes). Its size is governed by *how the resonances `jc'≈kc` land in `[0,1/14)`*:
- **Close speeds (`c'≈c`, i.e. CLUSTERS):** the small-`j` residues `m_j ≈ j(c'-c)` are all small, so many
  terms carry the near-maximal weight `(c+c')/14 − |m_j|` — the resonances **pile up near `0`** and `W` is
  LARGE (verified: `(99,101)` gives `W≈0.0051 ≈ 3.5×` the bulk `1/686`; `(50,99)`, `(23,45)` similar).
- **Far speeds:** the `m_j` spread and `W → bulk = (1/14)|bad_c∩bad_{c'}| ≈ 1/686`.

**Consequence (the honest negative).** The two-speed inclusion–exclusion refinement of the one-interval
bound (S292/HYP-6550) needs `W` *small*; it is not, for clusters. So **pairwise decorrelation FAILS near
`0` for close speeds** — the near-`0` equidistribution that `conc<7` needs is *intrinsically multi-speed*,
not a sum of pairwise terms. This is precisely why the milder one-interval cancellation is still not
elementary: clusters are exactly the close-speed regime where every low-order (pairwise) overlap is
correlated. See HYP-6570.

*Files: `04-computation/lrc14_pairwise_overlap_klein_S293.py`, `lrc14_windowed_overlap_klein_S294.py`
(+.out). HYP-6560, HYP-6570.*
