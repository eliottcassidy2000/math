---
id: THM-739
title: The pairwise coprime bad-arc overlap in EXACT closed form — for coprime speeds c,c', the two 1/14-bad sets overlap with |bad_c ∩ bad_{c'}| = 1/49 + (1/cc')·[B₂({(c'−c)/14}) − B₂({(c'+c)/14})], a B₂-Bernoulli correction to independence; hence |bad_c ∩ bad_{c'}| ≤ 1/49 + 1/(4cc') and → 1/49 as cc'→∞. This is the pairwise coprime-overlap bound (klein-S292): the two bad sets are independent up to an explicit O(1/cc') Farey-scale term
status: PROVED (klein-2026-07-13-S293). Rigorous 5-line Fourier derivation: expand both bad-arc indicators, integrate, only frequencies nc+mc'=0 survive; coprimality forces n=c'k, m=−ck; the k-sum is a cos/k² Bernoulli series = 2π²B₂. Verified numerically (NG=2²², |direct−formula| ≤ 8e-7 across 12 coprime pairs incl. (2,3),(90,101),(100,171); every overlap ≤ 1/49+1/(4cc')). The |B₂|≤1/6 range gives the difference ∈[−1/4,1/4], hence the ±1/(4cc') envelope.
source: klein-2026-07-13-S293
correction: >
  2026-09-06: the microscopic addendum omitted containment clipping and
  multiple resonances. The full-circle Bernoulli identity is unchanged.
  Its claimed uniform upper bound 1/49 is false; retain the stated error term.
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
1/49=(1/7)^2` as `cc'\to\infty`. This quantifies the earlier independence comparison. A uniform bound `≤1/49`
is false; the positive correction cannot be discarded. (Small coprime speeds carry a real correction — e.g. `(2,3)`: overlap `0.0476` vs `1/49=0.0204`;
`(90,101)`: `0.02041`.)

## Why it matters, and the two next steps

The one-interval leg (klein-S292/HYP-6550) bounds `|G(C)\cap[0,1/14)|` and gets stuck at threshold `6/49`
with a single speed because `\text{leading }6/7\approx1`. A **two-speed** inclusion–exclusion
`|G(\{c,c'\})| = 1-\tfrac27+|bad_c\cap bad_{c'}|` uses this identity to push the effective threshold down
toward `0.105`; the `k`-speed version (all `\binom{k}{2}` pairwise overlaps `\approx1/49`, higher-order
`\approx(1/7)^j`) is the route to the full multi-speed equidistribution — the milder cancellation.

**Two open extensions (same Fourier method):** (a) the **windowed** version
`|bad_c\cap bad_{c'}\cap[0,1/14)|` (multiply by `1_{[0,1/14)}` — a convolution of the Fourier coefficients,
still Bernoulli-shaped) is what the local one-interval bound literally needs; (b) for `gcd(c,c')=g>1`, substitution by the degree-g circle map gives exactly
the primitive-pair formula at `(c/g,c'/g)`. Its correction multiplier is
`g^2/(cc')`, with Bernoulli arguments `(c'±c)/(14g)`.

## Provenance in the fleet's B₂ line

This is a sibling of kps's **THM-732** exact edge-pair Bernoulli disc form and mac-mini's **THM-736** Farey
three-gap closed form: the same `B_2({·/14})` kernel governs the disc, the deep-well far peel, and now the
pairwise bad-overlap. The recurring object across the covering endgame is `B_2` evaluated at Farey points
`k/14`.

## Corrected microscopic resonance form — 2026-09-06

**PROVED correction to the klein-S294 addendum.** The original addendum
omitted the containment cap and restricted each tooth to a nearest partner.
At `(p,q)=(1,2)`, its origin formula gave `3/28`; the actual intersection
has length `1/14`. The old windowed formula also omitted clipping an overlap
arc to the window. Those formulas and their proposed asymptotic derivation
are SUPERSEDED. The full-circle Bernoulli proof above is unaffected.

Let `p<q` be positive coprime integers. Every overlapping pair of circle
teeth corresponds to exactly one signed integer `k` with
`|k|<(p+q)/14`. Its intersection length is

```text
ell_k = max(0, min(2p,p+q-14|k|))/(14pq).
```

The cap `2p` is the width of the smaller tooth in these units. Indeed the
two half-widths are `1/(14p),1/(14q)`, and their center distance is
`|k|/(pq)`; interval intersection is the smaller full width or the sum of
half-widths minus distance, whichever is smaller. Coprimality gives each
resonance once. Distinct tooth intersections have a positive gap, so the
number of open circle components is exactly

```text
J=2 ceil((p+q)/14)-1.
```

The complete measure is therefore

```text
mu = [p + sum_(k=1)^(ceil((p+q)/14)-1)
                  min(2p,p+q-14k)]/(7pq).
```

This agrees with the Bernoulli identity and, after primitive reduction,
works for arbitrary distinct speeds. For a window, use the length of the
intersection of **each actual overlap arc with that window**, including
partial and wrapping arcs. No unproved near-origin equidistribution follows
from the full-circle formula.

The historical window script evaluates literal masks numerically; its
computational body is unchanged, but its incorrect explanatory formula is
marked superseded. Its finite samples are not a universal window theorem.
The corrected exact geometry is independently checked in
[the extended grid audit](../../05-knowledge/results/third_20260906_grid_audit.md).
