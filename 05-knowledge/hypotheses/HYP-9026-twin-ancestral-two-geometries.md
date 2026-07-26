---
id: HYP-9026
title: "Two ancestral geometries on the twin ranks: harmonic staircase and entropic mediant tree"
status: >
  OPEN (structural law; both geometries FINITE-EXACT verified on the
  full census to centers 10^8 with hard checks). Inspired by the
  Farey/Stern-Brocot mediant construction (owner dispatch,
  arXiv:2011.06318 seed).
source: kind-pasteur-2026-07-26-S132
related:
  - THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry
  - THM-2443-twin-rank-parent-parity-margins-and-boundary-crossing-bridge
  - THM-2447-twin-gap-local-prime-harmonics-and-detrended-spectrum
  - HYP-9025-twin-gap-singular-series-partner-law
  - HYP-3003-summand-multiplicand-farey-basis-merge
script: 04-computation/twin_ancestral_geometries_kps_S132.py
output: 05-knowledge/results/twin_ancestral_geometries_kps_S132.out
---

# HYP-9026 -- the twins have a staircase and a tree

Typing per MISTAKE-268: `K = A002822`. The parent fibre `P(k)`
(THM-2422) admits many canonical selections; two extremes carry all
the structure, and their CONTRAST is the content.

## 1. The staircase (min partner; THM-2422's canonical pair)

`parent(k) = k - a*(k)`, `a*(k) = min partner`. Census facts:

- no orphans (`k >= 3` always has parents);
- depth is a constant-order fraction of the rank INDEX
  (`0.59 -> 0.29` across dyadic windows; max depth `118,503`);
- the partner histogram is the **first-hit transform of the
  pair-correlation harmonics**: top partners `5, 7, 30, 70` --
  `70 = 2 * 5 * 7` outranks every non-smooth smaller value because
  the distance-`a` correlation is boosted at both `p = 5` and
  `p = 7` (THM-2447 amplitudes). The Fourier reading is literal:
  a* is where the local constellation first comes "in phase".

## 2. The mediant tree (balanced partner; Stern-Brocot analog)

`parent_bal(k) = k - max{a <= k/2 : a, k-a in K}`. Census facts:

- a balanced parent EXISTS for every `k >= 3` (found within `26.6`
  candidates on average) -- each split is a windowed Seymour
  certificate `u + v = w`, `v <= k/2 < w`;
- the tree has genuinely logarithmic depth: max `28` vs
  `log2(k_max) = 24`; window means `log2 k + ~3.4` with the ratio
  to `log2 k` falling monotonically `1.378 -> 1.159`.

## 3. The law (OPEN)

1. **Mediant-depth law:** `depth_bal(k) = log2 k + O(log log k)`
   (or even `log2 k + O(1)` -- the census excess is ~3.4 and
   drifting down slowly). Any proof needs a Seymour-type statement
   in windows `[k/2 - H, k/2]` with `H = k^{o(1)}` -- i.e. a
   quantitative self-additivity, strictly between THM-2443's bridge
   and full Hardy-Littlewood.
2. **Staircase-depth law:** `depth_min(k) ~ index(k) / E[a*]` with
   `E[a*]` the singular-series-weighted first-hit mean, growing
   like `c log^2 k` -- so `depth_min(k)/index(k) -> 0` at the
   logarithmic rate visible in the census (`0.59 -> 0.29`).
3. **Contrast principle:** every canonical-parent selection
   interpolates between these two geometries; the selection's
   arithmetic content is measured by how far its depth profile sits
   from the entropic `log2 k` floor.

Cheapest decisive tests: (i) fit `depth_bal - log2 k` against
`log log k` on the census -- **DONE same session**: the window
excesses are `3.40, 3.54, 3.62, 3.64, 3.67, 3.65, 3.66, 3.66`
across `2^8..2^24`: a plateau at `~3.66`, NOT `log log` growth
(which would add ~0.9 over this range). The census favours the
sharp form `depth_bal(k) = log2 k + O(1)`, constant `~3.7`; (ii) prove
unconditionally that ANY additive-2-basis-with-density-`1/log^2`
sequence admits a balanced tree of depth `<= log2 k + O(width of the
needed window)`, reducing prediction 1 to a covering statement;
(iii) check whether the LRC owner-Farey rays (HYP-8871) admit the
same balanced-split certificate on their finite-state side -- the
mediant construction only needs the two-parent fibre, which the
Farey/Stern-Brocot tree supplies exactly (mediants ARE two-parent
sums in numerator and denominator).
