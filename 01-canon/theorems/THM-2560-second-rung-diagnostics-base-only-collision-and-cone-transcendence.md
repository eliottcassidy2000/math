---
id: THM-2560
title: "Second-rung diagnostics on the canonical typed row: base-only first collision (r = 3) and cone-transcendence of the positive drift"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT (two exact companions; the
  realization gate triple-verified by independent routes; the sidecar
  gated by full 13^3-table reconstruction; python/-O agree); independent
  hostile audit REQUESTED.  Both results are mechanism NEGATIVES with
  exact redirects; neither contradicts THM-2550's positive drift.
  SCOPE: canonical TYPED row THM-2309 (25) only (not an asserted scalar
  cover); no physical current; no row removed; LRC(14) OPEN.
source: opus-2026-07-27 (executing the two finite tests specified by the
  post-MISTAKE-281 realization/transport maps)
depends_on:
  - THM-2471 (collision index (24), fibre product (31), sidecar (44)-(47))
  - THM-2368 (rooted words (9), bare tensor (29), sidecar (34)-(37), Sec 8)
  - THM-2365 (drift dichotomy)
  - THM-2550 (the positive-drift/non-replica inputs, with MISTAKE-283 scope)
related:
  - THM-2512 (Section 5 bridge test; its u-slaved form needs r = 5)
  - MISTAKE-281 (common-base discipline, obeyed by both companions)
scripts:
  - 04-computation/lrc14_ancestry_realization_stage_opus_20260727.py
  - 04-computation/lrc14_phase_cone_sidecar_opus_20260727.py
outputs:
  - 05-knowledge/results/lrc14_ancestry_realization_stage_opus_20260727.out
  - 05-knowledge/results/lrc14_phase_cone_sidecar_opus_20260727.out
---

# THM-2560 -- what the second rung looks like, exactly

Row: `THM-2309 (25)`, `w = (1,14,27,40,53,66,13,2197,742586)`, owner
`j = 1`, word `sigma = {a}`.

## (A) The first collision is at depth 3 and the deep root is BASE-ONLY

The THM-2471 (24) collision index of the source-refined packet
(`e = 1_{E_1}`, `f = 1_Q P^2 1_{E_1}`, `c = c_1 = 13`) is

```text
r = 3:   I_1 = I_2 = 0 exactly,   I_3 = 9926558757352/109707098520974955,
```

each value confirmed by three independent exact routes; service mass
`169 I_3 = 9926558757352/649154429118195`. The collision depth
`d = 13^3 = 2197` COINCIDES with the middle blocker scale `c_2` (a new
typed-row fact: `r = nu_13(c_2)`). The deep-root sidecar (44)-(47) at
`C = c_3 = 2*13^5` gives `delta = 2197`, `h = C/d = 338 = 2*13^2` with
`13 | h`: the deep phase is **base-only** -- independent of the
collision root `u` AND the sheet `a`. Consequently the `r = 5` affine
invariant `theta = t - 2u` does not exist here, and THM-2512 Sec 5's
u-slaved toothpick bridge is NOT instantiable on this row's packet at
this clock. Redirects the files license: (i) a base-only bridge (pair
the collision colours against a y-side cut probe); (ii) other lawful
clocks/word strata where the collision index may hit the `r = 5`
window. (Neither the `r >= 6` structural negative nor `r = 5` occurred.)

## (B) The positive drift is CONE-TRANSCENDENT

The THM-2368 (37) phase/event sidecar of the bare tensor `H_E` was
built exactly: `1,489,966` chambers (`212,845` contributing; `19`
distinct rooted snapshots; `10,998` aggregated states), GATED by exact
reconstruction of the full `13^3` table (all `2197` cells) and by the
row-sum match to the stored successor numerators. All `168` nonzero
target-mode walk endpoints are nonzero (consistent with
`D_{H_E} > 0`, THM-2550). But the Sec 8 phase-cone certificate **fails
for every one of the 168 characters**: each walk's aggregated direction
classes (5,740-10,322 per character) span MORE than a half turn -- no
unit rotates any walk into a closed half-plane. The canonical drift is
therefore positive but NOT half-plane-explainable: every mode survives
only through incomplete cancellation of a full-rotation phase word.
First cancelling event (the exact diagnostic anchor): character
`(lam,mu) = (0,1)`, chamber left endpoint `y = 9653643/10396204`
(within `2.4e-6` of the owner-window edge `13/14`), first direction
leaving every open half-plane
`(-6,0,0,0,-1,-2,-2,-2,-2,-2,-1,0,0) in Z[zeta_13]`.

## Consequences for the uniform program

The THM-2368 (38) covering-row argument cannot be anchored on a plain
phase cone of the (37) sidecar: any uniform positivity mechanism must
(i) exploit the deep-coloured `a != 0` mode lines (untested here),
(ii) quotient by simultaneous-translation/factor-label symmetries
before the cone test, or (iii) break the rotating circulant word
BEFORE y-integration -- e.g. at the recorded first cancelling event on
the owner-window edge. Similarly, the ancestry/Boolean realization on
live rows must either go base-only or target row/clock/stratum
configurations with `r = 5`. Both companions obey MISTAKE-281's
common-base discipline (all pairings inside one integral / one chamber
partition; gates before interpretation).
