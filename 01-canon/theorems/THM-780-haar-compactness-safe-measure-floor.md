---
id: THM-780
title: Haar-compactness turns every strict lonely-runner margin into a uniform safe-measure floor
status: CLAIMED (namespace reservation; complete proof drafted, independent audit and integration in progress)
source: codex-2026-07-14-S10
depends_on:
  - LRC(<=13)
related: [THM-777, HYP-4452, HYP-4472, MISTAKE-117]
---

# THM-780 — Haar-compactness safe-measure floor

This file reserves THM-780 for the following general statement.

> Fix `d` and `0 <= alpha < beta <= 1/2`.  Suppose every admissible
> `d`-speed integer vector has a time at which every coordinate has circle
> distance at least `beta`.  Then there is a constant
> `c(d,alpha,beta)>0` such that every such vector has `alpha`-safe measure at
> least `c(d,alpha,beta)`.

The intended proof passes from each speed vector to Haar measure on its cyclic
subgroup of the `d`-torus, stabilizes every integer character relation along a
hypothetical zero-measure sequence, carries the `beta`-deep witnesses into the
annihilator limit subgroup, and contradicts Portmanteau on the open
`alpha`-safe box.  The LRC14 application is `d=12`, `beta=1/13`,
`alpha=1/14`, which would prove the qualitative global floor conjectured in
THM-777 (but not its sharp value `7/858`).

The status remains CLAIMED until the saturated-relation, Fourier convergence,
full-support, and admissibility quantifiers have been audited independently.
