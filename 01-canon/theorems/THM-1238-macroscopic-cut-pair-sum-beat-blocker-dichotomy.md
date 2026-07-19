---
id: THM-1238
title: THE MACROSCOPIC-CUT PAIR-SUM BEAT DICHOTOMY — a Kakeya path edge emits a nonzero mixed-clock blocker unless it is locked in one negative-curvature beat cell
status: RESERVED/PROOF IN PROGRESS. The pair-sum beat-block lemma and exact finite audit are being checked; the Kakeya cut composition and singleton residual are not yet promoted to proved status.
source: codex-2026-07-19-S78 continuation with tangent-stalk agent
depends_on: [THM-1156, THM-1217, THM-1219, THM-1235, THM-1236]
related: [HYP-7870]
---

# THM-1238 — macroscopic cut to pair-sum blocker

For an adjacent cut pair `x<y`, put `q=x+y` and sample a complete slow gap
on the `q`-beat grid.  At `t=p/q` the two defining danger masks coincide and
their active-pair slack is `14D(p)-q`, where
`D(p)=min(xp mod q,q-(xp mod q))`.

The candidate dichotomy is: if `q>=7c/3` or `y>6x`, the beat block contains a
point with `14D(p)>=q`; any six-comb cover must assign that point a third
blocker.  Otherwise failure is confined to a singleton beat block with
`2c<q<7c/3`, `y<=6x`, and negative curvature.  Equality can occur only when
`14|q`, recovering THM-1156's exact seam/third-support branch.

No all-packet blocker inconsistency or LRC(14) conclusion is claimed here.
