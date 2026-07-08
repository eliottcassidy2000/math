# The density floor was a circle-covering problem all along

*mac-mini-2026-07-07-S57. Written after THM-657 reframed the k=11,12,13 density-floor legs
as "do 13 arcs of length 1/7 cover the circle?" — the owner's "no diameter / Vitali" hint,
made literal.*

## The move

For a k-set `E`, the density-floor quantity is `mu_{1/7}(E) = P_x(maxgap{frac(e_i x)} > 1/7)`.
The whole fleet has read this through the **gap** functional `W = sum_i (g_i - 1/7)_+`
(mac-mini-S41's avoidance profile, kps's tent, klein's windows). One sign flip in
interpretation changes everything: `W` is not "excess gap," it is **uncovered measure**.
Give each phase `frac(e_i x)` the arc of length `1/7` *ahead* of it. A point is uncovered iff
its preceding-`1/7` arc holds no phase iff it sits more than `1/7` past its gap's start. So

> `W(x) = 1 - meas(union of the k arcs)`,  and  `mu = P(W > 0) = P(the arcs fail to cover).`

`mu_{1/7}(E)` is, exactly, the probability that `k = 13` arcs of length `1/7`, placed at the
Kronecker phases `{frac(e_i x)}`, fail to cover the circle. Total arc length `13/7 = 1.857 >
1`, so covering is possible; the question is how often the rotation `x` leaves a gap.

## Why this is the right frame

Three things fall out immediately that the gap-reading obscured:

1. **A free diameter-free bound.** `0 <= W <= 6/7` pointwise (one gap can be huge, the rest
   `>= 0`), so `mu >= (7/6) E[W]` by Markov — no diameter, no shape, no tent. The tent
   (THM-651) and the energy floor (THM-656) both DIE at `k >= 11` because their `E[F]`
   exceeds the toll; the window floor (THM-653) DECAYS as `146/(35 diam)` and falls below the
   bar at exactly `diam >= 76` (the k=13 residual). The covering bound has none of these
   failure modes — it is the correct tool for the tail the others can't reach.

2. **One lemma discharges three legs.** The most-efficient coverer is the configuration whose
   phases are most spread — the arithmetic progression `{jx}`, whose phases are the
   three-gap/Kronecker sequence, maximally anti-clustered. Any less-regular `E` covers worse,
   leaves gaps more often, so `mu(E) >= mu(block)`. That is *exactly* THM-530's "consecutive
   minimizes `mu`," and it now reads as an optimal-covering statement. And `mu(block)` clears
   every bar with room: `0.6263 / 0.5699 / 0.4425` at `k = 11/12/13` versus
   `0.3312 / 0.1993 / 0.0565` — margins **1.9x / 2.9x / 7.8x**. Prove one covering lemma,
   close three legs, diameter-free.

3. **It joins a literature.** Random covering of the circle by equal arcs is classical
   (Stevens 1939: `P(cover) = sum_j (-1)^j C(k,j)(1 - j l)^{k-1}`, giving iid baseline
   `mu_iid = 0.988` at k=13). Our arcs aren't iid — they are a *deterministic Kronecker
   sequence*, which covers BETTER than iid (block `0.4425 < 0.988`). "The AP is the tight LRC
   case" becomes "the AP is the optimal deterministic coverer" — a statement about which
   single-parameter arc family minimizes the gap probability, where tools exist (rearrangement,
   negative association of Kronecker phases, three-gap majorization) that never applied to the
   ad hoc tournament inequality.

## The honest wall, relocated

The reframing does not hand over the proof — it relocates the wall to a cleaner place.
Naive Bonferroni-3 on `E[W]` is useless (`-0.67`: 13 arcs of total length 1.86 overlap far
too much for inclusion-exclusion to truncate) — the same overlap that made mac-mini-S41's
truncated-lattice bounds a dead end and made the cherry-tree Bonferroni fragile. And a
subtlety worth its own line: **`mu` and `E[W]` have different minimizers** — the block
minimizes `mu` (it is rarely uncovered) but NOT `E[W]` (when it is uncovered, it is uncovered
by more; 30/150 random 13-sets have smaller `E[W]`). So the `mu`-route needs the genuine
extremal lemma, while the `E[W]`-route needs a family-independent floor `E[W] >= 0.0484`
(empirical min `0.11`). Two clean targets where there was one murky one.

## The pattern

This session's two moves rhyme with last session's. There, "average not sup" dissolved a
resonance wall by refusing to bound a sum by its largest term. Here, "uncovered not
excess-gap" dissolves the k>=11 vacuity by refusing to read a covering probability as a
first-moment tent. Both times the obstruction was in how we *wrote* the quantity, not in the
quantity. The tent normalization that goes vacuous at `k >= 11` (klein's "wrong-normalized")
is the tell: `E[F] = k(k-1) int f` grows like `k^2`, but coverage is a `k/7`-vs-`1` balance —
the tent was counting pairs when the problem was counting arcs. Read the object as what it is
— arcs covering a circle — and the diameter, which was only ever an artifact of the
pair-drift rate, disappears.

Follow the covering.
