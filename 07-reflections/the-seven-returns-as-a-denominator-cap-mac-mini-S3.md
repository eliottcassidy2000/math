# The seven returns as a denominator cap

**Session:** mac-mini-2026-06-18-S3. **Result:** HYP-2597 (universal good centers + μ_consec closed forms).

The task was to look at the *integer sequences* the LRC(14) floor throws off and reframe them
to simplify. The floor is μ(E) = meas{x : the points {frac(ex) : e∈E} leave a circular gap >
2/7}, and the open crux is that this is bounded below over all integer shapes E. I went hunting
for closed forms in the sequence μ_consec(k) — and what surfaced was not a formula so much as a
*place*.

## Where the good points live

For the consecutive cluster {0,1,…,k−1}, the set of good x is, for every k from 5 to 25,
**exactly five intervals** — neighborhoods of `0, 1/2, 1/3, 2/3`. Not "mostly"; exactly. No good
interval ever sits near 1/4, or 2/5, or any other rational. And the reason is a one-line fact
that holds for *every* integer shape, not just the AP:

> At `x = a/b`, the points `frac(e·a/b)` all land on the `(1/b)`-grid, so the largest gap is at
> least `1/b`. That gap exceeds `2/7` exactly when `1/b > 2/7`, i.e. **`b < 7/2`**, i.e.
> `b ∈ {1,2,3}`.

So `1/2, 1/3, 2/3` (and the integers) are good for *any* shape, unconditionally; `1/4` is not.
The threshold's denominator — the 14 of LRC(14), halved to the 7 of the `2/7` gap — comes back a
third time, now as the cap `b < 7/2` on which rationals are *forced* to be lonely. The same 7
that sets the band width sets the size of the universal skeleton of good times. Three of the
project's recurring sevens — the band `1/14`, the kernel period `7`, and now the center cap
`7/2` — are the same arithmetic seen at three scales.

## What simplified, and what didn't

Around each universal center the good interval has a width I could pin: near 0 (with its mirror
near 1) the total is exactly `10/(7(k−1))`, proved straight from the elementary near-0 lemma and
the reflection symmetry; near `1/2` it is `3/(7(k−2))`. So the consecutive floor obeys an
explicit, positive closed-form bound `μ_consec(k) ≥ 10/(7(k−1)) + 3/(7(k−2)) ~ 13/(7k)` — and the
three-gap literature, it turns out, has no closed form for this max-gap-exceeds-threshold
measure at all, so the decomposition is genuinely new.

But the simplification has an honest ceiling, and naming it is the point. The universal centers
are *fixed* — four points, not a recurring family. Their good intervals have width `~1/spread`.
So for a cluster of bounded spread (the AP, the perforated minimizers) they hand over a clean
positive floor; but as the spread grows the four widths shrink to nothing, and the measure that
keeps μ bounded must come from *somewhere else* — from the orbit, drifting fast, repeatedly
stumbling into a gap near rationals of denominator up to `k−1` that are *not* universal. Those
extra intervals are shape-dependent, governed by the continued fraction of the offsets, and they
are exactly the "colored discrepancy" that codex's lane is chasing and the spread-bound `B(k)`
that every angle reduces to.

So the sequence lane did what reframing usually does here: it made the *easy* part transparent
(the bounded-spread skeleton is five intervals at the b<7/2 centers, with closed-form widths) and
drew a sharp line at the *hard* part (the large-spread recurrence is irreducibly Diophantine). It
did not finalize the proof. What it added is a clean coordinate system: the floor is the
universal skeleton plus a Diophantine remainder, and the remainder is the whole game. The seven
tells you where to stand; it does not yet tell you the orbit always comes back.
