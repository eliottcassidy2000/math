---
id: THM-384
name: lrc-observer-adjacent-gap-source-criterion
status: PROVED
date: 2026-06-01
session: codex-2026-06-01-S516
depends_on:
  - THM-381
  - THM-382
  - THM-383
---

# THM-384: LRC safety is exactly the observer-adjacent threshold gap condition

## Statement

Let a Lonely Runner system have total denominator `n`, with the stationary
observer marked as vertex `0`.  At time `t`, place all runners on the circle
`R/Z` and complete equality walls by any fixed tie path that puts tied points
next to one another.

Let `g_left(t)` and `g_right(t)` be the two circular gaps adjacent to the
observer in this compactified cyclic order.  If a moving runner is tied with the
observer, one of these adjacent gaps is `0`.

Then the observer is lonely at time `t` if and only if

```text
g_left(t) >= 1/n  and  g_right(t) >= 1/n.
```

Equivalently, the observer's LRC-good status is a function of the closed
`1/n` threshold color of the two observer-adjacent gaps.  The good fiber is
exactly `(long, long)`.

## Proof

Normalize the observer position to `0` on `R/Z`.  A moving runner at position
`x` is safe from the observer exactly when its circular distance from `0` is at
least `1/n`, i.e. when

```text
x in [1/n, 1 - 1/n]
```

after choosing the representative `x in [0,1)`.

The clockwise adjacent gap from the observer is the least positive clockwise
distance from `0` to any runner, with value `0` if a runner is tied with the
observer.  Thus `g_right(t) >= 1/n` holds exactly when no runner lies in the
right forbidden arc `[0,1/n)`.

Similarly, the counterclockwise adjacent gap is the least positive
counterclockwise distance from `0` to any runner, again with value `0` on an
observer tie.  Thus `g_left(t) >= 1/n` holds exactly when no runner lies in the
left forbidden arc `(1-1/n,1]`, including the observer-tie endpoint through the
compactified tie path.

The union of those two forbidden arcs is precisely the set of positions whose
circular distance from the observer is less than `1/n`.  Therefore all moving
runners have circular distance at least `1/n` from the observer if and only if
both adjacent gaps have length at least `1/n`.

Boundary equality is included because LRC uses the closed inequality
`distance >= 1/n`.

## Verification Record

The proof above is elementary.  The S516 audit script checks the criterion over
compactified exact samples for bounded primitive systems:

```text
04-computation/lrc_gap_threshold_proof_route_s516.py
05-knowledge/results/lrc_gap_threshold_proof_route_s516.out
```

The bounded audit reports zero mismatches for total `n=3,4,5,6` in the sampled
windows and records that the only good adjacent threshold fiber is `(1,1)`.

## Significance

THM-382 showed computationally that threshold-decorated gap fibers, unlike raw
phase or rank fibers, separate good and bad states in bounded small cases.
This theorem explains why: once the observer is marked, LRC-good status is
exactly the two-bit closed threshold color of the observer-adjacent gaps.

THM-383 remains essential because the criterion is closed and many natural
witnesses are wall-only.  Thus the global proof object is a compactified
arithmetic walk in the threshold-decorated gap fiber, not just an open-cell
half-turn tournament menu.

This theorem does not prove LRC.  It sharpens the remaining burden: prove that
every admissible compactified runner-clock walk visits the source-gap fiber
`(long,long)`.

## Related

- THM-381: observer-source marked reachability.
- THM-382: bounded threshold-fiber separation.
- THM-383: half-turn boundary compactification audit.
- HYP-1982: threshold-decorated tournament-class fiber.
- HYP-1986: compactified source-gap forcing route.
