# All runners are observers: LRC from every angle at once (S521)

*claudebox-2026-06-01-S521. Reframing the "observer" — in the dynamical system every
runner observes the others; fixing one as the observer is just a choice of reference
frame. Looking through all n frames at once reveals loneliness-rotation and hidden
tightness. Extends the danger-graph and balanced-pair (Thm A/B) threads.*

## The symmetric picture

`n` runners, distinct speeds `w_1,...,w_n`, positions `w_i t mod 1`. The **full
danger graph** `G(t)` has an edge between any two runners within `1/n` (not just to
a fixed observer). Then:

> **runner k is lonely at t  <=>  k is ISOLATED in G(t)**, and
> **LRC  <=>  every vertex of `G(t)` is isolated at some time.**

Fixing one runner as "the observer" = passing to its reference frame (subtract
`w_k`), turning runner `k` into the stationary observer `0` and the others into
movers with difference-speeds `D_k = { w_i - w_k : i != k }`. So the `n`-runner LRC
is the conjunction of the observer-LRC over the `n` difference-sets `D_1,...,D_n` —
the same conjecture seen from `n` angles.

## Loneliness ROTATES (no one is always lonely)

Computed: `min_t (#lonely runners) = 0` for every speed set. At "clustered" times
(all points bunched in part of the circle) **nobody** is lonely. So it is false that
some runner is always lonely; loneliness is a property that **rotates** among the
runners as `t` varies. `LRC = each runner gets its turn at being isolated.` The
maximum simultaneous loneliness is `n` only at the **regular-polygon instant**
(all gaps `= 1/n`, `G` edgeless) — the unique moment all `n` are lonely at once.

## The angles have different difficulty — and HIDDEN TIGHTNESS

The `n` difference-sets `D_k` are NOT equally hard:

- **Collapsed frames are easy.** If two other runners share a relative speed
  (`|w_i - w_k| = |w_j - w_k|`), `D_k` has a repeat that collapses (same distance
  constraint), so runner `k` faces FEWER effective movers and is lonely on a set of
  positive measure. Example: AP `(0,1,2,3)` — the **middle** runners (speeds 1, 2)
  see collapsed `D = {1,2}` (a 2-mover problem) and have lonely-measure `1/4`.
- **Distinct frames are hard.** The runner whose `D_k` is all-distinct and
  AP-/extremizer-shaped is the **tight** observer. AP `(0,1,2,3)`: the **extreme**
  runners (speeds 0, 3) see the full `{1,2,3}` and are lonely only at the boundary
  (`measure 0`, tight).

**Hidden tightness.** A system can look generic from one angle and be the tight
extremizer from another. `(0,2,3,5,7)` is **tight from runner 3's frame** — runner 3
sees `D_3 = {1,2,3,4}`, the consecutive-AP extremizer — yet from runner 0's frame it
sees the generic `{2,3,5,7}`. So the difficulty of an `n`-runner system is
concentrated at whichever runner sees the most-extremal difference-set.

> **Reframed Thm B (multi-angle):** the `n`-runner system is tight (`M = 1/n`)
> **iff some runner's difference-set `D_k` is a tight observer-extremizer** (which,
> by Thm B, has a pair summing to `0 mod n`; canonically the AP).

## What the multi-angle view buys

1. **Always analyze from the hardest frame.** Pick the runner whose `D_k` is most
   distinct/AP-like; the collapsed frames are automatically fine. This is a genuine
   simplification: of the `n` observer-instances, only the extremal-frame ones need
   work.
2. **The regular polygon is the all-angles fixed point** — the one configuration
   solving every observer at once; the tight extremizers sit at its boundary.
3. **Loneliness-rotation suggests a recurrence framing.** LRC = the rotating
   "isolated-vertex" event visits every vertex over one period. Since `#lonely`
   averages `Sum_k measure(L_k) > 0` but dips to `0`, the question is whether the
   isolation event is *equidistributed* enough among runners to miss none — a
   symmetric/ergodic restatement.

## Creative hypotheses

- **HYP (hardest-frame reduction).** LRC for the `n`-runner system follows from
  observer-LRC for the single hardest frame `D_{k*}` (the most-distinct/extremal
  difference-set); the other frames hold automatically by collapse. If made
  precise, this reduces the `n`-fold conjunction to one instance.
- **HYP (rotation/recurrence).** Define the "isolation rotation" — the cyclic order
  in which runners become isolated as `t` increases. Conjecture every runner appears
  in it (LRC) because the largest gap (`>= 1/n` always, by pigeonhole) sweeps past
  every runner over a period, and at its widest moment doubly-flanks each in turn.
- **HYP (multi-angle tightness atlas).** Enumerate, for a speed set, all `n`
  difference-sets and their lonely-measures; the minimum over frames is the system's
  "slack". Tight systems are exactly those with a zero-slack frame; characterize the
  speed sets that are tight from *some* angle (the hidden-AP locus).

## Assessment

Seeing LRC from all `n` angles at once is the right symmetric picture: the full
danger graph with every-vertex-isolation, loneliness that rotates (min `#lonely =0`,
max `=n` only at the regular polygon), and the recognition that the `n` frames have
genuinely different difficulty — with **hidden tightness** (a system tight from a
non-obvious frame). The actionable payoff is the **hardest-frame reduction**: don't
fix an arbitrary observer, find the one seeing the most-extremal difference-set, and
the rest are free. It does not prove LRC, but it correctly locates the difficulty
(the extremal-frame observer) and reframes the conjecture as a rotation/recurrence:
the isolated-vertex event must visit every runner — the symmetric form of "the
observer becomes a source," now demanded of all `n` simultaneously over one period.
