# Good-Cut Buckets as a Merged Coordinate

**Session:** opus-2026-05-29-S13

Lean forced a useful sharpening of the good-cut story.

The informal theorem says "no tiling has exactly one good cut." The formal
proof does not naturally start with the cut. It starts with the witness tile.
If an upward tile `(hi,lo)` exists, it automatically contributes the two cuts
`lo+1` and `lo+2`. The singleton gap is therefore not a global connectivity
phenomenon first; it is an interval-union constraint in the staircase.

## The Formal Core

`TournamentH7.GoodCuts` now defines:

- `StTile.crossesCut t k` as `lo(t) < k ≤ hi(t)`.
- `StTiling.goodCuts b`.
- `StTiling.goodCutCount b`.

It proves, without project-specific axioms:

- bucket `0` is exactly the all-down tiling;
- bucket `1` is impossible;
- grid reflection preserves bucket size.

The reflection proof is the most revealing piece. Grid reflection sends tile
intervals `[lo+1, hi]` to reversed intervals, and cut `k` to cut `n-k`. So
the good-cut bucket is invariant under the same reflection that underlies the
merged/self-complementary story.

## The Surprise

The companion computation `goodcut_bucket_merged_s13.py` found that for
n=3..6, every merged tournament class is pure with respect to `goodCutCount`.
That is, no merged class contains tilings from two different good-cut buckets.

This is stronger than what Lean proved. Lean proved reflection-invariance.
The computation suggests iso-class invariance, and even merged-class
invariance:

> `g(τ)` may descend from tilings to vertices of `G_n/Z_2`.

If true, the good-cut count is not merely a base-path coordinate. It is a
coarse quotient invariant hiding in tiling language.

## New Picture

The merged metagraph may carry a natural "connectivity height"

`0, 2, 3, ..., n-1`

with the missing level `1` as a universal forbidden floor. This height is
different from H, score, and even-graph degree:

- `g=0` is the transitive class.
- `g=n-1` is exactly the strongly connected tiling bucket.
- intermediate `g` measures how many cut barriers have been pierced.

The bucket transition matrix under single-tile flips looks like a Morse graph
on this height, but with long jumps from the bottom because every tile spans
an interval of cuts. This may be the cleanest local model for why the
principal line begins with a large jump instead of a one-step gradient.

## Engineering Hint

For Tournament TDA, `g` is cheap: compute the union of upward-tile intervals.
If HYP-1764 holds, `g` is a quotient-stable feature rather than a labeling
artifact. Even if it eventually fails, its failure cases will measure exactly
where base-path geometry and merged iso-geometry stop agreeing.

## S95 Closure

THM-354 closes the main question: good-cut count is exactly

```text
g(T,path) = n - scc(T).
```

So `g` descends first to ordinary tournament isomorphism classes, and then to
merged classes because complement preserves strong-component count. The
coordinate was not merely compatible with the merged quotient; it was measuring
strong-component defect all along.
