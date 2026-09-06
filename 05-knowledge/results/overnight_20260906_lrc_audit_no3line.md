# Independent audit of the conditioned no-three-in-line skeleton

**Status: ACCEPT, analytical audit + independent FINITE-EXACT replay.**
I independently checked Sections 1--4 of
[overnight_20260906_no3line.md](overnight_20260906_no3line.md), including
all matching formulas, the sampling correction, and the complete `n=4`
distribution. No mathematical error was found in those claims.

## Matching identities and the actual random measure

A simple degree-two bipartite graph has `2n` edges and exactly `2n`
unordered adjacent edge pairs. An edge triple containing two such pairs
is a three-edge path, and there are exactly `2n` such paths, including
the four paths in an individual four-cycle. Three adjacent pairs would
be a triangle, excluded by bipartiteness. Inclusion-exclusion therefore
gives precisely

```text
m3=binom(2n,3)-2n(2n-2)+2n=2n(n-2)(2n-5)/3.
```

Independent uniform labels of the two shores act transitively on the
`binom(n,3)^2*3!` three-edge grid matchings. The inclusion probability of
one specified matching is the displayed `m3` divided by this denominator,
which simplifies to the coefficient in report equation (1). Every
non-axis Euclidean collinear triple is such a matching. The same argument
with four edges proves report equation (5); it does not include crossing
line-pair events and therefore gives no complete variance formula.

For a cycle of length `L`, the matching polynomial is
`lambda_+^L+lambda_-^L`, where `lambda_++lambda_-=1` and
`lambda_+ lambda_-=-t`. This follows from the path recurrence with the
closing-edge correction and is valid for the cycle lengths `L>=4` here.
Multiplying over components proves report equation (3). Since
`lambda_-/lambda_+=-t+O(t^2)`, a cycle of even length `2r` first contributes
its new count at degree `2r`, with coefficient one. The triangular recovery
of short-cycle counts is therefore valid.

For a single cycle on `2n` edges, the four-matching count is
`2n/(2n-4)*binom(2n-4,4)` when `n>=4`, which simplifies to
`n(n-3)(2n-7)(2n-5)/6`. Adding the degree-four correction `c4` proves report
equation (4). At `n=3` both the four-matching count and the asserted
polynomial vanish, so the boundary is correctly handled.

Each even cycle has exactly two alternating ordered edge colorings.
Choices on separate cycles are independent, so every uncolored board has
exactly `2^c` ordered decompositions into two disjoint permutations. This
is not an automorphism factor. Dividing by `2^c` compensates precisely
when grouping the ordered pairs into their uncolored union. Fixing the
first permutation leaves a derangement for the second; its cycles have
half the graph-cycle lengths. The weighted permutation-cycle exponential
formula gives the stated board-count series. No geometrical relabeling
invariance beyond the explicitly random shore labels is assumed.

## Independent small-grid census

[overnight_20260906_lrc_audit_no3line.py](../../04-computation/overnight_20260906_lrc_audit_no3line.py)
imports no repository code. It enumerates all `216` ordered disjoint
permutation pairs at `n=4`, groups their exact uncolored point sets into
`90` boards, and tests collinearity by literal determinants. This differs
from the producer's row/column-capacity recursion and line-bucketing path.

It independently reproduces:

| Cycle type | Boards | Triple mean | Triple variance | Quadruple mean | Zero probability |
|---|---:|---:|---:|---:|---:|
| `(4,4)` | 18 | 2 | `56/9` | `1/3` | `1/2` |
| `(8)` | 72 | 2 | `25/18` | `1/6` | `1/36` |

Every individual board has exactly its predicted `2^c` decompositions.
The resulting zero probabilities are `11/90` under uniform boards and
`5/27` under uniform ordered pairs. These verify the report's finite
mean/tail distinction and the sampling bias directly.

Reproduce with `python -B
04-computation/overnight_20260906_lrc_audit_no3line.py`; its matching output
is [overnight_20260906_lrc_audit_no3line.out](overnight_20260906_lrc_audit_no3line.out).
This audit does not independently replay the producer's larger `n=5,6`
censuses or make any asymptotic zero-event inference.
