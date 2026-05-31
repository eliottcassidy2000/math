# LRC Tournament 56 Bridge

Session: `codex-2026-05-31-S369`

The user's decomposition question was exactly pointed at the missing layer:

```text
56 - 12 = 44
56 - 14 = 42
```

At first this looked like loose numerology: `56` is both the S367 missed-cell
count and the number of unlabeled tournaments on six vertices; `12` is the
number of unlabeled tournaments on five vertices; `14` is the number of runners
in the current LRC frontier.

The exact ledger makes it sharper.

## The 12+44 Split

For tournaments on six vertices:

```text
T(6) = 56
self-converse classes = 12
chiral classes = 44 = 22 converse pairs
```

So `56 = 12 + 44` is not merely

```text
T(6) = T(5) + 44.
```

It is also

```text
six-vertex tournaments = self-converse core + chiral residue.
```

The fact that the self-converse layer at `n=6` has the same size as the entire
unlabeled layer at `n=5` feels like the real bridge.  It says the transition
from five to six vertices adds chirality as the main new mass.

The S367 LRC missed cells also have a mirror structure.  The eight alpha
stencils are four mirror pairs under

```text
bins -> 13 - bins.
```

That means the first comparison to six-vertex tournaments should not be a raw
class bijection.  It should compare the LRC mirror-pair quotient to the
self-converse/chiral split.

## The 14+42 Split

The S367 missed cells factor as

```text
56 = 8 alpha stencils * 7 odd shifts.
```

The eight stencils have four widths, each appearing with its mirror:

```text
1/728, 1/882, 1/1176, 1/1386.
```

The widest/outer mirror pair contributes

```text
2 stencils * 7 shifts = 14.
```

The remaining six interior stencils contribute

```text
6 stencils * 7 shifts = 42.
```

So the user's `56-14=42` is visible internally in the missed-cell ledger: the
fourteen-runner number removes one full mirror pair, leaving the six-stencil
interior.

Equivalently,

```text
8 stencils = 2 outer caps + 6 interior positions.
```

This is the cleanest direct contact with six-vertex tournaments: after the
outer mirror caps are removed, the LRC residue is literally a six-by-seven
interior grid.

On the tournament side, the exact six-vertex strong-connectivity split is

```text
strong = 35
non-strong = 21
```

so `42 = 2 * 21` is also twice the non-strong six-vertex layer.  This is
probably not accidental in repo language: `21` is the parity-forgotten half of
`42`, and non-strong tournaments are exactly the decomposable boundary layer.

## The 8 Gap

The old perspective conjecture says:

```text
sum of vertex orbits over n-classes predicts (n+1)-classes.
```

It works for `3 -> 4` and `4 -> 5`, but at `5 -> 6` the exact ledger gives:

```text
sum vertex orbits over T(5) = 48
T(6) = 56
gap = 8.
```

That `8` is exactly the number of LRC alpha stencils.  This may be the best
"what are we missing?" answer: five-vertex perspectives miss an eight-fold
stencil layer at the six-vertex boundary.  On the LRC side, the fourteen-runner
near-blocker is also missing exactly eight stencils, repeated across the seven
odd shifts.

## Working Hypothesis

The relevant comparison is:

```text
LRC:
  56 missed cells
  = 7 odd shifts * 8 stencils
  = 14 outer mirror-pair cells + 42 interior cells
  = 7 shifts * (2 outer caps + 6 interior positions)

Tournaments:
  56 six-vertex classes
  = 12 self-converse classes + 44 chiral classes
  = 35 strong + 21 non-strong
  = 48 five-perspectives + 8 missing perspective classes
```

The common missing object appears to be an eight-fold chirality/perspective
boundary.  It is not a class count alone; it is a failure of one-level
projection from five objects to see the mirror-paired six-object residue.

For the next search, I would test whether the eight LRC stencils can be
assigned to the eight "extra" six-vertex classes beyond five-perspective
prediction.  If that works, the LRC quotient proof may be able to borrow the
tournament repo's old self-converse/chiral grammar: prove that every nonzero
normalized vector either falls into the outer fourteen-runner mirror pair or
into the forty-two-cell interior boundary, and then discharge each layer with a
different symmetry.
