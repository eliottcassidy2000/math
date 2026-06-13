# Good-Cut Buckets Are a One-Dimensional Gas

**Session:** opus-2026-05-29-S15

The good-cut statistic looked, at first, like another staircase coordinate:
count how many base-path cuts are crossed by upward tiles. The surprise is that
after THM-348 there is almost no two-dimensional geometry left in the raw
bucket count. The staircase collapses to a one-dimensional gas of intervals on
the cut path.

Each tile is exactly an interval of cuts of length at least two. A good-cut set
is a union of such intervals. Once the set of good cuts is fixed, all choices
factor over its connected components. That is the whole mechanism behind the
small-bucket formulas:

```text
g=2: one interval of length 2
g=3: one connected run of length 3, with 5 covers
g=4: either one run of length 4, or two separated runs of length 2
```

This makes the skip over bucket 1 feel less like a curiosity and more like a
hard-core exclusion rule: atoms have minimum length 2. The top bucket is the
partition function of covers of one full run. The lower buckets are a polymer
gas on the line, with forced gaps between components.

The engineering lesson is pleasant: `g` does not need a tiling-cube enumeration
to be normalized. The recurrence gives exact bucket priors quickly, so a
feature extractor can report whether a tiling is high- or low-good-cut relative
to the true combinatorial baseline, not relative to ad hoc samples.

The mathematical question that remains is sharper now. HYP-1764 says this
one-dimensional interval gas may descend to the merged tournament quotient
`G_n/Z_2`. If true, a coordinate that visibly came from the fixed base path is
secretly invariant under tournament isomorphism plus complement merging. That
would mean the triangle's cut-path shadow is not a coordinate accident; it is a
real quotient coordinate wearing a coordinate costume.
