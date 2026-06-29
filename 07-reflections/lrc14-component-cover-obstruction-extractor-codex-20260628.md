# LRC14 Component-Cover Obstruction Extractor

The overlap-tax identity from HYP-3434 made the compression failure exact, but
it was still row-level.  HYP-3435 supplied branch-cell witnesses, but its
certificate was phrased as existence of a survivor.  HYP-3450 moves the problem
to the local component where a proof can actually bite.

For each component of `E_safe`, the exact audit asks whether branch `0` lives,
branch `1` lives, both live, or both die.  If a branch dies, it extracts a
minimal odd-bad cover.  If a branch lives, it records the endpoint rank of the
survivor window.  This is the missing structural data between scalar measure
and a proof certificate.

The strong finite signal is:

```text
rows_with_endpoint_rank_le_2_survivor = 135/135
dead_components_with_two_min_covers = 3770/3770
max dead paired cover rank = 6
```

So a counterexample is now much more constrained.  It cannot be just "too much
bad mass."  It would need every even-safe component to be simultaneously
covered by two odd-bad families, with explicit minimal blocker sets on each
component.  That object looks like a finite cut/flow certificate, not an
analytic tail estimate.

The graph route is now very concrete.  Put component vertices on one side,
odd-blocker labels on the other, split branch `0` and branch `1` as two colors,
and attach even endpoint gates as boundary vertices.  A failed row would be a
two-color saturation of every component.  The observed rank-`<=2` survivor in
every audited row suggests a bounded Menger cut, Green-current escape, or
algebraic-connectivity obstruction should exist.

The next session should try two things:

1. Prove the canonical `{1..11,13,84}` paired-cover obstruction by hand.
2. Build a graph Laplacian/conductance audit for the largest-rank example
   `random_covering_082`, where one dead component has paired cover rank `6`.

The proof slogan is now local:

```text
No primitive covering row can make every E_safe component dead_both.
```
