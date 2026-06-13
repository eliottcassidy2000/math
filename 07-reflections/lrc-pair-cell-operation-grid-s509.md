# LRC pair-cell operation grid after S509

S506b found that LRC loneliness is a vector of runner-level arc gauges, not one
tournament.  S509 adds the missing second layer: make the runner-pairs
themselves into vertices.

That changes the question.  A runner-level tournament asks "which runner is
more protected, exposed, central, or dangerous?"  A pair-cell tournament asks
"which distance cell is carrying the obstruction?"  This is closer to the
actual LRC geometry, because the forbidden intervals, half-turn walls, and
endpoint handoffs all live on pair or endpoint cells before they project down
to a runner score.

The clean new metric is `edge_danger_deficit`.  Its vertex value is

```text
max(0, 1/n - dist(i,j,t))
```

on each runner pair.  Completing the pair-cell comparison to a tournament and
reading its tie rate tracks the close-pair burden extremely well in the exact
small clocks: `tie_rate -> close_pair_count` has mean absolute Spearman
`0.986`.  This is a stronger local threshold metric than completed `H` on that
second-order tournament.

The static arithmetic criteria behaved exactly as they should have: they did
not become time-loneliness meters.  `edge_dyadic_row`, `edge_odd_core`,
`edge_same_odd_chain`, and `edge_product_sum_defect` are fixed for a speed set,
so within a clock movie their Spearman scores are zero.  Their job is not to
replace geometry.  Their job is to label the pair-cells by operation type.

That is where the A000568 analogy becomes precise.  A000568 counts unlabeled
tournaments, and its Burnside formula only sums odd cycle partitions.  In the
repo language, tournament structure survives on odd cores.  The user's `x+2`
and `x*2` grid is the natural-number version of the same split:

```text
horizontal: odd core -> odd core + 2
vertical:   odd core -> 2 * odd core -> 4 * odd core -> ...
```

LRC hard rows move on that grid.  The denominator is `2^r * odd_core`, and the
pair speed gaps inherit their own rows and odd cores.  This is why a scalar
gap number keeps losing information: it forgets which branch a close pair is
on.

Addition and multiplication also fit this picture.  The additive shadow is the
complete transitive order: `x -> z` whenever some positive `y` has `x+y=z`.
The multiplicative shadow is the sparse divisibility DAG.  Product-sum
equations are their critical pairs:

```text
P = product(F)
S = sum(F)
D = P - S
D ones repair the additive fold: D + S = P
```

For two nonunit factors at total arity `k`, this becomes

```text
(a-1)(b-1) = k-1.
```

S509 uses that equation as a pair-cell label.  It is not a loneliness meter by
itself, but it says which runner-pairs sit near an addition/multiplication
critical pair for the denominator.  For `k=14`, the best product-sum seed is
`1^11 + 2 2 5`; for `k=18`, it is `1^15 + 2 3 4`.  Those are static
arithmetic signatures of the hard first-even rows.

The updated proof-search object is:

```text
runner metric vector from S506b
+ pair-cell metric vector from S509
+ endpoint owner/pressure labels from S503-S504
```

The next corridor experiment should compute this combined vector on every
half-turn cell crossed by a positive lonely interval.  A real obstruction
should leave a nontrivial core after the pair-cell transitive tie-path
completion is quotiented away.
