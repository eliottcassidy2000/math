# LRC operation-grid arc criteria after S511

S509 taught the right caution: the `x+2`/`x*2` grid is real, but static
arithmetic labels do not move with time.  S511 was the next move: push those
labels back down onto runners by aggregating pair-cell danger through them.

The distinction is now sharp.

Pure branch labels are coordinates:

```text
dyadic_branch_pressure
odd_core_branch_count
same_odd_chain_degree
additive_shadow_degree
multiplicative_shadow_degree
product_sum_interface
```

They have zero within-clock Spearman signal because the speed set is fixed.
That is exactly what a coordinate label should do.  It says where a pair-cell
lives on the operation grid, not whether the current time is lonely.

The hybrid labels are moving gauges:

```text
dyadic_danger_curvature
odd_chain_danger
additive_danger_interface
multiplicative_danger_interface
product_sum_danger
```

These multiply the static operation label by the current pair-cell deficit

```text
max(0, 1/N - dist(i,j,t)).
```

That is the first convincing bridge between the user's operation-grid picture
and LRC loneliness.  The grid gives the branch coordinate; the deficit gives
the clock coordinate.

The scorecard makes this visible.  The taut but essential threshold gauge
`s506_lrc_close_sector` perfectly tracks close-pair count by tie rate.  After
that, the best new signals are hybrid operation gauges:

```text
additive_danger_interface  mean |rho| = 0.882
product_sum_danger         mean |rho| = 0.840
dyadic_danger_curvature    mean |rho| = 0.792
multiplicative_danger_interface  mean |rho| = 0.783
```

The lesson is not that addition beats multiplication or product-sum beats
dyadic height.  The lesson is that each operation reads a different shadow of
the same danger fiber.  Addition sees horizontal transport among odd-core
cells.  Multiplication sees vertical refinement along doubling/divisibility
branches.  Product-sum sees the interface where additive arity is paid for by
multiplicative factor structure:

```text
(a-1)(b-1)=N-1.
```

This is also the cleanest version of the A000568 analogy so far.  A000568 is
the base count of unlabelled complete orientations.  Its Burnside formula keeps
only odd cycle partitions, which is the quotient-level version of odd-core
survival.  LRC sits over that base with a fiber:

```text
endpoint labels
+ gap lengths
+ pair-cell operation-grid labels
+ operation-weighted current danger
```

So a counterexample is not just a bad tournament class.  It is a quotient walk
whose entire fiber remains endpoint-forbidden after every available operation
weighted repair is accounted for.

The current metric vector I would hand to the next LRC corridor script is:

```text
phase_H
lrc_close_tie_rate
incident_danger_origin_score
dyadic_danger_score_hist
additive_danger_score_hist
multiplicative_danger_score_hist
product_sum_danger_score_hist
strict pressure SCCs
endpoint-labelled origin bracket
```

That vector is deliberately not a single scalar.  The LRC seems to be asking
for a section through a labelled quotient space, and a scalar keeps throwing
away the labels that make the section visible.
