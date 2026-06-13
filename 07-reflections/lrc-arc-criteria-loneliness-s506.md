# LRC arc criteria after S506

The useful answer is not another single definition of an LRC tournament.  The
useful object is a small panel of arc criteria, each responsible for a different
piece of loneliness data.

`phase_half` remains the global spread gauge.  It is the only tested criterion
where `H` itself is a good small-clock feature: `H -> max_gap` has mean
absolute Spearman `0.790`, with `H~max_gap = -0.723`.  That agrees with
HYP-1971: low `H` detects unanchored semicircle bunching, while high `H`
belongs to the regularized side of the runner clock.

The best LRC-threshold gauge in the scorecard is not an H gauge.  The
`lrc_close_sector` switch puts strict arcs only when one runner lies inside
another runner's `1/n` danger sector; otherwise the tie path completes the
tournament.  Its useful statistic is the strict/tie shape: `tie_rate` tracks
`safe_gap_count` with mean absolute Spearman `0.881`.  This is the first clean
threshold-local metric from this pass.

The marked stationary score is the second useful channel.  `local_moat_min`
and `safe_deficit` mostly collapse to rank/tie-path tournaments if read through
completed `H`, but their marked origin score tracks origin clearance.  This is
the right way to use scalar vertex gauges: keep the marked vertex, do not
pretend the completed tournament's Hamiltonian paths are the signal.

The pressure gauges are proof-shape gauges.  `k2_relief_pressure`,
`threshold_deficit_pressure`, and `bracket_relief_pressure` should be read by
strict SCCs, directed triangles, and tie rates.  In the selected `n=14` and
`n=18` rows they still peel: strict SCC size stays `1` and strict pressure
triangles stay `0` in the sampled snapshots.  That supports the current
pressure-DAG certificate story rather than a pressure-core disproof story.

The proposed LRC tournament metric vector is:

```text
(phase_H, phase_score_width,
 origin_marked_score, safe_deficit_score_hist,
 origin_blocker_strict/tie_shape,
 pressure_largest_SCC, pressure_triangles,
 strict_tie_rate, edge_flip_rate)
```

This vector is the replacement for a single loneliness scalar.  `H` supplies
the cyclic spread coordinate.  Marked scores supply the stationary `1/n`
coordinate.  Pressure SCCs and triangles supply the counterexample-core
coordinate.

The next useful experiment is a corridor version: compute this vector on every
half-turn cell crossed by the positive lonely interval, not only at the
midpoint and endpoint.
