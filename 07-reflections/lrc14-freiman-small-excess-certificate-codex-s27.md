# LRC14 Freiman Small-Excess Certificate

**Source:** codex-2026-06-19-S27, HYP-2638 / T886.

The useful sharpening from this session is that the first additive-energy
pocket is no longer heuristic.  Freiman `3k-4` turns exact sumset excess
`exc(E)<=k-3` into a finite problem: after translation and dilation, every row
lives in `[0,2k-4]`, and the stronger hull check is `max(E)+1<=k+exc(E)`.

The exact enumeration passed for `k=8,9,10`:

```text
k=8:  225 small-excess rows, best positive L_y = 297/980
k=9:  620 small-excess rows, best positive L_y = 38681/79380
k=10: 1644 small-excess rows, best positive L_y = 3307/5880
```

All three best positive-excess rows are below the AP and below the relevant
cap.  The tight row is `k=9`, `E=(0,1,2,3,4,5,6,7,9)`, with cap gap
`39541/5675670`.  That is small enough that the near-AP finite table should
stay exact in any proof writeup.

The new guardrail is that excess is not a monotone scalar.  Later excess layers
can rise again.  The correct claim is finite certification by exact layer
maxima, not "larger excess implies smaller `L_y`."  This keeps the KPS S12
Freiman-dimension slogan honest: dimension is the right architecture, but
the first rigorous pocket is a finite theorem plus a table.

The incoming KPS S12 dimension-pocket note is compatible and useful, despite
the HYP-2637 namespace collision.  Its higher-level picture is:

```text
AP / dimension 1 -> exact top row
small excess -> finite 3k-4 pocket
higher-rank GAP -> dimension penalty still to prove
large doubling -> peel or independent-limit tail
```

So the next sharp obligation is not the near-AP small-excess side.  It is the
region `exc>k-3`, especially the HYP-2635 wide examples that have no
dissociated stranger.  Those should be routed through the relation-fiber split:
uncovered relation vertex gives a peel, while full bounded relation coverage
requires a signed higher-rank GAP/lattice tail bound.

Tournament Analysis used proof-obligation vertices:

```text
Freiman_3k4_pocket
> exact_excess_layer
> AP_invariance
> relation_fiber_cover
> higher_rank_GAP_tail
> raw_pair_sum_energy
> raw_runner_vertices
```

The challenged assumption remains important: ordinary pair sums and raw
runners are too thin.  They become useful only after the Freiman pocket has
preserved the exact sector predicate, and they do not retain the multiplicative
sign information needed by the reciprocal-tail side.
