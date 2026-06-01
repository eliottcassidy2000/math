# LRC Arc-Gauge Zoo S506

The prompt asked for a long session inventing criteria for assigning arcs in
Tournament Analysis of the LRC until the resulting tournament shapes, H values,
score sequences, or related invariants became useful loneliness metrics.

The main outcome is negative in the right way: there is no single tournament
that should be called the LRC loneliness tournament.  There is a useful bundle
of gauges.

## Gauges Tried

The script `04-computation/lrc_arc_gauge_zoo_s506.py` evaluates selected hard
rows for `n=14` and `n=18`:

```text
n14: initial, row-parent, gate, double-gate
n18: initial, row-parent, gate, double-gate
```

At the selected LRC boundary/gap time, it builds complete tournaments from:

- half-turn phase orientation,
- origin-danger and local-moat scalar ranks,
- close-pair and danger-band switches at thresholds `1/n` and `[1/n,2/n)`,
- quarter-turn lead switches,
- open-arc density and open-arc count,
- nearest and two-neighbor deletion-relief pressure,
- threshold-deficit relief,
- observer two-guard and observer relief ranks,
- kinetic escape and kinetic-close switches,
- endpoint cross-protection and endpoint-debt ranks.

For `n=14` the script computes exact Hamiltonian-path counts.  For `n=18` it
keeps score width, 3-cycles, SCCs, source/sink counts, and observer outdegree.
It also runs a deterministic n=14 time-sample correlation audit against the
actual observer margin `min ||v t|| / (1/n)`.

## What Became Useful

The useful gauges separate into roles.

The half-turn phase gauge is the HYP-1970 channel.  It measures global circular
spread.  On n=14 selected rows:

```text
initial       H=24104937
row-parent    H=22168229
gate          H=17826951
double-gate   H=17826951
```

In the n=14 time-sample audit, phase `Hratio` correlates positively with actual
observer margin (`+0.393`) and so do phase directed 3-cycles (`+0.416`).  This
does not make phase H an LRC endpoint meter, but it is a real spread meter.

Open-arc density is the best new bridge gauge.  It asks which direction between
a pair has lower runner density, rather than only whether the pair is separated
by a half-turn.  It distinguishes row-parent and gate distortion:

```text
n14 row-parent  H=2952757, c3=76
n14 gate        H=452385,  c3=60
```

Its c3 correlation with observer margin is `+0.403`.  The related
open-arc-count gauge mostly collapses back to the circular phase menu on the
selected rows, so the metric lengths, not just order counts, are doing work.

Close-pair, danger-band, and kinetic-close switches are local crowding alarms.
They are not monotone spread meters.  For example close-pair has observer-out
correlation `+0.506`, but H/c3 correlations are negative against observer
margin.  That is not a failure: the switch is telling us where local crowding
is active, not where the observer is globally lonely.

The relief gauges are pressure gauges.  Blocker-relief gives very low H on the
hard n=14 gates (`3601` then `161`) and smaller SCCs (`12` then `10`), while
two-neighbor relief keeps more structure (`183523` then `55701`) and stays
closer to the existing S492 pressure searches.  The metric to watch here is
not H alone; it is `largest_scc`, sources/sinks, and whether endpoint-private
debt can be peeled along the pressure order.

Observer-two-guard is a marked-channel gauge.  At selected hard rows it
collapses because the observer is already outside the strict danger threshold,
but across the time audit its observer outdegree tracks actual observer margin
with correlation `+0.506`.  This is the right kind of anchored scalar, but it
needs the other channels to see pairwise obstruction.

## What Failed Usefully

Scalar ranks are readable and bad as tournaments.  `origin_danger_rank`,
`local_moat_rank`, and `observer_relief_rank` usually give `H=1` and zero
3-cycles.  They are ordered ledgers, not Tournament Analysis switchboards.

The runner-level endpoint tournaments also collapse.  Both
`endpoint_cross_protect` and `endpoint_debt_rank` are transitive on every
selected row.  This is important: endpoint data is real, but the speed-level
projection is too coarse.  Endpoint protection should remain labelled by
owner, center, sign, and protector, or be lifted to endpoint-runner pair-cells
before tournament invariants are taken.

The velocity-only `kinetic_escape` gauge is inert on the hard selected rows.
Velocity becomes useful only when tied to a local danger condition, as in
`kinetic_close_switch`.

## Practical Metric

The next LRC Tournament Analysis report should not ask for one number.  It
should record a small vector:

```text
phase_H
open_arc_density score/c3/H
danger_band c3 and score width
close_pair or kinetic_close observer_out
two_neighbor_relief SCC/source/sink profile
observer_two_guard outdegree
labelled endpoint debt and private endpoint layer
```

This bundle has the right semantics:

- phase/open-arc gauges say whether the circle is globally spread or bunched;
- threshold switches say where the `1/n` crowding lives;
- relief gauges say whether the crowding peels like a DAG or closes into a core;
- observer gauges keep the marked LRC bracket visible;
- endpoint labels preserve the finite boundary certificate.

The new HYP-1975 packages this as: LRC loneliness is a gauge bundle, not a
single tournament.
