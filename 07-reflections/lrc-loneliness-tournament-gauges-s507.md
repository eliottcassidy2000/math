# LRC Loneliness Tournament Gauges S507

The user asked for a long session inventing arc definitions for Tournament
Analysis of the Lonely Runner Conjecture until the resulting tournament shapes,
`H` values, score sequences, or related fingerprints become useful loneliness
metrics.

The local starting point was S26/HYP-1968: the half-turn circular tournament is
a sharp meter for a `1/2` gap, but LRC needs `1/(k+1)`.  So I treated the
creative step as the gauge, not the invariant.  `H` can only become a loneliness
meter after the arcs remember the correct endpoint threshold.

During close-out, a concurrent S506b/HYP-1972 arc-criteria metric-vector session
landed upstream.  I kept that broader metric-vector work and renumbered this
benchmark to S507/HYP-1974.  The two threads are complementary: HYP-1972 builds
the wide vector of marked-origin/pressure criteria, while this session asks
which completed tournament fingerprints isolate safe times on small clock grids.

## Gauge Menu

I added `04-computation/lrc_loneliness_tournament_gauges_s507.py` and stored
`05-knowledge/results/lrc_loneliness_tournament_gauges_s507.out`.

The script tests ten gauges on six small speed families over `t=a/840`.
Vertices are the moving runners; the stationary runner enters through the
origin-distance and endpoint observables.

```text
half_turn
theta_close_switch
antipodal_open_switch
unsafe_dominance
safe_dominance
lrc_band_dominance
endpoint_wall_pressure
origin_bracket_flip
same_side_escape
relief_bottleneck
```

For each tournament it records `H`, score sequence, score width, score variance,
directed 3-cycles, SCC count, largest SCC, and the actual normalized loneliness
depth `(k+1) min_i ||s_i t||`.

The evaluation deliberately uses safe isolation, not just raw purity, because
safe times are rare.  A trivial "always unsafe" shadow looks good under raw
accuracy; it looks useless under safe isolation.

## What Worked

The best direct witness detectors are the status gauges:

```text
unsafe_dominance:
  H+score pure fraction 0.996
  safe isolation 0.845
  unsafe isolation 0.866

safe_dominance:
  same performance, reversed status order
```

These gauges are almost tautological, but that is not a defect.  They give the
repo a sanity-check tournament fingerprint for "is this time an LRC witness?"
The score sequence makes the zero-unsafe block visible, while half-turn ties
inside equal-status blocks keep some circular shape.

The best non-tautological shape gauges are:

```text
theta_close_switch:
  best numeric feature score_width, rho -0.266
  mixed-7 safe isolation 0.767

antipodal_open_switch:
  best numeric feature scc_count, rho -0.207
  mixed-7 safe isolation 0.433
```

These are less accurate globally than status dominance, but more interesting:
they do not directly ask whether each runner is safe.  They ask whether
LRC-scale pair distances or near-antipodal openings have changed the tournament
shape.  That makes them better candidates for corridor movies and proof-search
diagnostics.

## What Did Not Work As A Witness Meter

Pure half-turn remains the wrong threshold for LRC witness detection:

```text
half_turn safe isolation across all cases: 0.028
```

This agrees with S26.  Half-turn `H` is still important, but it is a spread
meter, not a `1/(k+1)` endpoint meter.

Endpoint-wall pressure and LRC band dominance also did not classify safe times
well by themselves.  Their value is different: they are boundary/pressure
diagnostics.  Endpoint-wall pressure should spike near
`||s_i t|| = 1/(k+1)`, and slack bands should feed endpoint-core or pressure-DAG
certificate searches rather than serve as final witness classifiers.

## Proposed Definition Stack

For future LRC Tournament Analysis, a good report should not say only "the
tournament has H=...".  It should say:

```text
pairwise observable:
  endpoint status, LRC-scale separation, antipodal opening, wall pressure, etc.

gauge:
  exact arc rule, including threshold and tie Hamiltonian path

fingerprint:
  H, score sequence, score width, c3, SCC shape, edge flips across corridor

target comparison:
  normalized loneliness depth and safe/unsafe isolation
```

The practical shortlist:

1. Use `unsafe_dominance` or `safe_dominance` as the witness sanity meter.
2. Use `theta_close_switch` for LRC-scale crowding shape.
3. Use `antipodal_open_switch` for opening/corridor shape.
4. Use `endpoint_wall_pressure` and `lrc_band_dominance` for boundary pressure,
   not standalone classification.
5. Keep `half_turn` as the `1/2` spread clock, especially for corridor context.

## New Hypothesis

I added HYP-1974:

```text
LRC loneliness metrics need endpoint-aware tournament gauges.
```

The useful object is not a single `H`; it is `(gauge, H, score sequence,
cycle/SCC fingerprints)` with the LRC threshold built into the gauge.
