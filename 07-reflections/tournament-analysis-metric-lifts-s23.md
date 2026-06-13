---
source: oracle-2026-06-01-S23
status: exploratory synthesis
tags:
  - tournament-analysis
  - metric-lifts
  - lonely-runner
  - hamiltonian-paths
  - wall-crossing
---

# Tournament Analysis as metric lifting

The working definition from this session:

```text
Tournament Analysis =
  choose a meaningful binary comparator for every pair,
  complete ties by a fixed Hamiltonian path,
  then study the tournament and its motion as the data changes.
```

The basketball example is the clean discrete model.  The five starters are the
vertices.  For each pair `i,j`, orient `i -> j` if `i` passed to `j` more than
`j` passed to `i`.  If the counts tie, use the label path
`1 -> 2 -> 3 -> 4 -> 5`.  That is exactly the repo pattern: a noisy or
degenerate relation is completed by a labelled Hamiltonian path.

The continuous model is runners or point clouds.  The raw data may be positions
on a circle, chord distances, cuboid `L_p` distances, simplex distance-profile
entropy, or oriented volumes.  The important question is not "what is the
metric?" by itself.  The important question is "what comparator turns that
metric into one bit for each pair?"

## Rankers versus analyzers

The computation makes a sharp split.

Rank lifts assign one scalar score to each object:

```text
distance to anchor
row-sum distance centrality
nearest-cell radius
distance-profile entropy
instantaneous receding speed
```

Then `i` beats `j` iff `s_i > s_j`.  These lifts are rankings.  When scores are
distinct, they are transitive tournaments.  In the S23 samples, score lifts were
transitive in all `316/316` continuous states.

Analyzer lifts decide each edge from pair-specific data:

```text
flux lift:       directed count i->j versus j->i
phase lift:      circular half-turn or chirality
lens lift:       compare i,j through shifted labels or third vertices
switch lift:     symmetric distance D_ij toggles a base-path edge
volume lift:     oriented area or determinant sign
dynamic lift:    edge flips across comparator walls
```

These preserve cycles.  In the S23 samples, edge-local and edge-switch lifts
were transitive in only `1/290` continuous states.

This is the main methodological point: a metric is not automatically useful for
tournament theory.  If it becomes a scalar ranking, it has thrown away most of
the cyclic information.  To get OCF/Hamiltonian-path structure, the comparator
must usually be edge-local or it must switch individual base-path edges.

## The switch lift

The user's "switch each of the N choose 2 things by a metric" is the switch
lift:

```text
start with fixed path 1 -> 2 -> ... -> n
compute a symmetric pair metric D_ij
if D_ij is in the chosen active band, keep i -> j
otherwise reverse it
```

The active band can be median threshold, high/low threshold, annulus, resonance
like `sin(k*pi*D_ij)`, or any other meaningful boolean function.  This directly
turns continuous geometry into a path through the tournament cube
`{0,1}^{binom(n,2)}`.

This is different from a score lift.  The edge `(i,j)` is not decided by
`s_i-s_j`; it is decided by the pair's own metric value.

## Geometry examples

On the circle, the half-turn phase lift and pivot-area lift produce exactly the
same tournament path.  That is not a coincidence: both are signs of the same
oriented `S^1` chord relation.

On the same circle, anchor arc distance and anchor chord distance also produce
the same tournament path, because chord length is monotone in circular distance
from the anchor.  This is a useful negative result: changing from arc to chord
does not necessarily create a new tournament lens.

On the cuboid samples, `L1`, `L2`, `Linf`, and entropy score lifts all remain
rankers, hence transitive.  But label-shift `L2`, median-switch `L2`,
`Linf`-resonance, `xy` area, and drift-volume are analyzers.  They generate
cycle-rich tournament paths with broad Hamiltonian path ranges.

## LRC connection

S22 showed that LRC loneliness is a marked two-neighbor bracket around the
stationary runner.  S23 adds that there are two different uses of metric data:

1. The bracket is a metric witness condition.
2. Tournament structure appears when the same positions are lifted by phase,
   lens, switch, or volume comparators.

For the active LRC witness rows, anchor-distance and local-cell score lifts are
transitive.  The phase, label-shift, and median-switch lifts remain highly
cyclic.  Example:

```text
n16 d=2: phase 162 cycles, label-shift 158, switch 161
n16 d=8: phase 168 cycles, label-shift 164, switch 156
```

So the right abstraction is not "distances imply a tournament" in one unique
way.  The abstraction is a stack:

```text
positions -> metrics -> comparators -> tournaments -> invariants -> wall crossings
```

The comparator is the creative step.

## What to try next

The next useful program should treat comparator families as first-class
objects.  For a fixed moving point system, scan:

```text
threshold switches
annulus switches
resonance switches
anchor-pair lenses
oriented volume lenses
nearest-neighbor lenses
```

Then cluster the resulting tournament paths by equality, flip distance, and
OCF statistics.  That would show which continuous metrics are genuinely new
Tournament Analysis lenses and which are the same path in disguise.
