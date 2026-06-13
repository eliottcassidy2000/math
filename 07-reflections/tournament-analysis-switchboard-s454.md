---
source: codex-2026-06-01-S454
status: exploratory synthesis
tags:
  - tournament-analysis
  - switchboards
  - metric-lifts
  - lonely-runner
  - basketball
  - wall-crossing
---

# Tournament Analysis Switchboards

The user phrase "switch each of the `N choose 2` things" is the right primitive.
S23 defined Tournament Analysis as:

```text
data -> metric/sensor -> comparator -> tournament -> invariants.
```

S454 adds the missing middle object:

```text
switchboard = all pairwise bit-channels c_ij(x), plus a fixed tie-break path.
```

The tournament at one time is a snapshot.  The switchboard is the moving
instrument panel that produces all snapshots.

## Rankers Summarize Objects

A rank lift gives each object one scalar `s_i` and orients

```text
i -> j iff s_i > s_j.
```

This is useful, but it is usually a total order in disguise.  Across the S454
samples, rank shadows were transitive in every state:

```text
172 / 172 rank states transitive.
```

That is not a bug.  It says the lift is answering an ordering question, not a
pairwise-relation question.

## Switchboards Summarize Relations

A switchboard gives each pair its own channel:

```text
pass_i_to_j - pass_j_to_i
D_ij - median(D)
sin(k*pi*D_ij)
oriented_area(i,j,anchor)
KL(P_j||P_i) - KL(P_i||P_j)
```

These are not scalar rankings.  They can preserve cycles because each edge has
its own reason for pointing one way.

Across the S454 analyzer samples:

```text
31 / 672 analyzer states transitive,
mean live edges per path = 22.95 out of 28 for the 8-vertex families.
```

That is the Tournament Analysis zone: many pair channels are actually moving,
and the motion is cyclic enough for Hamiltonian-path and OCF machinery to see.

## Basketball Is the Friendly Model

The basketball example is not just illustrative.  It is a clean human
switchboard.

Five starters give ten pair channels.  For each pair:

```text
pass-flux:          i -> j if i passes to j more than j passes to i
assist-rank:        i -> j if i has higher net outgoing pass balance
reciprocity-switch: total pair activity toggles the labelled base path
two-hop-lens:       i -> j if i creates more two-hop paths through teammates
pressure-switch:    pair imbalance toggles the base path
```

The same raw data produces different tournament paths.  This is the point:
there is no canonical tournament until we choose the comparator.  The fixed
Hamiltonian tie-break path is a real layer, analogous to the repo's labelled
base path in tiling models and lex completions of degenerate circular states.

## Continuous Switches

For runners on a circle, chord distance by itself is symmetric, so it cannot
orient an edge.  But it can switch an edge:

```text
keep i -> j from the base path if chord(i,j) is in the active band,
otherwise reverse it.
```

The active band can be a median threshold, annulus, resonance, approach rate,
or anything else that has mathematical meaning for the problem.  This turns a
moving geometric system into a path through the tournament cube.

The higher-dimensional examples confirm the same pattern:

```text
cuboid dominant-axis lens: a pair chooses its own coordinate judge
sphere normal lens:       a pair orients by the normal to its great-circle arc
simplex KL-skew lens:     information geometry supplies an asymmetric edge
```

The important move is not "use a metric."  The important move is to declare
which pairwise question the metric is answering.

## LRC Reading

The lonely-runner condition is still a marked two-neighbor bracket around the
stationary runner.  Tournament Analysis adds surrounding pairwise shadows.

In S454, the safe-distance switch is transitive on the initial skeletons but
cyclic on structured ladders:

```text
n14 initial:      safe-switch cycles 0
n14 seven-ladder: safe-switch cycles 93
n16 initial:      safe-switch cycles 0
n16 dyadic:       safe-switch cycles 139
```

This suggests a useful split:

```text
marked bracket:    witnesses loneliness or failure at the stationary vertex
pair switchboard:  records global crowding, antipodal debt, and cyclic pressure
```

The scalar LRC gap can be positive while the pairwise switchboard is already
very cyclic.  That is exactly the kind of hidden structure the repo keeps
finding: scalar safety does not mean relational simplicity.

## Metrics To Keep

Future Tournament Analysis runs should log:

```text
switchboard type
live edge count
wall count / total flips
maximum step flip
Hamming diameter
Hamiltonian-path range
cyclic-triple range
SCC range
exact tournament-path clusters
```

These are cheap but meaningful.  They tell us whether a comparator is merely a
ranker, whether two geometric metrics are the same lens in disguise, and where
the wall crossings concentrate.

The motto from this pass:

```text
rankers summarize objects;
switchboards summarize relations.
```

For this repo, that second line is the main line.
