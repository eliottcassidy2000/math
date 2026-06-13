---
source: codex-2026-06-01-S471
status: framework synthesis with computational probe
tags:
  - tournament-analysis
  - pairwise-metrics
  - lonely-runner
  - basketball
  - cuboid
  - simplex
---

# Tournament Analysis: Pairwise Data as a Binary Relation

This session names the repo's central method:

```text
Tournament Analysis =
  pairwise data + switching rule + tie Hamiltonian path
  -> tournament-valued observable
  -> tournament fingerprints and trajectories.
```

A tournament is a binary relation on `n` labelled objects.  A pairwise metric
or measurement gives one datum on each of the

```text
C(n,2)
```

edges.  Tournament Analysis asks how to turn each datum into an orientation,
then studies the resulting point in the tournament cube.

The important extra choice is the switch functional.  There is no single
correct switch.  Different switches expose different geometry.

## Tie-Breaks Are Base Hamiltonian Paths

The basketball example is perfect.  A lineup has labelled roles

```text
1,2,3,4,5
```

and a game gives directed pass counts `P_ij`.  Orient

```text
i -> j iff P_ij > P_ji.
```

If the counts tie, use the role path

```text
1 -> 2 -> 3 -> 4 -> 5
```

as the tie-break.

This is exactly a repo object.  The fixed role order is a Hamiltonian path.
The staircase/tiling model also starts from a base Hamiltonian path and then
flips selected pair orientations.  Basketball labels are not cosmetic; they
are the convenient path that turns a relation with ties into a tournament.

The script's synthetic pass examples show three regimes:

```text
point-hub:     transitive tournament, H=1
motion-ties:   one directed 3-cycle, H=3
inverted-big:  transitive tournament, H=1 with the big/wing side on top
```

The pass matrix is discrete, but the abstraction is identical to the
continuous runner case.

## Symmetric Metrics Need A Switch

A distance `d_ij=d_ji` does not orient an edge by itself.  There are several
natural ways to switch it:

### 1. Threshold cuboid switch

Fix a base path.  For `i<j`, set

```text
i -> j  if d_ij >= theta_ij,
j -> i  if d_ij < theta_ij.
```

This turns the symmetric metric into a vertex of the tournament hypercube
`{0,1}^{C(n,2)}`.  It is the cleanest interpretation of "each pair metric
switches one of the `C(n,2)` things."

For circle runners, a chord threshold can use

```text
theta = 2 sin(pi/n),
```

the chord form of the LRC distance `1/n`.

In the S471 computation, the initial five-runner chord-threshold trajectory
sampled at `t=q/60` visits:

```text
12 distinct tournaments,
H spectrum {1,5,9,11,13,15},
edge flip counts increasing with speed difference.
```

That is a literal tournament-cube trajectory driven by a continuous metric.

### 2. Semicircle switch

For points on a circle or sphere with a chosen oriented great circle, orient

```text
i -> j
```

when `j` lies in the forward open half from `i`.

This is the round-tournament lift already used in the LRC sessions.  It tends
to preserve cyclic structure: regular five-runner boundary data gives
`H=15` and five directed triangles.

### 3. Opening/closing switch

For moving points, use the derivative of a distance.  In the circle case, the
chord distance between runners can be opening or closing.  A useful rule is:

```text
if the chord is opening, faster endpoint wins;
if the chord is closing, slower endpoint wins.
```

This is not a static metric rule.  It is a dynamical relation: who is creating
or absorbing separation?

In the S471 audit, chord-opening and Fourier-opening often keep high-H cyclic
structure while centrality rules collapse to transitive order.

### 4. Isolation or heat switch

Compute a vertex scalar from pairwise distances:

```text
isolation_i = sum_j d_ij,
heat_i      = sum_j exp(-lambda d_ij^2).
```

Then orient more isolated over less isolated, or less crowded over more
crowded.  These switches are meaningful but frequently transitive, because
they reduce the metric to a one-dimensional ranking.  That is useful as a
warning: some Tournament Analysis choices throw away the same kind of
incidence data that scalar LRC methods lose.

### 5. Simplex stress switch

View the `n` objects as vertices of a weighted simplex.  For edge `{i,j}`,
delete that edge and compare the remaining incident stress:

```text
stress_i^(ij) = sum_{k notin {i,j}} d_ik.
```

Orient toward the endpoint with larger residual stress.  This is a
metric-simplex analogue of private-pivot thinking: the edge is decided by how
the two endpoints sit relative to the rest of the simplex.

### 6. Pressure switch

Use the two nearest neighbors.  For a pair `{i,j}`, measure deletion relief:

```text
relief_i(j) = nearest_i(after deleting j) - nearest_i.
```

Orient the irreplaceable blocker toward the runner it blocks.  This is the
mobile analogue of endpoint handoff protection.

## Runner Metrics Beyond The Circle

For LRC, the circle is the native object:

```text
p_i(t) = v_i t mod 1.
```

But Tournament Analysis can lift the same phases into richer spaces.

### Chords

Map each runner to the unit circle and use chord length

```text
d_ij = 2 sin(pi ||(v_i-v_j)t||).
```

This keeps LRC geometry but makes the metric Euclidean.

### Fourier sphere / torus

Map phase `x` to

```text
(cos 2pi x, sin 2pi x, cos 4pi x, sin 4pi x, ...)
```

and use Euclidean distance or distance derivative in that higher-dimensional
embedding.  This separates first-harmonic closeness from higher-harmonic
resonance.  In the audit, Fourier-opening changes `H` from the plain
semicircle rule while preserving the same broad cyclic score shape.

### Cuboid

Every thresholded pair distance is a bit.  All pair thresholds together place
the configuration in

```text
{0,1}^{C(n,2)}.
```

This is the direct cuboid/tournament-hypercube view.  The LRC chamber walk is
then a path through a cuboid of edge states.

### Simplex

The complete pairwise metric is an edge weighting of the simplex `K_n`.
Stress, cut load, and deleted-edge comparisons are simplex-native
tournamentization rules.  They are natural companions to the repo's existing
simplex/cuboid duality language.

## What The S471 Audit Shows

The script is:

```text
04-computation/tournament_analysis_framework_s471.py
```

It compares basketball pass matrices and runner configurations under nine
switch rules.

The central pattern is:

```text
centrality-like switches -> often transitive, H=1
cyclic/geometric switches -> many directed triangles, large H
threshold cuboid switches -> sparse H spectra and edge-flip trajectories
pressure switches -> peel/leaf information, often low H but not scalar-only
```

For initial `n=14` at `t=1/14`:

```text
semicircle:      H=24104937, 3cyc=112, score 6^7 7^7
chord-opening:   H=24104937, 3cyc=112, score 6^7 7^7
fourier-opening: H=24577317, 3cyc=112, score 6^7 7^7
centrality rules: H=1
```

For the `n=14` seven-ladder boundary:

```text
semicircle:      H=19622601, 3cyc=108
chord-opening:   H=21934885, 3cyc=110
fourier-opening: H=20727933, 3cyc=108
pressure:        H=5,        3cyc=2
```

So the choice of switch functional matters.  It is not just a visualization
choice; it selects which obstruction layer is visible.

## Reframing The Repo

The repo has been doing Tournament Analysis all along:

```text
tournament H        = binary relation -> Hamiltonian path count
tiling model        = base path + edge flips in a cuboid
OCF                 = cycle-packet compatibility graph
LRC endpoint work   = labelled protection relation
runner-distance work= continuous metric -> tournament trajectory
Zeckendorf bridge   = path/cycle independent-set normal forms
```

The unifying move is:

```text
Do not ask for one scalar summary.
Choose a meaningful binary relation, then study its tournament shadow.
```

The strongest future direction is to build a catalogue of switch functionals
and ask which ones preserve the residue that matters for each problem:

```text
passes       -> role path, hub/cycle offense, model-misspecification ties
circle LRC   -> bracket gaps, endpoint handoffs, pressure leaves
chords       -> metric cuboid flips and opening/closing dynamics
spheres      -> harmonic resonance tournaments
cuboids      -> threshold edge-bit trajectories
simplexes    -> stress and deleted-edge comparison tournaments
```

Tournament Analysis is therefore not a metaphor.  It is the repo's core
method for making pairwise structure computable, comparable, and portable
between discrete and continuous problems.
