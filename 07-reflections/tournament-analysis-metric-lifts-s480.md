# Tournament Analysis Metric Lifts S480

The user's phrase "Tournament Analysis" should become literal repo vocabulary.
A tournament is a binary relation on labelled things, but the world rarely hands
us binary relations directly.  It hands us pass counts, distances, angles,
nearest neighbors, forbidden intervals, endpoint invoices, velocities, scores,
and ties.  Tournament Analysis is the act of choosing the binary gauge that
turns that richer object into a tournament, then watching the tournament
invariants move as the original variables move.

The basic object is:

```text
(labels, pairwise/geometric/kinetic data, gauge, tie path) -> tournament
```

The tie path matters.  It is not a nuisance.  It is the same fixed Hamiltonian
path that keeps reappearing in the tiling model: a coordinate system for
ambiguity.

## Basketball Is Not A Toy Example

A five-player pass matrix is already a weighted directed relation.  For each
pair `{i,j}`, compare passes from `i` to `j` against passes from `j` to `i`.
If the counts tie, use the jersey order `1->2->3->4->5` as a Hamiltonian
tiebreak path.

S480's sample matrix has three tied pairs.  After completion it gives:

```text
score=(0,2,2,3,3), scc=2, c3=2, H=5
signature=1110111001
```

The human reading is nice: the pass matrix has too much information, and the
tournament skeleton asks a first discrete question: which pairs have directional
dominance?  The leftover weights are not discarded forever.  They become edge
margins, confidence, flip costs, or future gauges.

This is exactly the active-ranking/OCF instinct in a friendlier setting.  A
basketball team is a weighted comparison graph with a real tiebreak convention.

## Metric Switches Are The Clean Bridge Back To Tiling

The user suggested: for a metric collected pairwise on `N` things, decide a
metric by which to switch each of the `N choose 2` things.

That is a direct fixed-path tiling coordinate:

```text
base edge: lower label beats higher label
switch:    flip edge {i,j} when predicate(d(i,j)) is true
```

The predicate can be a threshold, an annulus, a shell, a sign of a derivative,
a comparison to a local quantile, or a learned classifier.  The important part
is that a symmetric metric becomes a tournament only after a gauge is chosen.
The gauge supplies the orientation.

S480 tested two simple versions:

```text
threshold switch: flip when d(i,j) > theta
annulus switch:   flip when lo < d(i,j) < hi
```

On the runner chord metric, these switch gauges were much richer than vertex
score gauges.  For speeds `(0,1,4,5,6,7,11)` sampled at `t=k/84`, chord row-sum,
nearest-radius, zero-distance, and x-projection stayed transitive.  But chord
threshold and chord annulus produced strongly connected cyclic tournaments with
Hamiltonian-path counts up to `157` and `159`.

That is the small answer of this session: symmetric distance data is not too
poor for tournaments.  It is poor only if we collapse it through a vertex score.
Pairwise switch gauges preserve the `N choose 2` edge-level life.

## Runner Gauges Split Into Two Families

For circular runners, S480 compared:

- `round-half`: orient by clockwise semicircle.
- `zero-distance`: compare each runner's distance from the stationary runner.
- `chord-row-sum`: compare each runner's total chord distance to all others.
- `nearest-radius`: compare each runner's nearest-neighbor radius.
- `chord-threshold`: flip base edges when chord distance exceeds `sqrt(2)`.
- `chord-annulus`: flip base edges when chord distance lies in a middle shell.
- `x-projection`: compare x-coordinate shadows.
- `pull-apart`: orient by instantaneous outward velocity from the pair midpoint.

At the S431 near-tight time `t=181/588`, the scalar/vertex-score gauges are
total orders:

```text
zero-distance, chord-row-sum, nearest-radius, x-projection:
score=(0,1,2,3,4,5,6), scc=7, c3=0, H=1
```

But the pairwise/circular/kinetic gauges retain cycle structure:

```text
round-half:      score=(2,2,2,3,4,4,4), scc=1, c3=11, H=105
pull-apart:      score=(2,2,2,3,4,4,4), scc=1, c3=11, H=105
chord-threshold: score=(2,2,3,3,3,3,5), scc=1, c3=11, H=111
chord-annulus:   score=(2,3,3,3,3,3,4), scc=1, c3=13, H=159
```

So there are two different meanings of "distance between runners":

1. Distance-to-a-mark creates a ranking around the stationary runner.
2. Pairwise metric switches create a cyclic relation among the moving runners.

LRC needs both.  The lonely condition is marked and local at zero; the possible
counterexample mechanism is global and pairwise.

## Sphere, Cube, Simplex

For a point cloud, projection and centroid gauges are usually vertex-score
gauges, so they tend to be transitive.  That is not a failure.  It says those
gauges are ranking gauges, not relation gauges.

Distance threshold and annulus gauges behave differently.  In the S480 six-point
cloud, median-distance and middle-annulus switches produced strongly connected
tournaments:

```text
median-distance switch: score=(2,2,2,3,3,3), scc=1, c3=8, H=41
middle-annulus switch:  score=(1,2,2,3,3,4), scc=1, c3=6, H=23
```

This suggests a useful rule of thumb:

```text
vertex-score gauge  -> order/ranking analysis
pairwise-switch gauge -> tournament/cycle analysis
```

The first is not lesser, but it answers a different question.

## LRC Reframe

The Lonely Runner Conjecture is a marked Tournament Analysis problem.

At time `t`, the configuration `{0} union {v_i t mod 1}` has:

- a marked vertex: the stationary runner at `0`;
- a circular order;
- a two-neighbor bracket around the mark;
- pairwise chord distances among all runners;
- collisions and antipodal ties;
- derivative information from velocities;
- endpoint-protection rows in the finite anti-Bohr model.

The scalar lonely gap is one projection.  S22 and S431 already showed this: a
configuration can be scalar-safe while its pairwise runner cloud is crowded or
has collision/antipodal structure.  S480 adds that chord threshold and annulus
gauges can expose cyclic tournament structure that score-sort gauges erase.

The proof-search question changes from:

```text
Is there a time with distance from zero at least 1/n?
```

to:

```text
What tournament movie is forced when every marked zero-bracket is blocked?
```

That question can see SCC defect, two-nearest cores, good-cut protection,
edge-switch surfaces, and endpoint rows in one language.

## New Questions

1. Which LRC gauges are invariant under scaling speeds, quotienting by gcd, and
   the usual denominator reductions?
2. Do near-counterexample ladders have a characteristic chord-threshold movie,
   perhaps many strongly connected frames but a stable marked-bracket deficit?
3. Can endpoint-protection hypergraphs be converted into tournament gauges by
   orienting rows according to shared protecting columns?
4. Is there a "gauge basis" for Tournament Analysis, analogous to choosing a
   basis of projections, from which other useful tournament movies can be
   reconstructed?
5. Can we define a completeness-defect statistic for partial or tied real-world
   metrics: how much of the tournament cube is genuinely decided by data versus
   by the Hamiltonian tiebreak path?

## Practical Method

Future Tournament Analysis computations should log:

- labels and base Hamiltonian path;
- raw observable type: directed weight, symmetric metric, marked distance,
  derivative, incidence row, projection coordinate;
- gauge predicate;
- tie convention;
- invariant movie: score sequence, SCC count, cyclic triples, Hamiltonian-path
  count, OCF/H where feasible;
- flip events and which pair caused each event;
- marked-vertex data when an LRC-style observer exists.

The repo has long been studying tournaments as finished discrete objects.  This
session reframes the primary process as studying how many different finished
tournaments can be functorially lifted out of one richer phenomenon.
