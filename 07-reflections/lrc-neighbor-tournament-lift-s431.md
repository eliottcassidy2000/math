---
source: codex-2026-05-31-S431
status: exploratory synthesis
tags:
  - lonely-runner
  - tournament-structure
  - nearest-neighbor
  - round-tournaments
  - distance-data
  - endpoint-protection
---

# LRC Neighbor Tournament Lift

The main reframing from S431 is:

```text
The lonely gap is the view from one vertex, not the configuration.
```

At time `t`, the LRC speed set gives the stationary runner plus the moving
positions

```text
0, v_1 t, ..., v_k t  mod 1.
```

The usual LRC scalar keeps only the nearest distance from `0`.  But the whole
configuration carries more structure:

1. a danger-rank tournament by distance from `0`;
2. a circular half-distance tournament, with ties/collisions making it
   incomplete;
3. a two-nearest-neighbor digraph;
4. a pairwise distance ledger, equivalently a difference-speed LRC shadow.

That is a much better match to the tournament side of the repo.  Tournaments
have taught us repeatedly that projected scalars are dangerous: `H`, score
sequence, or good-cut count matter, but the proof usually lives in residues,
private pivots, SCCs, and incidence layers.  LRC seems to be the same story.

## What The Computation Shows

`lrc_neighbor_tournament_lift_s431.py` audits seven representative speed sets
at exact witness times inherited from S356.

The cleanest lesson is that scalar equality can hide very different
configuration types.

```text
initial k=4      at t=1/5: regular 5-gon
sporadic n=5     at t=1/5: same regular 5-gon, relabelled
initial k=5      at t=1/6: regular 6-gon with antipodal half-tournament ties
sporadic n=6     at t=1/6: scalar-tight, but v3 and v9 collide
```

So the sporadic example is not mysterious in the scalar picture; it is a
quotient-polygon picture.  The scalar says "lonely."  The pairwise ledger says
"two moving runners are identified."  That is exactly the kind of residue the
tournament machinery is good at remembering.

The near-tight examples are even more suggestive:

```text
near-tight k=6:      gap/th = 1.071429, pair/th = 0.535714
n14 seven-ladder:    gap/th = 1.208333, pair/th = 0.169913
n14 exported-debt:   gap/th = 1.173972, pair/th = 0.053301
```

They are not near regular packings of all runners.  They are safe for the
stationary runner while the moving runners crowd each other.  That explains
why tiny visible LRC gaps can be misleading: the safe valley can be narrow
while the global circular configuration is extremely compressed elsewhere.

## The Two-Nearest Move

The user's suggestion to track two nearest neighbors feels exactly right.
Oracle S22 independently reached the sharper marked-bracket language: around
the stationary runner, the predecessor and successor gaps are the two local
quantities that matter.  S431 is complementary: it keeps the full pairwise
configuration and the round/incomplete tournament shadow around that bracket.

Tracking one nearest runner to `0` gives a lower envelope.  Tracking two
nearest runners gives a bracket.  At a positive-gap local maximum, if the two
nearest threats to `0` were on the same side, a small motion should improve the
minimum distance.  Thus a genuine local obstruction wants a left/right bracket:

```text
nearest left runner  < 0 <  nearest right runner
```

That bracket is a circular cut.  Now compare THM-354:

```text
tournament good cut:  a backward arc crosses a Hamiltonian-path cut
LRC zero bracket:     a protector interval crosses a circular zero cut
```

This suggests a new proof technology:

```text
Build a two-nearest protection graph.
Peel private/nonreciprocal contacts.
A counterexample must leave a nonempty arithmetic two-nearest core.
```

This should sit between the scalar endpoint-core program and the full
endpoint-protection hypergraph.  It keeps more geometry than THM-357 but is
still much smaller than the full pairwise arrangement.

## Tournament Methodology Translation

Here is the methodology shift I would keep for future LRC work.

```text
Old route:
  speed set -> forbidden interval union -> max scalar gap

Lifted route:
  speed set -> movie through round/incomplete tournament strata
            -> zero-bracket cuts
            -> two-nearest contact graph
            -> endpoint/protector incidence labels
            -> scalar gap as one projection
```

The circular half-distance tournament is not arbitrary.  It is a round
tournament when there are no antipodal ties or collisions.  That restriction
is useful: it means LRC is not asking for generic tournament theory, but for
tournament theory inside a rigid cyclic realization space.  The interesting
defects are precisely the ways that realization degenerates:

```text
collision       = quotient identification
antipodal tie   = missing/incomplete arc
crowding        = small difference-speed distance
zero bracket    = circular cut around the stationary vertex
private contact = peelable endpoint or nearest-neighbor leaf
```

## New Questions

1. Can the two-nearest graph be peeled in a way that always empties on known
   tight and near-tight examples, just like endpoint cores emptied in S359?
2. Does the half-distance round tournament at a safe time have a score defect
   or SCC condensation invariant predicting endpoint peel depth?
3. Can `pairwise_min_gap / threshold` become a danger filter, separating true
   global packing obstructions from scalar tiny-gap mirages?
4. Is every sporadic tight example a regular residue polygon after quotienting
   colliding moving runners?
5. Can THM-354 be ported to a circular bracket theorem: bad zero-cuts are
   exactly boundaries in a condensation of the two-nearest protection graph?

The useful slogan:

```text
LRC is not only a covering problem on the time circle.
It is a one-parameter walk through cyclic tournament configurations,
observed through the stationary vertex.
```
