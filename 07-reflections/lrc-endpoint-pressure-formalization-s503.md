# LRC Endpoint-Pressure Formalization S503

This session tightened the pressure-DAG story into a finite theorem package.
The endpoint core is the arithmetic/topological witness layer: endpoints are
rows, forbidden intervals are protectors, and a no-lonely open cover has to
leave every terminal endpoint strictly protected.

The pressure graph is the geometric Tournament Analysis layer.  Its vertices
are runners, and its arrows are chosen from a mobile pairwise gauge such as
nearest-neighbor deletion relief, two-neighbor relief, or threshold-deficit
relief.  Earlier sessions found that the hard `n=14` and `n=18` ladders keep
returning DAGs in these gauges.  S503 explains why that matters formally.

## Owner-Cycle Lemma

THM-379 proves the abstract engine:

```text
nonempty owner-compatible endpoint core
+ strict non-self protection
=> owner-protection digraph has a directed cycle.
```

The proof is deliberately small.  In a nonempty core, every active endpoint
owner has some strict protector owned by another active owner.  Therefore every
active owner has indegree at least one in the owner-protection graph.  A finite
digraph with indegree at least one at every vertex contains a directed cycle.

For LRC, the no-self part is automatic: a speed's own endpoints are boundary
points of its own forbidden intervals, not interior points.

## Certificate Trilemma

THM-380 packages the proof workflow:

```text
denominator sieve -> endpoint core -> pressure SCC
```

A full-open-cover counterexample at threshold `1/n` must pass all three tests:

1. it must be small-denominator sieve-complete, by THM-366/THM-369;
2. its terminal endpoint protection core must be nonempty, by THM-357/THM-359;
3. if the core's strict protection incidences are represented in the chosen
   pressure lift, the active pressure graph must contain a directed cycle.

Thus a proof search can now stop for three independent reasons:

```text
missing denominator -> rational lonely witness
empty endpoint core -> no all-protected open cover
pressure-realized core + pressure DAG -> no counterexample core
```

The remaining counterexample-shaped object is much narrower:

```text
sieve-complete
+ nonempty terminal endpoint core
+ labelled pressure SCC
```

or else a failure of pressure realization that points to a missing gauge.

## Relation To Tournament Analysis

This is a clean Tournament Analysis bridge.  Endpoint protection is not a
tournament by itself; it is an incidence hypergraph.  Pressure is not the LRC
proof by itself; it is a pairwise binary shadow.  The bridge is owner
projection.  Once every protected endpoint can choose a strict protecting
owner, the incidence core casts a directed owner graph.  Any pairwise pressure
gauge that contains those owner arrows inherits the cycle obstruction.

That makes pressure DAGs positive evidence.  They are not just failed SCC
hunts.  A pressure DAG says: either the endpoint core is already gone, or the
current pressure gauge is not yet realizing the core incidences.  The next
formal target is to prove realization for specific gauges, or to refine the
arc labels until every terminal endpoint protection arrow is represented.

## Computation

`04-computation/lrc_endpoint_pressure_formal_s503.py` is a sanity audit for
the finite graph engine.  It exhaustively checks loopless digraphs through
four vertices and functional protector selectors through eight active owners.
Every all-owned strict-protector shadow has a directed cycle.

The stored output is:

```text
05-knowledge/results/lrc_endpoint_pressure_formal_s503.out
```

## Next Proof Moves

1. For each pressure lift (`k1`, `k2`, deficit), define exactly when a strict
   endpoint protector produces a pressure edge.
2. Label pressure arcs by endpoint rows or endpoint-runner cells before
   projecting to runner-level pressure.
3. Prove that terminal endpoint cores are pressure-realized in one of these
   labelled lifts, or identify the missing gauge component.
4. Re-run `n=14` and `n=18` hard rows with core-incidence labels attached to
   every pressure DAG source/sink layer.
