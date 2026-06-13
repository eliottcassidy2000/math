# LRC Pressure DAG Audit S500

The phrase "pressure searches returning DAGs" is stronger than it first sounds.
It means the search is not finding a cyclic obstruction, but it may be finding
an induction order.

The pressure relation is the S470 two-neighbor relation:

```text
relief_i(j) = nearest_distance_i(after deleting j) - nearest_distance_i
```

Orient `j -> i` when deleting `j` helps `i` more than deleting `i` helps `j`.
In words: `j` is the more irreplaceable blocker of `i`.

If this directed graph has an SCC, then we have a pressure knot: every runner in
the core is mutually implicated in someone else's local moat.  That would be a
counterexample-shaped object.  If the graph is a DAG, the relation is telling us
the opposite: the dependency can be topologically sorted.

## S500 Result

S500 audited the difficult rows from the recent n=14/n=18 feedback loop:

```text
n14-d7:   1177/1177 pressure samples were DAGs
n14-d14:  2349/2349 pressure samples were DAGs
n18-d3:    100/100 pressure samples were DAGs
n18-d9:    233/233 pressure samples were DAGs
n18-d18:   425/425 pressure samples were DAGs
```

Total:

```text
0 cyclic pressure samples out of 4284
```

The maximal dependency heights were small but nontrivial:

```text
n14-d7:   height 5, longest path 4
n14-d14:  height 5, longest path 4
n18-d3:   height 4, longest path 3
n18-d9:   height 4, longest path 3
n18-d18:  height 5, longest path 4
```

So the hard ladders are not flat.  They have dependency cascades.  But the
cascades do not close.

## Why This Matters

A DAG has sources and sinks.  In the pressure language:

- sources are blockers that can be charged first;
- sinks are runners whose moat is most dependent on others;
- layers are a peel schedule;
- height is the number of dependency rounds before the structure empties.

This gives a new interpretation of failed disproof searches.  When a pressure
search returns a DAG, it has returned a partial proof object.  The right next
step is not to say "no cycle found" and move on.  The right next step is to
extract the topological order and match its sources to endpoint-private rows.

## Not Merely A Scalar Sort

S500 also checked a simple potential: do pressure edges just descend nearest
distance `d1`?  They often do, but not perfectly:

```text
n14 rows: about 69-70% d1-descent
n18 rows: about 72-76% d1-descent
```

That matters.  The DAG is not just the transitive order by one scalar.  It is a
relation-level object.  This fits the broader Tournament Analysis lesson:
rankers collapse; analyzers preserve structure.

## Proof Shape

The proof architecture suggested by S500 is:

```text
endpoint-private row + pressure source -> peel
repeat through topological layers
if no source exists -> pressure SCC
pressure SCC -> labelled arithmetic endpoint-cycle contradiction
```

The last line is where THM-365-style labelled endpoint cycles should enter.
The first line is where THM-359 endpoint-core peeling enters.

Thus pressure DAGs may be the bridge between endpoint incidence and mobile
runner geometry.

## Next Concrete Test

Run a transitive-reduction pass on the pressure DAGs.  The transitive reduction
should expose the minimal dependency arrows.  Then compare:

```text
pressure source/sink labels
endpoint private owner labels
denominator depth / product-tree depth
```

If the source labels consistently hit private endpoint rows, the DAG is already
half of a dual certificate.
