# The permutohedron is right, but not enough (S526)

The clean geometric picture is now:

```text
time t
  -> points {0, v_i t} on the circle
  -> circular order chamber of the affine braid arrangement
  -> adjacent-swap walk on permutohedral facets
  -> two observer-adjacent gap coordinates
  -> THM-384 source iff both are >= 1/14
```

That is the deep permutohedron geometry.  LRC at `n=14` is not a free walk on
tournament classes; it is a one-parameter rational sweep through this chamber
fan.  HYP-2000 was right to insist that the arcs are not independent.  The
permutohedron is the object that remembers the dependence.

But the naive proof attempt fails.  I wanted the local lemma:

```text
when the last blocker hands off to another blocker,
the sweep must cross the long-long source facet.
```

The S526 audit says no.  The hard rows have many source-avoiding handoffs:

```text
seven-ladder: raw_handoff_defects=142
S380 gate:    raw_handoff_defects=310
double gate:  raw_handoff_defects=600
```

So the permutohedron does not magically turn the n=14 proof into a pure
acyclicity statement.  The blocker-owner graph can have SCCs.  The free
metagraph obstruction from HYP-1997 survives, just in better coordinates.

The useful new invariant is labelled debt.  The hard ladders all keep the same
open source measure:

```text
seven-ladder -> S380 gate -> double gate: source_measure = 1/143
```

while the chamber count and raw handoff defects scale upward.  That is the
geometry of exported debt: the source corridor survives, but the bad handoffs
move to deeper quotient labels.

So the proof target is now sharper:

```text
No closed source-avoiding blocker circulation remains after
endpoint-private leaves and quotient-depth debt labels are charged.
```

This wants a Hall/Farkas certificate over labelled permutohedral handoffs.  The
permutohedron supplies the chamber/facet ledger; endpoint-debt work supplies
the labels that can make the ledger acyclic after peeling.

This dovetails with the parallel HYP-2002/HYP-2003 work: the same object is a
closed geodesic on the permutohedral torus piercing a central box, equivalently
a circle-covering problem by thirteen blocking arcs.  This S526 lane isolates
the extra labelled-handoff debt needed beyond those clean geometric
reformulations.
