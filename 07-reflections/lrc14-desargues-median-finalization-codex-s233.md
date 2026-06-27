# LRC14 Desargues-Median Finalization Lens

The useful part of "median graph" is not the graph class itself.  It is the
discipline it imposes on a proof quotient: if three routes through a packet
fiber are all legitimate, there should be a unique state that lies between
each pair.  That unique state is exactly what controlled forgetting has been
trying to protect under different names: sidecar, payload, owner strip, zeta,
closed-H1 support, period deck, deletion fiber, value origin.

The Desargues graph is the right warning object.  It is bipartite, cubic, and
incidence-shaped, with clean theta-like edge classes in the S233 audit.  Yet
it has `160` triples with no median center.  So a proof graph can look
beautifully organized and still fail the one property we need for final
assembly: three local certificates may have no common center after a quotient
forgets the wrong coordinate.

For LRC14 this suggests a stricter close-out test than "the residual count is
small."  Take a coarse fiber.  Let the vertices be proof states:

```text
topology state
owner state
period state
zeta state
automaton handoff
observer-cut orbit
certificate/debt state
```

Connect states that differ by one retained sidecar or one discharge decision.
Then test triples of routes.  A unique median center means the sidecars are
compatible enough to assemble.  An empty median is a Desargues defect: the
quotient has incidence structure but no consensus proof state.  Multiple
medians mean the vocabulary is too redundant or too coarse to pick a canonical
handoff.

This is a good finalization language because it absorbs the recent work:
HYP-3054 names the cut payload, HYP-3056 orbits it under visible
automorphisms, HYP-3057 tags where the small values came from, and HYP-3067
asks whether those retained coordinates medianize the proof graph.  Pair-good
decoys, residual capacitors, AP-tail clocks, automaton fibers, diagonal-layer
rectangle defects, and matrix observability columns become hyperplanes in the
same compatibility graph.

The next useful computation is not a bigger count.  It is a median-failure
table over HYP-2963 coarse fibers:

```text
route_triple
coarse_fiber
median_center_status
first_missing_sidecar
repair_or_debt
```

If that table empties except AP/GW and named THM-572/F7 debt, we would have a
much more theorem-shaped LRC14 closure target.
