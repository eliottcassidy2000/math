# LRC14 Component-Cover Conductance Router

HYP-3450 made the obstruction local.  HYP-3451 asks what the local obstruction
looks like as a graph.

The useful surprise is that the row with largest paired-cover rank is not the
row that feels closest to a counterexample.  `random_covering_082` has a
rank-`6` dead component, but it also has `106` low-rank escape components.  It
is a good stress test for blocker alphabet size, not the main proof bottleneck.

The bottleneck is still the AP-with-`84m` family.  Those rows have a connected
dead-cover projection and only four low-rank escapes.  This is exactly the
shape a hand proof can attack: the dead part is not fragmented noise; it is a
single connected obstruction graph whose boundary should be forced to leak at
the four survivor windows.

The graph-theoretic proof target is now:

```text
full two-colour blocker saturation of E_safe is impossible.
```

I would try the canonical row first.  Put a signed current on branch-coloured
blockers `B0:o` and `B1:o`, treat even endpoints as boundary gates, and show
that any current saturating all `22` dead components in `{1..11,13,84}` leaves
nonzero boundary at one of the four rank-`2` escape components.  If that works,
the AP-tail lift should preserve the same blocker alphabet while only refining
the even grid.

The information-theory lesson matches the earlier compression failures:
component count, dead mass, and max cover rank are all lossy.  The preserved
payload has to be the branch-coloured blocker incidence graph plus the
endpoint-rank escape labels.
