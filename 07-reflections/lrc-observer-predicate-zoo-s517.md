# LRC Observer Predicate Zoo (S517)

HYP-1981 is stronger than it first looks.  It does not merely say "source is a
nice metaphor for loneliness."  It gives a marked tournament coordinate whose
whole score stratification is arithmetic:

```text
observer indegree = number of runners inside the forbidden LRC caps.
```

So the source target has neighbors with exact meanings:

```text
source          zero blockers
almost-source   one blocker
k-score layer   k blockers
```

This is useful because a proof rarely jumps straight from arbitrary bad state
to source.  The natural near-target is the almost-source layer: one runner is
blocking, and one observer incident edge flip would make the witness.  The next
label is the side defect: is the blocker in the left forbidden cap, the right
forbidden cap, or tied with the observer on a wall?

The other creative predicate is observer 2-king.  In the marked tournament, the
observer is a 2-king when every blocker is beaten by some currently safe runner.
That is not LRC, and it should not be mistaken for LRC.  But it looks like a
repair or pressure predicate: the safe set can reach every current obstruction
in one runner-runner half-turn step.

The audit confirmed the hierarchy:

- Marked classes never mixed blocker counts in the bounded windows.
- Runner-subtournament classes did mix blocker counts.
- Initial systems stayed wall-only: open cells had one blocker, walls had zero.
- Wall-completed source menus are larger than open source menus, echoing
  THM-383's boundary warning.

The phrase I would keep is "source-distance fiber."  HYP-1981 names the target;
THM-385 names the distance from the target; HYP-1988 proposes the repair labels
that may explain why source-avoiding walks are impossible.
