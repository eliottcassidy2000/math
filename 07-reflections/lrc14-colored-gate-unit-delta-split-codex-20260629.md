# LRC14 Colored Gate Unit-Delta Split Reflection

The hidden structure behind HYP-3471 is that the local color alphabet is more
rigid than the headline says.  "Every dead row has an E/branch gate" still
sounds loose.  After splitting the bank into structured rows and random rows,
the minimum gate nearly freezes.

On the structured side the minimum gate is not just low-rank, and not just
E/branch.  It is always a unit-delta edge gate, and in fact it sits in a tiny
typed palette:

```text
B1:7|E:0
E:0|B1:5
E:0|B0:5
```

That is a real theorem-shaped packet.  The structured branch is not asking for
another coloring story or another AP84 sidecar.  It is asking for a direct
proof that the local boundary can only leak through those unit-delta gates.

The random side is also cleaner than it first looked.  The residual is not
"wild random geometry."  There is exactly one both-branch unit-delta row, and
the true non-unit residual is only `19` rows.  Even those rows never leave the
same local skeleton: they stay single-bad-edge and single-branch.  What
changes is only the adjacent cover delta, with the tiny menu

```text
(1,2), (2,1), (2,2), (1,3).
```

So the correct slogan is:

```text
HYP-3471 gives the boundary alphabet.
HYP-3472 shows that almost all minima are the unit-delta dialect.
The residual is a finite delta-sidecar accent, not a different language.
```

This matters because it changes what "prove the finite lemma" should mean.  A
bad formulation would still ask for an all-purpose colored gate theorem.
Better is:

1. prove the unit-delta packet directly on the structured branch;
2. isolate the lone both-branch row if it has a special symmetry;
3. route the `19`-row delta-sidecar packet through HYP-3451 or HYP-3455.

That is a much smaller frontier than "understand all colored survivor gates."
The controlled-forgetting lesson is the same one that keeps recurring in this
repo: quotient until the theorem predicate is almost rigid, then name the
small surviving packet honestly instead of pretending it has disappeared.
