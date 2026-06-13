# n=18 observer-source gate battlefield after S520

The instruction was to take THM-382 and THM-383 as the big inspiration, then
to use HYP-1981 as a major source too.  The natural move is to stack them:

```text
THM-382: threshold-decorated fibers are the separating LRC quotient.
THM-383: equality walls and tie completion must be compactified.
HYP-1981: an LRC witness is exactly an observer-source marked tournament.
THM-384/HYP-1986: locally, the source target is the compactified source-gap fiber.
HYP-1987: globally, reachable source targets form a tiny arc-confined A000568 menu.
THM-385/HYP-1988: observer score is exact blocker count around the source target.
S520: test the n=18 gate split as source-target survival and debt export.
```

The unit skeleton behaves exactly as hoped.  The unit points are:

```text
1/18, 5/18, 7/18, 11/18, 13/18, 17/18.
```

Only residue `0 mod 18` can cover them.  Without an `18`-multiple, all six
remain safe.  In HYP-1981 language, they are six observer-source boundary
targets.  So the no-gate branch is dead as a counterexample branch, just as in
the n=14 story.

The has-gate branch is where the debt goes.  The initial row is tight:

```text
initial: boundary_only, unit_safe=6, unprotected=6.
```

Inserting an `18`-gate kills the unit witnesses but makes a positive gap:

```text
8 -> 18:  gap/th=0.051136
17 -> 18: gap/th=0.033333
```

The S383 lpd ladder gives the cleanest shape:

```text
9-ladder:        gap/th=1/176, unprotected=176, first=11/162
18-gate ladder:  gap/th=1/352, unprotected=352, first=19/324
36-double gate:  gap/th=1/704, unprotected=704, first=37/648
```

This is delightfully rigid: double the gate scale, halve the normalized gap,
double the endpoint debt, and double the first descendant denominator.  It is
not an open-cover construction; it is an endpoint-debt export mechanism.

The targeted one-step repairs around the fragile coordinates

```text
6, 8, 9, 12, 16, 17
```

found only positive-gap rows.  The no-gate repairs keep the six unit witnesses
alive.  The gate repairs kill those witnesses but move the exposed endpoints
to denominators like `306`, `324`, `648`, and `1944`.

The branch tournament is transitive:

```text
H=1, c3=0.
```

So the sampled n=18 branch structure is not a cyclic repair menu.  It is a
debt ladder.  The next proof object should not be another speed search; it
should be an observer-source endpoint-debt certificate for the `9 -> 18 -> 36`
chain, split into the even gate and the `3^2` torsion layers.

That certificate should be phrased as a marked-walk statement: either the
clock reaches one of the threshold-decorated observer-source targets, or the
remaining endpoint debt contains a THM-380-style owner-compatible pressure
core.  S520 found the first side everywhere it looked and no cyclic repair
shape on the second side.
