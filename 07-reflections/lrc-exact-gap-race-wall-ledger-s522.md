# Exact gap races are wall ledgers (S522)

The S518 race model had the right shape but the wrong emphasis.  It asked
whether the right gap reaches `1/n` before the left gap falls below `1/n`,
approximating the race by the resetter speed and the nearest left runner.
S522 follows the actual wall sequence after each reset.

The finding is clean:

```text
closed-LL full audit == exact reset race reaches wall_LL
```

in the bounded `n=3..6` windows, with zero mismatches.  Every successful first
hit is a compactified wall hit, not an open-cell hit.  That is not a weakness;
it is exactly what THM-384/THM-387 predict.  To enter an open LL cell, the
trajectory first touches the threshold wall where the right gap becomes
`1/n` while the left gap is still long.

So the proof object is a reset-wall ledger:

1. A runner wraps through the observer.
2. The system starts an LS race.
3. Before the next reset or a left-gap loss, some threshold wall is hit with
   both adjacent gaps long.

Initial segments are the model case, not a nuisance.  They have zero open LL
measure through the audit, but the wall ledger still fires.  This makes the
compactification from THM-383 structural rather than cosmetic.

The Tournament Analysis over selected rows is transitive.  Gate rows with
open LL measure rank above the wall-only initial rows, but there is no cyclic
repair menu in the exact race strength vector.  The next step should be
arithmetical: classify reset walls by resetter residue, left-runway source,
and the threshold equation that closes the race.
